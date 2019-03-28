#define TMB_LIB_INIT R_init_fit
#include <TMB.hpp>
#include "pnorm4.hpp"
#include "Fcorr.hpp"
#include "MVMIX.hpp"

//Convert a matrix to it's proportions by row
template<class Type>
matrix<Type>  propM(matrix <Type > M){
  int A = M.cols();
  int Y = M.rows();

  matrix<Type> retM(Y,A);
  
  for(int y = 0;y < Y;y++){
    vector<Type> cy = M.row(y);
    Type tot = cy.sum();
    vector<Type> prop = cy/tot;
    retM.row(y) = prop;
  }
  
  return retM;
}

/*Calculate population weighted average F for some set of ages
 bAge is the beginning age to calculate average F over, eAge is the end Age.*/
template<class Type>
vector<Type> aveF(matrix<Type> F, matrix<Type> N, int bAge,int eAge){
  int Y = F.rows();
  int dist = (eAge-bAge)+1;
  matrix<Type> aveF = F.array()*N.array();
  //Get the block of the matrix starting at year 0 going to year Y, column of bAge extending the number of columns between eAge and bAge and take the sum rowwise.
  vector<Type> aveF_XY = aveF.block(0,bAge,Y,dist).rowwise().sum();
  vector<Type> tn = N.block(0,bAge,Y,dist).rowwise().sum();
  aveF_XY = aveF_XY/tn;
  return aveF_XY;
  
}


/*Calculate CRLs from catch proportions*/
template<class Type>
matrix<Type> makeCRLs(matrix <Type > cProps){
  int A = cProps.cols();
  int Y = cProps.rows();
  
  matrix<Type> crls(Y,A-1);

  for(int a = 0; a < A-1;a++){
    for(int y = 0; y < Y;y++){
      vector<Type> cprow = cProps.row(y);
      //get p_a+...+p_A
      Type denom = cprow.segment(a,A-a).sum();
      Type num = cProps(y,a);
      Type pi = num/denom;
      crls(y,a) = log(pi/(1-pi));
    }
  }

  return crls;
}

/*Beverton-Holt Recruitment function*/
template<class Type>
Type bHolt(Type SSBly,vector <Type> bhParm){
  Type ret = bhParm(0) + log(SSBly) -log(1.0+exp(bhParm(1))*SSBly);
  return ret;
}

/*Ricker recruitment function*/
template<class Type>
Type ricker(Type SSBly,vector <Type> rickParm){
  Type ret = rickParm(0) + log(SSBly) - exp(rickParm(1))*SSBly;
  return ret;
}

//' Smooth Hockey Stick Recruitment Function
template<class Type>
Type sHS(Type SSBly,vector <Type> segParm,Type gammaSq){
  //Alpha and delta must be postive.
  Type delta = exp(segParm(1));
  Type alpha = exp(segParm(0));
  
  Type gm = gammaSq/4.0;
  Type sq1 = sqrt(delta*delta + gm);
  Type sq2 = sqrt((SSBly-delta)*(SSBly-delta) +gm);
  Type ret = log(alpha*(SSBly+sq1-sq2));
  return ret;
}

  
template<class Type>
Type objective_function<Type>::operator() ()
{
  enum fleetType {
    survey = 0,
    Catch = 1,
    land = 2
  };

  enum recType {
    rw = 0,
    bh = 1,
    rick = 2,
    smoothHS = 3
  };
    
  
  // input data;  
  DATA_MATRIX(M);
  DATA_MATRIX(weight); 
  DATA_MATRIX(mat); 
  DATA_MATRIX(midy_weight);
  DATA_VECTOR(log_index);
  DATA_IVECTOR(i_zero); //A lot of this could be put in one IMATRIX?
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_IVECTOR(isurvey);
  DATA_IVECTOR(iq);
  DATA_VECTOR(fs);
  DATA_IVECTOR(ft);
  DATA_VECTOR(lowerMult); //Lower multiplier for censored bounds
  DATA_VECTOR(upperMult); //Upper multiplier for censored bounds
  DATA_IVECTOR(keyF);
  DATA_INTEGER(recflag);
  DATA_INTEGER(corflag);
  DATA_INTEGER(corflagCRL);
  DATA_SCALAR(gammaSq); //For smoothHS recruitment function, controls smoothness of breakpoint
  DATA_VECTOR(crlVec); //Vector of Continuation Ratio Logits for catch
  DATA_IVECTOR(aveFages); //Ages in model ages to calculate average F over...0 bAge, 1 for eAge

  
  DATA_VECTOR_INDICATOR(keep,log_index); //Used for one step predict, does not need to be read in from R
  DATA_VECTOR_INDICATOR(keepCRL,crlVec);
  
  int n = log_index.size();
  
  int A = mat.cols();
  int Y = mat.rows();
  vector<Type> index = exp(log_index);
  
  //define parameters;  
  PARAMETER(log_std_log_R);
  PARAMETER_VECTOR(log_qparm);
  PARAMETER_VECTOR(log_std_index); 
  PARAMETER_VECTOR(log_std_logF);          
  PARAMETER_VECTOR(log_std_CRL);
  PARAMETER_VECTOR(log_std_landings); //std. dev for landings, size 0 if fit_land=0, else 1
  PARAMETER(log_sdS);
  PARAMETER_VECTOR(rec_parm);
  PARAMETER_VECTOR(tRhoF);
  PARAMETER_VECTOR(tRhoCRL);

  //Random Effects
  PARAMETER_MATRIX(log_F);    
  PARAMETER_MATRIX(log_N);
  
  Type std_log_R = exp(log_std_log_R); 
  vector<Type> std_index = exp(log_std_index);
  vector<Type> std_logF = exp(log_std_logF); 
  vector<Type> std_CRL = exp(log_std_CRL);
  vector<Type> std_landings = exp(log_std_landings);
  Type sdS = exp(log_sdS);
    
  matrix<Type> F(Y,A);    
  matrix<Type> Z(Y,A);  
  matrix<Type> EC(Y,A);  
  matrix<Type> ECW(Y,A); 
  
  vector<Type> Elog_index(n); 
  vector<Type> std_index_vec(n);
     
  //**********  start the engine ***************;

  using namespace density;
  Type nll = 0.0;  
  
  //compute F & Z
  for(int a = 0;a < A;++a){ 
    for(int y = 0;y < Y;++y){
      //keyF really cuts down on number of random parameters
      F(y,a) = exp(log_F(y,keyF(a)));
      Z(y,a) = F(y,a) + M(y,a);
    }
  }

  matrix<Type> N = log_N.array().exp();

  //Because of the magic from making log_N a parameter you can calculate SSB whenever you want.
  matrix<Type> B_matrix = weight.array()*N.array();
  matrix<Type> SSB_matrix = mat.array()*B_matrix.array();
  vector<Type> biomass = B_matrix.rowwise().sum();
  vector<Type> ssb = SSB_matrix.rowwise().sum();
  vector<Type> log_biomass = log(biomass);
  vector<Type> log_ssb = log(ssb);       


  //crud for process residuals, could probably reuse sigmaFgen again but eh.
  array<Type> resN(A,Y-1);
  matrix<Type> nSig(A,A);
  for(int i = 0; i < A;i++){
    for(int j = 0; j < A;j++){
      if(i != j){
	nSig(i,j) = 0;
      }else{
	if(i == 0){
	  nSig(i,j) = std_log_R*std_log_R;
	}else{
	  nSig(i,j) = sdS*sdS;
	}
      }
    }
  }

  //This is the approach SAM takes to get joint sample process residuals
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovN(nSig);
  matrix<Type> LN = lltCovN.matrixL();
  matrix<Type> LinvN = LN.inverse();
  

  //Recruitment, adding recruitment curves much easier because ssb is already available and will be optimized later
  matrix<Type> mpredN(Y,A);
  Type predN;
  for(int y =1;y < Y;y++){
    switch(recflag){
    case rw:
      predN = log_N(y-1,0);
      break;
    case bh:
      predN = bHolt(ssb(y),rec_parm);
      break;
    case rick:
      predN = ricker(ssb(y),rec_parm);
      break;
    case smoothHS:
      predN = sHS(ssb(y),rec_parm,gammaSq);
      break;
    default:
      error("Not right rec type");
      break;
    }
    mpredN(y,0) = predN;
    nll -= dnorm(log_N(y,0),predN,std_log_R,true);
  }

      
  //Fill in N, add survival process
  for(int y = 1;y < Y;++y){
    for(int a = 1;a < A;++a){
      predN = log_N(y-1,a-1)-Z(y-1,a-1);
      if(a == A-1){//Plus group
	predN = logspace_add(predN,log_N(y-1,a)-Z(y-1,a));
      }
      mpredN(y,a) = predN;
      nll -= dnorm(log_N(y,a),predN,sdS,true);
    }
    vector<Type> Nrow = log_N.row(y);
    vector<Type> vpredN = mpredN.row(y);
    resN.col(y-1) = LinvN*(vector<Type>(Nrow-vpredN));
  }
  ADREPORT(resN);
  
  //F correlation matrix making, sigmaFgen lives in Fcorr.hpp
  matrix<Type> sigmaF = sigmaFgen(log_F.cols(),corflag,log_std_logF,tRhoF);
  
  array<Type> resF(log_F.cols(),Y-1);
  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovF(sigmaF);
  matrix<Type> LF = lltCovF.matrixL();
  matrix<Type> LinvF = LF.inverse(); 


  //Add the F random walk to the likelihood.
  density::MVNORM_t<Type> Fdens(sigmaF);
  for(int y = 1; y < Y;y++){
    resF.col(y-1) = LinvF*(vector<Type>(log_F.row(y)-log_F.row(y-1)));  
    nll += Fdens(log_F.row(y)-log_F.row(y-1));
  }
  ADREPORT(resF);

  //Expected catch, for making expected catch CRLs
  //Expected landings too
  for(int y = 0;y <Y;y++){
    for(int a = 0;a <A;a++){
      EC(y,a) = N(y,a)*(1.0-exp(-1.0*Z(y,a)))*F(y,a)/Z(y,a);
      ECW(y,a) = EC(y,a)*midy_weight(y,a);
    }
  }

  vector<Type> landings_pred = ECW.rowwise().sum();
  vector<Type> total_expected_catch = EC.rowwise().sum();
  matrix<Type> log_EC = EC.array().log();

  vector<Type> log_landings_pred = log(landings_pred);
  vector<Type> log_total_expected_catch = log(total_expected_catch);

  REPORT(log_landings_pred);
  REPORT(log_total_expected_catch);
  
  // Landings, Survey index predictions
  int ia,iy;
  for(int i = 0;i < n;++i){
    int fType = ft(i);
    ia = iage(i);
    iy = iyear(i);

    //Observation equations for survey and landings
    switch(fType){
    case survey:
      std_index_vec(i) = std_index(isurvey(i));
      Elog_index(i) = log_qparm(iq(i)) + log_N(iy,ia) - fs(i)*Z(iy,ia);
      if(i_zero(i) == 1){ //if zero use left-censor bound
	nll -= keep(i)*pnorm4(log_index(i),Elog_index(i),std_index_vec(i),true);
      }      
      else{
	nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
      }
      break;
    case land:
      Elog_index(i) = log_landings_pred(iy);
      std_index_vec(i) = std_landings(0);
      nll -= keep(i)*censored_bounds(log_index(i),Elog_index(i),std_index_vec(i),-log(lowerMult(iy)),log(upperMult(iy)));
      break;
    default:
      error("Fleet Type not implemented");
      break;
    }
    
  }
  
  matrix<Type> EPC = propM(EC);  //Expected Catch Proportions
  matrix<Type> ECRL = makeCRLs(EPC); //Expected Catch Continuation Ratio Logits

  matrix<Type> S = propM(F); //Selectivity F_{y,a}/(sum_{a}(F_{y,a}))

  REPORT(S);
  REPORT(EPC);
  REPORT(ECRL);

  //Correlation matrix for CRLs, might as well reuse sigmaFgen...
  matrix<Type> sigmaCRL = sigmaFgen(ECRL.cols(),corflagCRL,log_std_CRL,tRhoCRL);
  REPORT(sigmaCRL);

  //MV NORMAL MIXTURE distribution, second argument is the mixture probability, grabbed from SAM lives in MVMIX.hpp
  //MVNORM_t in TMB doesn't support OSA residuals yet, MVMIX_t from SAM does however there is a slight speed loss
  MVMIX_t<Type> CRLdens(sigmaCRL,0); //cause this is a MV MIXTURE duh. I could maybe actually try to use this part?
  for(int y = 0;y < Y;y++){
    //IT'S VERY VERY IMPORTANT CATCH CRLs ARE SORTED YEAR THEN AGE HERE.
    vector<Type> ECRLr = ECRL.row(y);
    nll += CRLdens(crlVec.segment(y*(A-1),(A-1))-ECRLr,keepCRL.segment(y*(A-1),(A-1)));
  }
      
      
    
  //pop size weighted ave F;  

  vector<Type> aveFXY = aveF(F,N,aveFages(0),aveFages(1));
  vector<Type> log_aveFXY = log(aveFXY);
  
  REPORT(std_logF);      
  REPORT(std_log_R);  
  REPORT(std_index);  
  
  REPORT(N);          
  REPORT(F);                  
  REPORT(Z);                  
  REPORT(B_matrix);             
  REPORT(SSB_matrix);                
  REPORT(biomass);                  
  REPORT(ssb);            
  REPORT(aveFXY);              
  
  REPORT(EC);                          
  REPORT(ECW);        
  
  REPORT(landings_pred); 
  
  REPORT(Elog_index);
  
  REPORT(log_F);                    
  
  REPORT(log_std_index); 
  REPORT(log_qparm);
  
  ADREPORT(log_landings_pred);     
  ADREPORT(log_biomass);           
  ADREPORT(log_ssb);   
  ADREPORT(log_aveFXY);  
  ADREPORT(log_qparm);
  
  ADREPORT(std_logF);      
  ADREPORT(std_log_R);
  ADREPORT(std_index);  
  
  return nll;
}

