#include <TMB.hpp>
#include "pnorm4.hpp"

/*Convert catch to proportion of catch at age*/
 template<class Type>
 vector<Type>  catchProp(vector <Type > catchy){
   Type tot = catchy.sum();
   vector<Type> prop = catchy/tot;
   return prop;
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

      
  

template<class Type>
Type objective_function<Type>::operator() ()
{
  enum censorType {
    noCensor = 0,
    basic = 1,
    stable = 2,
    bounds = 3 // implemented :)
  };

  enum fleetType {
    survey = 0,
    Catch = 1,
    land = 2
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
  DATA_INTEGER(index_censor); 
  DATA_INTEGER(use_pe);        
  DATA_INTEGER(use_cye);
  DATA_INTEGER(fit_land); //Fit landings or fit C@A?
  DATA_VECTOR(lowerMult); //Lower multiplier for censored bounds
  DATA_VECTOR(upperMult); //Upper multiplier for censored bounds
  DATA_INTEGER(use_cb); //option to use censored bounds on landings
  
  DATA_VECTOR_INDICATOR(keep,log_index); //Used for one step predict, does not need to be read in from R
  
  
  int n = log_index.size();  
  Type one = 1.0;
  Type zero = 0.0;
  
  int A = mat.cols();
  int Y = mat.rows();
  vector<Type> index = exp(log_index);
  
  //define parameters;  
  PARAMETER_VECTOR(log_No);  
  PARAMETER(log_Rec_mean);    
  PARAMETER(log_std_log_R);
  PARAMETER_VECTOR(log_qparm);
  PARAMETER_VECTOR(log_std_index); 
  PARAMETER_VECTOR(log_std_logF);          
  PARAMETER(log_std_pe);                  
  PARAMETER(log_std_cye);    
  PARAMETER(logit_ar_logF_year);        
  PARAMETER(logit_ar_logF_age);     
  PARAMETER(logit_ar_pe_year);    
  PARAMETER(logit_ar_pe_age);       
  PARAMETER(logit_ar_cye_year); 
  PARAMETER_VECTOR(log_std_log_C);
  PARAMETER_VECTOR(log_std_CRL);
  PARAMETER_VECTOR(log_std_landings); //std. dev for landings, size 0 if fit_land=0, else 1

  //Random Effects
  PARAMETER_VECTOR(log_Rec_dev);   
  PARAMETER_ARRAY(log_F);    
  PARAMETER_ARRAY(pe);      
  PARAMETER_MATRIX(cye);       
  
  Type std_log_R = exp(log_std_log_R); 
  Type std_pe = exp(log_std_pe);      
  Type std_cye = exp(log_std_cye);         
  vector<Type> std_index = exp(log_std_index);
  vector<Type> std_logF = exp(log_std_logF); 
  vector<Type> std_log_C = exp(log_std_log_C);
  vector<Type> std_CRL = exp(log_std_CRL);
  vector<Type> std_landings = exp(log_std_landings);
  vector<Type> llog_qparm = exp(log_qparm)/(one+exp(log_qparm));
  
  Type ar_logF_age = exp(logit_ar_logF_age)/(one + exp(logit_ar_logF_age));    
  Type ar_logF_year = exp(logit_ar_logF_year)/(one + exp(logit_ar_logF_year)); 
  Type ar_pe_year = exp(logit_ar_pe_year)/(one + exp(logit_ar_pe_year)); 
  Type ar_pe_age = exp(logit_ar_pe_age)/(one + exp(logit_ar_pe_age)); 
  Type ar_cye_year = exp(logit_ar_cye_year)/(one + exp(logit_ar_cye_year)); 
  
  matrix<Type> log_N(Y,A); 
  matrix<Type> N(Y,A);     
  matrix<Type> F(Y,A);    
  matrix<Type> Z(Y,A);  
  matrix<Type> EC(Y,A);  
  matrix<Type> log_EC(Y,A); 
  matrix<Type> ECW(Y,A); 
  matrix<Type> C_resid(Y,A); 
  matrix<Type> C_resid_std(Y,A);     
  matrix<Type> std_C_resid(Y,A);
  matrix<Type> B_matrix(Y,A);           
  matrix<Type> SSB_matrix(Y,A);
  matrix<Type> C(Y,A); //Yes a matrix of catch not read in
  matrix<Type> PC(Y,A);
  matrix<Type> EPC(Y,A);
  matrix<Type> CRL(Y,A-1);
  matrix<Type> ECRL(Y,A-1);
  matrix<Type> EPoC(Y,A);
  
  vector<Type> Elog_index(n); 
  vector<Type> resid_index(n); 
  vector<Type> std_resid_index(n); 
  vector<Type> std_index_vec(n);
   
  //**********  SD report objects ***************;
  
  vector<Type> biomass(Y); 
  vector<Type> log_biomass(Y);  
  vector<Type> ssb(Y); 
  vector<Type> log_ssb(Y); 
  vector<Type> aveF_46(Y);
  vector<Type> log_aveF_46(Y);    
  vector<Type> aveF_69(Y);
  vector<Type> log_aveF_69(Y);   
  
  //**********  start the engine ***************;
  
  using namespace density;
  Type nll = 0.0;  
  
  //compute F 
  for(int a = 0;a < A;++a){ 
    for(int y = 0;y < Y;++y){
      F(y,a) = exp(std_logF(a)*log_F(y,a)); 
      Z(y,a) = F(y,a) + M(y,a);
    }
  }
  
  // The cohort model;
  
  vector<Type> log_Rec = log_Rec_dev + log_Rec_mean;
  log_N(0,0) = log_Rec(0);
  for(int a = 1;a < A;++a){
    log_N(0,a) = log_No(a-1);
  }
  
  for(int y = 1;y < Y;++y){
    log_N(y,0) = log_Rec(y);  
    for(int a = 1;a < A;++a){
      log_N(y,a) = log_N(y-1,a-1) - Z(y-1,a-1) + std_pe*pe(y,a);
      if(a == A-1){//Plus group
	log_N(y,a) = logspace_add(log_N(y,a),log_N(y-1,a) - Z(y-1,a)); //log(exp(log x) + exp(log y))
      }
    }
  }
  N = log_N.array().exp();
  
  B_matrix = weight.array()*N.array();
  SSB_matrix = mat.array()*B_matrix.array(); 
  
  //  Catch, Landings, Survey index predictions
  //This was done entirely for the sake of oneStepPredict as you need it all in a vector
  //The easiest way I could think to do this was to prepare the catch/landings data in the same fashion as the survey data
  int ia,iy;
  vector<Type> landings_pred(Y);
  for(int i = 0;i < n;++i){
    int fType = ft(i);
    ia = iage(i);
    iy = iyear(i);

    //Observation equations
    switch(fType){
    case survey:
      std_index_vec(i) = std_index(isurvey(i));
      Elog_index(i) = log_qparm(iq(i)) + log_N(iy,ia) - fs(i)*Z(iy,ia);
      resid_index(i) = log_index(i) - Elog_index(i); 
      std_resid_index(i) = resid_index(i)/std_index_vec(i);
      break;
    case Catch:
      EC(iy,ia) = N(iy,ia)*(one - exp(-one*Z(iy,ia)))*F(iy,ia)/Z(iy,ia);
      ECW(iy,ia) = EC(iy,ia)*midy_weight(iy,ia);
      C(iy,ia) = index(i);
      log_EC(iy,ia) = log(EC(iy,ia));
      Elog_index(i) = log_EC(iy,ia);
      resid_index(i) = log_index(i) - log_EC(iy,ia) + cye(iy,ia);
      landings_pred(iy) += ECW(iy,ia);
      std_C_resid(iy,ia) = resid_index(i)/std_log_C(ia);
      std_index_vec(i) = std_log_C(ia);
      break;
    case land: // Must handle landings after catch!
      break;
    default:
      error("Fleet Type not implemented");
      break;
    }    
  }
  vector<Type> log_landings_pred = log(landings_pred);

  
  // Catch proportion stuff
  for(int y = 0;y < Y; y++){
    vector<Type> cy = C.row(y);
    vector<Type> cpy = catchProp(cy);
    PC.row(y) = cpy;

    vector<Type> ecy = EC.row(y);
    vector<Type> ecpy = catchProp(ecy);
    EPC.row(y) = ecpy;
    
  }
  
  CRL = makeCRLs(PC);
  ECRL = makeCRLs(EPC);
  
  REPORT(PC);
  REPORT(EPC);
  REPORT(CRL);
  REPORT(ECRL); 

  //Ensure landings are handled after catch by doing a seperate loop
  //DO NLL if fit_land == 1 TOO
  if(fit_land == 1){
    for(int i = 0;i < n;i++){
      int fType = ft(i);
      iy = iyear(i);
      ia = iage(i);
    
      if(fType == Catch){
	if(ia < A-1){ //CRL doesn't exist for A, so for ia = A-1 we don't want to add anything to nll
	  log_index(i) = CRL(iy,ia);
	  Elog_index(i) = ECRL(iy,ia);
	  std_index_vec(i) = std_CRL(ia);
	}
	  
      }
    
      if(fType == land){
	Elog_index(i) = log_landings_pred(iy);
	std_index_vec(i) = std_landings(0);
      }
    

      if(fType != land){
	if(fType != Catch){ 
	  if(i_zero(i) == 1){
	    switch(index_censor){
	    case noCensor:
	      nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
	      break;
	    case basic:
	      nll -= keep(i)*log(pnorm(log_index(i),Elog_index(i),std_index_vec(i)));
	      break;
	    case stable:
	      nll -= keep(i)*pnorm4(log_index(i),Elog_index(i),std_index_vec(i),true);
	      break;
	    default:
	      error("Invalid index_censor type");
	      break;
	    } 
	  }      
	  else{
	    nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
	  }
	}
	else{
	  if(ia < A-1){
	    nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
	  }//DONT DO ANYTHING FOR CATCH PROPORTIONS A
	
	}

      }else if(fType == land){
	switch(use_cb){
	case noCensor:
	  nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true); //Instead of NCAM's fit to the mean fit to obs.
	  break;
	case basic:
	  { //Forgot you can't normally declare vars in a switch, { } are workaround
	    Type upb = log_index(i) + log(upperMult(i));
	    Type lowb = log_index(i) + log(lowerMult(i));
	    Type upbp = pnorm(upb,Elog_index(i),std_index_vec(i));
	    Type lowbp = pnorm(lowb,Elog_index(i),std_index_vec(i));
	    nll -= keep(i)*log(upbp-lowbp);
	    break;
	  }
	case stable: //Not really any better if Z > 40 which easily happens, so for illustration purposes?
	  {
	    Type upb = log_index(i) + log(upperMult(i));
	    Type lowb = log_index(i) + log(lowerMult(i));
	    Type upbp = pnorm4(upb,Elog_index(i),std_index_vec(i),true);
	    Type lowbp = pnorm4(lowb,Elog_index(i),std_index_vec(i),true);
	    nll -= keep(i)*logspace_sub(upbp,lowbp);
	    break;
	  }
	case bounds: //The GOOD one
	  //Yeah I subtract lower bounds in the atomic stuff so lower = -log(lowerMult), bleh
	  //censored_bounds exists in pnorm4.hpp...
	  nll -= keep(i)*censored_bounds(log_index(i),Elog_index(i),std_index_vec(i),-log(lowerMult(i)),log(upperMult(i)));
	  break;
	default:
	  error("Unsupported use_cb type");
	  break;
	}

      }else{
	error("What fType is this?");
      }
    
    }
  }
  
  

  //Fit observations to nll here if not fitting landings
  else if(fit_land == 0){
  for(int i = 0; i < n;i++){
        int fType = ft(i);
	ia = iage(i);
	iy = iyear(i);

	

	//Fit everything but landings, if fitting landings
	if(fType != land){
	  if(i_zero(i) == 1){
	    switch(index_censor){
	    case noCensor:
	      nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
	      break;
	    case basic:
	      nll -= keep(i)*log(pnorm(log_index(i),Elog_index(i),std_index_vec(i)));
	      break;
	    case stable:
	      nll -= keep(i)*pnorm4(log_index(i),Elog_index(i),std_index_vec(i),true);
	      break;
	    default:
	      error("Invalid index_censor type");
	      break;
	    } 
	  }      
	  else{
	    nll -= keep(i)*dnorm(log_index(i),Elog_index(i),std_index_vec(i),true);
	  }
	}


	//"Fit" landings
	else if(fType == land){//do nothing for landings. Man is this stuff gross
	}

	else{
	  error("fType not supported in nll");
	}
	
  }
  }
  else{
    error("fit_land value wrong");
  }
  
  
  
  
  
  //recruitment
  //  nll += SCALE(AR1(phi_logR),std_log_R)(log_Rec_dev);
  nll -= dnorm(log_Rec_dev,zero, std_log_R, true).sum();
  
  //Log_F nll
  //RW on first age;
  vector<Type> del = log_F.col(0);
  nll -= dnorm(del(0),Type(-10.0),one, true);
  for(int y = 1;y < Y;++y){ 
    nll -= dnorm(del(y),del(y-1),one, true);
  }
  array<Type> log_F1(Y,A);
  for(int a = 1;a < A;++a){
    log_F1.col(a-1) = log_F.col(a);
  } 
  //year x age correlation on first+1:last ages;
  nll += SEPARABLE(AR1(ar_logF_age),AR1(ar_logF_year))(log_F1);
  
  //pe nll
  if(use_pe==1){nll += SEPARABLE(AR1(ar_pe_age),AR1(ar_pe_year))(pe);}
  
  //year effect nll;
  if(use_cye==1){
    for(int y = 0;y < Y;++y){
      vector<Type> del = vector<Type>(cye.row(y));
      nll += SCALE(AR1(ar_cye_year),std_cye)(del);
    }
  }       
  
  // create some REPORT output;  
  
  for(int y = 0;y < Y;++y){
    biomass(y) = zero;     
    ssb(y) = zero;
    for(int a = 0;a < A;++a){
      biomass(y) += B_matrix(y,a); 
      ssb(y) += SSB_matrix(y,a);}
  }          
  log_biomass = log(biomass);
  log_ssb = log(ssb);       
  
  //pop size weighted ave F;  
  
  Type tni;
  
  for(int y = 0;y < Y;++y){
    aveF_46(y) = zero; 
    tni = zero;
    for(int a = 2;a < 4;++a){
      aveF_46(y) += F(y,a)*N(y,a); 
      tni += N(y,a);
    }
    aveF_46(y) = aveF_46(y)/tni;  
    aveF_69(y) = zero; 
    tni = zero;
    for(int a = 4;a < 7;++a){
      aveF_69(y) += F(y,a)*N(y,a); 
      tni += N(y,a);
    }
    aveF_69(y) = aveF_69(y)/tni;
  }
  
  log_aveF_46 = log(aveF_46); 
  log_aveF_69 = log(aveF_69);
  
  REPORT(std_log_C); 
  REPORT(std_logF);      
  REPORT(std_cye);       
  REPORT(std_pe);    
  REPORT(std_log_R);  
  REPORT(std_index);  
  REPORT(ar_logF_age);           
  REPORT(ar_logF_year);  
  REPORT(ar_pe_year);  
  REPORT(ar_pe_age);  
  REPORT(ar_cye_year);
  
  REPORT(N);          
  REPORT(F);                  
  REPORT(Z);                  
  REPORT(B_matrix);             
  REPORT(SSB_matrix);                
  REPORT(biomass);                  
  REPORT(ssb);            
  REPORT(aveF_46);              
  REPORT(aveF_69); 
  
  REPORT(EC);                          
  REPORT(C_resid);             
  REPORT(std_C_resid);  
  REPORT(ECW);        
  REPORT(cye);         
  
  REPORT(landings_pred); 
  
  REPORT(Elog_index);                 
  REPORT(resid_index);             
  REPORT(std_resid_index); 
  
  REPORT(log_Rec);
  REPORT(log_Rec_dev);      
  REPORT(log_F);                    
  REPORT(pe);   
  
  REPORT(log_std_index); 
  REPORT(log_qparm);
  
  ADREPORT(log_landings_pred);     
  ADREPORT(log_biomass);           
  ADREPORT(log_ssb);   
  ADREPORT(log_aveF_46);  
  ADREPORT(log_aveF_69);  
  ADREPORT(log_Rec);   
  ADREPORT(log_qparm);
  
  ADREPORT(std_log_C); 
  ADREPORT(std_logF);      
  ADREPORT(std_cye);      
  ADREPORT(std_pe);      
  ADREPORT(std_log_R);
  ADREPORT(std_index);  
  ADREPORT(ar_logF_age);           
  ADREPORT(ar_logF_year);  
  ADREPORT(ar_pe_year);  
  ADREPORT(ar_pe_age);  
  ADREPORT(ar_cye_year);  
  
  return nll;
}

