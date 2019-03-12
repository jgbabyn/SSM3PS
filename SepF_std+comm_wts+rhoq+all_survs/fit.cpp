#include <TMB.hpp>
#include "pnorm4.hpp"

//Simplifed, in my opinion more readable, added stable versions of pnorm...

template<class Type>
Type objective_function<Type>::operator() ()
{
  enum censorType {
    noCensor = 0,
    basic = 1,
    stable = 2,
    bounds = 3 //Not implemented yet
  };
  
  
  // input data;  
  DATA_MATRIX(M);
  DATA_MATRIX(weight); 
  DATA_MATRIX(mat); 
  DATA_MATRIX(midy_weight);
  DATA_MATRIX(comm_wt);
  DATA_MATRIX(C);         
  DATA_VECTOR(landings); 
  DATA_VECTOR(index);
  DATA_VECTOR(log_OFF);
  DATA_VECTOR(diffs);
  DATA_IVECTOR(i_zero); 
  DATA_IMATRIX(C_zero);  
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_IVECTOR(isurvey); 
  DATA_IVECTOR(iq);
  DATA_IVECTOR(irho);
  DATA_VECTOR(fs);
  DATA_INTEGER(index_censor); 
  DATA_INTEGER(catch_censor); 
  DATA_INTEGER(use_pe);        
  DATA_INTEGER(use_cye);
  //DATA_SCALAR(i_detect); //Why read a vector, when a single value will do?
  //DATA_SCALAR(c_detect);
  
  int n = index.size();  
  Type one = 1.0;
  Type zero = 0.0;
  
  int A = mat.cols();
  int Y = mat.rows();
  matrix<Type> log_C = C.array().log(); //we can do this in a line inside too
  vector<Type> log_landings = log(landings);
  vector<Type> log_index = log(index);
  
  //define parameters;  
  PARAMETER_VECTOR(log_No);  
  PARAMETER(log_Rec_mean);    
  PARAMETER(log_std_log_R);
  PARAMETER_VECTOR(log_qparm);
  PARAMETER_VECTOR(log_rho);
  PARAMETER_VECTOR(log_std_index); 
  PARAMETER_VECTOR(log_std_rho); 
  //PARAMETER_VECTOR(log_std_ratio); 
  PARAMETER_VECTOR(log_std_logF);  
  PARAMETER_VECTOR(log_std_logF_new);  
  PARAMETER(log_std_pe);                  
  PARAMETER(log_std_cye);    
  PARAMETER(logit_ar_logF_year);        
  PARAMETER(logit_ar_logF_age);     
  PARAMETER(logit_ar_pe_year);    
  PARAMETER(logit_ar_pe_age);       
  PARAMETER(logit_ar_cye_year); 
  PARAMETER_VECTOR(log_std_log_C);             
  
  PARAMETER_VECTOR(log_Rec_dev);   
  PARAMETER_ARRAY(log_F);    
  PARAMETER_ARRAY(pe);      
  PARAMETER_MATRIX(cye);       
  
  Type std_log_R = exp(log_std_log_R); 
  Type std_pe = exp(log_std_pe);      
  Type std_cye = exp(log_std_cye);         
  vector<Type> std_index = exp(log_std_index);
  vector<Type> std_rho = exp(log_std_rho);
  //vector<Type> std_ratio = exp(log_std_ratio);
  vector<Type> std_logF = exp(log_std_logF); 
  vector<Type> std_logF_new = exp(log_std_logF_new); 
  vector<Type> std_log_C = exp(log_std_log_C);
  
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
  
  vector<Type> Elog_index(n);
  //vector<Type> Elog_ratio(n);
  vector<Type> resid_index(n); 
  vector<Type> resid_rho(n); 
  vector<Type> std_resid_index(n); 
  //vector<Type> std_resid_ratio(n); 
  vector<Type> log_q_vec = log_qparm(iq); 
  //vector<Type> new_log_q= log_qparm(iq); 
  vector<Type> std_index_vec = std_index(isurvey);
  //vector<Type> std_ratio_vec = std_ratio(isurvey);
  
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
     if(y<35){F(y,a) = exp(std_logF(a)*log_F(y,a));} 
     if(y>=35){F(y,a) = exp(std_logF_new(a)*log_F(y,a));}
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
    }
  }
  N = log_N.array().exp();
  
  B_matrix = weight.array()*N.array();
  SSB_matrix = mat.array()*B_matrix.array(); 
  
  // Catch model prediction from Baranov;
  
  vector<Type> landings_pred(Y);
  
  for(int y = 0;y < Y;++y){ 
    landings_pred(y)=zero;
    for(int a = 0;a < A;++a){
      EC(y,a) = N(y,a)*(one - exp(-one*Z(y,a)))*F(y,a)/Z(y,a);
      ECW(y,a) = EC(y,a)*comm_wt(y,a);
      log_EC(y,a) = log(EC(y,a));
      C_resid(y,a) = log_C(y,a) - log_EC(y,a) + cye(y,a);
      landings_pred(y) += ECW(y,a); 
      std_C_resid(y,a) = C_resid(y,a)/std_log_C(a);
    }}
  vector<Type> log_landings_pred = log(landings_pred);
  
  //  Survey index predictions, and residuals;
  

  //for(int i=0;i < n;++i){
  //  ia = iage(i);
  //  Elog_ratio(i) = log_index(i)-log_OFF(i);
  //  resid_ratio(i) = log_ratio(ia)-Elog_ratio(i);
  //  std_resid_ratio(i) = resid_ratio(i)/std_ratio_vec(i);
  //}
  
  //new_log_q<-log_q_vec*log_ratio;
  int ia,iy;   
  for(int i = 0;i < n;++i){
    ia = iage(i);
    iy = iyear(i);
    //log_rho(ia) = diffs(i);
    if(irho(i)==1){Elog_index(i) = log_rho(ia) + log_q_vec(i) + log(N(iy,ia)) - fs(i)*Z(iy,ia);}
    if(irho(i)==0){Elog_index(i) = log_q_vec(i) + log(N(iy,ia)) - fs(i)*Z(iy,ia);}
    resid_index(i) = log_index(i) - Elog_index(i); 
    std_resid_index(i) = resid_index(i)/std_index_vec(i);
  }
  
  for(int i = 0;i < n;++i){
    ia = iage(i);
    if(irho(i)==1){resid_rho(i) = diffs(i)-log_rho(ia);} 
    //std_resid_rho(i) = resid_rho(i)/std_rho_vec(i);
  }
  // End of model, now fit functions (nll - negative loglikelihoods);
  
  // Catch at age nll;
  
  for(int y = 0;y < Y;++y){
    for(int a = 0;a < A;++a){
      if(C_zero(y,a) == 1){
        switch(catch_censor){
        case noCensor:
          nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
          break;
        case basic:
          nll -= log(pnorm(std_C_resid(y,a)));
          break;
        case stable:
          nll -= pnorm4(std_C_resid(y,a));
          break;
        default:
          error("catch_censor type not implemented");
        break;
        }
      }
      else{
        nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
      }
      
      // if(catch_censor==0){
      // 	nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
      // }
      // if(catch_censor==1){  
      //   if(C_zero(y,a)==0){
      // 	  nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
      // 	}
      //   if(C_zero(y,a)==1){
      // 	  nll -= log(pnorm(std_C_resid(y,a)));
      // 	}
      // } 
    }
  }
  
  // Index nll;
  
  
 // for(int i =0; i < n;++i){
//    nll-=dnorm(resid_ratio(i), zero, std_ratio_vec(i), true);
//  }
  
  
  for(int i = 0;i < n;++i){   
    if(i_zero(i) == 1){
      switch(index_censor){
      case noCensor:
        nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);
        break;
      case basic:
        nll -= log(pnorm(std_resid_index(i)));
        break;
      case stable:
        nll -= pnorm4(std_resid_index(i));
        break;
      default:
        error("Invalid index_censor type");
      break;
      }	
      // if(index_censor==0){
      // 	nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);
      // }       
      // if(index_censor==1){
      // 	nll -= log(pnorm(std_resid_index(i)));
      // }
    }      
    else{
      nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);
        }
    }
  
  for(int i = 0;i < n;++i){ 
  ia = iage(i);
  if(irho(i)==1){nll-=dnorm(resid_rho(i),zero,std_rho(ia), true);}
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
  array<Type> log_F1(Y,A-1);
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
  REPORT(std_logF_new);
  REPORT(std_cye);       
  REPORT(std_pe);    
  REPORT(std_log_R);  
  REPORT(std_index);
  REPORT(std_rho);
  //REPORT(std_ratio);  
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
  
  REPORT(log_rho);
  REPORT(resid_rho); 
  //REPORT(Elog_ratio);                 
  //REPORT(resid_ratio);             
  //REPORT(std_resid_ratio); 
  
  REPORT(log_Rec);
  REPORT(log_Rec_dev);      
  REPORT(log_F);                    
  REPORT(pe);   
  
  REPORT(log_std_index);
  //REPORT(log_std_ratio);
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
  ADREPORT(std_logF_new);
  ADREPORT(std_cye);      
  ADREPORT(std_pe);      
  ADREPORT(std_log_R);
  ADREPORT(std_index);
  ADREPORT(log_rho);
  //ADREPORT(std_ratio);
  ADREPORT(ar_logF_age);           
  ADREPORT(ar_logF_year);  
  ADREPORT(ar_pe_year);  
  ADREPORT(ar_pe_age);  
  ADREPORT(ar_cye_year);  
  
  return nll;
}