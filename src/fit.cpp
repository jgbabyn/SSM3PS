#include <TMB.hpp>
#include "pnorm4.hpp"

//The Marie Kondo version

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // input data;  
  DATA_MATRIX(M);
  DATA_MATRIX(weight); 
  DATA_MATRIX(mat); 
  DATA_MATRIX(midy_weight);
  DATA_MATRIX(C);
  DATA_VECTOR(landings); 
  DATA_VECTOR(index);
  //DATA_IVECTOR(i_zero); //I reallly feel this two could just be replaced by a couple checks
  //DATA_IMATRIX(C_zero); //Well let's tidy up then! 
  DATA_SCALAR(i_detect); //The dection limit of the survey
  DATA_SCALAR(c_detect); //catch detection limit 
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_IVECTOR(isurvey); 
  DATA_IVECTOR(iq);
  DATA_VECTOR(fs);
  DATA_INTEGER(index_censor); 
  DATA_INTEGER(catch_censor); 
  DATA_INTEGER(use_pe);        
  DATA_INTEGER(use_cye); 

  //Let's reduced the risk of bombing by cleaning up a bit?
  int n = index.size();
  int A = M.cols();
  int Y = M.rows();
  matrix<Type> log_C = C.array().log();
  vector<Type> log_landings = landings.log();
  vector<Type> log_index = index.log();
  
  Type one = 1.0;
  Type zero = 0.0;
  
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
  
  Type ar_logF_age = exp(logit_ar_logF_age)/(one + exp(logit_ar_logF_age));    
  Type ar_logF_year = exp(logit_ar_logF_year)/(one + exp(logit_ar_logF_year)); 
  Type ar_pe_year = exp(logit_ar_pe_year)/(one + exp(logit_ar_pe_year)); 
  Type ar_pe_age = exp(logit_ar_pe_age)/(one + exp(logit_ar_pe_age)); 
  Type ar_cye_year = exp(logit_ar_cye_year)/(one + exp(logit_ar_cye_year)); 
  
  matrix<Type> log_N(Y,A); 
  matrix<Type> N(Y,A);     
  matrix<Type> F(Y-1,A);    
  matrix<Type> Z(Y-1,A);  
  matrix<Type> EC(Y-1,A);  
  matrix<Type> log_EC(Y-1,A); 
  matrix<Type> ECW(Y-1,A); 
  matrix<Type> C_resid(Y-1,A); 
  matrix<Type> C_resid_std(Y-1,A);     
  matrix<Type> std_C_resid(Y-1,A);
  matrix<Type> B_matrix(Y,A);           
  matrix<Type> SSB_matrix(Y,A); 
  
  vector<Type> Elog_index(n); 
  vector<Type> resid_index(n); 
  vector<Type> std_resid_index(n); 
  vector<Type> log_q_vec = log_qparm(iq); 
  vector<Type> std_index_vec = std_index(isurvey);
  
  //**********  SD report objects ***************;
  
  vector<Type> biomass(Y); 
  vector<Type> log_biomass(Y);  
  vector<Type> ssb(Y); 
  vector<Type> log_ssb(Y); 
  vector<Type> aveF_46(Y-1);
  vector<Type> log_aveF_46(Y-1);    
  vector<Type> aveF_69(Y-1);
  vector<Type> log_aveF_69(Y-1);   
  
  //**********  start the engine ***************;
  
  using namespace density;
  Type nll = 0.0;  
  
  //compute F 
  for(int a = 0;a < A;++a){ 
    for(int y = 0;y < Y-1;++y){
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
    for(int a = 1; a < A;++a){
      log_N(y,a) = log_N(y-1,a-1) - Z(y-1,a-1) + std_pe*pe(y,a);
    }
  }
  //Clean up an entire for loop!
  N = log_N.array().exp();
  B_matrix = weight.array()*N.array();
  SSB_matrix = mat.array()*B_matrix.array(); 
  
  // Catch model prediction from Baranov;
  
  vector<Type> landings_pred(Y-1);
  
  for(int y = 0;y < Y-1;++y){ 
    landings_pred(y)=zero;
    for(int a = 0;a < A;++a){
      EC(y,a) = N(y,a)*(one - exp(-one*Z(y,a)))*F(y,a)/Z(y,a);
      ECW(y,a) = EC(y,a)*midy_weight(y,a);
      log_EC(y,a) = log(EC(y,a));
      C_resid(y,a) = log_C(y,a) - log_EC(y,a) + cye(y,a);
      landings_pred(y) += ECW(y,a); 
      std_C_resid(y,a) = C_resid(y,a)/std_log_C(a);
    }}
  vector<Type> log_landings_pred = log(landings_pred);
  
  //  Survey index predictions, and residuals;
  
  int ia,iy;
  for(int i = 0;i < n;++i){
    ia = iage(i);
    iy = iyear(i);
    Elog_index(i) = log_q_vec(i) + log(N(iy,ia)) - fs(i)*Z(iy,ia);
    resid_index(i) = log_index(i) - Elog_index(i); 
    std_resid_index(i) = resid_index(i)/std_index_vec(i);
  }
  
  // End of model, now fit functions (nll - negative loglikelihoods);
  
  // Catch at age nll;
  
  for(int y = 0;y < Y-1;++y){
    for(int a = 0;a < A;++a){          
      if(catch_censor==0){
	nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
      }
      if(catch_censor==1){  
        if(C(y,a) >= c_detect){
	  nll -= dnorm(C_resid(y,a),zero,std_log_C(a),true);
	}
        else{
	  nll -= pnorm4(std_C_resid(y,a));
	}
      } 
    }
  }
  
  // Index nll;
  
  for(int i = 0;i < n;++i){   
    if(index(i) >= i_detect){
      if(index_censor==0){
	nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);
      }       
      if(index_censor==1){
	nll -= pnorm4(std_resid_index(i));
      } 
    }      
    else{
      nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);
    }
  }
  
  //recruitment
  //  nll += SCALE(AR1(phi_logR),std_log_R)(log_Rec_dev);
  nll -= dnorm(log_Rec_dev,zero, std_log_R, true).sum();
  
  //Log_F nll
  //RW on first age;
  vector<Type> del = log_F.col(0);
  nll -= dnorm(del(0),Type(-10.0),one, true);
  for(int y = 1;y < Y-1;++y){ 
    nll -= dnorm(del(y),del(y-1),one, true);
  }
  array<Type> log_F1(Y-1,A-1);
  for(int a = 1;a < A;++a){
    log_F1.col(a-1) = log_F.col(a);
  } 
  //year x age correlation on first+1:last ages;
  nll += SEPARABLE(AR1(ar_logF_age),AR1(ar_logF_year))(log_F1);
  
  //pe nll
  if(use_pe==1){
    nll += SEPARABLE(AR1(ar_pe_age),AR1(ar_pe_year))(pe);
  }
  
  //year effect nll;
  if(use_cye==1){
    for(int y = 0;y < Y-1;++y){
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

  int i,j;
  //average F 4-6 and 6-9 WHY INSIDE
  for(i = 0;i < Y-1;++i){
    aveF_46(i) = zero; 
    tni = zero;
    for(j = 2;j < 4;++j){
      aveF_46(i) += F(i,j)*N(i,j); 
      tni += N(i,j);
    }
    aveF_46(i) = aveF_46(i)/tni;  
    aveF_69(i) = zero; 
    tni = zero;
    for(j = 4;j < 7;++j){
      aveF_69(i) += F(i,j)*N(i,j); 
      tni += N(i,j);
    }
    aveF_69(i) = aveF_69(i)/tni;
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

