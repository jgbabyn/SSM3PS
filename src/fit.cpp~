#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // input data;  
  DATA_MATRIX(M);
  DATA_MATRIX(weight); 
  DATA_MATRIX(mat); 
  DATA_MATRIX(midy_weight);
  DATA_MATRIX(C);         
  DATA_MATRIX(log_C);  
  DATA_VECTOR(landings); 
  DATA_VECTOR(log_landings); 
  DATA_VECTOR(index);
  DATA_IVECTOR(i_zero); 
  DATA_IMATRIX(C_zero);  
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(iage);
  DATA_IVECTOR(isurvey); 
  DATA_IVECTOR(iq);
  DATA_VECTOR(fs);
  DATA_INTEGER(A);     
  DATA_INTEGER(Y);
  DATA_VECTOR(log_index);  
  DATA_INTEGER(index_censor); 
  DATA_INTEGER(catch_censor); 
  DATA_INTEGER(use_pe);        
  DATA_INTEGER(use_cye); 
  
  int n = index.size();  
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
  
  int i,j; 
  
  //compute F 
  for(j = 0;j < A;++j){ 
    for(i = 0;i < Y-1;++i){
      F(i,j) = exp(std_logF(j)*log_F(i,j)); 
      Z(i,j) = F(i,j) + M(i,j);
    }
  }
  
  // The cohort model;
  
  vector<Type> log_Rec = log_Rec_dev + log_Rec_mean;
  log_N(0,0) = log_Rec(0);
  for(i = 1;i < A;++i){log_N(0,i) = log_No(i-1);}
  
  for(i = 1;i < Y;++i){
    log_N(i,0) = log_Rec(i);  
    for(j = 1;j < A;++j){log_N(i,j) = log_N(i-1,j-1) - Z(i-1,j-1) + std_pe*pe(i,j);}
  }
  for(j = 0;j < A;++j){ 
    for(i = 0;i < Y;++i){
      N(i,j) = exp(log_N(i,j));
    }
  }
  B_matrix = weight.array()*N.array();
  SSB_matrix = mat.array()*B_matrix.array(); 
  
  // Catch model prediction from Baranov;
  
  vector<Type> landings_pred(Y-1);
  
  for(i = 0;i < Y-1;++i){ 
    landings_pred(i)=zero;
    for(j = 0;j < A;++j){
      EC(i,j) = N(i,j)*(one - exp(-one*Z(i,j)))*F(i,j)/Z(i,j);
      ECW(i,j) = EC(i,j)*midy_weight(i,j);
      log_EC(i,j) = log(EC(i,j));
      C_resid(i,j) = log_C(i,j) - log_EC(i,j) + cye(i,j);
      landings_pred(i) += ECW(i,j); 
      std_C_resid(i,j) = C_resid(i,j)/std_log_C(j);
    }}
  vector<Type> log_landings_pred = log(landings_pred);
  
  //  Survey index predictions, and residuals;
  
  int ia,iy;
  for(i = 0;i < n;++i){
    ia = iage(i);
    iy = iyear(i);
    Elog_index(i) = log_q_vec(i) + log(N(iy,ia)) - fs(i)*Z(iy,ia);
    resid_index(i) = log_index(i) - Elog_index(i); 
    std_resid_index(i) = resid_index(i)/std_index_vec(i);
  }
  
  // End of model, now fit functions (nll - negative loglikelihoods);
  
  // Catch at age nll;
  
  for(i = 0;i < Y-1;++i){
    for(j = 0;j < A;++j){          
      if(catch_censor==0){nll -= dnorm(C_resid(i,j),zero,std_log_C(j),true);}
      if(catch_censor==1){  
        if(C_zero(i,j)==0){nll -= dnorm(C_resid(i,j),zero,std_log_C(j),true);}
        if(C_zero(i,j)==1){nll -= log(pnorm(std_C_resid(i,j)));}
      } 
    }
  }
  
  // Index nll;
  
  for(i = 0;i < n;++i){   
    if(i_zero(i)==1){
      if(index_censor==0){nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);}       
      if(index_censor==1){nll -= log(pnorm(std_resid_index(i)));} 
    }      
    if(i_zero(i)==0){nll -= dnorm(resid_index(i),zero,std_index_vec(i),true);}
  }
  
  //recruitment
  //  nll += SCALE(AR1(phi_logR),std_log_R)(log_Rec_dev);
  nll -= dnorm(log_Rec_dev,zero, std_log_R, true).sum();
  
  //Log_F nll
  //RW on first age;
  vector<Type> del = log_F.col(0);
  nll -= dnorm(del(0),Type(-10.0),one, true);
  for(j = 1;j < Y-1;++j){ 
    nll -= dnorm(del(j),del(j-1),one, true);
  }
  array<Type> log_F1(Y-1,A-1);
  for(i = 1;i < A;++i){
    log_F1.col(i-1) = log_F.col(i);
  } 
  //year x age correlation on first+1:last ages;
  nll += SEPARABLE(AR1(ar_logF_age),AR1(ar_logF_year))(log_F1);
  
  //pe nll
  if(use_pe==1){nll += SEPARABLE(AR1(ar_pe_age),AR1(ar_pe_year))(pe);}
  
  //year effect nll;
  if(use_cye==1){
    for(j = 0;j < Y-1;++j){
      vector<Type> del = vector<Type>(cye.row(j));
      nll += SCALE(AR1(ar_cye_year),std_cye)(del);
    }
  }       
  
  // create some REPORT output;  
  
  for(i = 0;i < Y;++i){
    biomass(i) = zero;     
    ssb(i) = zero;
    for(j = 0;j < A;++j){
      biomass(i) += B_matrix(i,j); 
      ssb(i) += SSB_matrix(i,j);}
  }          
  log_biomass = log(biomass);
  log_ssb = log(ssb);       
  
  //pop size weighted ave F;  
  
  Type tni;
  
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

