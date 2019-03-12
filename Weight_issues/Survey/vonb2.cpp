#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{  
	DATA_INTEGER(nobs);
	DATA_INTEGER(nc);   
	DATA_INTEGER(na); 
	DATA_VECTOR(y);  
	DATA_VECTOR(sf);     
	DATA_VECTOR(cohort);
	DATA_VECTOR(age);
	DATA_IVECTOR(ic);  
	DATA_IVECTOR(iy); 
	DATA_IVECTOR(ia); 
	DATA_INTEGER(REk); 
	DATA_INTEGER(RELinf);
	DATA_INTEGER(REpo);
	
//Fixed effects
   PARAMETER(log_Linf); 
   PARAMETER(log_k_main);
   PARAMETER(logit_po); 
   PARAMETER(log_std_log_Linf);
   PARAMETER_VECTOR(log_std_log_k);    
   PARAMETER(log_std_lpo);  
   PARAMETER(log_std_me); 
   
//Random Effects;  
   PARAMETER_VECTOR(rlog_Linf);
   PARAMETER_MATRIX(rlog_k_dev);    
   PARAMETER_VECTOR(rlogit_po);
   
   Type zero = 0.0;
   Type one = 1.0;  
   int n=y.size(); 
   Type mZ = 0.0;

   Type std_log_Linf = exp(log_std_log_Linf); 
   vector<Type> std_log_k = exp(log_std_log_k);       
   Type std_lpo = exp(log_std_lpo); 
   Type std_me = exp(log_std_me); 

   vector<Type> logy = log(y);      
   vector<Type> lpo = logit_po + rlogit_po;
   vector<Type> po_c = exp(lpo)/(one + exp(lpo));
   vector<Type> log_Linf_c = log_Linf + rlog_Linf;
   vector<Type> Linf_c = exp(log_Linf_c); 
      
   vector<Type> pred(nobs);   
   vector<Type> log_pred(nobs);
   vector<Type> resid(nobs);  
   vector<Type> std_resid(nobs);
      
   matrix<Type> log_k(na,nc);
   matrix<Type> ka(na,nc);   
   for (int j = 0; j < nc; ++j){ 
     log_k(0,j) = log_k_main + rlog_k_dev(0,j); 
     ka(0,j) = exp(log_k(0,j)); 
     for (int i = 1; i < na; ++i){ 
       log_k(i,j) = log_k_main + rlog_k_dev(i,j);  
       ka(i,j) = ka(i-1,j) + exp(log_k(i,j));
     }
   }  

   Type nll = zero;
   
   vector<Type> ck(n);
   for (int i = 0; i < n; ++i){
     if(ia(i)==0){ck(i)=(one+sf(i))*ka(0,ic(i));}
     if(ia(i)>0){ck(i)=ka(ia(i)-1,ic(i))+(one+sf(i))*exp(log_k(ia(i),ic(i)));}
   } 
   
   log_pred = log_Linf_c(ic) + log(one - (one - po_c(ic))*exp(-ck)); 
   pred = exp(log_pred);
   resid = logy - log_pred;
   std_resid = resid/std_me;  
   
// nll for data;
   nll -= dnorm(resid, zero, std_me, true).sum();  

   if(REpo==1){ //Random effect on po= RW
     nll -= dnorm(rlogit_po(0), zero, std_lpo, true); 
	   for (int j = 1; j < nc; ++j) {  
		   nll -= dnorm(rlogit_po(j), rlogit_po(j-1), std_lpo, true);
     } 
   }    
   if(REk==1){  //random effect on Kdev = RW
//     for (int j = 0; j < nc; ++j){ 
//       for (int i = 0; i < na; ++i){  
//         nll -= dnorm(rlog_k_dev(i,j),zero, std_log_k, true);
//       }
//     }  
      
    for(int j = 0;j < nc;++j){
      for(int i = 0;i < na;++i){
       if((i==0)&(j == 0)){mZ = zero;} 
       if((i>0)&(j == 0)){mZ = rlog_k_dev(i-1,j);}
        if((i == 0)&(j > 0)){mZ = rlog_k_dev(i,j-1);}         
        if((i > 0)&(j > 0)){
          mZ = rlog_k_dev(i,j-1) + rlog_k_dev(i-1,j) - rlog_k_dev(i-1,j-1);
          }        
         nll -= dnorm(rlog_k_dev(i,j),mZ,std_log_k(i),true);
      }
    } 
   }    
   if(RELinf==1){
     nll -= dnorm(rlog_Linf(0), zero, std_log_Linf, true); 
	   for (int j = 1; j < nc; ++j) {  
		   nll -= dnorm(rlog_Linf(j), rlog_Linf(j-1), std_log_Linf, true);
     } 
   }

   REPORT(pred);
   REPORT(log_pred);
   REPORT(resid);
   REPORT(std_resid);
   REPORT(Linf_c);   
   REPORT(ka);       
   REPORT(log_k);     
   REPORT(rlog_k_dev);      
   REPORT(po_c);     
   REPORT(lpo);      
   REPORT(ck);      
   REPORT(rlogit_po); 
   REPORT(std_lpo); 
   REPORT(std_log_k);
   REPORT(std_log_Linf);
   
   ADREPORT(log_Linf_c);   
   ADREPORT(log_k);      
   ADREPORT(lpo);  
   ADREPORT(log_pred);

   return nll;
}















