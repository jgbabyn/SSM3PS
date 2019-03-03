#include <TMB.hpp>
#include "pnorm4.hpp"
#include "Fcorr.hpp"
#include "utilsFuns.hpp"



template<class Type>
Type objective_function<Type>::operator() ()
{

  enum fleetType{
		 Catch = 0,
		 survey = 1,
		 landings = 2
  };

  enum censorType{
		  none = 0,
		  basic = 1,
		  stable = 2
  };

  enum recruitType{
		   randomWalk = 0,
		   ricker = 1,
		   beverton = 2
  };

  //Handed by obsSetup
  DATA_VECTOR(obs); //indices
  DATA_IMATRIX(aux); //Ages, fleetType, etc.
  DATA_VECTOR(fs);
  DATA_INTEGER(plusGroup);
  DATA_IVECTOR(fleetTypes);
  DATA_IVECTOR(idx);

  //Blank keys created by keyGen 
  DATA_IVECTOR(keyF);
  DATA_IMATRIX(keyQ);
  DATA_IMATRIX(keySD);

  
  DATA_INTEGER(corrFFlag);
 
  
  DATA_INTEGER(recType);
  DATA_MATRIX(M);
  DATA_MATRIX(stock_wt);
  DATA_MATRIX(midy_wt);
  DATA_MATRIX(mat);
  DATA_INTEGER(minYear);
  DATA_INTEGER(minAge);

  PARAMETER_MATRIX(logN);
  PARAMETER_MATRIX(logF);
  PARAMETER(logSDrec);
  PARAMETER(logSDsur);
  PARAMETER_VECTOR(tRhoF);
  PARAMETER_VECTOR(logsdF);
  PARAMETER_VECTOR(logQ);
  PARAMETER_VECTOR(obsSD);
  

  vector<Type> ssb = ssbCalc(logN,stock_wt,mat);
  vector<Type> tsb = tsbCalc(logN,stock_wt);
  
  Type nll = 0.0;

  matrix<Type> logFF(logN.rows(),logN.cols());
  for(int y = 0; y < logN.rows();y++){
    for(int a = 0; a < logN.cols();a++){
      logFF(y,a) = logF(y,keyF(a));
    }
  }
  REPORT(logFF);

  
  //Recruitment
  Type sdRec = exp(logSDrec);
  Type sdSur = exp(logSDsur);
  Type predN;
  for(int y = 1; y < logN.rows(); y++){
    switch(recType){
    case randomWalk:
      predN = logN(y-1,0);
      break;
    case ricker:
      {
  	PARAMETER_VECTOR(rickpar);
  	predN = rickpar(0)+log(ssb(y))-exp(rickpar(1))*ssb(y);
  	break;
      }
    case beverton:
      {
  	PARAMETER_VECTOR(bhpar);
  	predN = bhpar(0)+log(ssb(y))-log(1.0+exp(bhpar(1))*ssb(y));
  	break;
      }
    default:
      error("Stock recruitment type not implemented");
      break;
    }
    nll -= dnorm(logN(y,0),predN,sdRec,true);
  }

  //Rest of N matrix
  for(int y = 1; y < logN.rows(); y++){
    for(int a = 1; a < logN.cols();a++){
      predN = logN(y-1,a-1)-exp(logFF(y-1,a-1))-M(y-1,a-1);
      if(a==(logN.cols()-1)){
      	if(plusGroup == true){
      	  predN = log(exp(predN) + exp(logN(y-1,a)-exp(logFF(y-1,a))-M(y-1,a)));
      	}
      }
      nll -= dnorm(logN(y,a),predN,sdSur,true);
    }
  }
        
  // F correlation
  matrix<Type> sigmaF = sigmaFgen(logF.cols(),corrFFlag,logsdF,tRhoF);
  density::MVNORM_t<Type> nldens(sigmaF);
  for(int y = 1; y < logF.rows(); y++){
    nll += nldens(logF.row(y)-logF.row(y-1));
  }
    
  //observation equations
  vector<Type> logPred(obs.size());
  matrix<Type> ECW(logN.rows(),logN.cols());
  Type Z;
  for(int i = 0; i < obs.size(); i++){
     int y = aux(i,0)-minYear;
     int a = aux(i,1);
     int f = aux(i,2)-1;
     int fType = fleetTypes(f);

     if(!isNAINT(a)){
       a = aux(i,1)-minAge;  //Very much very important to be here!    
       Z=exp(logFF(y,a))+M(y,a);
     }else{
       Z=0;
     }

    switch(fType){
    case Catch:
      logPred(i) = logN(y,a)-log(Z)+log(1-exp(-Z))+logFF(y,a);
      ECW(y,a) = exp(logPred(i))*midy_wt(y,a);
      break;
    case survey:
      logPred(i) = logQ(keyQ(f,a))+logN(y,a)-Z*fs(i);
      break;
    case landings:
      logPred(i) = log(ECW.row(y).sum()); //Gotta make sure catch is BEFORE landings in obs to do this
      break;
    default:
      error("Fleet type not implemented!");
    }

  }

  //Basic nlls for now, get something better later?
  for(int i = 0; i < logPred(i);i++){
       int a = aux(i,1)-minAge;
       int f = aux(i,2)-1;
       int fType = fleetTypes(f);

       switch(fType){
       case Catch:
  	 nll -= dnorm(log(obs(i)),logPred(i),obsSD(keySD(f,a)),true);
  	 break;
       case survey:
  	 nll -= dnorm(log(obs(i)),logPred(i),obsSD(keySD(f,a)),true);
  	 break;
       case landings:
  	 break; //Do nothing with landings for now
       default:
  	 error("How did you get down here?");
  	   break;
       }
       
  }
  
  


  ADREPORT(logPred);
  return nll;
}

