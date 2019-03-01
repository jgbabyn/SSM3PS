#include <TMB.hpp>
#include "pnorm4.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{

  enum fleetType{
		 catch = 0,
		 survey = 1
  };

  enum censorType{
		  none = 0,
		  basic = 1,
		  stable = 2
  };
  
  DATA_VECTOR(obs);
  DATA_INTEGER(censb);
  DATA_INTEGER(censd);
  DATA_FACTOR(fleetTypes);
  


  
  Type nll = 0.0;
  vector<Type> pred(obs.size());
  for(int i = 0; i < obs.size();i++){

    switch(fleetTypes(f)){
    case catch:
      if(censb)
      break;
    case survey:
      break;
    default:
      error("fleetType not implemented")
      }
      
  }
      
      
    
	
    return nll;
}

