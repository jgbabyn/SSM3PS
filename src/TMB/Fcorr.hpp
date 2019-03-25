template<class Type>
Type rhoTrans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2)*x))-Type(1);
}

VECTORIZE1_t(rhoTrans);


/*
  Generate correlation matrix for F's

  param nage number of ages
  param corrFlag Type of correlation structure
  param tRhoF vector of untransformed rhos for corr structure
 
*/
template<class Type>
matrix<Type> sigmaFgen(int nage,int corrFlag, vector<Type> logsdF,vector<Type> tRhoF){

  enum corrType{
		independent = 0,
		parallel = 1,
		compound = 2,
		ARone = 3,
		custom = 4
  };

  matrix<Type> sigma(nage,nage);
  vector<Type> sdF = exp(logsdF);
  vector<Type> rhoF = rhoTrans(tRhoF);

  sigma.setZero();

  switch(corrFlag){
  case independent:
    sigma.diagonal() = sdF*sdF;
    break;
  case parallel:
    sigma.diagonal() = sdF*sdF;
    for(int i=0; i < nage; i++){
      for(int j=0; j < i; j++){
	//Anders explained this is needed to get parallel to coverge.
	sigma(i,j) = 0.99999*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case compound:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = rhoF(0)*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case ARone:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = pow(rhoF(0),Type(i-j))*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case custom:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = rhoF(i)*rhoF(j)*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
    

  default:
    error("F correlation mode not supported");
      break;
  }

  return sigma;
}
