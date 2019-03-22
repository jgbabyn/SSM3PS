/*Code ripped from SAM's SAM/stockassessment/inst/include/define.hpp github file.
  I'm not going to attempt writing my own MVNORM with support for OSA residuals for a class project...
  Maybe for fun later...*/

template<class Type>
Type logspace_add_p (Type logx, Type logy, Type p) {
  return log((Type(1)-p)*exp(logy-logx)+p)+logx; // the order of x and y is taylored for this application 
}

template<class Type>
Type logdrobust(Type x, Type p){
  Type ld1=dnorm(x,Type(0.0),Type(1.0),true);
  if(p<Type(1.0e-16)){
    return ld1;
  }else{
    Type ld2=dt(x,Type(3),true);
    Type logres=logspace_add_p(ld2,ld1,p);
    return logres;
  }
}
VECTORIZE2_tt(logdrobust)

template <class Type>
class MVMIX_t{
  Type halfLogDetS;         
  Type p1;                  /*fraction t3*/
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  matrix<Type> cov(){return Sigma;}
  void setSigma(matrix<Type> Sigma_){
    Sigma = Sigma_;
    sd = sqrt(vector<Type>(Sigma.diagonal()));
    Eigen::LLT<Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> > llt(Sigma);
    L_Sigma = llt.matrixL();
    vector<Type> D=L_Sigma.diagonal();
    halfLogDetS = sum(log(D));
    inv_L_Sigma = atomic::matinv(L_Sigma);
  }
  void setSigma(matrix<Type> Sigma_, Type p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    vector<Type> z = inv_L_Sigma*x;
    return -sum(logdrobust(z,p1))+halfLogDetS;
  }
  Type operator()(vector<Type> x, vector<Type> keep){
    matrix<Type> S = Sigma;
    vector<Type> not_keep = Type(1.0) - keep;
    for(int i = 0; i < S.rows(); i++){
      for(int j = 0; j < S.cols(); j++){
	S(i,j) = S(i,j) * keep(i) * keep(j);
      }
      //S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(1)/M_PI),2); //(t(1))
      S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(2)/(M_PI*sqrt(Type(3)))),2);
    }
    return MVMIX_t<Type>(S,p1)(x * keep);
  }

  vector<Type> simulate() {
    int siz = Sigma.rows();
    vector<Type> x(siz);
    for(int i=0; i<siz; ++i){
      Type u = runif(0.0,1.0);
      if(u<p1){
        x(i) = rt(3.0);
      }else{
        x(i) = rnorm(0.0,1.0);
      }
    }
    x = L_Sigma*x;
    return x;
  }
};

template <class Type>
MVMIX_t<Type> MVMIX(matrix<Type> Sigma, Type p1){
  return MVMIX_t<Type>(Sigma,p1);
}
