#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' catchProp
//'
//' make a row prop
//'
//' Why doesn't this work
//'
//' @param catchy a row to proportionize
//' @export
// [[Rcpp::export]]
Eigen::VectorXd catchProp(Eigen::Map<Eigen::VectorXd > catchy){
  double tot = catchy.sum();
  Eigen::VectorXd prop = catchy/tot;
  return prop;
}

//'makeCRLs
//'
//' Turn matrix of proportions into CRLs
//'
//' @param cProps matrix of proportions
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd makeCRLs(Eigen::Map<Eigen::MatrixXd> cProps){
  int A = cProps.cols();
  int Y = cProps.rows();

  Eigen::MatrixXd crls(Y,A-1);

  for(int a = 0;a < A-1;a++){
    for(int y = 0;y < Y;y++){
      Eigen::VectorXd cprow = cProps.row(y);
      double denom = cprow.segment(a,A-a).sum();
      double num = cProps(y,a);
      double pi = num/denom;
      crls(y,a) = pi/(1-pi);
    }
  }
  return crls;
}


