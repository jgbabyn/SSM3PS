template<class Type>
vector<Type> ssbCalc(matrix<Type> logN, matrix<Type> weight, matrix<Type> mat){
  vector<Type> ret(logN.rows());
  ret.setZero();
  for(int y = 0; y < logN.rows();y++){
    for(int a = 0; a < logN.cols();a++){
      ret(y) += exp(logN(y,a))*weight(y,a)*mat(y,a);
    }
  }
  return ret;
}

template<class Type>
vector<Type> tsbCalc(matrix<Type> logN,matrix<Type> weight){
  vector<Type> ret(logN.rows());
  ret.setZero();
  for(int y = 0; y < logN.rows();y++){
    for(int a = 0; a < logN.cols();a++){
      ret(y) += exp(logN(y,a))*weight(y,a);
    }
  }
  return ret;
}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}
