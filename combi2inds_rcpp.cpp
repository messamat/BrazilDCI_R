#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

DataFrame combi2inds(const Rcpp::CharacterVector inputVector){
  const int len = inputVector.size();
  const int retLen = len * (len-1) / 2;
  Rcpp::IntegerVector outputVector1(retLen);
  Rcpp::IntegerVector outputVector2(retLen);
  int indexSkip;
  for (int i = 0; i < len; ++i){
    indexSkip = len * i - ((i+1) * i)/2;
    for (int j = 0; j < len-1-i; ++j){
      outputVector1(indexSkip+j) = i+1;
      outputVector2(indexSkip+j) = i+j+1+1;
    }
  }
  return(Rcpp::DataFrame::create(Rcpp::Named("xid") = outputVector1,
                                 Rcpp::Named("yid") = outputVector2));
};
