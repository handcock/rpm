#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(name = "rowMax")]]
NumericVector rowMax(const NumericMatrix &stat) {
  NumericVector ret = NumericVector(stat.nrow());
  
 for (unsigned int i = 0u; i < static_cast<unsigned int>(stat.nrow()); ++i) {
   ret[i] = max(stat(i, _));
 }
  
  return ret;
}
