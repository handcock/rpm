#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(name = "PropMax")]]
IntegerMatrix PropMax(const NumericMatrix &stat) {
  NumericVector ret = NumericVector(stat.nrow());
  IntegerMatrix Prop = IntegerMatrix(stat.nrow(),stat.ncol());
  
 for (unsigned int i = 0u; i < static_cast<unsigned int>(stat.nrow()); ++i) {
   ret(i) = max(stat(i,_));
 }
 for (unsigned int i = 0u; i < static_cast<unsigned int>(stat.nrow()); ++i) {
   if(ret(i) < 0.0) ret(i)=0.0;
 }
 for (unsigned int i = 0u; i < static_cast<unsigned int>(stat.nrow()); ++i) {
  for (unsigned int j = 0u; j < static_cast<unsigned int>(stat.ncol()); ++j) {
   if(stat(i,j) >= ret(i)-1.0e-14){
    Prop(i,j) = 1;
   }else{
    Prop(i,j) = 0;
   }
  }
 }
  
 return Prop;
}
