#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
NumericMatrix augpmfWM(NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, double gw, double gm,
                       NumericVector pmfW_N, NumericVector pmfM_N, double gwN, double gmN){
  int NumGammaW=pmfW.size();
  int NumGammaM=pmfM.size();
  int j,k;

  NumericMatrix apmf(NumGammaW+1,NumGammaM+1);
  apmf(NumGammaW,NumGammaM)=0.0;
  for (j = 0; j < NumGammaW; ++j) {
   for (k = 0; k < NumGammaM; ++k) {
    apmf(j,k) = pmf(j,k)*exp(gmN + gwN - gm - gw)*pmfW_N(j)*pmfM_N(k)/(pmfW(j)*pmfM(k));
   }
  }

  for (j = 0; j < NumGammaW; ++j) {
    apmf(j,NumGammaM) = pmf(j,NumGammaM)*exp(gwN - gw)*pmfW_N(j)/pmfW(j);
  }
  for (k = 0; k < NumGammaM; ++k) {
    apmf(NumGammaW,k) = pmf(NumGammaW,k)*exp(gmN - gm)*pmfM_N(k)/pmfM(k);
  }
  return apmf;
}
