#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
double loglik(NumericVector beta, NumericVector GammaW, NumericVector GammaM, NumericVector Sd, NumericVector Xd, NumericVector Zd, IntegerVector Sdim, IntegerVector Xdim, IntegerVector Zdim, NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, NumericMatrix counts, double gw, double gm, int constraints){
  arma::cube S(Sd.begin(), Sdim[0], Sdim[1], Sdim[2]);
  arma::cube X(Xd.begin(), Xdim[0], Xdim[1], Xdim[2]);
  arma::cube Z(Zd.begin(), Zdim[0], Zdim[1], Zdim[2]);
  int Sdm=Sdim(2);
  int Wdm=Xdim(2);
  int Mdm=Zdim(2);
  int NumGammaW=GammaW.size();
  int NumGammaM=GammaM.size();
  int i,j,k;
  double Ustar, Vstar, Wstar;
  double llik;

  llik=0.0;
  for (j = 0; j < NumGammaW; ++j) {
   for (k = 0; k < NumGammaM; ++k) {
    Wstar = 0.0;
    for (i = 0; i < Sdm; ++i) {
     Wstar += beta(i)*S(j,k,i);
    }
    Ustar = 0.0;
    for (i = 0; i < Wdm; ++i) {
     Ustar += beta(i+Sdm)*X(j,k,i);
    }
    Vstar = 0.0;
    for (i = 0; i < Mdm; ++i) {
     Vstar += beta(i+Sdm+Wdm)*Z(k,j,i);
    }
    Wstar += Ustar + Vstar;
    llik += 0.5*counts(j,k)*(Wstar + GammaW(j)+ GammaM(k) + gm + gw + log(pmfW(j)*pmfM(k))-log((1.0+exp(GammaW(j)))*(1.0+exp(GammaM(k))))) ;
   }
  }

  for (j = 0; j < NumGammaW; ++j) {
    llik += counts(j,NumGammaM)*(GammaW(j) + gw + log(pmfW(j)) - log(1.0+exp(GammaW(j))));
  }
  for (k = 0; k < NumGammaM; ++k) {
    llik += counts(NumGammaW,k)*(GammaM(k) + gm + log(pmfM(k)) - log(1.0+exp(GammaM(k))));
  }
  if(std::isnan(llik) | !std::isfinite(llik)){
    llik= 10000.0;
  }
  return llik;
}
