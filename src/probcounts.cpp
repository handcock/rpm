#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
double probcounts(arma::vec beta, arma::vec Gamma, arma::cube S, arma::cube X, arma::cube Z, arma::vec pmfW, arma::vec pmfM, arma::mat counts, double gw, double gm){
  int Sdm=S.n_slices;
  int Wdm=X.n_slices;
  int Mdm=Z.n_slices;
  int NumGammaW=pmfW.size();
  int NumGammaM=pmfM.size();
  int i,j,k;
  double Ustar, Vstar, Wstar;
  double llik;

  arma::vec GammaW(NumGammaW);
  arma::vec GammaM(NumGammaM);
  GammaW = Gamma.head(NumGammaW);
  GammaM = Gamma.tail(NumGammaM);

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
    llik += 0.5*counts(j,k)*(log(2.0) + Wstar + GammaW(j)+ GammaM(k) + gm + gw + log(pmfW(j)*pmfM(k))-log((1.0+exp(GammaW(j)))*(1.0+exp(GammaM(k))))) ;
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
