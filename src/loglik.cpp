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
//  for (i = 0; i < Sdm; ++i) {
//   Vstar += beta(i)*S(k,j,i);
//  }
    for (i = 0; i < Mdm; ++i) {
     Vstar += beta(i+Sdm+Wdm)*Z(k,j,i);
    }
    Wstar += Ustar + Vstar;
    llik += counts(j,k)*(log(2.0) + Wstar + GammaW(j)+ GammaM(k) + gm + gw + log(pmfW(j)*pmfM(k))-log((1.0+exp(GammaW(j)))*(1.0+exp(GammaM(k))))) ;
   }
  }

  for (j = 0; j < NumGammaW; ++j) {
    llik += counts(j,NumGammaM)*(GammaW(j) + gw + log(pmfW(j)) - log(1.0+exp(GammaW(j))));
//Rprintf("j %d counts %f\n",j,counts(j,NumGammaM));
  }
  for (k = 0; k < NumGammaM; ++k) {
    llik += counts(NumGammaW,k)*(GammaM(k) + gm + log(pmfM(k)) - log(1.0+exp(GammaM(k))));
//Rprintf("k %d counts %f\n",k,counts(NumGammaW,k));
  }
  if(std::isnan(llik) || !std::isfinite(llik)){
  Rprintf("bad lik %f GammaM(0) %f  GammaM(1) %f\n",llik, GammaM(0), GammaW(0));
    llik= 10000.0;
  }
//if(std::isnan(llik) || !std::isfinite(llik) || (llik < -1000000000.0)){
//Rprintf("lik %f\n",llik);
  return llik;
}
