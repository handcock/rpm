#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
NumericVector gloglik(NumericVector beta, NumericVector GammaW, NumericVector GammaM, NumericVector Sd, NumericVector Xd, NumericVector Zd, IntegerVector Sdim, IntegerVector Xdim, IntegerVector Zdim, NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, NumericMatrix counts, double gw, double gm, int constraints){
  arma::cube S(Sd.begin(), Sdim[0], Sdim[1], Sdim[2]);
  arma::cube X(Xd.begin(), Xdim[0], Xdim[1], Xdim[2]);
  arma::cube Z(Zd.begin(), Zdim[0], Zdim[1], Zdim[2]);
  int Sdm=Sdim(2);
  int Wdm=Xdim(2);
  int Mdm=Zdim(2);
  int NumBeta=beta.size();
  int NumGammaW=GammaW.size();
  int NumGammaM=GammaM.size();
  int i,j,k,l;

  NumericVector gf(NumBeta+NumGammaW+NumGammaM);
  for (l = 0; l < NumBeta+NumGammaW+NumGammaM; ++l) {
    gf(l)=0.0;
  }

  // \partial beta
  for (j = 0; j < NumGammaW; ++j) {
   for (k = 0; k < NumGammaM; ++k) {
    for (i = 0; i < Sdm; ++i) {
     gf(i) += counts(j,k)*S(j,k,i);
    }
    for (i = 0; i < Wdm; ++i) {
     gf(i+Sdm) += counts(j,k)*X(j,k,i);
    }
    for (i = 0; i < Mdm; ++i) {
     gf(i+Sdm+Wdm) += counts(j,k)*Z(k,j,i);
    }
   }
  }

  // \partial g(x,*)
  for (j = 0; j < NumGammaW; ++j) {
    gf(j + Mdm+Sdm+Wdm) = counts(j,NumGammaM);
    for (l = 0; l < NumGammaM; ++l) {
      gf(j + Mdm+Sdm+Wdm) += counts(j,l);
    }
    gf(j + Mdm+Sdm+Wdm) *= 1.0 / (1.0+exp(GammaW(j)));
  }
  // \partial g(*,z)
  for (k = 0; k < NumGammaM; ++k) {
   gf(k + NumGammaW+Mdm+Sdm+Wdm) = counts(NumGammaW,k);
   for (l = 0; l < NumGammaW; ++l) {
    gf(k + NumGammaW+Mdm+Sdm+Wdm) += counts(l,k);
   }
   gf(k + NumGammaW+Mdm+Sdm+Wdm) *= 1.0 / (1.0+exp(GammaM(k)));
  }
  return gf;
}
