#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
NumericMatrix hloglik(NumericVector beta, NumericVector GammaW, NumericVector GammaM, NumericVector Sd, NumericVector Xd, NumericVector Zd, IntegerVector Sdim, IntegerVector Xdim, IntegerVector Zdim, NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, NumericMatrix counts, double gw, double gm, int constraints){
  arma::cube S(Sd.begin(), Sdim[0], Sdim[1], Sdim[2]);
  arma::cube X(Xd.begin(), Xdim[0], Xdim[1], Xdim[2]);
  arma::cube Z(Zd.begin(), Zdim[0], Zdim[1], Zdim[2]);
  int dm=Sdim(2)+Xdim(2)+Zdim(2);
  int NumBeta=beta.size();
  int NumGammaW=GammaW.size();
  int NumGammaM=GammaM.size();
  int j,k,l;

  NumericMatrix gf(NumBeta+NumGammaW+NumGammaM+1,NumBeta+NumGammaW+NumGammaM+1);
  for (l = 0; l < NumBeta+NumGammaW+NumGammaM+1; ++l) {
   for (k = 0; k < NumBeta+NumGammaW+NumGammaM+1; ++k) {
    gf(l,k)=0.0;
   }
  }

  for (j = 0; j < NumGammaW; ++j) {
    gf(j + dm,j + dm) = -counts(j,NumGammaM);
    for (l = 0; l < NumGammaM; ++l) {
      gf(j + dm,j + dm) += -0.5*counts(j,l);
    }
    gf(j + dm,j + dm) *= ( exp(GammaW(j)) / ((1.0+exp(GammaW(j)))*(1.0+exp(GammaW(j)))) );
  }
  for (k = 0; k < NumGammaM; ++k) {
   gf(k + NumGammaW+dm, k + NumGammaW+dm) = -counts(NumGammaW,k);
   for (l = 0; l < NumGammaW; ++l) {
    gf(k + NumGammaW+dm, k + NumGammaW+dm) += -0.5*counts(l,k);
   }
   gf(k + NumGammaW+dm, k + NumGammaW+dm) *= ( exp(GammaM(k)) / ((1.0+exp(GammaM(k)))*(1.0+exp(GammaM(k)))) );
  }
  for (j = 0; j < NumGammaW; ++j) {
    gf(NumBeta+NumGammaW+NumGammaM,NumBeta+NumGammaW+NumGammaM) -= counts(NumGammaW,j)*exp(gw-2.0*gm);
    for (l = 0; l < NumGammaM; ++l) {
      gf(NumBeta+NumGammaW+NumGammaM,NumBeta+NumGammaW+NumGammaM) -= 0.5*counts(j,l)*exp(gw-2.0*gm);
    }
  }
  return gf;
}
