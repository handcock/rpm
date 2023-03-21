#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "PSeqcond.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat augpmfnew(arma::vec beta, arma::vec GammaW, arma::vec GammaM, arma::cube S, arma::cube X, arma::cube Z, arma::vec pmfW, arma::vec pmfM, double gw, double gm) {
  unsigned int it=0u;
//unsigned int NumBeta=beta.size();
  unsigned int NumGammaW=GammaW.size();
  unsigned int NumGammaM=GammaM.size();
  unsigned int NumGamma=NumGammaW+NumGammaM;

  unsigned int i,j,k;
  double Ustar, Vstar, Wstar;

  arma::mat delta=arma::eye(NumGamma,NumGamma);
  arma::vec Gamma(NumGamma);
  arma::vec GammaOrig(NumGamma);
  arma::vec GammaDiff(NumGamma);
  arma::vec eqcond(NumGamma);
  arma::mat J(NumGamma,NumGamma);
  GammaOrig = join_cols(GammaW,GammaM);
  Gamma = join_cols(GammaW,GammaM);
  delta = 1E-6 * delta;

//Next are Newton-Raphson iterations to find Gamma as the root of eqcond
//given beta
  eqcond = PSeqcond(beta, Gamma, S, X, Z, pmfW, pmfM, gw, gm);
  J = PSgeqcond(beta, Gamma, S, X, Z, pmfW, pmfM, gw, gm)+delta;
  GammaDiff = arma::inv(J)*eqcond;
  Gamma = Gamma - GammaDiff;
  while (dot(GammaDiff,GammaDiff) > 1E-9 && it < 40u) {
   eqcond = PSeqcond(beta, Gamma, S, X, Z, pmfW, pmfM, gw, gm);
   J = PSgeqcond(beta, Gamma, S, X, Z, pmfW, pmfM, gw, gm)+delta;
// Gamma = Gamma - arma::inv(J)*eqcond;
   GammaDiff = arma::inv(J)*eqcond;
   Gamma = Gamma - GammaDiff;
   if(std::isnan(Gamma(0))){
     Gamma = GammaOrig;
   }
   it++;
  }

  GammaW = Gamma.head(NumGammaW);
  GammaM = Gamma.tail(NumGammaM);

  arma::mat apmf(NumGammaW+1,NumGammaM+1);
  apmf(NumGammaW,NumGammaM)=-100.0;
  for (j = 0; j < NumGammaW; ++j) {
   for (k = 0; k < NumGammaM; ++k) {
    Wstar = 0.0;
    for (i = 0; i < S.n_slices; ++i) {
     Wstar += beta(i)*S(j,k,i);
    }
    Ustar = 0.0;
    for (i = 0; i < X.n_slices; ++i) {
     Ustar += beta(i+S.n_slices)*X(j,k,i);
    }
    Vstar = 0.0;
//  for (i = 0; i < S.n_slices; ++i) {
//   Vstar += beta(i)*S(k,j,i);
//  }
    for (i = 0; i < Z.n_slices; ++i) {
     Vstar += beta(i+S.n_slices+X.n_slices)*Z(k,j,i);
    }
    Wstar += Ustar + Vstar;
    apmf(j,k) = log(2.0)+(Wstar + GammaW(j)+ GammaM(k) + gm + gw + log(pmfW(j)*pmfM(k))-log((1.0+exp(GammaW(j)))*(1.0+exp(GammaM(k))))) ;
//  apmf(j,k) = (Wstar + GammaW(j)+ GammaM(k) + gm + gw + log(pmfW(j)*pmfM(k))-log((1.0+exp(GammaW(j)))*(1.0+exp(GammaM(k))))) ;
   }
  }

  for (j = 0; j < NumGammaW; ++j) {
    apmf(j,NumGammaM)=(GammaW(j) + gw + log(pmfW(j)) - log(1.0+exp(GammaW(j))));
  }
  for (k = 0; k < NumGammaM; ++k) {
    apmf(NumGammaW,k)=(GammaM(k) + gm + log(pmfM(k)) - log(1.0+exp(GammaM(k))));
  }
  return apmf;
}
