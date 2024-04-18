#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "PSeqcond.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::vec auxGamma(arma::vec beta, arma::vec GammaW, arma::vec GammaM, arma::cube S, arma::cube X, arma::cube Z, arma::vec pmfW, arma::vec
pmfM, double gw, double gm) {
  unsigned int it=0u;
//unsigned int NumBeta=beta.size();
  unsigned int NumGammaW=GammaW.size();
  unsigned int NumGammaM=GammaM.size();
  unsigned int NumGamma=NumGammaW+NumGammaM;

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
//Rprintf("it %d eqcond(0) %f eqcond(1) %f eqcond(2) %f eqcond(3) %f\n", it, eqcond(0), eqcond(1), eqcond(2), eqcond(3) );
   if(std::isnan(Gamma(0))){
     Gamma = GammaOrig;
   }
   it++;
  }

  return Gamma;
}
