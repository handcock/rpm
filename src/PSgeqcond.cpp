#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "PSgeqcond")]]
arma::mat PSgeqcond(arma::vec beta, arma::vec Gamma, arma::cube S, arma::cube X, arma::cube Z, arma::vec pmfW, arma::vec pmfM, double gw, double gm) {
  int Sdm=S.n_slices;
  int Wdm=X.n_slices;
  int Mdm=Z.n_slices;
  int NumGamma=Gamma.size();
  int NumGammaW=pmfW.size();
  int NumGammaM=pmfM.size();
  int i,j,k,l;
  double Ustar, Vstar, Wstar;
  double pfps;

  arma::vec GammaW(NumGammaW);
  arma::vec GammaM(NumGammaM);
  GammaW = Gamma.subvec(0, NumGammaW-1);
  GammaM = Gamma.subvec(NumGammaW, NumGamma-1);

// rows are constraints, columns parameters
  arma::mat gf(NumGamma,NumGamma);
  gf.zeros();

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
    pfps = exp(Wstar+GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)));
    gf(j,NumGammaW+k) = pfps / (1.0 + exp(GammaM(k)));
   }
   gf(j,j)= exp(-GammaW(j));
  }

// Now g(*,z)
  for (k = 0; k < NumGammaM; ++k) {
   for (j = 0; j < NumGammaW; ++j) {
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
    pfps = exp(Wstar+GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j)));
    gf(k+NumGammaW,j) = pfps / (1.0 + exp(GammaW(j)));
   }
   gf(k+NumGammaW,NumGammaW+k)= exp(-GammaM(k));
  }
  for (l = 0; l < NumGamma; ++l) {
   for (i = 0; i < NumGamma; ++i) {
    if(std::isnan(gf(l,i)) | !std::isfinite(gf(l,i))){gf(l,i) = 100.0;}
   }
  }
  return gf;
}
