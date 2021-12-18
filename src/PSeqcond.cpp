#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "PSeqcond")]]
arma::vec PSeqcond(arma::vec beta, arma::vec Gamma, arma::cube S, arma::cube X, arma::cube Z, arma::vec pmfW, arma::vec pmfM, double gw, double gm) {
  int Sdm=S.n_slices;
  int Wdm=X.n_slices;
  int Mdm=Z.n_slices;
  int NumGamma=Gamma.size();
  int NumGammaW=pmfW.size();
  int NumGammaM=pmfM.size();
  int i, j, k;
  double Ustar, Vstar, Wstar;
  double pfps;
  
  arma::vec pfp(NumGamma);
   
  arma::vec GammaW(NumGammaW);
  arma::vec GammaM(NumGammaM);
  GammaW = Gamma.head(NumGammaW);
  GammaM = Gamma.tail(NumGammaM);

  for (j = 0; j < NumGammaW; ++j) {
   pfps=0.0;
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
    pfps += exp(Wstar+GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)));
   }
   pfp(j)=pfps-(exp(-GammaW(j)));
  }

  for (k = 0; k < NumGammaM; ++k) {
   pfps=0.0;
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
    pfps += exp(Wstar+GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j)));
   }
   pfp(k+NumGammaW)=pfps-(exp(-GammaM(k)));
  }
  for (k = 0; k < NumGamma; ++k) {
   if(std::isnan(pfp(k)) | !std::isfinite(pfp(k))){pfp(k) = -10.0;}
  }
  return pfp;
}
