#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
NumericMatrix geqcond(NumericVector beta, NumericVector GammaW, NumericVector GammaM, NumericVector Sd, NumericVector Xd, NumericVector Zd, IntegerVector Sdim, IntegerVector Xdim, IntegerVector Zdim, NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, NumericMatrix counts, double gw, double gm, int constraints){
  arma::cube S(Sd.begin(), Sdim[0], Sdim[1], Sdim[2]);
  arma::cube X(Xd.begin(), Xdim[0], Xdim[1], Xdim[2]);
  arma::cube Z(Zd.begin(), Zdim[0], Zdim[1], Zdim[2]);
  int Sdm=Sdim(2);
  int Wdm=Xdim(2);
  int Mdm=Zdim(2);
  int NumBeta=beta.size();
  int NumGammaW=GammaW.size();
  int NumGammaM=GammaM.size();
  int NumGamma=NumGammaW+NumGammaM;
  int i,j,k,l;
  double Ustar, Vstar, Wstar;
  double pfps;

// rows are constraints, columns parameters
  NumericMatrix gf(NumGamma+NumGammaM*(constraints==0),NumBeta+NumGamma+1);
  for (l = 0; l < NumGamma+NumGammaM*(constraints==0); ++l) {
   for (i = 0; i < NumBeta+NumGamma+1; ++i) {
    gf(l,i)=0.0;
   }
  }

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
    gf(j,NumBeta+NumGammaW+k) = pfps / (1.0 + exp(GammaM(k)));
    for (i = 0; i < Sdm; ++i) {
     gf(j,i) += S(j,k,i)*pfps;
    }
    for (i = 0; i < Wdm; ++i) {
     gf(j,i+Sdm) += X(j,k,i)*pfps;
    }
    for (i = 0; i < Mdm; ++i) {
     gf(j,i+Sdm+Wdm) += Z(k,j,i)*pfps;
    }
    gf(j,NumBeta+NumGamma) += -pfps*exp(gw-gm);
   }
   gf(j,NumBeta+j)= exp(-GammaW(j));
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
    gf(k+NumGammaW,NumBeta+j) = pfps / (1.0 + exp(GammaW(j)));
    for (i = 0; i < Sdm; ++i) {
     gf(k+NumGammaW,i) += S(j,k,i)*pfps;
    }
    for (i = 0; i < Wdm; ++i) {
     gf(k+NumGammaW,i+Sdm) += X(j,k,i)*pfps;
    }
    for (i = 0; i < Mdm; ++i) {
     gf(k+NumGammaW,i+Sdm+Wdm) += Z(k,j,i)*pfps;
    }
    gf(k+NumGammaW,NumBeta+NumGamma) += pfps;
   }
   gf(k+NumGammaW,NumBeta+NumGammaW+k)= exp(-GammaM(k));
  }
//
  if(constraints==0){
   for (l = 0; l < NumGammaM; ++l) {
    pfps = 1.0 + exp(GammaM(l));
    gf(l+NumGamma,NumBeta+NumGammaW+l) = - exp(GammaM(l)+gm)*pmfM(l) / (pfps*pfps);
   }
  }
  return gf;
}
