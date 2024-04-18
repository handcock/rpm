#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp attribute tag require to create interface to R.
//' @export
// [[Rcpp::export]]
NumericVector eqcond(NumericVector beta, NumericVector GammaW, NumericVector GammaM, NumericVector Sd, NumericVector Xd, NumericVector Zd, IntegerVector Sdim, IntegerVector Xdim, IntegerVector Zdim, NumericVector pmfW, NumericVector pmfM, NumericMatrix pmf, NumericMatrix counts, double gw, double gm, int constraints){
  arma::cube S(Sd.begin(), Sdim[0], Sdim[1], Sdim[2]);
  arma::cube X(Xd.begin(), Xdim[0], Xdim[1], Xdim[2]);
  arma::cube Z(Zd.begin(), Zdim[0], Zdim[1], Zdim[2]);
  int Sdm=Sdim(2);
  int Wdm=Xdim(2);
  int Mdm=Zdim(2);
  int NumGammaW=GammaW.size();
  int NumGammaM=GammaM.size();
  int NumGamma=NumGammaW+NumGammaM;
  int i, j, k;
  double Ustar, Vstar, Wstar;
  double pfps;
//double checksum;
  
//Rprintf("constraints %d \n",constraints,NumGamma+NumGammaM*(constraints==0));
  NumericVector pfp(NumGamma+NumGammaM*(constraints==0));
   
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
//  for (i = 0; i < Sdm; ++i) {
//   Vstar += beta(i)*S(k,j,i);
//  }
    for (i = 0; i < Mdm; ++i) {
     Vstar += beta(i+Sdm+Wdm)*Z(k,j,i);
    }
    Wstar += Ustar + Vstar;
    pfps += exp(Wstar+GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)));
   }
   pfp(j)=pfps-(exp(-GammaW(j)));
// pfp(j)=log(pfps)+GammaW(j);
// pfp(j)=pfps*(exp(GammaW(j)))-1.0;
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
//  for (i = 0; i < Sdm; ++i) {
//   Vstar += beta(i)*S(k,j,i);
//  }
    for (i = 0; i < Mdm; ++i) {
     Vstar += beta(i+Sdm+Wdm)*Z(k,j,i);
    }
    Wstar += Ustar + Vstar;
    pfps += exp(Wstar+GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j)));
   }
   pfp(k+NumGammaW)=pfps-(exp(-GammaM(k)));
// pfp(k+NumGammaW)=log(pfps)+GammaM(k);
// pfp(k+NumGammaW)=pfps*(exp(GammaM(k)))-1.0;
  }
//pfps=0.0;
//for (j = 0; j < NumGammaW; ++j) {
// pfps=pmf(j,NumGammaM)-(exp(GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j))));
//}
//Rprintf("W sum %f\n",pfps);
//pfp(NumGamma)=pfps;
//pfps=0.0;
//for (k = 0; k < NumGammaM; ++k) {
// pfps=pmf(NumGammaW,k)-(exp(GammaM(k)+gw)*pmfM(k)*0.5 / (1.0 + exp(GammaM(k))));
//}
//pfp(NumGamma+1)=pfps;
//for (k = 0; k < NumGammaM; ++k) {
// pfp(k+NumGamma+NumGammaW)=pmf(NumGammaW,k)-(exp(GammaM(k)+gw)*pmfM(k)*0.5 / (1.0 + exp(GammaM(k))));
//}
//pfps=2.0*exp(gw)-1.0;
//for (j = 0; j < NumGammaW; ++j) {
// pfps-=exp(GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j)));
//}
//for (k = 0; k < NumGammaM; ++k) {
// pfps+=exp(GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)));
//}
//pfp(NumGamma)=pfps;

// Next is the (# single W in x)/(#people)
//for (j = 0; j < NumGammaW; ++j) {
// pfp(j+NumGamma)=pmf(j,NumGammaM)/(exp(GammaW(j)+gw)*pmfW(j) / (1.0 + exp(GammaW(j)))) - 1.0;
//}

  if(constraints==0){
   for (k = 0; k < NumGammaM; ++k) {
//  pfp(k+NumGamma)=pmf(NumGammaW,k)/(exp(GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)))) - 1.0;
    pfp(k+NumGamma)=pmf(NumGammaW,k) - (exp(GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k))));
   }
  }
// for (k = 0; k < NumGammaM; ++k) {
//Rprintf("k %d pmf(NumGammaW,k) %f LHS %f\n",k,pmf(NumGammaW,k),(exp(GammaM(k)+gm)*pmfM(k) / (1.0 + exp(GammaM(k)))));
// }
//for (j = 0; j < NumGammaW; ++j) {
//  Rprintf("j %d pmf %f p %f diff %f\n",j,pmf(j,NumGammaM),exp(GammaW(j)+gw)*pmfW(j)*0.5 / (1.0 + exp(GammaW(j))),pfp(j+NumGamma));
//}
//for (k = 0; k < (NumGammaM-1); ++k) {
// pfp(k+NumGamma+NumGammaW-1)=pmf(NumGammaW,k)-(exp(GammaM(k)+gm)*pmfM(k)*0.5 / (1.0 + exp(GammaM(k))));
//}
//for (j = 0; j < (NumGamma+NumGammaM); ++j) {
//  checksum=pfp(j);
//  if(std::isnan(checksum) || !std::isfinite(checksum)){
//   if(std::signbit(checksum)){
//    pfp(j)= -100.0-1.0;
//   }else{
//    pfp(j)=  100.0-1.0;
//   }
//  }else{
//   if(checksum < -100.0){
//    pfp(j)= -100.0-1.0;
//   }else{
//    if(checksum > 100.0){
//     pfp(j)= 100.0-1.0;
//    }
//   }
//  }
//}
//Rprintf("checksum %f\n",fsum);
  return pfp;
}
