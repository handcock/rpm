#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "rowwise_imax_idx.h"

using idx_type = arma::uword;
using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "GSi_TU")]]
arma::urowvec GSi_TU(arma::imat U, arma::imat V) {
  arma::urowvec m(U.n_rows);
  arma::urowvec w(U.n_cols);
  arma::irowvec cmax(U.n_cols);
  // arma::imat W(U.n_rows,U.n_cols); // Total utility for each pairing
  int reject, nit, nmax;

  U = U+V.t(); // Now the total utility of each pairing

  int nw = U.n_rows;
  int nm = U.n_cols;
  nmax = nw*nm;
  reject = 1;
  nit = 0;
  while(reject && nit < nmax){
   nit++;
   reject=0;

   m = rowwise_imax_idx(U); // index of man with maximum utility for woman
   w.zeros();
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
     cmax(j) = -100000;
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
     if(U(i,m(i)) >= 0){ // best man better than single
      if(U(m(i),i) > cmax(m(i))){ // man m(i) has higher utility for i than his current woman
        w(m(i)) = i; // best woman for man
        cmax(m(i)) = U(m(i),i); // utility of best woman for man
      }
      m(i) = m(i) + 1;
     }else{
      m(i)=0;
     }
   }
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
    if(cmax(j) > -10000 && U(w(j),j) > 0){
       cmax(j)=U(w(j),j);
       w(j) = w(j) + 1;
    }else{
       cmax(j)=0;
       w(j) = 0;
    }
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
     if(m(i) > 0 && U(i,m(i)-1) < cmax(m(i)-1)){
       U(i,m(i)-1) = -10000;
       reject++;
     }
   }
  }
  m = rowwise_imax_idx(U); // index of max of each row
  for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
   if(U(i,m(i)) < 0){
     m(i)=0;
   }else{
     m(i)=m(i)+1;
   }
  }
  return m; 
}

