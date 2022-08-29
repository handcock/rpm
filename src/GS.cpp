#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "rowwise_max_idx.h"

using idx_type = arma::uword;
using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "GS")]]
arma::urowvec GS(arma::mat U, arma::mat V) {
  arma::urowvec m(U.n_rows);
  arma::urowvec w(U.n_cols);
  arma::vec cmax(U.n_cols);
  int reject, nit, nmax;

  int nw = U.n_rows;
  int nm = U.n_cols;
  nmax = nw*nm;
  
  reject = 1;
  nit = 0;
  while(reject && nit < nmax){
   nit++;
   reject=0;

   m = rowwise_max_idx(U); // index of man with maximum utility for woman
   w.zeros();
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
     cmax(j) = -100000.0;
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
     if(U(i,m(i)) >= 0.0){ // best man better than single
      if(V(m(i),i) > cmax(m(i))){ // man m(i) has higher utility for i than his current woman
        w(m(i)) = i; // best woman for man
        cmax(m(i)) = V(m(i),i); // utility of best woman for man
      }
      m(i) = m(i) + 1;
     }else{
      m(i)=0;
     }
   }
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
    if(cmax(j) > -10000.00 && V(j,w(j)) > 0.0){
       cmax(j)=V(j,w(j));
       w(j) = w(j) + 1;
    }else{
       cmax(j)=0.0;
       w(j) = 0;
    }
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
     if(m(i) > 0 && V(m(i)-1,i) < (cmax(m(i)-1)-1.0e-14)){
       U(i,m(i)-1) = -10000.0;
       reject++;
     }
   }
  }
  m = rowwise_max_idx(U); // index of max of each row
  for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
   if(U(i,m(i)) < 0.0){
     m(i)=0;
   }else{
     m(i)=m(i)+1;
   }
  }
//Rprintf("GS: nit %d of nmax %d, reject %d \n",nit, nmax, reject);
  return m; 
}
