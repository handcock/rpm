#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "rowwise_imax_idx.h"

using idx_type = arma::uword;
using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "GSi_NTU_O")]]
arma::urowvec GSi_NTU_O(arma::imat U, arma::imat V, arma::imat O, arma::imat Ua, int nmax) {
  arma::urowvec m(U.n_rows);
  arma::urowvec w(U.n_cols);
  arma::irowvec cmax(U.n_cols);
  unsigned int reject, nit, ji;

  int nw = U.n_rows;
  int nm = U.n_cols;
//nmax = nw*nm;
  reject = 1;
  nit = 0;
  while(reject && nit < nmax){
   nit++;
   reject=0;

 //m = rowwise_imax_idx(U); // index of man with maximum utility for woman
   w.zeros();
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
     cmax(j) = -1000000;
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
//Rprintf("GSi: Onrow %d Oncol %d nw %d i %d O(i,0) %d \n",O.n_rows,O.n_cols,nw,i, O(i,0));
     m(i) = O(i,0)-1;
     if(U(i,m(i)) >= 0){ // best man better than single
      if(V(m(i),i) > cmax(m(i))){ // man m(i) has higher utility for i than his current woman
        w(m(i)) = i; // best woman for man
        cmax(m(i)) = V(m(i),i); // utility of best woman for man
      }
      m(i)++;
     }else{
      m(i)=0;
     }
   }
   for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
    if(cmax(j) > -100000 && V(j,w(j)) > 0){
       cmax(j)=V(j,w(j));
       w(j)++;
    }else{
       cmax(j)=0;
       w(j) = 0;
    }
   }
   for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
     if(m(i) > 0 && V(m(i)-1,i) < cmax(m(i)-1)){
//Rprintf("made: i %d m(i) %d \n",i, m(i));
       U(i,m(i)-1) = -100000;
// manually update U and m

    // ji = Ua(i,m(i)-1)-1;
       ji = 0;

// Rprintf("ji: %d O(i,ji)-1 %d m(i) %d\n", ji, O(i,ji), m(i));

   //  while(O(i,ji) != m(i) ){Rprintf("ji %d O(i,ji) %d m(i) %d\n",ji,O(i,ji),m(i)); ji--;}

    // while(O(i,ji) != m(i) ){Rprintf("ji %d O(i,ji) %d m(i) %d\n",ji,O(i,ji),m(i)); ji++;}
       while(O(i,ji) != m(i) ){ji++;}
    // while(O(i,ji) != m(i) ){ji--;}
//Rprintf("done: i %d m(i) %d ji %d O(i,ji) %d\n",i, m(i),ji,O(i,ji));
       for (unsigned int k = ji; k < static_cast<unsigned int>(nm-3); ++k) {
//Rprintf("GSi: k %d O(i,k) %d \n",k, O(i,k));
         O(i, k) = O(i, k+1);
       }
       reject++;
     }
   }
  }
//Rprintf("GSi: nit %d of nmax %d, reject %d \n",nit, nmax, reject);
  m = rowwise_imax_idx(U); // index of man with maximum utility for woman
  for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
   if(U(i,m(i)) < 0){
     m(i)=0;
   }else{
     m(i)++;
   }
  }
//Rprintf("GSi: nit %d of nmax %d, reject %d \n",nit, nmax, reject);
  return m; 
}

