#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "rowwise_max_idx.h"

using idx_type = arma::uword;
using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "GS_TU")]]
arma::urowvec GS_TU(arma::mat U, arma::mat V) {
  arma::urowvec m(U.n_rows);
  arma::urowvec w(U.n_cols);
  arma::rowvec cmax(U.n_cols);
  arma::mat W(U.n_rows,U.n_cols); // Total utility for each pairing
  int reject, nit, nmax;
  
  int nw = U.n_rows;
  int nm = U.n_cols;
  nmax = nw*nm;
  W = U+V.t();
  
  reject = 1;
  nit = 0;
  while(reject && nit < nmax){
    nit++;
    reject=0;
    
    m = rowwise_max_idx(W); // index of man who has max joint utility with each woman
    w.zeros(); // all men begin with no partner
    for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
      cmax(j) = -100000.0; // the max utility each man can get in this round
    }
    for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
      if(W(i,m(i)) >= 0.0){ // best man better than single
        if(W(i,m(i)) > cmax(m(i))){ // man m(i) has higher utility for i than his current match ()
          w(m(i)) = i; // update man's preferred woman to i
          cmax(m(i)) = W(i,m(i)); // update the max utility a man could get from his proposers in this round
        }
        m(i) = m(i) + 1; 
      }else{
        m(i)=0;
      }
    }
    for (unsigned int j = 0u; j < static_cast<unsigned int>(nm); ++j) {
      if(cmax(j) > -10000.00 && W(w(j),j) > 0.0){ // If man j was proposed to and he prefers his woman to being single
        cmax(j)=W(w(j),j); //update man j's utility
        w(j) = w(j) + 1; //update man j's partner
      }else{
        cmax(j)=0.0; // man j's utility is 0
        w(j) = 0; // man j remains single
      }
    }
    for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
      if(m(i) > 0 && W(i,m(i)-1) < (cmax(m(i)-1)-1.0e-14)){// If woman i proposes to man j and man's utility for her is less than the maximum utility he could have
        W(i,m(i)-1) = -10000.0; // woman is rejected and she won't propose to man again
        reject++;
      }
    }
  }
  m = rowwise_max_idx(W); // index of max of each row
  for (unsigned int i = 0u; i < static_cast<unsigned int>(nw); ++i) {
    if(W(i,m(i)) < 0.0){
      m(i)=0;
    }else{
      m(i)=m(i)+1;
    }
  }
  return m; 
}
