#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "stablearma")]]
unsigned int stablearma(arma::mat U, arma::mat V, arma::ivec m, arma::ivec w) {
  
// woman's perspective
 for (unsigned int i = 0u; i < static_cast<unsigned int>(U.n_rows); ++i) {
   if(m(i)==0){
//   i is single
     for (unsigned int j = 0u; j < static_cast<unsigned int>(V.n_rows); ++j) {
       if(U(i,j+1) > U(i,0) && V(j,i+1) > V(j,w(j)) ){ return(0);}
     }
   }else{
//   i is partnered with m(i), so w(m(i))=i
     if( U(i,0) > U(i,m(i)) ){return(0);}
     for (unsigned int j = 0u; j < static_cast<unsigned int>(V.n_rows); ++j) {
       if( (j != ((unsigned)(m(i))-1u)) && U(i,j+1) > U(i,m(i)) && V(j,i+1) > V(j,w(j)) ){return(0);}
     }
   }
 }
 return(1);
}
