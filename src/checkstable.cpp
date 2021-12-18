#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export(name = "checkstable")]]
int checkstable(const NumericMatrix &U, const NumericMatrix &V, const IntegerVector &m, const IntegerVector &w ) {
  
 for (int i = 0u; i < static_cast<int>(U.nrow()); ++i) {
   if(m(i)==0){
//   i is single
     for (int j = 0u; j < static_cast<int>(V.nrow()); ++j) {
       if(U(i,j+1) > U(i,0) && V(j,i+1) > V(j,w(j)) ){ return 0;}
     }
   }else{
//   i is partnered with m(i), so w(m(i))=i
     if( U(i,0) > U(i,m(i)) ){return(0);}
     for (int j = 0u; j < static_cast<int>(V.nrow()); ++j) {
       if(((j+1) != m(i)) && U(i,j+1) > U(i,m(i)) && V(j,i+1) > V(j,w(j)) ){return 0;}
     }
   }
 }
//  for (unsigned int j = 0u; j < static_cast<unsigned int>(V.nrow()); ++j) {
//    if(w(j)==0){
// //   j is single
//      for (unsigned int i = 0u; i < static_cast<unsigned int>(U.nrow()); ++i) {
//        if(V(j,i+1) > V(j,0) && U(i,j+1) > U(i,m(i))){return 0;}
//      }
//    }else{
// //   j is partnered with w(j), so m(w(j))=j
//      if( V(j,0) > V(j,w(j)) ){return(0);}
//      for (unsigned int i = 0u; i < static_cast<unsigned int>(U.nrow()); ++i) {
//        if(((i+1) != w(j)) && V(j,i+1) > V(j,w(j)) && U(i,j+1) > U(i,m(i)) ){return 0;}
//      }
//    }
//  }
 return 1;
}
