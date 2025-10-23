#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>

using idx_type = arma::uword;
using namespace std;
using namespace arma;

//' Row-wise argmax (0-based column indices)
//'
//' Returns the column index of the maximum element in each row of an integer matrix.
//' Uses Armadillo's `index_max()` to avoid deprecated `.max(uword&)`.
//'
//' @param A Integer matrix.
//' @return A row vector of 0-based column indices of the maxima.
//' @export
// [[Rcpp::export(name = "rowwise_imax_idx")]]
arma::urowvec rowwise_imax_idx(const arma::imat& A) {
    arma::urowvec res(A.n_rows);
    for (idx_type i = 0; i < A.n_rows; ++i) {
        res(i) = A.row(i).index_max();
    }
    return res;
}
