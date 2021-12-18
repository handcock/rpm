#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rmath.h>

using idx_type = arma::uword;
using namespace std;
using namespace arma;

// @export
// [[Rcpp::export(name = "rowwise_imax_idx")]]
arma::urowvec rowwise_imax_idx(const arma::imat A) {
    std::vector<idx_type> res;
    for (idx_type i = 0; i != A.n_rows; ++i) {
        idx_type col_idx;
        A.row(i).max(col_idx);
        res.push_back(col_idx);
    }
    return res;
}
