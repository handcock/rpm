#' Modeling of Revealed Preferences Matchings
#' 
#' An integrated set of tools to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using a large-population likelihood based on ideas from Menzel (2015).
#' 
#' For a complete list of the functions, use \code{library(help="rpm")} or
#' read the rest of the manual.
#' 
#' When publishing results obtained using this package the original authors are
#' to be cited as:
#' 
#' Mark S. Handcock, Ryan M. Admiraal, Fiona C. Yeung, Heide M. Jackson, Michael S. Rendall and Shuchi Goyal (2022) \pkg{rpm}: Modeling of Revealed Preferences Matchings
#' R package, Los Angeles, CA.  Version 0.70, \url{https://github.com/handcock/rpm}.
#' 
#' All programs derived from this package must cite it. For complete citation
#' information, use\cr \code{citation(package="rpm")}.
#' 
#' For details on how to construct data for input to \code{rpm()} see the documentation:
#' 
#' \code{help(fauxmatching)}
#'
#' For information on the current terms that can be used in formulas for \code{rpm()} see the documentation:
#' 
#' \code{help("rpm-terms")}
#' @name rpm-package
#' @docType package
#' @author Mark S. Handcock <handcock@stat.ucla.edu>
#' @import Rcpp
#' @import abind
#' @import methods
#' @import coda
#' @import foreach
#' @importFrom Rcpp evalCpp
#' @importFrom doRNG "%dorng%" registerDoRNG
#' @importFrom foreach "foreach"
#' @import ggplot2
#' @import graphics
#' @import dplyr
#' @importFrom future plan
#' @importFrom utils capture.output
#' @useDynLib rpm
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' \donttest{
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled",sampling_design="stock-flow")
#' summary(fit)
#' }
#' # For details on how to construct data for input:
#' help(fauxmatching)
#' # For information on the current terms that can be used in formulas:
#' help("rpm-terms")
#' 
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.18637/jss.v024.i07} 
#' Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' 
NULL
globalVariables(c("gender","cluster","multisession","multicore","Mean"))
