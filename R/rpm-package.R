#' Modeling of Revealed Preferences Matchings
#' 
#' An integrated set of tools to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' The "rpm" packages is part of
#' the "statnet" suite of packages for the analysis of social network 
#' data.  For a list of functions type:
#' help(package='rpm')
#' 
#' For a complete list of the functions, use \code{library(help="rpm")} or
#' read the rest of the manual.
#' 
#' When publishing results obtained using this package the original authors are
#' to be cited as:
#' 
#' Handcock, Mark S. (2020) \pkg{rpm}: Modeling of Revealed Preferences Matchings
#' R package, Los Angeles, CA.  Version 0.51, \url{http://statnet.org}.
#' 
#' All programs derived from this package must cite it. For complete citation
#' information, use\cr \code{citation(package="rpm")}.
#' 
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
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled",sampling_design="stock-flow")
#' summary(fit)
#' 
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' 
NULL
globalVariables(c("gender","cluster","multisession","multicore","Mean"))
