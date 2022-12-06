#' A \code{\link{logLik}} method for [`rpm`] fits.
#' 
#' A function to return the log-likelihood associated with an
#' \code{\link[=rpm.object]{rpm}} fit
#' 
#' @param object An \code{\link[=rpm.object]{rpm}} fit, returned by
#'   \code{\link{rpm}}.
#' @param \dots Other arguments to the likelihood functions.
#' @return 
#'   a \code{\link{logLik}} object.
#' @seealso \code{\link{logLik}}, \code{\link{logLikNull}}
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' \donttest{
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled",sampling_design="stock-flow")
#' logLik(fit)
#' }
#' 
#' @export
logLik.rpm<-function(object, ...){

  llk<- object$loglik
  
  class(llk)<-"logLik"
  attr(llk,"df")<-(object$NumBeta+1)
  attr(llk,"nobs")<- object$nobs

  llk
}

#' Calculate the null model likelihood
#'
#' @param object a fitted model.
#' @param ... further arguments to lower-level functions.
#' 
#' \code{logLikNull} computes, when possible the log-probability of
#' the data under the null model (reference distribution).
#' 
#' @return
#' \code{logLikNull} returns an object of type \code{\link{logLik}} if it is
#' able to compute the null model probability, and \code{NA} otherwise.
#' @export
logLikNull <- function(object, ...) UseMethod("logLikNull")

#' @describeIn logLikNull A method for [`rpm`] fits to compute the null likelihood (that is, relative to the constant only model).
#' @export
logLikNull.rpm <- function(object, ...){

  llk<- object$loglik.null
  
  class(llk)<-"logLik"
  attr(llk,"df")<-1
  attr(llk,"nobs")<- object$nobs

  llk
}
