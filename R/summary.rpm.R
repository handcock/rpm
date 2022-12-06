#' Summarizing rpm Model Fits
#'
#' [base::summary()] method for [rpm()] fits.
#'
#' @order 1
#' @param object an object of class \code{\link{rpm}}, usually, a result of a call to
#'   [rpm()].
#' @param digits significant digits for the coefficients;default=
#'   max(3,getOption("digits")-3), but the hard-coded value is 5.
#' @param correlation logical; if `TRUE`, the correlation matrix of the
#'   estimated parameters is returned and printed.
#' @param covariance logical; if `TRUE`, the covariance matrix of the estimated
#'   parameters is returned and printed.
#' @param include.single logical; if `TRUE`, include in the summary table the 
#'   coefficients of the log-odds of being single for each category of women and men.
#' @param ... For [summary.rpm()] additional arguments are passed to
#'   [logLik.rpm()]. For [print.summary.rpm()], to [stats::printCoefmat()].
#'
#' @details [summary.rpm()] tries to be smart about formatting the
#' coefficients, standard errors, etc.
#'
#' @return The function [summary.rpm()] computes and returns a list of summary
#'   statistics of the fitted [rpm()] model given in `object`. Note that for
#'   backwards compatibility, it returns the coefficient table.
#'
#'   The returned object is a list of class "summary.rpm" with the following
#'   elements:
#'   
#' \item{formula}{ERGM model formula}
#' \item{digits}{the 'digits' inputted to <summary.rpm> or the default value (despite the fact the digits will be 5)}
#' \item{correlation, covariance}{whether to print correlation/covariance matrices of the estimated parameters}
#' \item{iterations}{object$iterations}
#' \item{control}{the [control.rpm()] object used}
#' \item{samplesize}{MCMC sample size}
#' \item{message}{optional message on the validity of the standard error estimates}
#' \item{aic.null,bic.null}{values of AIC and BIC for the null model}
#' \item{aic, bic}{values of AIC and BIC}
#' \item{coefficients}{data frames with model parameters and associated statistics}
#' \item{asycov}{asymptotic covariance matrix}
#' \item{asyse}{asymptotic standard error matrix}
#' \item{offset, drop, estimate, iterations, mle.lik, null.lik}{
#' see documentation of the object returned by [rpm()]
#' }
#'
#' @seealso The model fitting function [rpm()], [print.rpm()], and
#'   [base::summary()]. Function [stats::coef()] will extract the data frame of
#'   coefficients with standard errors, t-statistics and p-values.
#'
#'
#'
#'
#' @keywords regression models
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
#' @importFrom MASS ginv
#' @method summary rpm
#' @export
summary.rpm <- function (object, ..., 
                        digits = max(3, getOption("digits") - 3),
                        correlation=FALSE, covariance=FALSE,
                        include.single=TRUE)
{
  control <- object$control
  
  nobs<- object$nobs
  df <- length(object$coefficients)

  devtext <- "Deviance:"
  mle.lik<-try(stats::logLik(object,...), silent=TRUE)
  null.lik<-try(logLikNull(object,...), silent=TRUE)

  if(include.single){
   if(length(object$coefficients) < length(object$solution)) object$coefficients <- c(object$coefficients,object$LOGODDS_SW,object$LOGODDS_SM)
   object$covar <- object$ext.covar
  }
  if(is.null(object$hessian) && is.null(object$covar)){
    object$covar <- diag(NA, nrow=length(object$coefficients))
  }
  
  if(is.null(object$covar)){
    asycov <- .catchToList(MASS::ginv(-object$hessian))
    if(!is.null(asycov$error)){
      asycov <- diag(1/diag(-object$hessian))
    }else{
      asycov <- asycov$value
    }
  }else{
    asycov <- object$covar
  }
  asycov <- asycov[seq_along(object$coefficients),seq_along(object$coefficients)]
  colnames(asycov) <- rownames(asycov) <- names(object$coefficients)
  
  asyse <- diag(asycov)
  asyse[asyse<0|is.infinite(object$coefficients)] <- NA
  asyse <- sqrt(asyse)
  asyse <- matrix(asyse, ncol=length(asyse))
  colnames(asyse) <- colnames(asycov)
  
  ans <- list(formula=object$formula,
              digits=digits, correlation=correlation,
              covariance=covariance,
              iterations=object$iterations,
              control=object$control)
  
  ans$null.lik.0 <- is.na(null.lik)
  
  rdf <- nobs - df
  tval <- object$coefficients / asyse
  pval <- rep(NA,length=length(tval))
  pval[!is.na(tval) & !is.nan(tval)] <- 2 * stats::pt(q=abs(tval[!is.na(tval) & !is.nan(tval)]), df=rdf, lower.tail=FALSE)
  
  count <- 1
  templist <- NULL
  while (count <= length(names(object$coefficients)))
  {
    templist <- append(templist,c(object$coefficients[count],
                                  asyse[count],pval[count]))
    count <- count+1
  }
  
  tempmatrix <- matrix(templist, ncol=3,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "p-value")
  rownames(tempmatrix) <- names(object$coefficients)
  
  
  ans$devtable <- c("",apply(cbind(paste(format(c("    Null", "Residual"), width = 8), devtext), 
                                   format(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik), digits = digits), " on",
                                   format(c(nobs-1, rdf), digits = digits)," degrees of freedom\n"), 
                             1, paste, collapse = " "),"\n")
  
  ans$aic_null <- stats::AIC(null.lik)
  ans$bic_null <- stats::BIC(null.lik)

  ans$aic <- stats::AIC(mle.lik)
  ans$bic <- stats::BIC(mle.lik)
  
  ans$coefficients <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.rpm"
  ans
}

#' @rdname summary.rpm
#' @order 2
#' 
#' @param x object of class `summary.rpm` returned by [summary.rpm()].
#' @param digits significant digits for coefficients. The default is max(3, getOption("digits")-3).
#' @param correlation logical whether the correlation matrix of the estimated parameters
#'                should be printed (T or F); default=FALSE
#' @param covariance logical whether the covariance matrix of the estimated parameters
#'                should be printed (T or F); default=FALSE
#' @param signif.stars whether to print dots and stars to signify
#'   statistical significance. See [print.summary.lm()].
#' @param eps.Pvalue \eqn{p}-values below this level will be printed
#'   as "<`eps.Pvalue`".
#' @param
#'   print.formula,print.fitinfo,print.coefmat,print.message,print.deviances,print.drop,print.header
#'   which components of the fit summary to print.
#'   
#' @details The default printout of the summary object contains the
#'   call, number of iterations used, null and residual deviances, and
#'   the values of AIC and BIC.
#'   The coefficient table contains the following
#'   columns:
#'   
#'   - `Estimate`, `Std. Error` - parameter estimates and their standard errors
#'   - `z value`, `Pr(>|z|)` - z-test and p-values
#'    
#' @method print summary.rpm
#' @export
print.summary.rpm <- function (x, 
                              digits = max(3, getOption("digits") - 3),
                              correlation=FALSE, covariance=FALSE,
                              signif.stars= getOption("show.signif.stars"),
                              eps.Pvalue=0.0001, print.header=TRUE, print.formula=FALSE, print.fitinfo=TRUE, print.coefmat=TRUE, print.message=TRUE, print.deviances=TRUE, print.drop=TRUE, ...){
  if(missing(digits)) digits <- x$digits
  
  control <- x$control
  if(print.header){
    cat("\n==========================\n")
    cat("Summary of RPM model fit\n")
    cat("==========================\n\n")
  }
  
  if(print.formula){
    cat("Formula:   ")
    print(x$formula)
    cat("\n")
  }
  
  if(print.fitinfo){
    if (!is.null(x$iterations)) {
      cat("Iterations: ", x$iterations, "\n")
    }
  }
  
  if(print.coefmat){
    stats::printCoefmat(x$coefficients, digits=digits, signif.stars=signif.stars,
                        P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
                        eps.Pvalue=eps.Pvalue, ...)
  }
  
  if(print.message){
    if(!is.null(x$message)){ 
      cat(x$message)
    }
    cat("\n")
  }
  
  if(print.deviances){
    if(!is.null(x$devtable)){
      cat(x$devtable)
      
      if(x$null.lik.0) cat("Note that the null model likelihood and deviance are defined to be 0.\n\n")
      
      cat(paste("null AIC:", format(x$aic_null, digits = digits), "  ", 
                "null BIC:", format(x$bic_null, digits = digits), "\n", sep=" "))
      cat(paste("     AIC:", format(x$aic, digits = digits), "  ", 
                "     BIC:", format(x$bic, digits = digits), "  ",
                "(Smaller is better.)", "\n", sep=" "))
    } 
  }
  
  if((missing(covariance)&x$covariance)|covariance == TRUE){
    cat("Asymptotic covariance matrix:\n")
    print(x$asycov)
  }
  
  if((missing(correlation)&x$correlation)|correlation == TRUE){
    cat("\nAsymptotic correlation matrix:\n")
    asycor <- x$asycov / crossprod(x$asyse)
    dimnames(asycor) <- dimnames(x$asycov)
    print(asycor)
  }
  
  invisible(x)
}

#' @method print summary.rpm
#' @export
show.summary.rpm <- print.summary.rpm
