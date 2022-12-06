#' ANOVA for rpm Fits
#' 
#' Compute an analysis of variance table for one or more rpm fits.
#' 
#' Specifying a single object gives a sequential analysis of variance table for
#' that fit.  That is, the reductions in the residual sum of squares as each
#' term of the formula is added in turn are given in the rows of a table, plus
#' the residual sum of squares.
#' 
#' The table will contain F statistics (and P values) comparing the mean square
#' for the row to the residual mean square.
#' 
#' If more than one object is specified, the table has a row for the residual
#' degrees of freedom and sum of squares for each model.  For all but the first
#' model, the change in degrees of freedom and sum of squares is also given.
#' (This only make statistical sense if the models are nested.)  It is
#' conventional to list the models from smallest to largest, but this is up to
#' the user.
#' 
#' Optionally the table can include test statistics.  Normally the F statistic
#' is most appropriate, which compares the mean square for a row to the
#' residual sum of squares for the largest model considered.  If \code{scale}
#' is specified chi-squared tests can be used. Mallows' \eqn{C_p}{Cp} statistic
#' is the residual sum of squares plus twice the estimate of
#' \eqn{\sigma^2}{sigma^2} times the residual degrees of freedom.
#' 
#' If any of the objects do not have estimated log-likelihoods, produces an
#' error, unless \code{eval.loglik=TRUE}.
#' 
#' @aliases anova.rpm
#' @param object,... objects of class \code{\link{rpm}}, usually, a result of a
#' call to \code{\link{rpm}}.
#' @return An object of class \code{"anova"} inheriting from class
#' \code{"data.frame"}.
#' @section Warning: The comparison between two or more models will only be
#' valid if they are fitted to the same dataset. This may be a problem if there
#' are missing values.
#' @seealso The model fitting function \code{\link{rpm}}, \code{\link{anova}},
#' \code{\link{logLik.rpm}} for adding the log-likelihood to an existing
#' \code{\link[=rpm.object]{rpm}} object.
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
#' anova(fit)
#' }
#' @export
anova.rpm <- function (object, ...) 
{
  if (length(list(object, ...)) > 1)
    return(anova_rpmlist(object, ...))

  logl <- try(stats::logLik(object), silent=TRUE)
  if(inherits(logl,"try-error"))
    stop("The log-likelihood was not estimated for this fit.")

  n <- object$nobs
  df <- (object$NumBeta+1)
  Rdf <- n - df
  logl.null <- if(is.null(object$loglik.null)) 0 else object$loglik.null

  df <- c(1, df)
  Rdf <- c(n-1, Rdf)
  logl <- c(logl.null, logl)
  pv <- stats::pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)
  table <- data.frame(c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- paste(deparse(stats::formula(object)), collapse = "\n")
  colnames(table) <- c("Df", "Deviance", "Resid. Df", "Resid. Dev", 
                       "Pr(>|Chisq|)")
    rownames(table) <- c("NULL", "Model 1:")
  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1), ": ", variables, sep = "", 
                   collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
