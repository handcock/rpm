#' @rdname anova.rpm
#' @export
anova_rpmlist <- function (object, ...) 
{
  objects <- list(object, ...)
  nmodels <- length(objects)
  if (nmodels == 1) 
    return(anova.rpm(object))
  n <- object$nobs
  logl <- df <- Rdf <- rep(0, nmodels)
  logl.null <- if(is.null(objects[[1]][["loglik.null"]])) 0 else objects[[1]][["loglik.null"]]
  for (i in 1:nmodels) {
    n <- objects[[i]]$nobs
    df[i] <- (objects[[i]]$NumBeta+1) 
    Rdf[i] <- n - df[i]
    logl[i] <- stats::logLik(objects[[i]])
  }
  df <- c(1, df)
  Rdf <- c(n-1, Rdf)
  logl <- c(logl.null, logl)
  pv <- stats::pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)

  table <- data.frame(c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- lapply(objects, function(x) paste(deparse(stats::formula(x)), 
                                                 collapse = "\n"))
  colnames(table) <- c("Df","Deviance", "Resid. Df",
                              "Resid. Dev", "Pr(>|Chisq|)")
  rownames(table) <- c("NULL", 1:nmodels)

  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
