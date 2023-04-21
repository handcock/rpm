#' Fit a Revealed Preference Matchings Model
#' 
#' \code{\link{rpm}} estimates the parameters of a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using an large-population likelihood based on ideas from Dagsvik (2000), Menzel (2015) and Goyal et al (2023).
#' 
#' @aliases rpm.object
#' @param formula formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param Xdata data.frame for women. Each row is a woman, each column is a variable on that women
#' or her partnerships. It must contain the women's ID variable (see \code{Xid}) and
#' a variable with the ID of the women's partner. If the women is single the men's ID should be NA.
#' @param Zdata data.frame for men. Each row is a man, each column is a variable on that men
#' It must contain the men's ID variable (see \code{Zid}).
#' @param Xid string The name of the variable in \code{Xdata} containing the IDs of the women.
#' @param Zid string The name of the variable in \code{Zdata} containing the IDs of the men.
#' @param pair_id string The name of the variable in \code{Xdata} containing the ID of the men paired with the women in
#' \code{Xid}. If the women is not paired it must be NA.
#' @param X_w string The name of the variable in \code{Xdata} containing the individual weight of the women.
#' If this is NULL then it is assumed the sample is unweighted from a population with 2000 women in it.
#' @param Z_w string The name of the variable in \code{Zdata} containing the individual weight of the man
#' If this is NULL then it is assumed the sample is unweighted from a population with 2000 men in it.
#' @param pair_w string The name of the variable in \code{Xdata} containing the pair weight of that women.
#' If the women is not paired it should be NA.
#' If this is NULL then it is computed from the individual weights using the \code{sampling_design}.
#' Note that the pair weights currently do not play a role in the estimation. They do in the quasi-likelihood version of the code.
##' If this is NULL then it is assumed the sample is unweighted from a population with 2000 men in it.
#' @param sampled string The name of the logical variable in \code{Xdata} and \code{Zdata} containing the 
#' indicator that the person was sampled directly (as distinct from being included as the match of a directly sampled
#' person. All single people are directly sampled. 
#' @param sampling_design string; The name of the sampling protocol used to select the survey data. Valid values are
#' \code{"stock-flow"} (default) (individuals are sampled, data contains both
#' singles and couples);
#' \code{"stock-stock"} (households are sampled, each household can be a single or a couple);
#' \code{"census"} (the sample is a census of the population of people).
## (only couples are included in the data).
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including data summary statistics.
#' @details The pairings are determined by the \code{pair_id} variable in \code{Xdata}. 
#' If that variable is NA then the women is
#' assumed to be single. If men are listed in \code{Zdata} and are not partnered then they are assumed single.
#' Weights are specified by three optional variables in \code{Xdata}.
#' \itemize{
#' \item{X_w}{: This is character string of the name of the weight variable for women. The sum of the weights should be the
#' number of women in the population.}
#' \item{Z_w}{: This is character string of the name of the weight variable for men. The sum of the weights should be the
#' number of men in the population.}
#' \item{pair_w}{: This is character string of the name of the weight variable for pairs.}
#' }
#' @return \code{\link{rpm}} returns an object of class \code{\link{rpm.object}}
#' that is a list consisting of the following elements: 
#' \item{coef}{The maximum psuedo-likelihood estimate of \eqn{\theta}, the vector of
#' coefficients for the model parameters. This includes the model \eqn{\beta} and the model \eqn{\Gamma}.}
#' \item{coefficients}{The bias-corrected bootstrap estimate of \eqn{\theta}, the vector of
#' coefficients for the model parameters. This includes the model \eqn{\beta} and the model \eqn{\Gamma}.}
#' \item{loglik}{The value of the maximized log-likelihood.}
#' \item{exitflag}{integer value with the status of the optimization (4 is success as 
#' \code{xtol_rel} or \code{xtol_abs} was reached). Other codes are 1 = generic success; 2 = optimization stopped because 
#' \code{ftol_rel} or \code{ftol_abs} was reached; 3 = optimization stopped 
#' because \code{stopval} was reached; 4 = optimization stopped because \code{xtol_rel} or \code{xtol_abs} was reached; 
#' 5 = optimization stopped because
#' \code{maxeval} was reached; 6 = optimization stopped because \code{maxtime} was reached.}
#' \item{call}{the call that was made to \code{nloptr}.}
#' \item{x0}{vector with starting values for the optimization.}
#' \item{message}{more informative message with the status of the optimization.}
#' \item{iterations}{number of iterations that were executed.}
#' \item{objective}{value if the objective function in the solution.}
#' \item{solution}{optimal value of the controls.}
#' \item{version}{version of NLopt that was used.}
#' \item{covar}{Approximate covariance matrix of the estimates.}
#' \item{eq}{Values from the equality constraints. Larger values indicate non-convergence.}
#' \item{sample}{A matrix with the number of rows the MCMC sample size and the number of rows the number of parameters.}
#' @seealso control.rpm, summary.rpm, print.rpm
#' @references
#'
#' Goyal, Shuchi; Handcock, Mark S.; Jackson, Heide M.; Rendall, Michael S. and Yeung, Fiona C. (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#'
#' Dagsvik, John K. (2000) \emph{Aggregation in Matching Markets} \emph{International Economic Review}, Vol. 41, 27-57.
#' JSTOR: https://www.jstor.org/stable/2648822, \doi{10.1111/1468-2354.00054}
#'
#' Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941. \doi{10.3982/ECTA12299}
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
#' summary(fit)
#' }
#' @importFrom doRNG registerDoRNG
#' @importFrom MASS cov.mcd
#' @importFrom coda mcmc mcmc.list
#' @export rpm
rpm <- function(formula, Xdata, Zdata,
    Xid=NULL, Zid=NULL, pair_id=NULL,
    X_w=NULL, Z_w=NULL, pair_w=NULL,
    sampled=NULL,
    sampling_design="stock-flow",
    control=control.rpm(), verbose=FALSE){

    if(!is.null(control$seed)) set.seed(control$seed)

    if(is.null(Xdata) | !is.data.frame(Xdata)){
      stop("Xdata should be a data frame containing the women's information.")
    }
    if(is.null(Zdata) | !is.data.frame(Zdata)){
      stop("Zdata should be a data frame containing the men's information.")
    }

    if(is.null(sampled) & sampling_design != "census"){
      stop("The variable name in Xdata and Zdata of the (directly) sampled indicator must be specified as the character string 'sampled'.")
    }

    if(!is.null(sampled) & sampling_design=="census"){
      if(is.null(Xdata[[sampled]])){ Xdata[[sampled]] <- rep(TRUE, nrow(Xdata)) }
      if(is.null(Zdata[[sampled]])){ Zdata[[sampled]] <- rep(TRUE, nrow(Zdata)) }
    }

#   if(is.null(sampled) & sampling_design!="census"){
#     if(is.null(Xdata[[sampled]]) | is.null(Zdata[[sampled]])){
#       warning("The data does not have a 'sampled' variable. It will assumed to be a census.")
#       sampling_design="census"
#     }
#   }

    if(missing(Xid)){
      warning("The variable name in Xdata of the women's IDs is missing and is presumed to be 'Xid'.")
      Xid <- "Xid"
      if(is.null(Xdata[["Xid"]])){
        stop("The variable name in Xdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(missing(Zid)){
      warning("The variable name in Zdata of the men's IDs is missing and is presumed to be 'Zid'.")
      Zid <- "Zid"
      if(is.null(Zdata[["Zid"]])){
        stop("The variable name in Zdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(missing(pair_id)){
      warning("The variable name in Xdata of the men's IDs of the women's partners has not been specified. It is presumed to be 'pair_id'.")
      pair_id <- "pair_id"
      if(is.null(Xdata[[pair_id]])){
        stop("The variable 'pair_id' in Xdata of the men's IDs of the women's partners is NULL.")
      }
    }
    if(missing(sampled) & sampling_design != "census"){
      stop("The variable name in Xdata and Zdata of the (directly) sampled indicator must be specified as the character string 'sampled'.")
    }
    if(is.null(X_w) & sampling_design != "census"){
      warning("The variable name in Xdata of the women's weights has not been specified. It is presumed to be 'X_w'.")
      X_w <- "X_w"
      if(is.null(Xdata[["X_w"]])){
    #   This is the unweighted case
        warning("The variable name in Xdata of the women's weights has not been specified. It is specified by the character string 'X_w'. We will presume the sample is unweighted from a population with 2000 women in it.")
        nw <- 2000
        Xdata[["X_w"]]=rep(nw/nrow(Xdata),length=nrow(Xdata))
      }
    }

    if(is.null(Z_w) & sampling_design != "census"){
      warning("The variable name in Zdata of the men's weights has not been specified. It is presumed to be 'Z_w'.")
      Z_w <- "Z_w"
      if(is.null(Zdata[["Z_w"]])){
      # This is the unweighted case
        warning("The variable name in Zdata of the men's weights has not been specified. It is specified by the character string 'Z_w'. We will presume the sample is unweighted from a population with 2000 men in it.")
        nm <- 2000
        Zdata[["Z_w"]]=rep(nm/nrow(Zdata),length=nrow(Zdata))
      }
    }

    if(is.null(X_w)){ X_w <- "X_w" }
    if(is.null(Z_w)){ Z_w <- "Z_w" }
    
    if(sampling_design == "stock-flow" && all(Xdata[,sampled]) && all(Zdata[,sampled])){
      warning("The sampling design is specified as stock-flow, but all men and women are designated as directly sampled.")
    }

    if(sampling_design!="census" && any(is.na(Xdata[Xdata[,sampled],X_w]))){
      warning("Some of the sampled women have NA as a weight. These have been converted to 0 weight. If this is not correct, please specify the correct weights.")
      Xdata[Xdata[,sampled] & is.na(Xdata[,X_w]),X_w] <- 0
    }
    if(sampling_design!="census" && any(is.na(Zdata[Zdata[,sampled],Z_w]))){
      warning("Some of the sampled men have NA as a weight. These have been converted to 0 weight. If this is not correct, please specify the correct weights.")
      Zdata[Zdata[,sampled] & is.na(Zdata[,Z_w]),Z_w] <- 0
    }

    fit <- rpm_MLPLE(formula=formula, Xdata=Xdata, Zdata=Zdata,
                    Xid=Xid, Zid=Zid, pair_id=pair_id,
                    X_w=X_w, Z_w=Z_w, pair_w=pair_w,
                    sampled=sampled,
                    sampling_design=sampling_design,
                    control=control, verbose=verbose)

    fit$aic = 2*fit$NumBeta-2*fit$loglik
    fit$bic = log(fit$nobs)*fit$NumBeta-2*fit$loglik

    return(fit)
}
