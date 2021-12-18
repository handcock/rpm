#' Compute the population distribution of pairs and singles from a Revealed Preference Matchings Model
#' 
#' \code{\link{rpmpopulationpmf}} computes the probability mass function for a population of the pairs and singles
#' from a Revealed Preference Matchings Model based on arbitary availability distribution and 
#' preferences. It is typically based on the estimate from a \code{rpm()} call.
#' 
#' The function \code{\link{rpm}} is used to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' 
#' The model represents the dyadic utility functions as deterministic linear utility functions of
#' dyadic variables. These utility functions are functions of observed characteristics of the women
#' and men.
#' These functions are entered as terms in the function call
#' to \code{\link{rpm}}. This function simulates from such a model.
#'
#' @param object list; an object of class\code{rpm} that is typically the result of a call to \code{rpm()}.
#' @param N integer; The total population size. This must be set. The number of women and men are derived from the
#' (weighted) data.
#' @param num_women integer; (Optional) The number of women in the population.
#' @param num_men integer; (Optional) The number of men in the population.
#' @param pmfW vector; (Optional) The population proportions of the numbers of women of each type. 
#' This should be compatible with the type in the object.
#' @param pmfM vector; (Optional) The population proportions of the numbers of men of each type. 
#' This should be compatible with the type in the object.
#' @param verbose  logical; Should verbose messages be printed out.
#' @return A list of data.frame, each a simulation from the the population. 
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' a <- rpmpopulationpmf(fit)
#' @references Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @export rpmpopulationpmf
rpmpopulationpmf <- function(object, N = 2000, num_women=NULL, num_men=NULL, pmfW=NULL, pmfM=NULL, verbose = FALSE) 
{
      if(!("rpm" %in% class(object))){stop("'rpmpopulationpmf' only works on objects of class 'rpm'.")}

      if (object$sampling_design != "census") {
       # Presumes individual weights
       N_w = sum(object$Xdata[object$Xdata[,object$sampled] & is.na(object$Xdata[,object$pair_id]),object$X_w])
       N_w = N_w + sum(object$Zdata[object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id]),object$Z_w]) # number of women
       N_m = sum(object$Zdata[object$Zdata[,object$sampled] & is.na(object$Zdata[,object$pair_id]),object$Z_w])
       N_m = N_m + sum(object$Xdata[object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]),object$X_w]) # number of men
      }else{
        N_w = nrow(object$Xdata)
        N_m = nrow(object$Zdata)
        if(is.null(object$sampled)){
          object$sampled <- "sampled"
          object$Xdata$sampled <- rep(TRUE,nrow(object$Xdata))
          object$Zdata$sampled <- rep(TRUE,nrow(object$Zdata))
        }else{
          object$Xdata[[object$sampled]] <- rep(TRUE,nrow(object$Xdata))
          object$Zdata[[object$sampled]] <- rep(TRUE,nrow(object$Zdata))
        }
      }
      if(missing(N)){
       N = N_w + N_m
       gw = log(N_w/N) # to ensure exp(gw)+exp(gm) = 1
       gm = log(N_m/N) # to ensure exp(gw)+exp(gm) = 1
      }

      if(!missing(num_women) & !missing(num_men)){
        N <- num_men + num_women
        gw = log(num_women/N)
        gm = log(num_men/N)
      }else{
        if(missing(num_women)){
          if(!missing(num_men)){
           if(num_men >= N){stop("The number of people in the population should be at least as large as the number of men.")}
            num_women = N-num_men
            gw = log(num_women/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(num_men/N) # to ensure exp(gw)+exp(gm) = 1
          }else{
           # missing num_women, num_men
           if(missing(N)){
            num_women <- N_w
            num_men <- N_m
            N <- N_w + N_m
           }else{
            num_women <- N*N_w/(N_w+N_m)
            num_men   <- N*N_m/(N_w+N_m)
            gw = log(num_women/N)
            gm = log(num_men/N)
           }
          }
        }else{ # !missing(num_women) & missing(num_men)
          if(num_women >= N){stop("The number of people in the population should be at least as large as the number of women.")}
  	  num_men = N-num_women
          gw = log(num_women/N) # to ensure exp(gw)+exp(gm) = 2
          gm = log(num_men/N) # to ensure exp(gw)+exp(gm) = 2
        }
      }

      if(!missing(pmfW)){pmfWN <- pmfW}else{pmfWN <- object$pmfW}
      if(!missing(pmfM)){pmfMN <- pmfM}else{pmfMN <- object$pmfM}

      num_women = round(N*exp(gw))
      num_men   = round(N*exp(gm))
      numX <- nrow(object$pmf)-1
      numZ <- ncol(object$pmf)-1
      NumBeta <- object$NumBeta
      NumGamma <- object$NumGamma
      NumGammaW <- object$NumGammaW
      NumGammaM <- object$NumGammaM
      th_hat <- object$coefficients

      hat_gw <- object$coefficients[object$NumBeta+object$NumGammaW+object$NumGammaM+1]
      hat_gm <- log(1-exp(hat_gw))
      pmf_target <- exp(augpmfnew(object$coefficients[1:object$NumBeta],
            GammaW=object$coefficients[object$NumBeta+(1:object$NumGammaW)],
            GammaM=object$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
            object$Sd, object$Xd, object$Zd,
            pmfWN/sum(pmfWN), pmfMN/sum(pmfMN), gw=hat_gw, gm=hat_gm))
      pmf_target[nrow(pmf_target),ncol(pmf_target)] <- 0
      if (object$sampling_design %in% c("census", "stock-stock")) { 
        pmf_target[-nrow(pmf_target), -ncol(pmf_target)] <- 2*pmf_target[-nrow(pmf_target), -ncol(pmf_target)]
        pmf_target <- pmf_target/sum(pmf_target)
      }
      out <- list(pmf=pmf_target)

      out$PMF_SW <- pmf_target[-nrow(pmf_target),ncol(pmf_target)] / (apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],1,sum)*0.5+pmf_target[-nrow(pmf_target),ncol(pmf_target)])
      out$PMF_SM <- pmf_target[nrow(pmf_target),-ncol(pmf_target)] / (apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],2,sum)*0.5+pmf_target[nrow(pmf_target),-ncol(pmf_target)])
      out$PMF_PW <- apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],1,sum)*0.5 / (apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],1,sum)*0.5+pmf_target[-nrow(pmf_target),ncol(pmf_target)])
      out$PMF_PM <- apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],2,sum)*0.5 / (apply(pmf_target[-nrow(pmf_target),-ncol(pmf_target)],2,sum)*0.5+pmf_target[nrow(pmf_target),-ncol(pmf_target)])
      names(out$PMF_SW) <- names(object$PMF_SW)
      names(out$PMF_SM) <- names(object$PMF_SM)
      names(out$PMF_PW) <- names(object$PMF_PW)
      names(out$PMF_PM) <- names(object$PMF_PM)
      if(object$control$logodds_single){
        b <- sum(out$PMF_SW*pmfWN)
        out$LOGODDS_SW <- th_hat[(NumBeta+1):(NumBeta+NumGammaW)] - log(b/(1-b))
        b <- sum(out$PMF_SM*pmfMN)
        out$LOGODDS_SM <- th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)] - log(b/(1-b))
      }else{
        out$LOGODDS_SW <- log(out$PMF_SW/(1-out$PMF_SW))
        out$LOGODDS_SM <- log(out$PMF_SM/(1-out$PMF_SM))
      }
      names(out$LOGODDS_SW) <- names(object$LOGODDS_SW)
      names(out$LOGODDS_SM) <- names(object$LOGODDS_SM)
      return(out)
}
