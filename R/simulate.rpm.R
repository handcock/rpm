#' Simulate a sample of pairs and singles from a Revealed Preference Matchings Model
#' 
#' \code{\link{simulate.rpm}} simulates a population of the pairs and singles
#' from a Revealed Preference Matchings Model. It is typically based on the estimate from a \code{rpm()} call.
#' 
#' The function \code{\link{rpm}} is used to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using a large-population likelihood based on ideas from Menzel (2015).
#' 
#' The model represents the dyadic utility functions as deterministic linear utility functions of
#' dyadic variables. These utility functions are functions of observed characteristics of the women
#' and men.
#' These functions are entered as terms in the function call
#' to \code{\link{rpm}}. This function simulates from such a model.
#'
#' @param object list; an object of class\code{rpm} that is typically the result of a call to \code{rpm()}.
#' @param nsim Number of matchings to be randomly drawn from the given
#' model on the set of all matchings / singles.
#' @param seed integer; (Optional) random number seed.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @param N integer; The total population size. This must be set. The number of women and men are derived from the
#' (weighted) data.
#' @param num_women integer; (Optional) The number of women in the population.
#' @param num_men integer; (Optional) The number of men in the population.
#' @param pmfW vector; (Optional) The population proportions of the numbers of women of each type. 
#' This should be compatible with the type in the object.
#' @param pmfM vector; (Optional) The population proportions of the numbers of men of each type. 
#' This should be compatible with the type in the object.
#' @param large.population logical; If TRUE a large population approximation is used to generate
#' the matchings (rather than the individual level generation of utilities). This is much faster and
#' uses a lot less memory. It is TRUE by default. If used, a sample is drawn rather than the
#' population being returned. The sample size is controlled by \code{num_sampled}.
#' @param num_sampled integer; The size of the sample to be drawn. For "stock-stock" sampling this is the 
#' number of sampled households. For "stock-flow" it is the number of sampled people.
#' For "census" it is the total population size, N. If NULL the size is the same as
#' the passed fitted object (that is, the original data), although this is only a guess and it should be explicitly set.
#' @param bootstrap logical; If TRUE the original population is sampled from. If FALSE the
#' population underlying the fitted model is sampled from.
#' @param sampling_design string; The name of the sampling protocol used to select the survey data. Valid values are
#' \code{"stock-flow"} (individuals are sampled, data contains both
#' singles and couples);
#' \code{"stock-stock"} (households are sampled, each household can be a single or a couple);
#' \code{"census"} (the sample is a census of the population of people).
#' The final option, the default, is NULL whereby the design is taken from the passed object.
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @param verbose  logical; Should verbose messages be printed out.
#' @return A list of data.frame, each a simulation from the the population. 
#' @keywords models
#' @examples
#' library(rpm)
#' \donttest{
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' a <- simulate(fit)
#' }
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#' Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @name simulate.rpm
#' @importFrom stats simulate
#' @export
simulate.rpm <- function(object, nsim=1, seed = NULL, ..., N = NULL, num_women=NULL, num_men=NULL, pmfW=NULL, pmfM=NULL, 
                         large.population=TRUE, num_sampled = NULL, bootstrap=FALSE, sampling_design=NULL,
                         control=control.rpm(), verbose = FALSE) 
{
      if(is.null(object$Xdata) | is.null(object$Zdata)){
        stop('The passed rpm fitted object has a null Xdata or Zdata component. Did you specify "save.data=FALSE" in the fit?')
      }

      if(!is.null(seed)) set.seed(seed)
      if(!is.null(sampling_design)){
        sampling_design <- match.arg(sampling_design, c("stock-stock","stock-flow","census"))
      }else{
        sampling_design <- object$sampling_design
      }

      if(!is.null(N) & !is.null(num_women) & !is.null(num_men)){
       if(abs(N - (num_women + num_men)) > 1){
        stop('You have passed N, num_women and num_men and N != (num_women + num_men). These values are inconsistent.')
       }
      }

      rescale <- FALSE
      # The number of people in the population
      if(object$sampling_design == "stock-stock"){
        N_w = sum(object$Xdata[object$Xdata[,object$sampled] & is.na(object$Xdata[,object$pair_id]),object$X_w])
        N_w = N_w + sum(object$Zdata[object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id]),object$Z_w]) # The population size
        N_m = sum(object$Zdata[object$Zdata[,object$sampled] & is.na(object$Zdata[,object$pair_id]),object$Z_w])
        N_m = N_m + sum(object$Xdata[object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]),object$X_w]) # The population size
      }
      if(object$sampling_design == "stock-flow"){
        N_w = sum(object$Xdata[object$Xdata[,object$sampled],object$X_w])
        N_w = N_w + sum(object$Zdata[object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id]),object$Z_w])
        N_m = sum(object$Zdata[object$Zdata[,object$sampled],object$Z_w])
        N_m = N_m + sum(object$Xdata[object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]),object$X_w])
      }
      if(object$sampling_design == "census"){
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
        if(is.null(object$X_w)){
          object$X_w <- "X_w"
          object$Xdata$X_w <- rep(1,nrow(object$Xdata))
        }else{
          object$Xdata[[object$X_w]] <- rep(1,nrow(object$Xdata))
        }
        if(is.null(object$Z_w)){
          object$Z_w <- "Z_w"
          object$Zdata$Z_w <- rep(1,nrow(object$Zdata))
        }else{
          object$Zdata[[object$Z_w]] <- rep(1,nrow(object$Zdata))
        }
      }
      if(is.null(N)){
       Nt = N_w + N_m
       gw = log(N_w/Nt) # to ensure exp(gw)+exp(gm) = 1
       gm = log(N_m/Nt) # to ensure exp(gw)+exp(gm) = 1
      }else{
       Nt <- N
      }

      if(!is.null(num_women) & !is.null(num_men)){
        rescale <- TRUE
        N <- num_men + num_women
        gw = log(num_women/N)
        gm = log(num_men/N)
      }else{
        if(is.null(num_women)){
          if(!is.null(num_men)){
           rescale <- TRUE
           if(num_men >= Nt){stop("The number of people in the population should be at least as large as the number of men.")}
            num_women = Nt-num_men
          }else{
           # missing num_women, num_men
           if(is.null(N)){
            num_women <- N_w
            num_men <- N_m
           }else{
            rescale <- TRUE
            num_women <- N*N_w/(N_w+N_m)
            num_men   <- N*N_m/(N_w+N_m)
           }
          }
        }else{ # !missing(num_women) & missing(num_men)
          rescale <- TRUE
          if(num_women >= Nt){stop("The number of people in the population should be at least as large as the number of women.")}
  	  num_men = Nt-num_women
        }
      }

      N = num_women + num_men
      gw = log(num_women/N) # to ensure exp(gw)+exp(gm) = 1
      gm = log(num_men/N)   # to ensure exp(gw)+exp(gm) = 1
      
      if(is.null(num_sampled)){
        num_sampled = object$nobs
      }

      if(num_sampled >  N){
        stop('The sample size requested is greater than the population size.')
      }

      if(!large.population &  N > 20000){
        warning('To use the small population exact method you need a population size of at most 20000. Setting the population size to 5000.')
        N = 5000
        rescale <- TRUE
      }

      if(!is.null(pmfW)){object$pmfW_N <- round(num_women*pmfW); rescale <- TRUE}else{object$pmfW_N <- round(num_women*object$pmfW)}
      if(!is.null(pmfM)){object$pmfM_N <- round(  num_men*pmfM); rescale <- TRUE}else{object$pmfM_N <- round(  num_men*object$pmfM)}

      if(!is.logical(bootstrap) | length(bootstrap) != 1){
        stop('The option "bootstrap" should be logical and of length 1.')
      }

      num_women = round(N*exp(gw))
      num_men   = round(N*exp(gm))
      numX <- nrow(object$pmf)-1
      numZ <- ncol(object$pmf)-1

      if(control$ncores > 1){
        doFuture::registerDoFuture()
        if(Sys.info()[["sysname"]] == "Windows"){
          future::plan(multisession, workers=control$ncores)  ## on MS Windows
        }else{
          future::plan(multisession, workers=control$ncores)     ## on Linux, Solaris, and macOS
        }
        ### initialize parallel random number streams
        if (!is.null(seed) | !is.null(control$seed)) {
          if(!is.null(seed)){
            doRNG::registerDoRNG(seed)
          }else{
            doRNG::registerDoRNG(control$seed)
          }
        }
      }

      if(large.population | bootstrap){
        # Use the large population aggregation
        # Or large.population==FALSE and bootstrap==TRUE
        # bootstrap
        if(rescale){
          if(!bootstrap){
           # Note that this recomputes the Gammas based on beta, gw, gm, pmfW, pmfM
            pmf_target <- exp(augpmfnew(beta=object$coefficients[1:object$NumBeta],
                  GammaW=object$coefficients[object$NumBeta+(1:object$NumGammaW)],
                  GammaM=object$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                  S=object$Sd, X=object$Xd, Z=object$Zd,
                  pmfW=object$pmfW_N/sum(object$pmfW_N), pmfM=object$pmfM_N/sum(object$pmfM_N), gw=gw, gm=gm))
          }else{
            pmf_target <- augpmfWM(
                  pmfW=object$pmfW, pmfM=object$pmfM, pmf=object$pmf, gw=object$gw, gm=object$gm,
                  pmfW_N=object$pmfW_N/sum(object$pmfW_N), pmfM_N=object$pmfM_N/sum(object$pmfM_N), gwN=gw, gmN=gm)
          }
          pmf_target[nrow(pmf_target),ncol(pmf_target)] <- 0
        }else{
          if(!bootstrap){
           #pmf_target <- object$pmf_est
            pmf_target <- exp(augpmfnew(beta=object$coefficients[1:object$NumBeta],
                  GammaW=object$coefficients[object$NumBeta+(1:object$NumGammaW)],
                  GammaM=object$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                  S=object$Sd, X=object$Xd, Z=object$Zd,
                  pmfW=object$pmfW_N/sum(object$pmfW_N), pmfM=object$pmfM_N/sum(object$pmfM_N), gw=gw, gm=gm))
          }else{
            pmf_target <- object$pmf
          }
        }

        # Store the original weights
        object$Xdata$X_w_raw <- object$Xdata[,object$X_w]
        object$Zdata$Z_w_raw <- object$Zdata[,object$Z_w]

        rescale <- !bootstrap | rescale
        X_w_rel <- rep(1,nrow(object$Xdata))
        Z_w_rel <- rep(1,nrow(object$Zdata))
        X_w_new <- object$Xdata[,object$X_w]
        Z_w_new <- object$Zdata[,object$Z_w]
        if(rescale){
          paired_and_sampled_W <- object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]) 
          paired_and_sampled_M <- object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id])
          single_and_sampled_W <- object$Xdata[,object$sampled] & is.na(object$Xdata[,object$pair_id]) 
          single_and_sampled_M <- object$Zdata[,object$sampled] & is.na(object$Zdata[,object$pair_id])
          # Rescale the weights of $pmf to $pmf_est
          for( j in 1:numX){
            a <- paired_and_sampled_W & object$Xdata$Xtype==j
            b <- match(object$Xdata[a,object$pair_id], object$Zdata[,object$Zid])
            if(any(is.na(b))){stop("Non matched W pairs")}
            w_rel <- object$Xdata[a,object$X_w]
            w     <- object$Xdata[a,object$X_w]
            for( k in 1:numZ){
              d <- object$Zdata$Ztype[b] == k
              if(any(d)){
               w[d] <- w[d] * pmf_target[j,k] / object$pmf[j,k]
               w_rel[d]    <- pmf_target[j,k] / object$pmf[j,k]
              }
            }
            X_w_rel[a] <- w_rel
            X_w_new[a] <- w
          }
          for( k in 1:numZ){
            a <- paired_and_sampled_M & object$Zdata$Ztype==k
            b <- match(object$Zdata[a,object$pair_id], object$Xdata[,object$Xid])
            w_rel <- object$Zdata[a,object$Z_w]
            w     <- object$Zdata[a,object$Z_w]
            for( j in 1:numX){
              d <- object$Xdata$Xtype[b] == j
              if(any(d)){
               w_rel[d]    <- pmf_target[j,k] / object$pmf[j,k]
               w[d] <- w[d] * pmf_target[j,k] / object$pmf[j,k]
              }
            }
            Z_w_rel[a] <- w_rel
            Z_w_new[a] <- w
          }
          # Now do singles
          w     <- object$Xdata[single_and_sampled_W,object$X_w]
          w_rel <- object$Xdata[single_and_sampled_W,object$X_w]
          for( j in 1:numX){
            a <- object$Xdata$Xtype[single_and_sampled_W]==j
            if(any(a)){
             w[a] <- w[a] * pmf_target[j,ncol(object$pmf)] / object$pmf[j,ncol(object$pmf)]
             w_rel[a]    <- pmf_target[j,ncol(object$pmf)] / object$pmf[j,ncol(object$pmf)]
            }
          }
          X_w_rel[single_and_sampled_W] <- w_rel
          X_w_new[single_and_sampled_W] <- w
#
          w     <- object$Zdata[single_and_sampled_M,object$Z_w]
          w_rel <- object$Zdata[single_and_sampled_M,object$Z_w]
          for( k in 1:numZ){
            a <- object$Zdata$Ztype[single_and_sampled_M]==k
            if(any(a)){
             w[a] <- w[a] * pmf_target[nrow(object$pmf),k] / object$pmf[nrow(object$pmf),k]
             w_rel[a]    <- pmf_target[nrow(object$pmf),k] / object$pmf[nrow(object$pmf),k]
            }
          }
          Z_w_rel[single_and_sampled_M] <- w_rel
          Z_w_new[single_and_sampled_M] <- w
#
         object$Xdata[,object$X_w] <- X_w_new
         object$Zdata[,object$Z_w] <- Z_w_new
#        Do the actual sampling
        }
# Old approach
        if(object$sampling_design %in% c("stock-stock", "census")){
          paired_W <- !is.na(object$Xdata[,object$pair_id]) 
          paired_M <- !is.na(object$Zdata[,object$pair_id])
          M_paired_to_sampled_W <- match(object$Xdata[paired_W,object$pair_id], object$Zdata[,object$Zid])
          I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_M), replace=TRUE)
          a <- X_w_rel
          a[paired_W] <- 0.0
          a[paired_W][ I] <- X_w_rel[paired_W][ I]
          X_w_rel <- a
          object$Xdata[paired_W,object$sampled] <- FALSE
          object$Xdata[paired_W,object$sampled][ I] <- TRUE
          a <- Z_w_rel
          a[paired_M] <- 0.0
          a[M_paired_to_sampled_W][!I] <- Z_w_rel[M_paired_to_sampled_W][!I]
          Z_w_rel <- a
          object$Zdata[paired_M,object$sampled] <- FALSE
          object$Zdata[M_paired_to_sampled_W,object$sampled][!I] <- TRUE
        }
  object$Xdata[,"X_w_rel"] <- X_w_new
  object$Zdata[,"Z_w_rel"] <- Z_w_new

        rpm.simulate.large.population.worker <- function(pmf_target, object, X_w_rel, Z_w_rel, num_sampled, num_women, num_men,numX,numZ){
          cts <- -pmf_target
          ntries <- 0
          while( any(cts < -0.000000001) & ntries < 1000){
           for(i in 1:numX){ 
             if(object$pmfW_N[i] > 0 & any(pmf_target[i,] > 0)){
               a <- stats::rmultinom(n=1,size=object$pmfW_N[i],prob=c(sum(pmf_target[i,1:numZ]),pmf_target[i,(numZ+1)]))
               cts[i,1:numZ] <- stats::rmultinom(n=1,size=a[1],prob=pmf_target[i,1:numZ])
               cts[i,numZ+1] <- a[2]
             }else{
               cts[i,] <- 0
             }
           }
           cts[numX+1,1:numZ] <- object$pmfM_N-apply(cts[-(numX+1),-(numZ+1),drop=FALSE],2,sum)
           ntries <- ntries + 1
          }
#         M-H
          cts.orig <- cts
          num_MH <- 1000*numX*numZ
          Z <- matrix(0,nrow=numX+1,ncol=numZ+1)
          pX1 <- sample(1:numX,size=num_MH,replace=TRUE)
          pX2 <- sample(1:numX,size=num_MH,replace=TRUE)
          pZ1 <- sample(1:numZ,size=num_MH,replace=TRUE)
          pZ2 <- sample(1:numZ,size=num_MH,replace=TRUE)
          alpha <- pmf_target[pX2+(pZ1-1)*(numX+1)]*pmf_target[pX1+(pZ2-1)*(numX+1)]/(pmf_target[pX1+(pZ1-1)*(numX+1)]*pmf_target[pX2+(pZ2-1)*(numX+1)])
          rtest0 <- stats::runif(num_MH)
          rtest1 <- stats::runif(num_MH)
          for(i in 1:num_MH){ 
            if(rtest0[i] > 0.5){
             if(cts[pX1[i],pZ1[i]] > 0 & cts[pX2[i],pZ2[i]] > 0){
              alphai <- alpha[i]*cts[pX1[i],pZ1[i]]*cts[pX2[i],pZ2[i]]/((cts[pX2[i],pZ1[i]]+1)*(cts[pX1[i],pZ2[i]] + 1))
              accept <- (alphai > 1) | (alphai > rtest1[i])
              if(accept){
               cts[pX1[i],pZ1[i]] <- cts[pX1[i],pZ1[i]] - 1 
               cts[pX2[i],pZ1[i]] <- cts[pX2[i],pZ1[i]] + 1 
               cts[pX2[i],pZ2[i]] <- cts[pX2[i],pZ2[i]] - 1 
               cts[pX1[i],pZ2[i]] <- cts[pX1[i],pZ2[i]] + 1 
              }
             }
            }else{
             if(cts[pX2[i],pZ1[i]] > 0 & cts[pX1[i],pZ2[i]] > 0){
              alphai <- (1/alpha[i])*cts[pX2[i],pZ1[i]]*cts[pX1[i],pZ2[i]]/((cts[pX1[i],pZ1[i]]+1)*(cts[pX2[i],pZ2[i]] + 1))
              accept <- (alphai > 1) | (alphai > rtest1[i])
              if(accept){
               cts[pX1[i],pZ1[i]] <- cts[pX1[i],pZ1[i]] + 1 
               cts[pX2[i],pZ1[i]] <- cts[pX2[i],pZ1[i]] - 1 
               cts[pX2[i],pZ2[i]] <- cts[pX2[i],pZ2[i]] + 1 
               cts[pX1[i],pZ2[i]] <- cts[pX1[i],pZ2[i]] - 1 
              }
             }
            }
          }

          num_jk <- cts
#         So far we have a census, like microsimulate. Now sub-sample
#         to get a sample
          if(num_sampled < sum(cts) ){
            a <- sample.int(sum(cts),size=num_sampled)
            rs <- factor(rep(row(cts),cts)[a],1:(numX+1))
            cs <- factor(rep(col(cts),cts)[a],1:(numZ+1))
            num_jk = stats::xtabs(~ rs + cs)
          }

          I <- NULL
          Is <- NULL
          for( j in 1:numX){
            for( k in 1:numZ){
              selX <- object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]) & object$Xdata$Xtype==j
              M_paired_to_sampled_W <- match(object$Xdata[selX,object$pair_id], object$Zdata[,object$Zid])
              selX[selX] <- object$Zdata$Ztype[M_paired_to_sampled_W] == k 
              selZ <- object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id]) & object$Zdata$Ztype==k
              W_paired_to_sampled_M <- match(object$Zdata[selZ,object$pair_id], object$Xdata[,object$Xid])
              selZ[selZ] <- object$Xdata$Xtype[W_paired_to_sampled_M] == j 
              if(num_jk[j,k]> 0 & sum(selX)+sum(selZ) > 0){
                I <- c(I, sample(
                        c(object$Xdata[selX,object$Xid], object$Zdata[selZ,object$Zid]),
                        replace=TRUE,size=num_jk[j,k],
                        prob=c(X_w_rel[selX], Z_w_rel[selZ])
                                )
                      )
              }
            }
#           resample Xdata singles with Xtype==j
            a <- object$Xdata[object$Xdata[,object$sampled] & object$Xdata$Xtype==j & is.na(object$Xdata[,object$pair_id]),object$Xid]
            if(length(a)>0){
              Is <- c(Is, sample(a,
                      replace=TRUE,size=num_jk[j,numZ+1],
                      prob=X_w_rel[object$Xdata[,object$sampled] & object$Xdata$Xtype==j & is.na(object$Xdata[,object$pair_id])]
                    ) )
            }
          }
#         resample Zdata singles with Ztype==k
          for( k in 1:numZ){
            a <- object$Zdata[object$Zdata[,object$sampled] & object$Zdata$Ztype==k & is.na(object$Zdata[,object$pair_id]),object$Zid]
            if(length(a)>0){
              Is <- c(Is, sample(a,
                      replace=TRUE,size=num_jk[numX+1,k],
                      prob=Z_w_rel[object$Zdata[,object$sampled] & object$Zdata$Ztype==k & is.na(object$Zdata[,object$pair_id])]
                    ) )
            }
          }
          Xmatch <- match(Is,object$Xdata[,object$Xid])
          Zmatch <- match(Is,object$Zdata[,object$Zid])
          XdataSingle <- object$Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataSingle <- object$Zdata[Zmatch[!is.na(Zmatch)],]
#        As the sampling is with replacement, create new people for the
#        duplicated cases by creating unique IDs for them (and their partners).
          Idups <- Is[duplicated(Is)]
          Iunique <- unique(Idups)
          for( l in Iunique ){
           a <- XdataSingle[,object$Xid] == l
           if(any(a)){
            # fix pid
            XXid <- XdataSingle[a,object$Xid]
            XdataSingle[a,object$Xid] = paste(XXid,seq_along(XXid),sep="_W_")
           }
           a <- ZdataSingle[,object$Zid] == l
           if(any(a)){
            ZZid <- ZdataSingle[a,object$Zid]
            ZdataSingle[a,object$Zid] = paste(ZZid,seq_along(ZZid),sep="_M_")
           }
          }
#         sampled and partnered
          Xmatch <- match(I,object$Xdata[,object$Xid])
          Zmatch <- match(I,object$Zdata[,object$Zid])
          XdataSampledandPartnered <- object$Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataSampledandPartnered <- object$Zdata[Zmatch[!is.na(Zmatch)],]
#         Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataSampledandPartnered[,object$pair_id])
          M_paired_to_sampled_W <- match(XdataSampledandPartnered[paired_and_sampled_W,object$pair_id], object$Zdata[,object$Zid])
          ZdataM_paired_to_sampled_W <- object$Zdata[M_paired_to_sampled_W,]
          paired_and_sampled_M <- !is.na(ZdataSampledandPartnered[,object$pair_id])
          W_paired_to_sampled_M <- match(ZdataSampledandPartnered[paired_and_sampled_M,object$pair_id], object$Xdata[,object$Xid])
          XdataW_paired_to_sampled_M <- object$Xdata[W_paired_to_sampled_M,]
#         Assigned sampled
          if(nrow(XdataSingle)>0) XdataSingle[,object$sampled] <- TRUE
          if(nrow(ZdataSingle)>0) ZdataSingle[,object$sampled] <- TRUE
          if(nrow(XdataSampledandPartnered)>0) XdataSampledandPartnered[,object$sampled] <- TRUE
          if(nrow(ZdataSampledandPartnered)>0) ZdataSampledandPartnered[,object$sampled] <- TRUE

#        As the sampling is with replacement, create new people for the
#        duplicated cases by creating unique IDs for them (and their partners).
          Idups <- I[duplicated(I)]
          Iunique <- unique(Idups)
          for( l in Iunique ){
           a <- XdataSampledandPartnered[,object$Xid] == l
           if(any(a)){
            # fix pid
            XXid <- XdataSampledandPartnered[a,object$Xid]
            XdataSampledandPartnered[a,object$Xid] = paste(XXid,seq_along(XXid),sep="_W_")
            b <- a & !is.na(XdataSampledandPartnered[,object$pair_id])
            if(any(b)){
             # create more pair_ids
             XXid <- XdataSampledandPartnered[b,object$pair_id]
             XdataSampledandPartnered[b,object$pair_id] = paste(XXid,seq_along(XXid),sep="_M_")
              for( m in unique(XXid) ){
               d <- ZdataM_paired_to_sampled_W[,object$Zid] == m
               if(any(d)){
                ZZid <- ZdataM_paired_to_sampled_W[d,object$Zid]
                ZdataM_paired_to_sampled_W[d,object$Zid] = paste(ZZid,seq_along(ZZid),sep="_M_")
                ZZid <- ZdataM_paired_to_sampled_W[d,object$pair_id]
                ZdataM_paired_to_sampled_W[d,object$pair_id] = paste(ZZid,seq_along(ZZid),sep="_W_")
               }
              }
            }
           }
           a <- ZdataSampledandPartnered[,object$Zid] == l
           if(any(a)){
            ZZid <- ZdataSampledandPartnered[a,object$Zid]
            ZdataSampledandPartnered[a,object$Zid] = paste(ZZid,seq_along(ZZid),sep="_M_")
            b <- a & !is.na(ZdataSampledandPartnered[,object$pair_id])
            if(any(b)){
             ZZid <- ZdataSampledandPartnered[b,object$pair_id]
             ZdataSampledandPartnered[b,object$pair_id] = paste(ZZid,seq_along(ZZid),sep="_W_")
              for( m in ZZid ){
               d <- XdataW_paired_to_sampled_M[,object$Xid] == m
               if(any(d)){
                XXid <- XdataW_paired_to_sampled_M[d,object$Xid]
                XdataW_paired_to_sampled_M[d,object$Xid] = paste(XXid,seq_along(XXid),sep="_W_")
                XXid <- XdataW_paired_to_sampled_M[d,object$pair_id]
                XdataW_paired_to_sampled_M[d,object$pair_id] = paste(XXid,seq_along(XXid),sep="_M_")
               }
              }
            }
           }
          }
           
#         Merge the de-duplicated sampled and paired-to-sampled
          XdataS <- rbind(XdataSingle, XdataSampledandPartnered)
          ZdataS <- rbind(ZdataSingle, ZdataSampledandPartnered)
          XdataS[,object$X_w] <- XdataS[,"X_w_rel"]
          ZdataS[,object$Z_w] <- ZdataS[,"Z_w_rel"]

          if(sampling_design == "stock-flow"){
            XdataS[!is.na(XdataS[,object$pair_id]),object$X_w] <- 2*XdataS[!is.na(XdataS[,object$pair_id]),object$X_w]
            ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w] <- 2*ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w]
            ZdataS[,object$Z_w] <- 20*ZdataS[,object$Z_w]
            if(nrow(XdataW_paired_to_sampled_M)>0) XdataW_paired_to_sampled_M[,object$sampled] <- FALSE
            if(nrow(ZdataM_paired_to_sampled_W)>0) ZdataM_paired_to_sampled_W[,object$sampled] <- FALSE
          } else { #if (sampling_design == "stock-stock"){
            XdataS[!is.na(XdataS[,object$pair_id]),object$X_w] <- 2*XdataS[!is.na(XdataS[,object$pair_id]),object$X_w]
            ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w] <- 2*ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w]
            if(nrow(XdataW_paired_to_sampled_M)>0) XdataW_paired_to_sampled_M[,object$sampled] <- TRUE
            if(nrow(ZdataM_paired_to_sampled_W)>0) ZdataM_paired_to_sampled_W[,object$sampled] <- TRUE
          }

          w <- XdataS[,object$X_w]
          for( j in 1:numX){
            a <- XdataS$Xtype==j
            if(any(a)){
             w[a] <- w[a] * object$pmfW_N[j] / sum(w[a])
            }
          }
          XdataS[,object$X_w] <- w

          w <- ZdataS[,object$Z_w]
          for( k in 1:numZ){
            a <- ZdataS$Ztype==k
            if(any(a)){
             w[a] <- w[a] * object$pmfM_N[k] / sum(w[a])
            }
          }
          ZdataS[,object$Z_w] <- w

          Xdata <- rbind(XdataS, XdataW_paired_to_sampled_M)
          Zdata <- rbind(ZdataS, ZdataM_paired_to_sampled_W)
    
          Xdata[Xdata[,object$sampled],object$X_w] <- num_women * Xdata[Xdata[,object$sampled],object$X_w] / sum(Xdata[Xdata[,object$sampled],object$X_w])
          Zdata[Zdata[,object$sampled],object$Z_w] <-   num_men * Zdata[Zdata[,object$sampled],object$Z_w] / sum(Zdata[Zdata[,object$sampled],object$Z_w])
          Xdata$X_w_raw <- NULL
          Zdata$Z_w_raw <- NULL
#
          list(Xdata=Xdata[sample.int(nrow(Xdata)),],Zdata=Zdata[sample.int(nrow(Zdata)),])
        }

        if(control$ncores > 1){
          if(verbose) message(sprintf("Starting large population parallel simulation using %d cores.",control$ncores))
          out.list <-
            foreach::foreach (i=1:nsim, .packages=c('rpm')
            ) %dorng% {
            rpm.simulate.large.population.worker(pmf_target, object, X_w_rel, Z_w_rel, num_sampled, num_women, num_men,numX,numZ)
            }
        }else{
          out.list <- vector(nsim, mode="list")
          for( i in 1:nsim ){
            out.list[[i]] <- rpm.simulate.large.population.worker(pmf_target, object, X_w_rel, Z_w_rel, num_sampled, num_women, num_men,numX,numZ)
          }
        }
  
      }else{
        # Use the direct (small population) method
        # These are the categories of women
        Ws <- rep(seq_along(object$pmfW_N),round(object$pmfW_N))
        if(length(Ws) > num_women){
         Ws <- Ws[-sample.int(length(Ws),size=length(Ws)-num_women,replace=FALSE)]
        }
        if(length(Ws) < num_women){
         Ws <- c(Ws,sample.int(length(object$pmfW_N),size=num_women-length(Ws),prob=object$pmfW_N,replace=TRUE))
        }
        Ws <- Ws[sample.int(length(Ws))]
        Ms <- rep(seq_along(object$pmfM_N),round(object$pmfM_N))
        if(length(Ms) > num_men){
         Ms <- Ms[-sample.int(length(Ms),size=length(Ms)-num_men,replace=FALSE)]
        }
        if(length(Ms) < num_men){
         Ms <- c(Ms,sample.int(length(object$pmfM_N),size=num_men-length(Ms),prob=object$pmfM_N,replace=TRUE))
        }
        Ms <- Ms[sample.int(length(Ms))]
        #
        NumBetaS <- dim(object$Sd)[3]
        beta_S <- object$coefficients[1:NumBetaS]

        if(object$NumBetaW>0){
          beta_w <- object$coefficients[NumBetaS+(1:object$NumBetaW)]
        }else{
          beta_w <- NULL
        }
        if(object$NumBetaM>0){
          beta_m <- object$coefficients[NumBetaS+object$NumBetaW + (1:object$NumBetaM)]
        }else{
          beta_m <- NULL
        }

        # create utility matrices
        U_star = matrix(0, nrow=length(Ws), ncol = length(Ms))
        V_star = matrix(0, nrow=length(Ms), ncol = length(Ws))
        for (ii in 1:dim(object$Sd)[3]) {
          U_star = U_star + object$Sd[Ws,Ms,ii] * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_w)){
         for (ii in 1:dim(object$Xd)[3]) {
          U_star = U_star + object$Xd[Ws,Ms,ii] * beta_w[ii]
         }
        }
        for (ii in 1:dim(object$Sd)[3]) {
          V_star = V_star + t(object$Sd[Ws,Ms,ii]) * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_m)){
         for (ii in 1:dim(object$Zd)[3]) {
          V_star = V_star + object$Zd[Ms,Ws,ii] * beta_m[ii]
         }
        }
  
        # adjust for outside option (remain single)
        Jw <- sqrt(num_men)
        Jm <- sqrt(num_women)

        rpm.simulate.small.population.worker <- function(i, num_women, num_men, Jw, Jm,
            U_star, V_star, Ws, Ms, Xid, Zid, Xu, Zu,
            X_w, Z_w, pair_id, sampling_design){
          eta  <- -log(-log(matrix(stats::runif(num_women * num_men), nrow=num_women)))
          zeta <- -log(-log(matrix(stats::runif(num_women * num_men), nrow=num_men)))
  
          eta0  <- -log(-log(stats::runif(num_women)))+log(Jw)
          zeta0 <- -log(-log(stats::runif(num_men)))+log(Jm)
    
          ############## temp to match one-to-one#################
          U <- cbind(eta0,  U_star + eta)
          V <- cbind(zeta0, V_star + zeta)
  
          ########################################################
  
          # generate the matching (W-optimal)
          # uses Menzel's GS which allows remaining single
          ############## temp to match one-to-one#################
          mu = Gale_Shapley(U, V, return.data.frame=TRUE)
          mu=data.frame(mu, Xtype=c(Ws, Ms))
          colnames(mu)[c(1,3)] <- c(Xid, pair_id)
          Xdata <- subset(mu, gender=="F")
          Xdata <- data.frame(Xdata,as.data.frame(Xu[Ws,-1,drop=FALSE]))
          Zdata <- subset(mu, gender=="M")
          colnames(Zdata)[match(Xid,colnames(Zdata))] <- Zid
          colnames(Zdata)[match("Xtype",colnames(Zdata))] <- "Ztype"
          Zdata <- data.frame(Zdata,as.data.frame(Zu[Ms,-1,drop=FALSE]))

          X_w_rel <- rep(1,nrow(Xdata))
          Z_w_rel <- rep(1,nrow(Zdata))
          if(sampling_design %in% c("stock-stock", "census")){
            paired_W <- !is.na(Xdata[,pair_id])
            paired_M <- !is.na(Zdata[,pair_id])
            M_paired_to_sampled_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
            I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_M), replace=TRUE)
            X_w_rel[paired_W] <- 0.0
            X_w_rel[paired_W][ I] <- 1# num_sampled / (nrow(Xdata)+nrow(Zdata))
            Z_w_rel[paired_M] <- 0.0
            Z_w_rel[M_paired_to_sampled_W][!I] <- 1.0
          }
          if(sampling_design %in% c("stock-flow")){
            paired_W <- !is.na(Xdata[,pair_id])
            M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
            I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
            X_w_rel[paired_W] <- 0.0
            X_w_rel[paired_W][ I] <- 2
            Z_w_rel[M_paired_to_W] <- 0.0
            Z_w_rel[M_paired_to_W][!I] <- 2
          }
          if(num_sampled > sum(c(X_w_rel,Z_w_rel) > 0)){
           stop(sprintf("The number of households to be sampled, num_sampled=%d, is greater than the number of households in the population, %d.",
                num_sampled, sum(c(X_w_rel,Z_w_rel) > 0)))
          }
          I <- sample(c(Xdata[,Xid], Zdata[,Zid]),
                        prob=c(X_w_rel,Z_w_rel),
                      replace=FALSE,size=num_sampled)
          Xmatch <- match(I,Xdata[,Xid])
          Xmatch <- Xmatch[!is.na(Xmatch)]
          Zmatch <- match(I,Zdata[,Zid])
          Zmatch <- Zmatch[!is.na(Zmatch)]
          XdataS <- Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataS <- Zdata[Zmatch[!is.na(Zmatch)],]
  #       Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataS[,pair_id])
          M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,pair_id], Zdata[,Zid])
          paired_and_sampled_M <- !is.na(ZdataS[,pair_id])
          W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,pair_id], Xdata[,Xid])
          XdataP <- Xdata[W_paired_to_sampled_M,]
          ZdataP <- Zdata[M_paired_to_sampled_W,]

          num_sampled_pairs <- nrow(XdataP)

          if (sampling_design == "stock-flow") {
            XdataP$sampled <- rep(FALSE,nrow(XdataP))
            ZdataP$sampled <- rep(FALSE,nrow(ZdataP))
           }else{
            # If not "stock-flow" then all are sampled
            XdataP$sampled <- rep(TRUE,nrow(XdataP))
            ZdataP$sampled <- rep(TRUE,nrow(ZdataP))
           }

          XdataS$sampled <- rep(TRUE,nrow(XdataS))
          ZdataS$sampled <- rep(TRUE,nrow(ZdataS))

          Xdata <- rbind(XdataP, XdataS)
          Zdata <- rbind(ZdataP, ZdataS)
  
          a <- rep(0, nrow(Xdata))
          a[Xdata$sampled] <- num_women/sum(Xdata$sampled)
          Xdata[[X_w]] <- a
          a <- rep(0, nrow(Zdata))
          a[Zdata$sampled] <- num_men/sum(Zdata$sampled)
          Zdata[[Z_w]] <- a

          # random permute to add randomness
          list(Xdata=Xdata[sample.int(nrow(Xdata)),],
               Zdata=Zdata[sample.int(nrow(Zdata)),])
        }
        if(control$ncores > 1){
          if(verbose) message(sprintf("Starting small population parallel simulation using %d cores.",control$ncores))
          out.list <-
            foreach::foreach (i=1:nsim, .packages=c('rpm')
            ) %dorng% {
            rpm.simulate.small.population.worker(i, length(Ws), length(Ms), Jw, Jm, U_star, V_star, Ws, Ms, object$Xid, object$Zid, object$Xu, object$Zu, object$X_w, object$Z_w, object$pair_id, sampling_design)
            }
        }else{
          out.list <- vector(nsim, mode="list")
          for( i in 1:nsim ){
            out.list[[i]] <- rpm.simulate.small.population.worker(i, length(Ws), length(Ms), Jw, Jm, U_star, V_star, Ws, Ms, object$Xid, object$Zid, object$Xu, object$Zu, object$X_w, object$Z_w, object$pair_id, sampling_design)
          }
        }
      }
      if(nsim > 1)
	{
		names(out.list) <- 1 : nsim
	}
	else
	{
		out.list <- out.list[[1]]
	}
	return(out.list)
}
