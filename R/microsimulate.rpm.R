#' Micro simulate a population from a Revealed Preference Matchings Model
#' 
#' \code{\link{microsimulate}} simulates a population of the pairs and singles
#' from a Revealed Preference Matchings Model. It is typically based on the estimate from a \code{rpm()} call.
#' 
#' The function requites the numbers of women of each type and the number of men of each type to be specified. 
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
#' to \code{\link{rpm}}. This function simulates a population from such a model.
#'
#' @param object list; an object of class\code{rpm} that is typically the result of a call to \code{rpm()}.
#' @param nsim Number of matchings to be randomly drawn from the given
#' model on the set of all matchings / singles.
#' @param seed integer; (Optional) random number seed.
#' @param pmfW_N vector; The population count of the number of women of each type. 
#' This should be compatible with the type in the object.
#' @param pmfM_N vector; The population count of the number of men of each type. 
#' This should be compatible with the type in the object.
#' @param large.population logical; If TRUE a large population approximation is used to generate
#' the matchings (rather than the individual level generation of utilities). This is much faster and
#' uses a lot less memory. It is TRUE by default. If used, a sample is drawn rather than the
#' population being returned. The sample size is controlled by \code{pmfW_N} and \code{pmfM_N}.
#' @param bootstrap logical; If TRUE the original population is sampled from. If FALSE the
#' population underlying the fitted model is sampled from.
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @param counts.only logical; If TRUE only the matrices of counts and the PMF of the population of households is returned.
#' If FALSE 
#' It is FALSE by default. 
#' @param verbose  logical; Should verbose messages be printed out.
#' @return A list of lists, each a simulation from the population. Each of the simulation lists contains
#' components \code{population} being a list with components \code{Xdata} and \code{Zdata} (for use with \code{rpm()}).
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' num_women = fit$N*exp(fit$gw)
#' num_men   = fit$N*exp(fit$gm)
#' pmfW_N <- round(fit$pmfW * num_women)
#' pmfM_N <- round(fit$pmfM * num_men)
#' a <- microsimulate(fit, pmfW_N=pmfW_N, pmfM_N=pmfM_N)
#' @references Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @export
microsimulate <- function(object, nsim=1, seed = NULL, pmfW_N=NULL, pmfM_N=NULL, large.population=TRUE, bootstrap=FALSE, control=control.rpm(), counts.only=FALSE, verbose = FALSE) 
{
      if(!is.null(seed)) set.seed(seed)

      if(missing(pmfW_N)){
        stop('pmfW_N, the population count of the number of women of each type is missing and must be specified.')
      }
      if(missing(pmfM_N)){
        stop('pmfM_N, the population count of the number of men of each type is missing and must be specified.')
      }

      if(any(abs(round(pmfW_N)-pmfW_N)>0.0000001)){
        warning("pmfW_N should all be natural numbers. They have been rounded to ensure this.")
        pmfM_N <- round(pmfM_N)
      }
      if(any(abs(round(pmfM_N)-pmfM_N)>0.0000001)){
        warning("pmfM_N should all be natural numbers. They have been rounded to ensure this.")
        pmfW_N <- round(pmfW_N)
      }

      num_women = round(sum(pmfW_N))
      num_men = round(sum(pmfM_N))
      N = num_women + num_men

      rescale <- TRUE
      pmfW <- pmfW_N / num_women 
      pmfM <- pmfM_N / num_men 

      gw = log(num_women/N) # to ensure exp(gw)+exp(gm) = 1
      gm = log(num_men/N) # to ensure exp(gw)+exp(gm) = 1

      if(!is.logical(bootstrap) | length(bootstrap) != 1){
        stop('The option "bootstrap" should be logical and of length 1.')
      }

      numX <- nrow(object$pmf)-1
      numZ <- ncol(object$pmf)-1

      num_sampled <- sum(pmfW_N) + sum(pmfM_N)

      if(control$ncores > 1){
        doFuture::registerDoFuture()
        cl <- parallel::makeCluster(control$ncores)
        future::plan("cluster", workers = cl)
        if(Sys.info()[["sysname"]] == "Windows"){
          future::plan("multisession")  ## on MS Windows
        }else{
          future::plan("multisession")  ## on Linux, Solaris, and macOS
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
        # bootstrap
        if(rescale){
          if(!bootstrap){
            hat_gw <- object$coefficients[object$NumBeta+object$NumGammaW+object$NumGammaM+1]
            hat_gm <- log(1-exp(hat_gw))
            pmf_target <- exp(augpmfnew(object$coefficients[1:object$NumBeta],
                GammaW=object$coefficients[object$NumBeta+(1:object$NumGammaW)],
                GammaM=object$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                object$Sd, object$Xd, object$Zd,
                pmfW, pmfM, gw=hat_gw, gm=hat_gm))
            pmf_target[nrow(pmf_target),ncol(pmf_target)] <- 0
            pmf_target[-nrow(pmf_target), -ncol(pmf_target)] <- 2*pmf_target[-nrow(pmf_target), -ncol(pmf_target)]
            pmf_target <- pmf_target/sum(pmf_target)
          }else{
            pmf_target <- augpmfWM(
                pmfW=object$pmfW, pmfM=object$pmfM, pmf=object$pmf, gw=object$gw, gm=object$gm,
                pmfW_N=pmfW_N/sum(pmfW_N), pmfM_N=pmfM_N/sum(pmfM_N), gwN=gw, gmN=gm)
            pmf_target[nrow(pmf_target),ncol(pmf_target)] <- 0
          }
        }else{
          if(!bootstrap){
            pmf_target <- object$pmf_est
          }else{
            pmf_target <- object$pmf
          }
        }

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
        # Use the large population aggregation
        # Or large.population==FALSE and bootstrap==TRUE
        # Store the original weights
        object$Xdata$X_w_raw <- object$Xdata[,object$X_w]
        object$Zdata$Z_w_raw <- object$Zdata[,object$Z_w]

        if(object$sampling_design == "census"){
          # For the census, sample pairs only via women only (to avoid duplicates)
          paired_M <- !is.na(object$Zdata[,object$pair_id])
          object$Zdata[paired_M,object$sampled] <- FALSE
        }

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
            w <- object$Xdata[a,object$X_w]
            for( k in 1:numZ){
              d <- object$Zdata$Ztype[b] == k
              if(any(d)){
               w[d] <- w[d] * pmf_target[j,k] / object$pmf[j,k]
               w_rel[d] <- pmf_target[j,k] / object$pmf[j,k]
              }
            }
            X_w_rel[a] <- w_rel
            X_w_new[a] <- w
          }
          for( k in 1:numZ){
            a <- paired_and_sampled_M & object$Zdata$Ztype==k
            b <- match(object$Zdata[a,object$pair_id], object$Xdata[,object$Xid])
            w_rel <- object$Zdata[a,object$Z_w]
            w <- object$Zdata[a,object$Z_w]
            for( j in 1:numX){
              d <- object$Xdata$Xtype[b] == j
              if(any(d)){
               w_rel[d] <- pmf_target[j,k] / object$pmf[j,k]
               w[d] <- w[d] * pmf_target[j,k] / object$pmf[j,k]
              }
            }
            Z_w_rel[a] <- w_rel
            Z_w_new[a] <- w
          }
          # Now do singles
          w <- object$Xdata[single_and_sampled_W,object$X_w]
          w_rel <- object$Xdata[single_and_sampled_W,object$X_w]
          for( j in 1:numX){
            a <- object$Xdata$Xtype[single_and_sampled_W]==j
            if(any(a)){
             w[a] <- w[a] * pmf_target[j,ncol(object$pmf)] / object$pmf[j,ncol(object$pmf)]
             w_rel[a] <- pmf_target[j,ncol(object$pmf)] / object$pmf[j,ncol(object$pmf)]
            }
          }
          X_w_rel[single_and_sampled_W] <- w_rel
          X_w_new[single_and_sampled_W] <- w
#
          w <- object$Zdata[single_and_sampled_M,object$Z_w]
          w_rel <- object$Zdata[single_and_sampled_M,object$Z_w]
          for( k in 1:numZ){
            a <- object$Zdata$Ztype[single_and_sampled_M]==k
            if(any(a)){
             w[a] <- w[a] * pmf_target[nrow(object$pmf),k] / object$pmf[nrow(object$pmf),k]
             w_rel[a] <- pmf_target[nrow(object$pmf),k] / object$pmf[nrow(object$pmf),k]
            }
          }
          Z_w_rel[single_and_sampled_M] <- w_rel
          Z_w_new[single_and_sampled_M] <- w
#
         object$Xdata[,object$X_w] <- X_w_new
         object$Zdata[,object$Z_w] <- Z_w_new
#        Do the actual sampling
        }
        X_w_rel <- rep(1,nrow(object$Xdata))
        Z_w_rel <- rep(1,nrow(object$Zdata))
        if(object$sampling_design == "stock-stock"){
          paired_W <- !is.na(object$Xdata[,object$pair_id]) 
          M_paired_to_W <- match(object$Xdata[paired_W,object$pair_id], object$Zdata[,object$Zid])
          I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
          X_w_rel[paired_W] <- 0.0
          X_w_rel[paired_W][ I] <- 1
          Z_w_rel[M_paired_to_W] <- 0.0
          Z_w_rel[M_paired_to_W][!I] <- 1
          I <- sample(c(object$Xdata[,object$Xid],
                        object$Zdata[,object$Zid]),
                        prob=c(X_w_rel,Z_w_rel),
                      replace=FALSE,size=num_sampled)
        }
        if(object$sampling_design == "census"){
          paired_W <- !is.na(object$Xdata[,object$pair_id]) 
          M_paired_to_W <- match(object$Xdata[paired_W,object$pair_id], object$Zdata[,object$Zid])
          I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
          I <- c(object$Xdata[paired_W,object$Xid][I],object$Zdata[M_paired_to_W,object$Zid][!I])
          I <- c(I,sample(c(object$Xdata[!paired_W,object$Xid], object$Zdata[is.na(object$Zdata[,object$pair_id]),object$Zid]),
               replace=TRUE,size=num_sampled-sum(paired_W)))
        }
        if(object$sampling_design == "stock-flow"){
          X_w_rel[!object$Xdata[,object$sampled]] <- 0.0
          Z_w_rel[!object$Zdata[,object$sampled]] <- 0.0
          I <- sample(c(object$Xdata[,object$Xid],
                        object$Zdata[,object$Zid]),
                        prob=c(X_w_rel,Z_w_rel),
                      replace=FALSE,size=num_sampled)
        }

        object$Xdata$X_w_rel <- object$Xdata[,object$X_w]
        object$Zdata$Z_w_rel <- object$Zdata[,object$Z_w]
        object$Xdata[object$Xdata[,object$sampled],"X_w_rel"] <- X_w_new[object$Xdata[,object$sampled]]
        object$Zdata[object$Zdata[,object$sampled],"Z_w_rel"] <- Z_w_new[object$Zdata[,object$sampled]]

        rpm.simulate.large.population.worker <- function(object, coefficients, N, X_w_rel, Z_w_rel, pmfW_N, pmfM_N, pmfW, pmfM, gm, gw, numX, numZ, bootstrap){

          if(!bootstrap){
            hat_gw <- coefficients[object$NumBeta+object$NumGammaW+object$NumGammaM+1]
            hat_gm <- log(1-exp(hat_gw))
            pmf_target_boot <- exp(augpmfnew(coefficients[1:object$NumBeta],
                GammaW=coefficients[object$NumBeta+(1:object$NumGammaW)],
                GammaM=coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                object$Sd, object$Xd, object$Zd,
                pmfW, pmfM, gw=hat_gw, gm=hat_gm))
            pmf_target_boot[nrow(pmf_target_boot),ncol(pmf_target_boot)] <- 0
          }else{
            pmf_target_boot <- augpmfWM(
                pmfW=object$pmfW, pmfM=object$pmfM, pmf=object$pmf, gw=object$gw, gm=object$gm,
                pmfW_N=pmfW_N/sum(pmfW_N), pmfM_N=pmfM_N/sum(pmfM_N), gwN=gw, gmN=gm)
            pmf_target_boot[nrow(pmf_target_boot),ncol(pmf_target_boot)] <- 0
          }

          cts <- -pmf_target_boot
#
          ntries <- 0
          while( any(cts < -0.000000001) & ntries < 1000){
           for(i in 1:numX){ 
             if(pmfW_N[i] > 0 & any(pmf_target_boot[i,] > 0)){
               a <- stats::rmultinom(n=1,size=pmfW_N[i],prob=c(sum(pmf_target_boot[i,1:numZ]),pmf_target_boot[i,(numZ+1)]))
               cts[i,1:numZ] <- stats::rmultinom(n=1,size=a[1],prob=pmf_target_boot[i,1:numZ])
               cts[i,numZ+1] <- a[2]
             }else{
               cts[i,] <- 0
             }
           }
           cts[numX+1,1:numZ] <- pmfM_N-apply(cts[-(numX+1),-(numZ+1),drop=FALSE],2,sum)
           ntries <- ntries + 1
          }
#         M-H
          cts.orig <- cts
#ptm <- proc.time()
          num_MH <- 1000*numX*numZ
          Z <- matrix(0,nrow=numX+1,ncol=numZ+1)
          pX1 <- sample(1:numX,size=num_MH,replace=TRUE)
          pX2 <- sample(1:numX,size=num_MH,replace=TRUE)
          pZ1 <- sample(1:numZ,size=num_MH,replace=TRUE)
          pZ2 <- sample(1:numZ,size=num_MH,replace=TRUE)
          alpha <- pmf_target_boot[pX2+(pZ1-1)*(numX+1)]*pmf_target_boot[pX1+(pZ2-1)*(numX+1)]/(pmf_target_boot[pX1+(pZ1-1)*(numX+1)]*pmf_target_boot[pX2+(pZ2-1)*(numX+1)])
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

          if(counts.only){
            return( list(pmf_target=pmf_target_boot, counts=cts, cts.orig=cts.orig) )
          }
          num_jk <- cts
#         End count estimation
          I <- NULL
          Is <- NULL
          for( j in 1:numX){
            for( k in 1:numZ){
              selX <- object$Xdata[,object$sampled] & !is.na(object$Xdata[,object$pair_id]) & object$Xdata$Xtype==j
#             select Xtype==j
              M_paired_to_sampled_W <- match(object$Xdata[selX,object$pair_id], object$Zdata[,object$Zid])
#             select X paired Xtype==j & Ztype==k
              selX[selX] <- object$Zdata$Ztype[M_paired_to_sampled_W] == k 
              selZ <- object$Zdata[,object$sampled] & !is.na(object$Zdata[,object$pair_id]) & object$Zdata$Ztype==k
              W_paired_to_sampled_M <- match(object$Zdata[selZ,object$pair_id], object$Xdata[,object$Xid])
#             select Z paired Xtype==j & Ztype==k
              selZ[selZ] <- object$Xdata$Xtype[W_paired_to_sampled_M] == j 
              if(num_jk[j,k]> 0 & sum(selX)+sum(selZ) > 0){
#             resample Xdata and Zdata pairs with Xtype==j & Ztype==k
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
#         partnered
          Xmatch <- match(I,object$Xdata[,object$Xid])
          Zmatch <- match(I,object$Zdata[,object$Zid])
          XdataS <- object$Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataS <- object$Zdata[Zmatch[!is.na(Zmatch)],]
#         Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataS[,object$pair_id])
          M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,object$pair_id], object$Zdata[,object$Zid])
          ZdataW <- object$Zdata[M_paired_to_sampled_W,]
          paired_and_sampled_M <- !is.na(ZdataS[,object$pair_id])
          W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,object$pair_id], object$Xdata[,object$Xid])
          XdataM <- object$Xdata[W_paired_to_sampled_M,]
#        As the sampling is with replacement, create new people for the
#        duplicated cases by creating unique IDs for them (and their partners).
          Idups <- I[duplicated(I)]
          Iunique <- unique(Idups)
          for( l in Iunique ){
           a <- XdataS[,object$Xid] == l
           if(any(a)){
            # fix pid
            XXid <- XdataS[a,object$Xid]
            XdataS[a,object$Xid] = paste(XXid,seq_along(XXid),sep="_W_")
            b <- a & !is.na(XdataS[,object$pair_id])
            if(any(b)){
             # create more pair_ids
             XXid <- XdataS[b,object$pair_id]
             XdataS[b,object$pair_id] = paste(XXid,seq_along(XXid),sep="_M_")
              for( m in unique(XXid) ){
               d <- ZdataW[,object$Zid] == m
               if(any(d)){
                ZZid <- ZdataW[d,object$Zid]
                ZdataW[d,object$Zid] = paste(ZZid,seq_along(ZZid),sep="_M_")
                ZZid <- ZdataW[d,object$pair_id]
                ZdataW[d,object$pair_id] = paste(ZZid,seq_along(ZZid),sep="_W_")
               }
              }
            }
           }
           a <- ZdataS[,object$Zid] == l
           if(any(a)){
            ZZid <- ZdataS[a,object$Zid]
            ZdataS[a,object$Zid] = paste(ZZid,seq_along(ZZid),sep="_M_")
            b <- a & !is.na(ZdataS[,object$pair_id])
            if(any(b)){
             ZZid <- ZdataS[b,object$pair_id]
             ZdataS[b,object$pair_id] = paste(ZZid,seq_along(ZZid),sep="_W_")
              for( m in ZZid ){
               d <- XdataM[,object$Xid] == m
               if(any(d)){
                XXid <- XdataM[d,object$Xid]
                XdataM[d,object$Xid] = paste(XXid,seq_along(XXid),sep="_W_")
                XXid <- XdataM[d,object$pair_id]
                XdataM[d,object$pair_id] = paste(XXid,seq_along(XXid),sep="_M_")
               }
              }
            }
           }
          }
           
#         Merge the de-duplicated sampled and paired-to-sampled
          Xdata <- rbind(XdataSingle, XdataS, XdataM)
          Zdata <- rbind(ZdataSingle, ZdataS, ZdataW)
          Xdata <- Xdata[!duplicated(Xdata[,object$Xid]),]
          Zdata <- Zdata[!duplicated(Zdata[,object$Zid]),]
#
          Xdata[,object$X_w] <- NULL
          Zdata[,object$Z_w] <- NULL
          Xdata[,object$X_w_raw] <- NULL
          Zdata[,object$Z_w_raw] <- NULL
          Xdata[,object$sampled] <- rep(TRUE,nrow(Xdata))
          Zdata[,object$sampled] <- rep(TRUE,nrow(Zdata))

          list(population=list(Xdata=Xdata,Zdata=Zdata), pmf_target=pmf_target_boot, counts=cts, cts.orig=cts.orig)
        }
# End rpm.simulate.large.population.worker

        if(is.null(object$boot.coef)){
          boot.coefficients <- matrix(rep(object$coefficients, nsim),byrow=TRUE,nrow=nsim)
        }else{
          if(nsim <= nrow(object$boot.coef)){
            boot.coefficients <- object$boot.coef[
              sample.int(nrow(object$boot.coef), size=nsim, replace=FALSE),,drop=FALSE]
          }else{
            warning("The number of boostrap samples in the passed fit is less than the requested 'nsim'. The effective size of 'nsim' is capped at the number of bootstrap samples. To increase this, refit the 'rpm' model with the 'control.rpm(nbootstrap=nsim)' and rerun this command.")
            boot.coefficients <- rbind(object$boot.coef,
             object$boot.coef[sample.int(nrow(object$boot.coef),
                              size=nsim-nrow(object$boot.coef), replace=TRUE),] )
          }
        }
   #    boot.coefficients <- matrix(rep(object$coefficients, nsim),byrow=TRUE,nrow=nsim)
        if(control$ncores > 1 & nsim > 1){
          if(verbose) message(sprintf("Starting parallel simulation using %d cores.",control$ncores))
          out.list <-
           foreach::foreach (i=1:nsim, .packages=c('rpm')
           ) %dorng% {
           rpm.simulate.large.population.worker(object, boot.coefficients[i,], N, X_w_rel, Z_w_rel, pmfW_N, pmfM_N, pmfW, pmfM, gm, gw, numX, numZ, bootstrap)
           }
          parallel::stopCluster(cl)
        }else{
          out.list <- vector(nsim, mode="list")
          for( i in 1:nsim ){
           out.list[[i]] <- rpm.simulate.large.population.worker(object, boot.coefficients[i,], N, X_w_rel, Z_w_rel, pmfW_N, pmfM_N, pmfW, pmfM, gm, gw, numX, numZ, bootstrap)
          }
        }
  
      }else{
        # Use the direct (small population) method
        # These are the categories of women
        Ws <- rep(x=seq_along(pmfW_N),times=round(pmfW_N))
        if(length(Ws) > num_women){
         Ws <- Ws[-sample.int(length(Ws),size=length(Ws)-num_women,replace=FALSE)]
        }
        if(length(Ws) < num_women){
         Ws <- c(Ws,sample.int(length(pmfW_N),size=num_women-length(Ws),prob=pmfW_N,replace=TRUE))
        }
        Ws <- Ws[sample.int(length(Ws))]
        Ms <- rep(x=seq_along(pmfM_N),times=round(pmfM_N))
        if(length(Ms) > num_men){
         Ms <- Ms[-sample.int(length(Ms),size=length(Ms)-num_men,replace=FALSE)]
        }
        if(length(Ms) < num_men){
         Ms <- c(Ms,sample.int(length(pmfM_N),size=num_men-length(Ms),prob=pmfM_N,replace=TRUE))
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

        rpm.simulate.small.population.worker <- function(object, num_women, num_men, Jw, Jm, U_star, V_star, Ws, Ms){
          eta  <- -log(-log(matrix(stats::runif(num_women * num_men), num_women)))
          zeta <- -log(-log(matrix(stats::runif(num_women * num_men), num_men)))
  
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
          mu=data.frame(mu, type=c(Ws, Ms), sampled=rep(TRUE,nrow(mu)))
          colnames(mu)[c(1,3)] <- c(object$Xid, object$pair_id)
          Xdata <- subset(mu, gender=="F")
          Xdata <- data.frame(Xdata,as.data.frame(object$Xu[Ws,-1,drop=FALSE]))
          Zdata <- subset(mu, gender=="M")
          Zdata <- data.frame(Zdata,as.data.frame(object$Zu[Ms,-1,drop=FALSE]))
          colnames(Zdata)[match(object$Xid,colnames(Zdata))] <- object$Zid
          # random permute to add randomness
          list(population=list(Xdata=Xdata[sample.int(nrow(Xdata)),],
                               Zdata=Zdata[sample.int(nrow(Zdata)),]) )
        }
        if(control$ncores > 1 & nsim > 1){
          if(verbose) message(sprintf("Starting parallel simulation using %d cores.",control$ncores))
          out.list <-
           foreach::foreach (i=1:nsim, .packages=c('rpm')
           ) %dorng% {
           rpm.simulate.small.population.worker(object, num_women, num_men, Jw, Jm, U_star, V_star, Ws, Ms)
           }
          parallel::stopCluster(cl)
        }else{
          out.list <- vector(nsim, mode="list")
          for( i in 1:nsim ){
           out.list[[i]] <- rpm.simulate.small.population.worker(object, num_women, num_men, Jw, Jm, U_star, V_star, Ws, Ms)
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
