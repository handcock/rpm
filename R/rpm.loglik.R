    loglikfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      gw <- theta[(NumBeta+NumGammaW+NumGammaM)+1]
      gm <- log(1-exp(gw))
      -loglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    } 
    eqfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      gw <- theta[(NumBeta+NumGammaW+NumGammaM)+1]
      gm <- log(1-exp(gw))
      eqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
#   ggloglikfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw,
#   gm, N, sampling, constraints){
#     nloptr::nl.grad(theta, loglikfun,Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
#             pmfW=pmfW,pmfM=pmfM,pmf=pmf,counts=counts, gw=gw, gm=gm, N=N, sampling=sampling)
#   }
#   gloglikfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
##    numDeriv::grad(func=loglikfun,x=theta,Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
##            pmfW=pmfW,pmfM=pmfM,pmf=pmf,counts=counts, gw=gw, gm=gm, N=N, sampling=sampling, constraints=constraints)
#     nloptr::nl.grad(theta, loglikfun,Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
#             pmfW=pmfW,pmfM=pmfM,pmf=pmf,counts=counts, gw=gw, gm=gm, N=N, sampling=sampling, constraints=constraints)
#   }
    gloglikfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      gw <- theta[(NumBeta+NumGammaW+NumGammaM)+1]
      gm <- log(1-exp(gw))
      -gloglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
#   jeqfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N,
#   sampling, constraints){
#     nloptr::nl.jacobian(theta, eqfun,Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
#                 pmfW=pmfW,pmfM=pmfM, pmf=pmf,counts=counts,gw=gw,gm=gm,N=N,sampling=sampling, constraints=constraints)
#   }
    jeqfun_default <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      gw <- theta[(NumBeta+NumGammaW+NumGammaM)+1]
      gm <- log(1-exp(gw))
 #    d <- numDeriv::jacobian(func=eqfun,x=theta,Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
 #                pmfW=pmfW,pmfM=pmfM, pmf=pmf,counts=counts,gw=gw,gm=gm,N=N,sampling=sampling)
      geqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
