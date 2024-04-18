    loglikfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      -loglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    } 
    eqfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      eqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
    gloglikfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      -gloglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
    jeqfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      jeqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
#
#   These are the intercept only versions
#
    loglikfun_intercept <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      GammaW <- theta[(1:NumGammaW)]
      GammaM <- theta[(NumGammaW)+(1:NumGammaM)]
      -loglik(0, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    } 
    eqfun_intercept <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      GammaW <- theta[(1:NumGammaW)]
      GammaM <- theta[(NumGammaW)+(1:NumGammaM)]
      eqcond(0, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
    gloglikfun_intercept <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      GammaW <- theta[(1:NumGammaW)]
      GammaM <- theta[(NumGammaW)+(1:NumGammaM)]
      -gloglik(0, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)[-1]
    }
    jeqfun_intercept <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      GammaW <- theta[(1:NumGammaW)]
      GammaM <- theta[(NumGammaW)+(1:NumGammaM)]
      jeqcond(0, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)[,-1]
    }
