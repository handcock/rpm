    loglikfun_nog <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      -loglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    } 
    eqfun_nog <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      eqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
    gloglikfun_nog <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      -gloglik_nog(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
    jeqfun_nog <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, N, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      jeqcond_nog(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
