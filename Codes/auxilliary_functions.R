mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}
impulsdtrf <- function(B,smat,nstep)
  ### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
  ### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
  ###                 impulse responses" it should have the property that crossprod(t(smat))=sigma,
  ###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
  ###                 way to get such a smat is to set smat=t(chol(sigma)).  To get the smat
  ###                 corresponding to a different ordering, use
  ###                 smat = t(chol(P %*% Sigma %*% t(P)) %*% P), where P is a permutation matrix.
  ###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
  ###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
  ###                 positive powers of L in A), smat=A_0^(-1).
  ###                 In general, though, it is not required that smat be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
{
  ##-----debug--------
  ##browser()
  ##------------------
  neq <- dim(B)[1]
  nvar <- dim(B)[2]
  lags <- dim(B)[3]
  dimnB <- dimnames(B)
  if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
  response <- array(0,dim=c(neq,nvar,nstep+lags-1));
  response[ , , lags] <- smat
  response <- aperm(response, c(1,3,2))
  irhs <- 1:(lags*nvar)
  ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol=neq)
  B <- B[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
  B <- matrix(B,nrow=nvar)
  for (it in 1:(nstep-1)) {
    #browser()
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar
    ilhs <- ilhs + nvar
  }
  ## for (it in 2:nstep)
  ##       {
  ##         for (ilag in 1:min(lags,it-1))
  ##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
  ##       }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  response <- aperm(response[ , -(1:(lags-1)), ], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
  dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
  ## dimnames(response)[2] <- dimnames(smat)[1]
  ## dimnames(response)[1] <- dimnames(B)[2]
  return(response)
}
get.pred <- function(zt=t(X[nrow(X),,drop=F]), laglen=laglen, m=m, K=K, At=At, SIGMA=SIGMA, fhorz=horz.fcst, Yraw.full=Y, Xraw.full=X, sl.joint=c(1,2,3)){
  get.comp <- get_companion(t(At),c(m,0,laglen))
  Mm <- get.comp$MM
  Jm <- get.comp$Jm
  
  z1 <- zt
  Mean00 <- zt
  Sigma00 <- matrix(0,K,K)
  Sigma00[1:m,1:m] <- SIGMA
  
  LPS.store <- matrix(NA,fhorz)
  LPS.marg.store <- matrix(NA,length(sl.joint),fhorz)
  pred.store <- matrix(NA,length(sl.joint),fhorz)
  for (ih in 1:fhorz){
    #Predict E(Y(t+n|t)) n=1,..,fhorz and compute cholesky factor of forecast error variance-covariance
    z1 <- Mm%*%z1
    cholSig <- try(t(chol(Sigma00[1:m,1:m])),silent=TRUE)
    if (is(cholSig, "try-error")){
      yf <- try(rmvnorm(1,z1[1:M],Sigma00[1:m,1:m]),silent=TRUE)
      if (is(yf,"try-error")) yf <- z1[1:m] #should be the option of last resort (happens in approx. 1 out of 100,000 cases)
    }else{
      yf <- z1[1:m]+cholSig%*%rnorm(m,0,1)
    }
    Mean00 <- Mm%*%Mean00
    Sigma00 <- Mm%*%Sigma00%*%t(Mm)+Jm%*%SIGMA%*%t(Jm)
    meanyt <- Mean00[1:m]
    varyt <- Sigma00[1:m,1:m]
    # sl.forecast <- nrow(Xraw.full)-hold.out+h.out-1+ih
    # if (sl.forecast>nrow(Xraw.full)) next
    # LPS.store[ih] <- mvtnorm::dmvnorm(Yraw.full[sl.forecast,sl.joint],meanyt[sl.joint],varyt[sl.joint, sl.joint],log=TRUE)
    # LPS.marg.store[,ih] <- dnorm(Yraw.full[sl.forecast,sl.joint], meanyt[sl.joint], sqrt(diag(varyt[sl.joint, sl.joint])),log=TRUE)
    pred.store[,ih] <- yf[sl.joint]
  }
  
  pred.list <- list(yhat=pred.store, LPS=LPS.store, LPS.marg=LPS.marg.store)
  return(pred.list)
}


norm_vec <- function(x) sqrt(sum(x^2))

sparsify <- function(X,b,K.max,nu){
  b.sparse <- matrix(0,K.max,1)
  for (jj in seq_len(K.max)){
    mu.jj <- 1/abs(b[jj,1])^nu
    norm.j <- norm_vec(X[,jj])^2
    if ((abs(b[jj,1])*norm.j) < mu.jj){
      b.sparse[jj,1] <- 0
    }else{
      b.sparse[jj,1] <- sign(b[jj,1])* 1/norm.j * ((abs(b[jj,1])*norm.j)-mu.jj)
    }
  }
  return(b.sparse)
}

get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  Jm <- matrix(0,nkk,nn)
  MM <- matrix(0,nkk,nkk)
  Jm[1:nn,1:nn] <- diag(nn)
  if (nd==1){
    MM <- rbind(t(Beta_), cbind(diag((nl-1)*nn),matrix(0,(nl-1)*nn,nn+nd)),c(matrix(0,nd,nn*nl),nd))
  }else{
    MM <- rbind(t(Beta_), cbind(diag((nl-1)*nn),matrix(0,(nl-1)*nn,nn+nd)))
  }
  return(list(MM=MM,Jm=Jm))
}
get.a.theta <- function(a.theta,theta,lambda2,scale1){
  
  a_prop <- exp(rnorm(1,0,scale1))*a.theta
  post_a_prop <- atau_post(atau=a_prop,thetas = theta,lambda2 = lambda2)
  post_a_old <- atau_post(atau=a.theta,thetas =theta,lambda2 = lambda2)
  post.diff <- post_a_prop-post_a_old+log(a_prop)-log(a.theta)
  post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
  
  if (post.diff> log(runif(1,0,1))){
    a.theta <- a_prop
    #accept <- accept+1
  }
  return(a.theta)
}
atau_post <- function(atau=a_tau,lambda2=lambda2_tau,thetas=theta,k=length(theta),rat=1/.1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)  
  return(logpost)
}


