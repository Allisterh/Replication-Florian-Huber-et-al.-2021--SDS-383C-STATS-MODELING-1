
Rcpp::sourceCpp('threshold_functions.cpp')
#sim.data <- sim_piecewise()
sparse.SAVS <- function(Yraw, Xraw, h = 0, nsave=1000,nburn=1000,thin=1,p=0,a_i=0.001,b_i=0.001,cons=0,nu=0.3,prior.mean=c(1,0),sv=FALSE,store.draws=TRUE,cons.model=FALSE,shrinkage="DL",store.last=FALSE,  MCMC.sparse = TRUE, sample.nu = FALSE){
  
#  Yraw=Y.i; Xraw=X.i; h = 0; nsave=1000; nburn=1000; thin=1; p=0; a_i=0.001; b_i=0.001; cons=0; nu=0.3;prior.mean=c(11,1);sv=TRUE; store.draws=TRUE;cons.model=FALSE;shrinkage="DL";store.last=FALSE;  MCMC.sparse = TRUE; sample.nu = FALSE
  
  #Some standard values for prior hyperparameters that are borrowed from the literature. a.tau refers to the parameter that controls the excess kurtosis of the NG prior while the remaining
  #things are related to the MNIG prior and closely resemble the values outlined in Malsiner-Walli and Wagner (2015)
  a.tau <- 0.1
  scale0.ssvs <- 0.000025 #add in function argument
  scale1.ssvs <- 1
  d.0 <- 5
  d.1 <- 4
  
  #store.last <- FALSE
  #--------------s-----------------------------Load packages -----------------------------------------------------------------#
  require(GIGrvg)
  require(mvtnorm)
  require(stochvol)
  require(MASS)
  # -------------------------------------------Some settings-----------------------------------------------------------------#
  #Number of saved draws and burn-ins
  X.full <- Xraw
  ntot <- nsave+nburn
  #-----------------------------------------Data preliminaries---------------------------------------------------------------#
  if (p>=1){
    if (cons==1){
      Xraw <- cbind(mlag(Xraw,p),1)
    }else{
      Xraw <- cbind(mlag(Xraw,p))
    }
  }else{
    if (cons==1) Xraw <- cbind(Xraw,1)  
  }
  
  Y <- Yraw[(p+1):(nrow(Yraw)-h),]
  X <- Xraw[(p+1):(nrow(Xraw)-h),]
  T <- nrow(X)
  M <- ncol(Y)
  K <- ncol(X)
  #-----------------------------------------Initialize Gibbs sampler---------------------------------------------------------#
  A_OLS <- A_draw <- ginv(crossprod(X))%*%crossprod(X,Y)
  Em <- Em.str <-  Y-X%*%A_OLS
  SIGMA <- SIGMA_OLS <- crossprod(Y-X%*%A_OLS)/(T-K)
  #chol_SIG <- t(chol(SIGMA_OLS))
  D <- diag(diag(SIGMA))
  V_prior <- diag(K)*10^2
  V_priorinv <- diag(1/diag(V_prior))
  
  A_prior <- matrix(0,2*K,1)
  A_prior[prior.mean[1],] <- prior.mean[2]
  a_prior <- as.vector(A_prior)
  
  theta <- matrix(10,K,1)
  omega_mat <- matrix(0.01,K,1)
  sqrtomega <- diag(K)*0.01
  omega.sd <- omega.prior <- matrix(10,K,1)
  
  
  phi0 <- psi0  <- matrix(1/(2*(K)), 2*K, 1)
  T_i <- matrix(0,2*K,1)
  a.DL <- 1/(2*(K))
  
  if (shrinkage=="SSVS"){
    indicators <- matrix(1,2*K,1)
  }
  
  #SV for unrestricted and restricted model
  startpara.sparse <- startpara <- list(mu = 0, phi = 0.9, sigma = 0.1, nu = Inf, rho = 0, beta = NA,latent0 = 0)
 

 # hv <- svdraw$latent
  para <- list(mu=-3,phi=.9,sigma=.2)
  H.sparse <- H <- matrix(-3,T,1)
  eta <- list()
  zeta.list <- list()
  #MH preliminaries
  accept <- 0
  accept2 <- 0
  scale1 <- .43
  scale2 <- .45
  #------------------------------------------Create storage matrices----------------------------------------------------------#
  if (store.draws){
    if (!store.last){
      A_store <- array(NA,c(round(thin*nsave),T,K,1))
      ALPHA_store <- array(NA,c(round(thin*nsave),T,K,1))
      hv_store <- array(NA,c(round(thin*nsave),T,1))
    }else{
      A_store <- array(NA,c(round(thin*nsave),K,1))
      ALPHA_store <- array(NA,c(round(thin*nsave),K,1))
      hv_store <- array(NA,c(round(thin*nsave),1))
    }
    A0_store <- matrix(NA,round(thin*nsave),K)
    coefsparse_store <- matrix(NA,round(thin*nsave),2*K)
    
    theta_store <- array(NA,c(round(thin*nsave),K,1))
    omega.sd_store <- array(NA,c(round(thin*nsave),K,1))
    omega.prior_store <- array(NA,c(round(thin*nsave),K,1))
    pars_store <- array(NA,c(round(thin*nsave),3,1))
    hvsparse_store <- array(NA,c(round(thin*nsave),T,1))
    
    kappa_store <- matrix(NA,round(thin*nsave),1)
    Lik_store <- matrix(NA,round(thin*nsave),1)
    nu_store <- matrix(NA,round(thin*nsave),1)
  }
  a_store <- matrix(NA,round(thin*nsave),p+1)
  #--------------------------------------------Start big Gibbs loop-----------------------------------------------------------#
  pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
  save.set <- round(seq(nburn,ntot,length.out=round(thin*nsave)))
  save.count <- 0
  start <- Sys.time()
  for (irep in 1:ntot){
    #Step I: Sample regression coefficients
    scale.cons <- ifelse(cons.model,1e-8,1)
    Y.i <- Y-(X%*%A_draw)#*exp(-0.5*H[,mm])
    X.i <- X%*%diag(as.numeric(omega_mat))#*exp(-0.5*H[,mm])
    #Draw for the time-varying part ALPHA
    if (cons.model){
      ALPHA1 <- matrix(0,ncol(X.i),T)
    }else{
      ALPHA1 <- KF_fast(t(as.matrix(Y.i)), X.i,as.matrix(exp(H)),t(matrix(scale.cons,K,T)),K, 1, T, matrix(0,K,1), diag(K)*1e-10)
     }
    
    G <- X*t(ALPHA1)
    #Create a new X matrix X_new = (X__,G)
    Xnew <- cbind(X,G)*as.numeric(exp(-H/2))
    Ynew <- Y*as.numeric(exp(-H/2))
    V.prior.inv <- diag(1/c(theta,omega.prior))
    V_post <- try(solve((crossprod(Xnew)) + V.prior.inv),silent=TRUE)
    if (is(V_post,"try-error")) V_post <- ginv((crossprod(Xnew) + V.prior.inv))
    A_mean <- V_post %*% (crossprod(Xnew, Ynew)+V.prior.inv%*%A_prior)
    
    A.draw.i <- try(A_mean+t(chol(V_post))%*%rnorm(2*K,0,1),silent=T)
    if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_mean,V_post)
    
    
    A_draw <- A.draw.i[1:K]
    diag(sqrtomega) <- A.draw.i[(K+1):(2*K)]
    # #Introduce non-identified element
    # uu <- runif(1,0,1)
    # if (uu>0.5){
    #   sqrtomega <- -sqrtomega
    #   ALPHA1 <- -ALPHA1
    # }
    omega_mat <- diag(sqrtomega)
 
    ALPHA <- t((as.numeric(A_draw)+diag(sqrtomega)*ALPHA1))
    
    if (MCMC.sparse){
    #Do sparsification within MCMC
    X.tilda <- cbind(X,G)*as.numeric(exp(-H/2))
    b.tilda <- matrix(c(as.numeric(A_draw),diag(sqrtomega)))
    grid.nu <- seq(0.3,4,length.out=25)
    if (sample.nu){
      grid.lik.nu <- matrix(NA,length(grid.nu))
      for (ii in 1:length(grid.nu)){
        b.sparse.i <- sparsify(X.tilda,b.tilda,2*K,grid.nu[[ii]])  
        Lik.a <- sum(dnorm(Y,X.tilda%*%b.sparse.i,exp(H.sparse/2),log=TRUE))  
        grid.lik.nu[ii,] <- Lik.a
      }
      probs.i <- exp(grid.lik.nu-max(grid.lik.nu))/sum(exp(grid.lik.nu-max(grid.lik.nu)))
      nu <- sample(grid.nu,1,prob=probs.i)
    }
    b.sparse <- sparsify(X.tilda,b.tilda,2*K,nu)  
    A.sparse <- t((as.numeric(b.sparse[1:K])+(b.sparse[(K+1):(2*K)])*ALPHA1))
    }
    Em.sparse <- matrix(0,T,1)
    for (t in 1:T){
      Em.str[t,] <- Y[t] - X[t,]%*%(as.numeric(A_draw)+diag(sqrtomega)*ALPHA1)[,t]
      if (MCMC.sparse){
        Em.sparse[t,] <- Y[t]-X[t,]%*%A.sparse[t,]
      } 
    }
    
    if (shrinkage=="DL"){
      #This block samples the individual and lag-specific shrinkage parameters for the constant parameter part
      theta.full <- c(A_draw,omega_mat)
      omega.full <- c(theta,omega.prior)
      lambda <- a.DL*(2*K)
      
      kappa.draw <- GIGrvg::rgig(n=1,lambda=lambda - (2*K),2*sum(abs(theta.full)/phi0)+1e-20,1)
      
      for (jj in 1:(2*K)){
        mu0 <- phi0[jj,1]*kappa.draw/abs(theta.full[jj])
        psi0[jj,1] <- 1 / GIGrvg::rgig(1, lambda = -.5, 1, 1 / mu0^2)
        T_i[jj,1] <-  GIGrvg::rgig(n = 1, lambda = a.DL - 1, 2 * abs(theta.full[jj])+1e-20, 1)
      }
      phi0 <- T_i/sum(T_i)
      
      # Draw from proposal density
      prop.sd <- 1e-1
      a.star <- truncnorm::rtruncnorm(1, a=1/(2*K), b=1/2, mean=a.DL, sd=prop.sd)
      
      # Ratio q(a.star | a)/q(a | a.star)
      q.ratio <- (pnorm((1-a.DL)/prop.sd)-pnorm((1/(2*K)-a.DL)/prop.sd))/(pnorm((1-a.star)/prop.sd)-pnorm((1/(2*K)-a.star)/prop.sd))
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      pi.ratio <- 2^(2*K*(a.DL-a.star))*(gamma(a.DL)/gamma(a.star))^(2*K) * prod(psi0^(a.star-a.DL))
      
      # Acceptance probability
      alpha <- min(1, (q.ratio*pi.ratio) )
      
      # Accept/reject algorithm
        if ( runif(1) < alpha){
          a.DL <- a.star
        }
      theta  <- matrix(as.vector(psi0[1:K]) * phi0[1:K]^2 * (kappa.draw)^2, K, 1) #1 VAR COEFF 0 SQRT OMEGA PART
      omega.prior <- matrix(as.vector(psi0[(K+1):(2*K)]) * phi0[(K+1):(2*K)]^2 * (kappa.draw)^2, K, 1) #1 VAR COEFF 0 SQRT OMEGA PART
    }else if (shrinkage=="NG" || shrinkage=="LASSO"){
      if (shrinkage=="LASSO") a.tau <- 1
      theta.full <- c(A_draw[1:(K)],omega_mat[1:(K)])
      omega.full <- c(theta,omega.prior)
      
      kappa_j <- rgamma(1,a_i+a.tau*2*K,b_i+a.tau/2*sum(omega.full))
      for (jj in 1:(2*K)){
       if (jj<=K){
         theta[jj,] <- rgig(n=1,lambda=a.tau-0.5, theta.full[[jj]]^2+1e-10,a.tau*kappa_j)
       }else{
         omega.prior[jj-K,] <-  rgig(n=1,lambda=a.tau-0.5, theta.full[[jj]]^2+1e-10,a.tau*kappa_j)
       }
      }
      
      if (shrinkage !="LASSO"){
        
        a.tau <- get.a.theta(a.tau,omega.full,kappa_j,scale1)
        
      }
      
    }else if  (shrinkage=="SSVS"){
      theta.full <- c(A_draw[1:(K)],omega_mat[1:(K)])

      
      psi0 <- 1/rgamma(2*K,d.0+1/2, d.1+.5*theta.full^2/(indicators*scale1.ssvs+(1-indicators)*scale0.ssvs))
      
      w <- rbeta(1, 1+sum(indicators==1),1+2*K-sum(indicators==1)) #implies uniform prior on weights
      
      
      for (jj in 1:(2*K)){
        p0 <- dnorm(theta.full[jj],0,sqrt(scale0.ssvs*psi0[[jj]]))*(1-w)#+1e-45
        p1 <- dnorm(theta.full[jj],0,sqrt(scale1.ssvs*psi0[[jj]]))*w#+1e-45
        prob.jj <- p1/(p0+p1)
        
        if (prob.jj>runif(1)){
          indicators[jj,1] <- 1
          if (jj<=K){
            theta[jj,] <- scale1.ssvs*psi0[[jj]]
          }else{
            omega.prior[jj-K,] <-  scale1.ssvs*psi0[[jj]]
          }
        }else{
          indicators[jj,1] <- 0
          if (jj<=K){
            theta[jj,] <- scale0.ssvs*psi0[[jj]]
          }else{
            omega.prior[jj-K,] <-  scale0.ssvs*psi0[[jj]]
          }
        }    
      }
    }else if (shrinkage=="horse"){
      #move upwards later on
      if (irep ==1){
        tau.hs <- 1
        lambda.hs <- matrix(1,2*K,1)
        nu.hs <- matrix(1,2*K,1)
        zeta.hs <- 1
      }
      theta.full <- c(A_draw[1:(K)],omega_mat[1:(K)])
      lambda.hs <- invgamma::rinvgamma(2*K,shape=1, rate  =1/nu.hs+theta.full^2/(2*tau.hs))#1/rgamma(2*K,1,1/nu.hs+theta.full^2/(2*tau.hs))
      tau.hs <- invgamma::rinvgamma(1,shape=(2*K+1)/2, rate=1/zeta.hs+sum(theta.full^2/lambda.hs)/2)
      nu.hs <- invgamma::rinvgamma(2*K, shape=1, rate=1+1/lambda.hs)
      zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
      
      theta[1:K,1] <- lambda.hs[1:K]*tau.hs
      omega.prior[1:K,1] <- lambda.hs[(K+1):(2*K)]*tau.hs
    }else if (shrinkage=="RW"){
      theta <- matrix(1e-10,K,1)
      omega.prior <- matrix(1e-10,K,1)
    }
    if (shrinkage != "RW"){
    theta[theta<1e-10] <- 1e-10
    omega.prior[omega.prior<1e-10] <- 1e-10
    }
    #Step IV: Sample covariance parameters
    
    if (sv){
      #Sample stochastic volatilities
      res <- svsample_fast_cpp(Em.str, startpara = startpara, startlatent = H)
      startpara[c("mu","phi","sigma")] <- as.list(res$para[, c("mu", "phi", "sigma")])
      H <- drop(res$latent)
      
      res.sparse <- svsample_fast_cpp(Em.sparse, startpara = startpara.sparse, startlatent = H.sparse)
      startpara.sparse[c("mu","phi","sigma")] <- as.list(res.sparse$para[, c("mu", "phi", "sigma")])
      H.sparse <- drop(res.sparse$latent)
      
    }else{
      S_1 <- a_i+T/2
      S_2 <- b_i+crossprod(Em.str)/2
      
      sig_eta <- 1/rgamma(1,S_1,S_2)
      H[,1] <- log(sig_eta)
      
      H.sparse[,1] <- log(1/rgamma(1,a_i+T/2,b_i+crossprod(Em.sparse)/2))
    }
    
    H[H< -14] <- -14 #Offset in case we overfit heavily
    
    #Step VIII: Store draws after burn-in/ Compute forecasts/ Impulse responses etc.
    if (irep %in% save.set){
      save.count <- save.count+1
      if (store.draws){
        if (!store.last){
          if (MCMC.sparse)  A_store[save.count,,,] <- A.sparse else A_store[save.count,,,] <- t(ALPHA1)
          ALPHA_store[save.count,,,] <- (ALPHA)
          hv_store[save.count,,] <- (H)
        }else{
          if (MCMC.sparse)  A_store[save.count,,] <- A.sparse[T,] else A_store[save.count,,] <- t(ALPHA1[,T])
          ALPHA_store[save.count,,] <- (ALPHA[T,])
          hv_store[save.count,] <- (H[T,])
        }
        coefsparse_store[save.count,] <- b.sparse
        A0_store[save.count,] <-  A_draw
        theta_store[save.count,,] <- matrix(theta,K,1)
        omega.sd_store[save.count,,] <- omega_mat
        omega.prior_store[save.count,,] <- omega.prior
        if (shrinkage=="DL") kappa_store[save.count,] <- a.DL
        
        nu_store[save.count,] <- nu
        if (MCMC.sparse){
          hvsparse_store[save.count,,] <- H.sparse
        } 
        if (sv){
          pars_store[save.count,,] <- unlist(startpara[c("mu","phi","sigma")])
        } 
      }
    }
       setTxtProgressBar(pb, irep)
    
  }
  end <- Sys.time()
  print(end-start)   
  
#  A_mean <- apply(A_store,c(2,3,4),mean)
  A0_mean <- apply(A0_store,2,mean)
  omega_mean <- apply(omega.sd_store,c(2,3),mean)
  
  #H.mean <- apply(hv_store,c(2,3),mean)
  #G.mean <- X*A_mean[,,1]
  #Create a new X matrix X_new = (X__,G)
  #X.mean <- cbind(X,G.mean)
  #b.full <- matrix(c(A0_mean,omega_mean))
  
  if (sv){
    if (!store.last){
      H_post <- apply(exp(hv_store),c(2,3),mean)
      pars_post <- apply(pars_store,c(2,3),mean)
    }else{
      H_post <- apply(exp(hv_store),c(1),mean)
      pars_post <- apply(pars_store,c(2,3),mean)
    }
  } 
  
  if (!MCMC.sparse){
  
  norm_vec <- function(x) sqrt(sum(x^2))
  
  b.sparse <- matrix(0,2*K,1)
    for (jj in seq_len(2*K)){
      mu.jj <- 1/abs(b.full[jj,1])^nu
      norm.j <- norm_vec(X.mean[,jj])^2
      if ((abs(b.full[jj,1])*norm.j) < mu.jj){
        b.sparse[jj,1] <- 0
      }else{
        b.sparse[jj,1] <- sign(b.full[jj,1])* 1/norm.j * ((abs(b.full[jj,1])*norm.j)-mu.jj)
      }
    }

  
    ALPHA <- t((as.numeric(A0_mean)+diag(sqrtomega)*t(A_mean[,,1])))
    
    ret.object <- list(A.non=ALPHA_store,A.thrsh=b.sparse,A.non.thrsh=c(A0_mean,omega_mean),hv=hv_store)

    
  }else{
   # A.sparse <- apply(A_store[,,,1],c(2,3),median)
    ret.object <- list(A.non=ALPHA_store,A.thrsh=A_store, coef.state.thrsh=coefsparse_store, coef.state.non=c(A0_mean,omega_mean), hv=hv_store,nu=nu_store)
  
  }
# browser()
  return(ret.object)
}

