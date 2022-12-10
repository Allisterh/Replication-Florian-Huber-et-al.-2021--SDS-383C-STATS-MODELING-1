require(bvarsv)
source("auxilliary_functions.R")
source("ng_SAVS.R")

priors <- c("horse") #Set the prior (horse=horseshoe, LASSO=LASSO, SSVS=MNIG, NG=NG, DL= DL)
nu.grid <- c(2) #Sets the penalty term for SAVS ( 1 / abs(beta)^nu)
nsave <- 500
nburn <- 500
horz.fcst <- 4 #select forecast horzion
h <- 4
sl.marg <- seq(1,3) #singles out variables in the set of focus variables
nhor <- 20 #Selects the impulse response horizon

grid.full <- expand.grid(priors,nu.grid) #in case one would like to compute the model over a grid of different priors and nu (add above and write a loop with loop index run that goes from 1:nrow(grid.full))
run <- 1
sl.grid <- grid.full[run,] #This selects a configuration from the grid

nu <- sl.grid[[2]]
prior.choice <- as.character(sl.grid[[1]])
#Model preliminaries
laglen <- 2 #sets the number of lags
load("datasetMNg.RData")
size <- "small" #select the model size (small, medium, large)

X.full <- apply(X.list[["large"]],2,function(x) (x-mean(x))/sd(x)) #Standardize the data prior to estimation to avoid scaling issues with the different priors
if (size=="large"){
  sl.X <- seq(1,ncol(X.full))
}else if (size=="medium"){
  sl.X <- c("GDPC1", "PCECC96","FPIx","CE16OV", "CES0600000007", "GDPCTPI", "CES0600000008","FEDFUNDS")
}else{
  sl.X <- c("GDPC1","GDPCTPI","FEDFUNDS")
}
Yraw <- X.full[,sl.X] #slct the model size selected above

#Set up data for a large BVAR
Xraw <- mlag(Yraw,laglen) #Lags the columns in Yraw 
Y <- Yraw[(laglen+1):(nrow(Yraw)-h),]
X <- Xraw[(laglen+1):(nrow(Xraw)-h),]

#Get dimensions of the dataset
m <- ncol(Y)
T <- nrow(X)
K <- ncol(X)

#Creates matrices that store the posterior distributions
A.matrix <- array(NA,c(nsave,T,m,K)) #Posterior distribution of coefficients across equations
A0.matrix <- array(0,c(nsave,T,m,m)) #Posterior of coefficients on contemporenous y's
H.matrix <- array(NA,c(nsave,m,T)) #Stores the log-volatilities

#This is embarassingly parallel and one could use snowfall etc. here; Code for the Sun Grid Engine (SGE) is available upon request such that each equation is estimated on a separate cluster-CPU.
for (i in seq_len(m)){
  #Augment X matrix with contemporenous values of the preceeding i-1 endogenous quantities
  sl.i <- 1:(i-1)
  if (i==1) sl.i <- NULL else sl.i <- 1:(i-1)
  X.i <- cbind(Y[,sl.i],X)
  Y.i <- Y[,i,drop=F]
  #Estimates a state space model using SAVS
  sim.1 <- sparse.SAVS(Y.i,X.i,h=0,p=0,nsave=nsave,nburn=nburn,sv=TRUE,prior.mean=c(0,0),nu=nu,shrinkage=prior.choice)
  
  A0.matrix[,,i,sl.i] <- -1*sim.1$A.thrsh[,,sl.i,1]
  if (i==1){
    A.matrix[,,i,] <- sim.1$A.thrsh[,,,1]
  }else{
    A.matrix[,,i,]<- sim.1$A.thrsh[,,-sl.i,1]
  }
  H.matrix[,i,] <- sim.1$hv
}


A.red.array <- array(NA,c(nsave,T,K,m)) #Posterior distribution of reduced-form coefficients
Sig.red.array <- array(NA,c(nsave,T,m,m)) #Posterior distribution of error variance-covariance
pred.array <- array(0,c(nsave,length(sl.marg),horz.fcst)); dimnames(pred.array) <- list(NULL,c(colnames(Y)[sl.marg]),NULL)
#Creates a big IRF object that stores the IRFs for each point in time. The second dimension refers to the responses to shocks to the third dimension.
IRF.array <- array(NA, c(nsave, m, m, nhor, T)); dimnames(IRF.array) <- list(NULL, colnames(Y), colnames(Y), NULL)
for (irep in 1:nsave){
  for (t in 1:T){
    A0 <- A0.matrix[irep,t,,]
    diag(A0) <- 1
    A0inv <- t(solve(A0))
    At <- A0inv%*%(A.matrix[irep,t,,])
    A.red.array[irep,t,,] <- t(At)
    Sig.red <- A0inv%*%diag(exp(H.matrix[irep,,t]))%*%t(A0inv)
    Sig.red.array[irep,t,,] <- A0inv%*%diag(exp(H.matrix[irep,,t]))%*%t(A0inv)
    
    #This part computes impulse responses based on using Chris Sim's impulsdtrf function
    #Step I: Construct an array that includes the VAR coefs in time t
    PHI_array <- array(0,c(m,m,laglen))
    Att <- t(At)
    for (ss in 1:laglen){
      PHI_array[,,ss] <- t(Att[((ss-1)*m+1):(ss*m),])
    }
    
    
    IRF.array[irep,,,,t] <- impulsdtrf(PHI_array, t(chol(Sig.red)), nhor)
  }
  #Do predictions out-of-sample under the assumption that the states are set equal to their expected values (see Koop et al. 2009 for details)
  pred.sparse <- get.pred(zt=t(Xraw[nrow(Xraw)-h+1,,drop=F]), laglen=laglen, m=m, K=K, At=At, SIGMA=Sig.red, fhorz=horz.fcst, Yraw.full=Yraw, Xraw.full=Xraw, sl.joint=sl.marg)
  pred.array[irep,,] <- pred.sparse$yhat
}

fcst.median <- apply(pred.array,c(2,3),median) #Posterior median of predictive distribution for Y(T+j) (j=1,...,horz.fcst)
Sig.mean <- apply(Sig.red.array, c(2,3,4),median) #Posterior median of Sigma(t)
A.mean <- apply(A.red.array, c(2,3,4),median) #Posterior mean of A(t) (VAR coefs)

PIP.mean <- apply((A.matrix[,1,,] != 0)*1,c(2,3),mean) #PIPs across equations
PIP.mean.sigma <- apply((A0.matrix[,1,,] !=0)*1, c(2,3),mean)
diag(PIP.mean.sigma) <- 1
colnames(PIP.mean.sigma) <- rownames(PIP.mean.sigma) <- colnames(Yraw)

#Computes posterior quantiles of the IRFs
IRF.median <- apply(IRF.array,c(2,3,4,5), median)
IRF.low <- apply(IRF.array,c(2,3,4,5), quantile, 0.16)
IRF.high <- apply(IRF.array,c(2,3,4,5), quantile, 0.84)

