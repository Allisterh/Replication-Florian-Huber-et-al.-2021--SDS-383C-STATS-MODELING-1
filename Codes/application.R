library(bvarsv)
library(coda)
library(shrinkTVP)
library(parallel)

source("auxilliary_functions.R")
source("ng_SAVS.R")

data("usmacro.update")

lags <- usmacro.update[1:(nrow(usmacro.update) - 1), ]
colnames(lags) <- paste0(colnames(lags), "_lag")
us_data <- data.frame(inf = usmacro.update[2:nrow(usmacro.update), "inf"], + lags)
View(us_data)


us_res <- shrinkTVP(inf ~ inf_lag + une_lag + tbi_lag, us_data[1:236,], niter = 60000, nburn = 10000, nthin = 10,display_progress = TRUE)


summary(us_res, showprior = FALSE)
par(mar=c(2,4,2,2))
par(mfcol=c(4,1))
plot(us_res$beta$beta_Intercept,probs = c(0.025,0.975),ylab=expression(paste(beta,"_Intercept")))
plot(us_res$beta$beta_inf_lag,probs = c(0.025,0.975),ylab=expression(paste(beta,"_Inflation Lag")))
plot(us_res$beta$beta_une_lag,probs = c(0.025,0.975),ylab=expression(paste(beta,"_Unemployment Lag")))
plot(us_res$beta$beta_tbi_lag,probs = c(0.025,0.975),ylab=expression(paste(beta,"_Treasury Bill Lag")))




## Simulation Graph
sim_data_1<-function(K,t,sparsity=0.3)
{
  sigma_error<-0.1^2 
  X<-sapply(1:K,function(i) runif(t,-1,1))
  mean_beta0<-rep(0,K)
  var_beta0<-diag(0.1^2,nrow=K,ncol=K)
  beta_0<-mvtnorm::rmvnorm(1,mean = mean_beta0,var_beta0)
  v<-rnorm(K,0,0.1)
  beta<-matrix(NA,ncol=K,nrow=t)
  beta_tilda<-matrix(NA,ncol=K,nrow=t)
  beta[1,]<-beta_0
  for (i in 2:t) {
    beta[i,]<-beta[(i-1),] + mvtnorm::rmvnorm(1,mean=rep(0,K),diag(abs(v^2),nrow=K,ncol=K))
  }
  index_sparsity<-sample(1:(1*(K-1)),size=round(sparsity*K,digits = 0),replace=FALSE)
  beta[,index_sparsity]<-0
  beta[,K]<- 0.8
  beta_tilda=t(sapply(1:t, function(i)(beta[i,]-beta_0)/abs(v)))
  alpha<-c(beta_0,abs(v))
  yt<-X %*% alpha[1:K] + (beta_tilda*X) %*% alpha[(K+1):(2*K)] + rnorm(t,0,sqrt(sigma_error))
  
  return(list("yt"=yt,"X"=X,"true_beta"=beta,"index_sparsity"=index_sparsity))
}

data_simulated<-sim_data_1(15,250,sparsity=0.7)

sim_res <- shrinkTVP(data_simulated$yt ~ data_simulated$X, niter = 5000, nburn = 1000, nthin = 10,display_progress = TRUE)

library(latex2exp)
par(mar=c(2,4,2,2))
par(mfrow=c(1,3))
plot(sim_res$beta$`beta_data_simulated$X1`,probs = c(0.025,0.975),main=TeX(r'($\beta_{jt} = 0$)'),ylab=TeX(r'($\beta_{jt} $)'))
plot(sim_res$beta$`beta_data_simulated$X20`,probs = c(0.025,0.975),main=TeX(r'($\beta_{jt} = 0.8$)'),ylab=TeX(r'($\beta_{jt} $)'))
plot(sim_res$beta$`beta_data_simulated$X9`,probs = c(0.025,0.975),main=TeX(r'($\beta_{jt}$ time-varying)'),ylab=TeX(r'($\beta_{jt} $)'))


