source("auxilliary_functions.R")
source("ng_SAVS.R")

sim_data<-function(K,t,sparsity=0.3) #simulate data
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
  index_sparsity<-sample(1:(1*K),size=round(sparsity*K,digits = 0),replace=FALSE)
  beta[,index_sparsity]<-0
  beta_tilda=t(sapply(1:t, function(i)(beta[i,]-beta_0)/abs(v)))
  alpha<-c(beta_0,abs(v))
  yt<-X %*% alpha[1:K] + (beta_tilda*X) %*% alpha[(K+1):(2*K)] + rnorm(t,0,sqrt(sigma_error))
  
  return(list("yt"=yt,"X"=X,"true_beta"=beta,"index_sparsity"=index_sparsity))
}

prior<-"NG"

l=1
j=1
mean_mae_sparse<-matrix(NA,nrow=3,ncol=3) #store mean absolute error
mean_mae_nonsparse<-matrix(NA,nrow=3,ncol=3) #store mean absolute error for non-sparse estimates
mean_hit_rate<-matrix(NA,nrow=3,ncol=3) #find hit rates
for (K in c(5,15,30)) 
{ 
  j=1
  for(sparsity in c(0.3,0.7,0.9))
  {
    mae_sparse<-c()
    mae_nonsparse<-c()
    hit_rate<-c()
    for (i in 1:100) 
    {
      data_simulated<-sim_data(K,250,sparsity)
      
      output_sim<-sparse.SAVS(Yraw=as.matrix(data_simulated$yt),Xraw=as.matrix(data_simulated$X),nsave=10000,nburn=1000,thin=1,shrinkage =prior,cons=0)
      
      mae_sparse[i]<-mean(sapply(1:K, function(i) abs(apply(output_sim$A.thrsh[,,i,][seq(1,10000,by=5),],2,quantile,0.5,na.rm=TRUE) - data_simulated$true_beta[,i])))
      
      mae_nonsparse[i]<-mean(sapply(1:K, function(i) abs(apply(output_sim$A.non[,,i,][seq(1,10000,by=5),],2,quantile,0.5,na.rm=TRUE) - data_simulated$true_beta[,i])))
      
      hit_rate[i]<-(sum(sapply(1:length(data_simulated$index_sparsity),function(i) mean(round(apply(output_sim$A.thrsh[,,data_simulated$index_sparsity[i],],2,quantile,0.5),2)) == 0))/length(data_simulated$index_sparsity))*100
    }
    mean_mae_sparse[l,j]<-mean(mae_sparse)
    mean_mae_nonsparse[l,j]<-mean(mae_nonsparse)
    mean_hit_rate[l,j]<-mean(hit_rate)
    j=j+1
  }
  l=l+1
}

mega_output_NG<- list(mean_mae_nonsparse,mean_mae_sparse,mean_hit_rate)

save(mega_output_NG,file="NG.RData")