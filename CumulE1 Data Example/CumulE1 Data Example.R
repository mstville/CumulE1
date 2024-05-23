#################################################################
#################################################################
## Cumulative-E1 Data Example: 
## Cumulative effect of third-trimester GWG on 37-week EFW
#################################################################
#################################################################

library(mvtnorm)
library(msm)
library(stats)
library(stringr)
library(missMDA)

########################################
# Install R package from Github
library(devtools)
devtools::install_github('mstville/CumulE1')
library(CumulE1)

gest.ages=27:37
m=length(gest.ages)

data=read.csv("EFW37_Weight_Dataset.csv",header=T)  ## Set working directory appropriately to access data

head(data)
summary(data)
data$Education=as.factor(data$Education);summary(data$Education)
data$Race=as.factor(data$Race);summary(data$Race)
data$InfSex=as.factor(data$InfSex);summary(data$InfSex)
data$Parity=as.factor(data$Parity);summary(data$Parity)
data$GDM=as.factor(data$GDM);summary(data$GDM)
data$HTN=as.factor(data$HTN);summary(data$HTN)


pats=unique(data$ID)
n=length(pats);n

################################################
## Create GWG Exposure matrix
w.id=which(str_detect(names(data),"Weight")==TRUE)
W.dat=data[,w.id]

Z=matrix(-99,nrow=n,ncol=m)
for (t in 1:m){
  id=which(names(W.dat)==paste0("Weight",gest.ages[t]))
  Wt=W.dat[,c(id-1,id)]
  Zt=apply(Wt,1,diff)
  Z[,t]=Zt
}
colnames(Z)=paste0("GWG",gest.ages)
head(Z)
data=cbind.data.frame(data,as.data.frame(Z));head(data)

################################################
## Create covariate matrix
X=model.matrix(~Race+Age+Education+InfSex+Cotinine+GDM+HTN+Parity+PreBMI,data=data)
colnames(X)[c(1,15)]=c("Intercept","Parity2")
p=ncol(X);p

################################################
## Outcome: 37-week EFW
Y=data$EFW  


#-------------------------------------------------------------------------------------
###############################################################################################
## MODEL ESTIMATION
set.seed(221)

G=5000 ## MCMC samples

# initial values
a_e=b_e=0.01
sig2_e=1

sig2_b=10000
beta=rep(0,p)

sig2_th=10000
theta=0       ### CUMULATIVE EFFECT PARAMETER

gamma=rep(1,m)

eta=rep(0,m)

a_phi=b_phi=1
phi=-log(0.05)/(m-1)
temp_corr_info=temporal_corr_fun(m,phi)

alpha=theta*gamma


sig2.e.mcmc=rep(-99,G)
beta.mcmc=matrix(-99,nrow=G,ncol=p)
theta.mcmc=rep(-99,G)
gamma.mcmc=matrix(-99,nrow=G,ncol=m)
eta.mcmc=matrix(-99,nrow=G,ncol=m)
phi.mcmc=rep(-99,G)
alpha.mcmc=matrix(-99,nrow=G,ncol=m)

acctot_phi_trans = 0

#####################################
### GIBBS SAMPLER STARTS HERE
for (g in 1:G){
  ###########################
  # Update residual variance
  sig2_e = sig2_epsilon_update(a_e,b_e,Y,X,Z,beta,alpha)
  sig2.e.mcmc[g]=sig2_e
  
  w=rep(1/sig2_e, n)
  
  ###########################
  # Update beta
  beta = beta_update(X,Z,sig2_b,w,Y,alpha)
  beta=as.vector(beta)
  beta.mcmc[g,]=beta
  
  ###########################
  # Update theta
  theta = theta_update(Y,X,Z,beta,gamma,sig2_th,sig2_e)
  theta.mcmc[g]=theta
  
  ###########################
  # Update gamma and alpha
  gamma=gamma_update(Y,X,Z,beta,w,gamma,alpha,theta,eta)
  gamma=as.vector(gamma)
  for (t in 1:m){alpha[t]=theta*gamma[t]}
  gamma.mcmc[g,]=gamma
  alpha.mcmc[g,]=alpha
  
  ###########################
  # Sample latent variables
  # from truncated normal
  gamma_star = gamma_star_update(gamma,eta)
  gamma_star=as.vector(gamma_star)
  
  ###########################
  # Update eta
  eta = eta_update(gamma_star,corr_inv = temp_corr_info$temporal_corr_inv)
  eta=as.vector(eta)
  eta.mcmc[g,]=eta
  
  ###########################
  # Update phi
  phi.up = phi_update(phi,eta,temp_corr_info,a_phi,b_phi,metrop_var_phi_trans = 1,acctot_phi_trans)
  phi = phi.up$phi
  acctot_phi_trans = phi.up$acctot_phi_trans
  temp_corr_info = phi.up$temp_corr_info
  
  phi.mcmc[g]=phi
  
  print(g)
}
### GIBBS SAMPLER ENDS HERE
#####################################

burn=2500:G

gamma=apply(gamma.mcmc[burn,],2,mean);gamma

theta=mean(theta.mcmc[burn]);theta
theta.ci=quantile(theta.mcmc[burn],probs=c(.025,.975));theta.ci
