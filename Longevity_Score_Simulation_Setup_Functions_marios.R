#Function for generating complete data
library(EnvStats)
library(flexmix)
library(extraDistr)
library(spatstat)
library(numDeriv)
library(emdbook)
library(ReIns)


### determine the right variance##############
# gamma = rgamma(10000,1/0.1,1/0.1)
# u = rnorm(10000,0,0.3) # random effect 
# s = exp(u) # standard deviation
# plot(density(gamma), col = "red", ylim = c(0,2), xlim = c(0,2), main = "Selected Standard Deviation")
# lines(density(s), col = "blue")
# legend("topright", legend = c("Gamma Distribution", "exp(u)"), col = c("red", "blue"), lty = 1)
# 


#n = vector with number of families per family size e.g. c(100,100)
#kvec = vector with possible family sizes e.g. c(4,8)
#sigma = specification for standard deviation of random-effect

sigma = 0.3

Percentile.data_u <- function(n, kvec, sigma){
  alpha = 0.0690573905
  beta = 0.0004263967
  u<-rnorm(100000,mean=0,sd=sigma)
  t = rgompertz(100000, alpha*exp(u), beta)
  pt<-quantile(t,prob=seq(0.01,1,by=0.01))
  k.list<-lapply(1:length(n),function(i)rep(kvec[i],n[i]))
  k<-c(k.list[[1]])
  for (i in 2:length(n)){
    k<-c(k,k.list[[i]])
  }
  # print(k)
  k<-sample(k,length(k))
  # print(table(k))
  nn<-sum(table(k)*sort(unique(k)))
  # print(nn)
  u<- rep(rnorm(length(k), mean = 0, sd = sigma),k) # random effects
  y = numeric()
  perc = numeric()
  for (i in 1:nn){
    y[i] = rgompertz(1, alpha*exp(u[i]), beta)
  }
  perc = ecdf(y)(y)
  # perc<-which(min(abs(y[i]-pt))==abs(y[i]-pt))/100
  idi<-rep(1:nn)
  id<-rep(1:length(k),k)
  kk<-rep(k,k)
  mat.random<-data.frame(idi=idi, id=id, y=y, perc = perc, kk=kk, u=u)
  mat.random
}
# head(Percentile.data_u(c(100,100),c(8,8), 0.3))

######## sigma = 0.3


##### Plot for the distibution of the percentiles ############
# perc_sigma_0.3 = as.data.frame(Percentile.data_u(c(250,250),c(8,8), 0.3))
# # plot(perc_sigma_0.1$u,perc_sigma_0.1$perc)
# u_q1<-quantile(perc_sigma_0.3$u)[2]
# u_q2<-quantile(perc_sigma_0.3$u)[3]
# u_q3<-quantile(perc_sigma_0.3$u)[4]
# plot(density(perc_sigma_0.3$perc[perc_sigma_0.3$u <=u_q1]), col = "darkgreen",main = "Percentile Distribution",xlim = c(-0.3,1.2))
# lines(density(perc_sigma_0.3$perc[perc_sigma_0.3$u>u_q1 & perc_sigma_0.3$u<u_q2]), col = "blue")
# lines(density(perc_sigma_0.3$perc[perc_sigma_0.3$u>u_q2 & perc_sigma_0.3$u<u_q3]), col = "cyan")
# lines(density(perc_sigma_0.3$perc[perc_sigma_0.3$u>u_q3]), col = "green")
# legend("topleft", c("1st quartile of u", "2nd quartile of u", "3rd quartile of u", "4th quartile of u"), lty=1, cex=0.8, col=c("darkgreen", "blue","cyan", "green"))
# 

# idi: Individual ID
# id: Family ID
# y: random gompertz value
# perc: survival percentile
# u: random effect
# kk: number of family members

###### bimodal residuals for Beta distribution

Percentile.data_umix<-function(n, kvec, beta0 ,phi, sigma){
  k.list<-lapply(1:length(n),function(i)rep(kvec[i],n[i]))
  k<-c(k.list[[1]])
  for (i in 2:length(n)){
    k<-c(k,k.list[[i]])
  }
  
  k<-sample(k,length(k))
  nn<-sum(table(k)*sort(unique(k)))
  u<-rep(rnormMix(length(k), mean1 = 0, sd1 = 0.3, mean2 = 1, sd2 = sigma, p.mix = 0.7),k) # change random normal values using a mixture
  # p.mix defines the strength you want to give in the second part of the mixture
  mu = exp(beta0 + u)/(1+exp(beta0 + u))
  perc = numeric()
  for (i in 1:length(mu)){
    perc[i] = rbeta(1, mu[i] * exp(phi), (1-mu[i]) * exp(phi))
  }
  
  idi<-rep(1:nn)
  id<-rep(1:length(k),k)
  kk<-rep(k,k)
  mat.random<-data.frame(idi=idi,id=id,perc=perc,kk=kk,u=u)
  mat.random
}

#Function for generating data with right-censoring

Percentile.data.cens<-function(n,kvec,sigma, censoring = F){
  alpha = 0.0004263967
  beta = 0.0690573905
  u<-rnorm(100000,mean=0,sd=sigma)
  t = rgompertz(100000, alpha*exp(u), beta)
  pt<-quantile(t,prob=seq(0.01,1,by=0.01))
  
  k.list<-lapply(1:length(n),function(i)rep(kvec[i],n[i]))
  k<-c(k.list[[1]])
  for (i in 2:length(n)){
    k<-c(k,k.list[[i]])
  }
  
  k <- sample(k,length(k))
  nn <- sum(table(k)*sort(unique(k)))
  u <- rep(rnorm(length(k), mean = 0, sd = sigma),k) ### change random normal values using a mixture
 
  #t are the ages
  fun_ecdf<-ecdf(t) #Survival function for t
  my_ecdf<-fun_ecdf(t)
  Ft<-data.frame(t,my_ecdf)
  # plot(Ft$t,1-Ft$my_ecdf)
  
  surv_prob = numeric()
  pred_age = numeric()
  # surv_cond = numeric()
  perc = numeric()
  y = numeric()
  c = numeric()
  d = numeric()
  for(i in 1:100){
    index<-which(round(Ft$t,0)==i)
    surv_prob[i]<-1-mean(unique(round(Ft[index,2],1)))
    surv_cond<-(1-Ft[,2])/surv_prob[i] #only meaningful from 60
    
    ind <- which(round(surv_cond,1) == 0.5)
    pred_age[i] = round(mean(round(Ft$t[ind],0)))
    pred_age[i] = ifelse(is.nan(pred_age[i]), 97, pred_age[i])
  }

  for (i in 1:nn){
    y[i] = rgompertz(1, alpha*exp(u[i]), beta) # random ages
    if (censoring == F){  # No censoring
      c[i] <- 100
    }else{
      c[i] <- runif(1, 40,100)  # assume uniform max at 110.
    }
    d[i] <- ifelse(y[i] < c[i], 1, 0)   # indicator for censoring: 1 if survival time observed; 0 if censored (survival time remains unobserved)
    y[i] <- ifelse(d[i] == 0, pred_age[round(y[i],0)], min(y[i], c[i])) # survival time is minimum of the survival and censoring time; survival time is not observed if censoring time is earlier
    # y[i] = min(y[i], c[i])
    #Transform age at death/censoring to percentile: go to the lifetable py previously calculated and check the closest percentile to age y[i]
    # perc[i]<-which(min(abs(y[i]-pt))==abs(y[i]-pt))/100
  }
  perc <- ecdf(y)(y)

  idi<-rep(1:nn)
  id<-rep(1:length(k),k)
  kk<-rep(k,k)
  # mat.random<-data.frame(idi=idi,id=id,y=y,perc = perc,cens=cens,death=death,kk=kk,u=u)
  mat.random<-data.frame(idi=idi,id=id,y=y,perc = perc, d=d, kk=kk,u=u)
  mat.random
}
prop = Percentile.data.cens(c(200,200),c(2,2),0.3,censoring = T)
# head(prop)
# hist(prop$y)
# table(prop$d)


##### Data for beta regression scores with censoring
Percentile.data.cens.beta<-function(n,kvec,sigma, censoring = F){
  alpha = 0.0004263967
  beta = 0.0690573905
  u<-rnorm(100000,mean=0,sd=sigma)
  t = rgompertz(100000, alpha*exp(u), beta)
  pt<-quantile(t,prob=seq(0.01,1,by=0.01))
  
  k.list<-lapply(1:length(n),function(i)rep(kvec[i],n[i]))
  k<-c(k.list[[1]])
  for (i in 2:length(n)){
    k<-c(k,k.list[[i]])
  }
  
  k <- sample(k,length(k))
  nn <- sum(table(k)*sort(unique(k)))
  u <- rep(rnorm(length(k), mean = 0, sd = sigma),k) ### change random normal values using a mixture
  perc = numeric()
  y = numeric()
  c = numeric()
  d = numeric()
  for (i in 1:nn){
    y[i] = rgompertz(1, alpha*exp(u[i]), beta) # random ages
    if (censoring == F){  # No censoring
      c[i] <- 100
    }else{
      c[i] <- runif(1, 40,100)  # assume uniform max at 110.
    }
    d[i] <- ifelse(y[i] < c[i], 1, 0)   # indicator for censoring: 1 if survival time observed; 0 if censored (survival time remains unobserved)
    y[i] <- min(y[i], c[i]) # survival time is minimum of the survival and censoring time; survival time is not observed if censoring time is earlier
    #Transform age at death/censoring to percentile: go to the lifetable py previously calculated and check the closest percentile to age y[i]
    # perc[i]<-which(min(abs(y[i]-pt))==abs(y[i]-pt))/100
  }
  perc <- ecdf(y)(y)
  
  #plots for specific cohort given the age but need to use the above for loop 
  # plot(Ft$t,Surv_cond60)
  # abline(h=1)
  # abline(h=0.5)
  
  idi<-rep(1:nn)
  id<-rep(1:length(k),k)
  kk<-rep(k,k)
  # mat.random<-data.frame(idi=idi,id=id,y=y,perc = perc,cens=cens,death=death,kk=kk,u=u)
  mat.random<-data.frame(idi=idi,id=id,y=y,perc = perc, d=d, kk=kk,u=u)
  mat.random
}

prop_beta = Percentile.data.cens.beta(c(200,200),c(2,2),0.3,censoring = T)
# head(prop_beta)


# idi: Individual ID
# id: Family ID
# y: random gompertz value
# perc: Observed Survival Percentile
# u: random effect
# kk: number of family members
# cens: survival percentile at right-censoring point
#death: death indicator variable


#####Functions for calculation of mLRC

dbin.person = function(u, y, beta){
  p = exp(beta+u)/(1+exp(beta+u)) 
  fy = p^y * (1-p)^(1-y)
  prod(fy)
}

### example
# dbin.person(1.0559919, 0.8517372, 1)


dbeta.mixed = function(par,Y,X){
  #Fixed dispersion parameter
  #marginal likelihood per person
  loglik = c()
  N=length(Y)
  for(i in 1:N)
    loglik[i] = log(gauss.hermite(dbin.person, mu=0, sd=exp(par[2]),  y=Y[X==i], beta=par[1]))
  -sum(loglik)
}



# Y = as.numeric(data$perc > TH_sig_expmin1)
# X = data$id
# par=c(1,1)
# dbeta.mixed(par=par,X=X,Y=Y)

###########


predict.fam = function(par, famnr, nsamples){
  U = rnorm(nsamples, 0, exp(par[2]))
  quantiles = sapply(U, function(u){
    exp(par[1]+u)/(1+exp(par[1]+u)) 
  })
  
  quantiles.ordered = quantiles[order(quantiles)]
  U.ordered = U[order(quantiles)]
  
  y=Y[X==famnr]
  
  # print(table(y))
  dens.y = sapply(U.ordered, function(u) dbin.person(u, y=y, beta=par[1]))
  dens.y = dens.y / sum(dens.y)
  
  sum(quantiles.ordered * dens.y)
}

par = c(1,1)

#Functions for calculation of Beta Regression Scores
#######################################################
# Functions
#######################################################
#Percentile: Y
#Death status: Death
#B: Start follow-up
#Family ID: Z
#Covariates: X
################################################
###################################

#### pos: c(beta, phi, s)

loglik.sum = function(par, pos, Y, B, Death, Z){
  Zunique = unique(Z)
  Q = sapply(Zunique, function(q)
    log(gauss.hermite(likcluster, mu=0, sd=exp(par[pos$s]),  Y=Y[Z==q], B=B[Z==q],  Death=Death[Z==q], par=par, pos=pos))
  )
  -sum(Q)
}


betareg.mixed = function(Y,Death,B,Z){
  pos = list(beta=1, phi = 2, s = 3)
  par = rep(0, pos$s)
  par[pos$phi] = 1 
  
  est = nlminb(par, loglik.sum, Y=Y, Death=Death, B=B, Z=Z, pos=pos, control=list(trace=TRUE))$par
  H = hessian(loglik.sum, est, pos=pos, Y=Y, B=B,Death=Death, Z=Z)
  SE = diag(solve(H))
  
  #beta: effects covariates
  #phi: precision parameter
  #s: st. dev. random effects
  list(beta=est[pos$beta], phi=est[pos$phi], s=est[pos$s], pos=pos, estimates=data.frame(est=est, se=SE))
} 



likcluster = function(ranef, par, pos, Y,B, Death, prod=TRUE){
  xb = par[pos$beta] + ranef
  mu = exp(xb)/(1+exp(xb))
  phi = exp(par[pos$phi])
  
  alpha = mu*phi
  beta = (1-mu)*phi
  
  #If you are still alive, we require P(Y>y|p,q); this can be obtained from the 1-cum density function of the beta distribution
  if(prod){
    prod(
      pbeta(Y, alpha, beta)^B *
        (dbeta(Y, alpha, beta)^Death * pbeta(Y, alpha, beta, lower.tail=FALSE)^(1-Death))^(1-B))
    
  }else{
    pbeta(Y, alpha, beta)^B *
      (dbeta(Y, alpha, beta)^Death * pbeta(Y, alpha, beta, lower.tail=FALSE)^(1-Death))^(1-B)  
  }
}


#Predict for a certain family
predict.ranef = function(estimates, Yclus, Deathclus, Bclus, nsamples){
  U = rnorm(nsamples, 0, exp(estimates$s))
  
  #Estimate their conditional likelihood (given their random effect)
  dens.y = sapply(U, function(u){
    likcluster(ranef=u, par=c(estimates$beta, estimates$phi, estimates$s), pos=estimates$pos, Y=Yclus,Death=Deathclus, B=Bclus)
  })
  
  dens.y = dens.y / sum(dens.y)
  sum(U * dens.y)
}


#Process for determining Thresholds when sigma = 0.2 + TH_expmin1

#Specify parameters for simulating a large dataset with
set.seed(1234)
Sig = 0.3
kvec=c(4,4)
n=c(250,250)
sigma=0.3

#Commented code creates a large life-table from which imputed expected survival
#percentiles can be drawn. This is necessary when dealing with right-censored
#observations in the LRC and mLRC scores

data_sig_expmin1<-as.data.frame(Percentile.data_u(n,kvec,sigma))
# data_sig_expmin1
rml_sig_expmin1<-c()
perc_sort_sig_expmin1<-sort(data_sig_expmin1$perc)
for (i in 1: length(perc_sort_sig_expmin1)){
  rml_sig_expmin1[i]<-mean(perc_sort_sig_expmin1[i:length(perc_sort_sig_expmin1)])
}
lifetable_sig_expmin1<-data.frame(perc=perc_sort_sig_expmin1,expectedy=rml_sig_expmin1)
TH_sig_expmin1 <- quantile(lifetable_sig_expmin1$perc,0.9)
TH_sig_expmin1
save(lifetable_sig_expmin1, file = 'lifetable_sig_expmin1.rdata')
#CORRECT BETA_TH Threshold
# y at 90th percentile becomes threshold










