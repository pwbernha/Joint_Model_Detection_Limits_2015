##############Application Code#####################


####################Necessary Packages######################
library(mnormt)	 	#generate/evaluate (multivariate) normal distributions
library(corpcor)		#check positive definiteness of matrices
library(lme4)		#fit linear mixed models
library(MCMCpack)		#generate samples for calculating expectations using MCMC
library(lmec)		#fit linear mixed models with censored responses


#############Reading in the Data (and getting in appropriate form)#############
#Reading in  GenIMS longitudinal data
GenIMS <- read.csv("C://Users/Paul/Documents/Data.csv")

#Reading in GenIMS single time point data
GenIMSInd <- read.csv("C://Users/Paul/Documents/Data2.csv")

#Renumbering patient ids in long. data from 1:2320
t <- 1
s <- GenIMS[1,1]
for(i in 1:length(GenIMS[,1])){
  if(GenIMS[i,1]!=s) {s <- GenIMS[i,1]
  t <- t +1
  GenIMS[i,1] <- t}
  if(GenIMS[i,1]==s) GenIMS[i,1] <- t
}

#Renumbering pateint ids in cross. data to match long. data
GenIMSInd[,1] <- 1:2320



#Removing patients who did not actually have CAP & who were immediately discharged
remove <- which(GenIMSInd[,68]==0)
remove2 <- GenIMS[which(GenIMS[,66]=="OUTPATIENT"),1]
for(i in remove) GenIMS <- GenIMS[-which(GenIMS[,1]==i),]
GenIMS <- GenIMS[-which(GenIMS[,66]=="OUTPATIENT"),]

#removing those without all three biomarkers observed
for(i in c(451,595,924,932,1188,1312,1402,2109,2152,2207)) GenIMS <- GenIMS[-which(GenIMS[,1]==i),] 

#removing sames patients from cross. data
GenIMSInd <- GenIMSInd[-unique(remove2),]
GenIMSInd <- GenIMSInd[-which(GenIMSInd[,68]==0),]
for(i in c(451,595,924,932,1188,1312,1402,2109,2152,2207)) GenIMSInd <- GenIMSInd[-which(GenIMSInd[,1]==i),]
for(i in c(142,468,844,1249,1333,1718,1757,2057,2088,2191)) GenIMSInd <- GenIMSInd[-which(GenIMSInd[,1]==i),]

#Keeping only columns of interest for analysis
GenIMS <- GenIMS[,c(1,2,22,23,24)]

#renumbering patient ids from 1:1875 for long. data
t <- 1
s <- GenIMS[1,1]
for(i in 1:length(GenIMS[,1])){
  if(GenIMS[i,1]!=s) {s <-GenIMS[i,1]
  t <- t +1
  GenIMS[i,1] <- t}
  if(GenIMS[i,1]==s) GenIMS[i,1] <- t
}

#renumbering for cross. data
GenIMSInd[,1] <- 1:1875 

#Selecting TNF, IL6, and IL10 variables
X1 <- na.omit(GenIMS[,c(1,2,3)]) #TNF
X2 <- na.omit(GenIMS[,c(1,2,4)]) #IL6
X3 <- na.omit(GenIMS[,c(1,2,5)]) #IL10

#adding in a censoring variable, initially set to 1
X1 <- cbind(X1, 1)
X2 <- cbind(X2, 1)
X3 <- cbind(X3, 1)

#setting censoring variable to 0 when DL observed
X1[which(X1[,3]==-4),4] <- 0
X2[which(X2[,3]==-2 | X2[,3]==-5 ),4] <- 0
X3[which(X3[,3]==-5),4] <- 0

#taking log transformation of data
X1[,3] <- log(abs(X1[,3]))
X2[,3] <- log(abs(X2[,3]))
X3[,3] <- log(abs(X3[,3]))

#yet again, renumbering id, this time by covariates (X1 here)
t <- 1
s <- X1[1,1]
for(i in 1:length(X1[,1])){
  if(X1[i,1]!=s) {s <-X1[i,1]
  t <- t +1
  X1[i,1] <- t}
  if(X1[i,1]==s) X1[i,1] <- t
}

#renumbering for X2
t <- 1
s <- X2[1,1]
for(i in 1:length(X2[,1])){
  if(X2[i,1]!=s) {s <-X2[i,1]
  t <- t +1
  X2[i,1] <- t}
  if(X2[i,1]==s) X2[i,1] <- t
}

#renumbering for X3
t <- 1
s <- X3[1,1]
for(i in 1:length(X3[,1])){
  if(X3[i,1]!=s) {s <-X3[i,1]
  t <- t +1
  X3[i,1] <- t}
  if(X3[i,1]==s) X3[i,1] <- t
}

#Defining initial observation as baseline time (0) instead of 1
X1[,2] <- X1[,2]-1
X2[,2] <- X2[,2]-1
X3[,2] <- X3[,2]-1

#setting race for non-whites to  0
GenIMSInd[which(GenIMSInd[,5]>1),5] <- 0

#creating matrix of baseline covariates
Z <- cbind(1, GenIMSInd[,4]-1,GenIMSInd[,5],GenIMSInd[,69])

#defining response as 0 for those individuals who died and 1 for those who survived 90 days
Y <- as.numeric(GenIMSInd[,18])
for(i in 1:1875){
  if(Y[i]==90) Y[i] <-1 else Y[i] <- 0
}

#sample size
n <- 1875

#Gaussian quadrature points
u <-c(1.4806867,1.1904181,1.0813192,1.0328036,1.0185664,1.0328036,1.0813192,1.1904181,1.4806867)
v <- c(4.5127459, 3.205429, 2.076848, 1.0232557,0,-1.023256, -2.076848,-3.205429,-4.512746)


#Getting Initial Values for EM algorithm
X21 <- cbind(1,X1[,2])
X22 <- cbind(1,X2[,2])
X23 <- cbind(1,X3[,2])

Dhat <- matrix(0,6,6)
Bhat <- rep(0,6)
tauhat2 <- rep(0,3)
b1s <- b2s <- b3s <- matrix(0,2,1875)
vals1 <- lmec(yL=X1[,3],cens=c(1-X1[,4]),X=X21,Z=X21,cluster=X1[,1],epsstop=0.1)
Bhat[1:2] <- vals1$beta
tauhat2[1] <- vals1$sigma^2
Dhat[1:2,1:2] <- vals1$Psi
b1s <- vals1$bi
vals2 <- lmec(yL=X2[,3][1:5192],cens=c(1-X2[,4][1:5192]),X=X22[1:5192,],Z=X22[1:5192,],cluster=X2[,1][1:5192],epsstop=50)
vals2b <- lmec(yL=X2[,3][5193:9081],cens=c(1-X2[,4][5193:9081]),X=X22[5193:9081,],Z=X22[5193:9081,],cluster=(X2[,1][5193:9081]-1023),epsstop=50)
Bhat[3:4] <- vals2$beta*5192/9081+ vals2b$beta*3889/9081
tauhat2[2] <- vals2$sigma^2*5192/9081+vals2b$sigma^2*3889/9081
Dhat[3:4,3:4] <- vals2$Psi*5192/9081+vals2b$Psi*3889/9081
b2s <- cbind(vals2$bi,vals2b$bi)
vals3 <- lmec(yL=X3[,3][1:5194],cens=c(1-X3[,4][1:5194]),X=X23[1:5194,],Z=X23[1:5194,],cluster=X3[,1][1:5194],epsstop=50)
vals3b <- lmec(yL=X3[,3][5195:9079],cens=c(1-X3[,4][5195:9079]),X=X23[5195:9079,],Z=X23[5195:9079,],cluster=(X3[,1][5195:9079]-1023),epsstop=50)
Bhat[5:6] <- vals3$beta*5194/9079+vals3b$beta*3885/9079
tauhat2[3] <- vals3$sigma^2*5194/9079+vals3b$sigma^2*3885/9079
Dhat[5:6,5:6] <- vals3$Psi*5194/9079+	vals3b$Psi*3885/9079
b3s <- cbind(vals3$bi,vals3b$bi)
coefs <- coef(brglm(Y ~ Z[,2]+Z[,3]+Z[,4]+b1s[1,]+b1s[2,]+b2s[1,]+b2s[2,]+b3s[1,]+b3s[2,],family="binomial"))

sumt <- 0
for(i in 1:n){
  Xi1 <- X1[X1[,1]==i,]
  Xi2 <- X2[X2[,1]==i,]
  Xi3 <- X3[X3[,1]==i,]
  if((sum(Xi1[,4])==0 & sum(Xi2[,4])==0) | (sum(Xi1[,4])==0 & sum(Xi3[,4])==0) | (sum(Xi2[,4])==0 & sum(Xi3[,4])==0)) sumt <- sumt+1
}

Pars <- list(coefs[1:4],coefs[5:10],tauhat2,Bhat,Dhat)

#Running EM Algorithm Iterations
ParsMatrix <- rbind(rep(1,55),c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])))
Iter <- 1
while(max(abs((ParsMatrix[(Iter+1),1:10]-ParsMatrix[(Iter),1:10])/ParsMatrix[(Iter),1:10]))>0.05 & Iter < 50 & STOP==FALSE){
  print(MCMC)
  Pars <- EM(Y,X1,X2,X3,Z,Pars[[1]],Pars[[2]],Pars[[3]],Pars[[4]],Pars[[5]],10000)
  ParsMatrix <- rbind(ParsMatrix,c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])))
  Iter <- Iter + 1
  if(Iter[k]>=10) if(max(abs((colMeans(ParsMatrix[(Iter[k]-3):(Iter[k]+1),])-colMeans(ParsMatrix[(Iter[k]-8):(Iter[k]-4),]))/colMeans(ParsMatrix[(Iter[k]-8):(Iter[k]-4),])))<epsilon) STOP <- TRUE
}

if(STOP==FALSE) Final[k,] <- c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])) else Final[k,] <- colMeans(ParsMatrix[(Iter[k]-3):(Iter[k]+1),])		
Variances <- VarEM(Y,X1,X2,X3,Z,as.vector(Final[k,1:4]),as.vector(Final[k,5:10]),as.vector(Final[k,11:13]),as.vector(Final[k,14:19]),matrix(Final[k,20:56],6,6),2000,3000)
SE[k,] <- sqrt(diag(Val[[1]]))	#Variance estimator based on information matrix
SE2[k,] <- sqrt(diag(Val[[2]]))	#Sandwich variance estimator




#######################Functions#########################

#Function called in an nlm or optim (maximization function built in R) to find posterior mode and variance of f(b|data;thetahat)
Maxb <- function(Y,X12,X13,X14,X22,X23,X24,X32,X33,X34,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat){
  like <- function(bv){
    llike <- -log(dbinom(Y,1,(exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat)/(1+exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat))))*dmnorm(c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6]),Bhat,Dhat))-sum(log((X14*dnorm(X13,(bv[1] + bv[2]*X12),sqrt(tauhat2[1]))+(1-X14)*pnorm(X13,(bv[1] + bv[2]*X12),sqrt(tauhat2[1])))))-sum(log((X24*dnorm(X23,(bv[3] + bv[4]*X22),sqrt(tauhat2[2]))+(1-X24)*pnorm(X23,(bv[3] + bv[4]*X22),sqrt(tauhat2[2])))))-sum(log((X34*dnorm(X33,(bv[5] + bv[6]*X32),sqrt(tauhat2[3]))+(1-X34)*pnorm(X33,(bv[5] + bv[6]*X32),sqrt(tauhat2[3])))))	
    return(llike)}
  like
}

#Function (whose argument is proportional to the likelihood of f(b|data;thetahat)) called in MCMCmetrop1R for obtaining samples from f(b|data;thetahat)

Sampler <- function(bv,Y,X12,X13,X14,X22,X23,X24,X32,X33,X34,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat){
  log(dbinom(Y,1,(exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat)/(1+exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat))))*dmnorm(c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6]),Bhat,Dhat))+sum(log((X14*dnorm(X13,(bv[1] + bv[2]*X12),sqrt(tauhat2[1]))+(1-X14)*pnorm(X13,(bv[1] + bv[2]*X12),sqrt(tauhat2[1])))))+sum(log((X24*dnorm(X23,(bv[3] + bv[4]*X22),sqrt(tauhat2[2]))+(1-X24)*pnorm(X23,(bv[3] + bv[4]*X22),sqrt(tauhat2[2])))))+sum(log((X34*dnorm(X33,(bv[5] + bv[6]*X32),sqrt(tauhat2[3]))+(1-X34)*pnorm(X33,(bv[5] + bv[6]*X32),sqrt(tauhat2[3])))))	
}

#EM algorithm update function: takes data and current parameter estimates and updates the parameter estimates using one EM step
#(see simulation code for more details)
EM <- function(Y,X1,X2,X3,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat,n.MC){
  
  #Printing current parameter values
  print(c(betaAhat,betaBhat,tauhat2,Bhat,Dhat))
  
  #Initializing vectors and matrices
  Dhats <- dScoreBetas <- dScoretau2s <- NULL
  Bhati <-  ScoreBetaB <- matrix(0,n,6)
  ScoreBetaA <- matrix(0,n,4)
  Scoretau2 <- matrix(0,n,3)
  
  for(i in 1:n){
    print(i)
    
    #Defining X matrix for ith individual
    Xi1 <- X1[(X1[,1]==i),]
    Xi2 <- X2[(X2[,1]==i),]
    Xi3 <- X3[(X3[,1]==i),]
    Ni1 <- dim(Xi1)[1]
    Ni2 <- dim(Xi2)[1]
    Ni3 <- dim(Xi3)[1]
    
    #Initializing individual specific matrices
    Dhati <-  matrix(0,6,6)
    dScoretau2i <- matrix(0,3,3)
    dScoreBetai <- matrix(0,10,10)
    
    if(sum(X1[(X1[,1]==i),4],X2[(X2[,1]==i),4],X3[(X3[,1]==i),4])/(Ni1+Ni2+Ni3)>=0.2){
      
      #Finding mode and variance of random effects
      inits <- coef(lm(Xi1[,3] ~ Xi1[,2]))
      inits <- c(inits,coef(lm(Xi2[,3] ~ Xi2[,2]))) 
      inits <- c(inits,coef(lm(Xi3[,3] ~ Xi3[,2])))
      tryCatch(f <-nlm(Maxb(Y[i],Xi1[,2],Xi1[,3],Xi1[,4],Xi2[,2],Xi2[,3],Xi2[,4],Xi3[,2],Xi3[,3],Xi3[,4],Z[i,],betaAhat,betaBhat,tauhat2,Bhat,Dhat),Bhat,hessian=TRUE), error=function(...) f <<- optim(inits,Maxb(Y[i],Xi1[,2],Xi1[,3],Xi1[,4],Xi2[,2],Xi2[,3],Xi2[,4],Xi3[,2],Xi3[,3],Xi3[,4],Z[i,],betaAhat,betaBhat,tauhat2,Bhat,Dhat),method="BFGS",hessian=TRUE))
      bhat <- f$estimate 
      if(is.positive.definite(solve(f$hessian))==TRUE) sigmahat <- solve(f$hessian) else sigmahat <- make.positive.definite(solve(f$hessian))
      if(is.positive.definite(solve(f$hessian))==FALSE & is.positive.definite(-solve(f$hessian))==TRUE) sigmahat <- -solve(f$hessian)
      sigmahat.5 <- t(chol(sigmahat))
      
      #Defining some terms to make input easier
      mu <- betaAhat%*%Z[i,] + betaBhat%*%bhat
      vari <- betaBhat%*%sigmahat%*%betaBhat
      mu21 <- vari21 <- rep(0, Ni1)
      mu22 <- vari22 <- rep(0, Ni2)
      mu23 <- vari23 <- rep(0, Ni3)
      for(j in 1:Ni1) mu21[j] <- c(1,Xi1[j,2])%*%bhat[1:2]
      for(j in 1:Ni2) mu22[j] <- c(1,Xi2[j,2])%*%bhat[3:4]
      for(j in 1:Ni3) mu23[j] <- c(1,Xi3[j,2])%*%bhat[5:6]
      for(j in 1:Ni1) vari21[j] <- c(1,Xi1[j,2])%*%sigmahat[1:2,1:2]%*%c(1,Xi1[j,2])
      for(j in 1:Ni2) vari22[j] <- c(1,Xi2[j,2])%*%sigmahat[3:4,3:4]%*%c(1,Xi2[j,2])
      for(j in 1:Ni3) vari23[j] <- c(1,Xi3[j,2])%*%sigmahat[5:6,5:6]%*%c(1,Xi3[j,2])
      Pi <- exp(mu+sqrt(vari)*v)/(1+exp(mu+sqrt(vari)*v))
      
      #Calculating computings for random effects mean and matrix estimates
      for(j in 1:6) Bhati[i,j] <- sum(u/sqrt(2*pi)*exp(-v^2/2)*(bhat[j]+sqrt(sigmahat[j,j])*v))
      for(j in 1:6) Dhati[j,j] <- sum(u/sqrt(2*pi)*exp(-v^2/2)*(bhat[j]+sqrt(sigmahat[j,j])*v-Bhat[j])^2)
      for(p in 1:5) for(q in (p+1):6) for(j in 1:length(u)) for(k in 1:length(u)) Dhati[p,q] <- Dhati[p,q] + u[j]*u[k]*(bhat[p]+sigmahat.5[p,c(p,q)]%*%c(v[j],v[k])-Bhat[p])*(bhat[q]+sigmahat.5[q,c(p,q)]%*%c(v[j],v[k])-Bhat[q])/(2*pi)*exp(-v[j]^2/2-v[k]^2/2)
      for(p in 1:5) for(q in (p+1):6) Dhati[q,p] <- Dhati[p,q]
      Dhats <- c(Dhats,list(Dhati))
      
      #Finding Scores and derivatives of Scores for each 
      ScoreBetaA[i,1] <- Y[i]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi)
      ScoreBetaA[i,2] <- Y[i]*Z[i,2]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*Z[i,2])
      ScoreBetaA[i,3] <- Y[i]*Z[i,3]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*Z[i,3])
      ScoreBetaA[i,4] <- Y[i]*Z[i,4]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*Z[i,4])
      for(j in 1:6) ScoreBetaB[i,j] <- Y[i]*Bhati[i,j]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*(bhat[j]+v*betaBhat%*%sigmahat[,j]/(sqrt(vari))))
      dScoreBetai[1,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi-Pi^2))
      dScoreBetai[2,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]^2-Pi^2*Z[i,2]^2))
      dScoreBetai[3,3] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,3]^2-Pi^2*Z[i,3]^2))
      dScoreBetai[4,4] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,4]^2-Pi^2*Z[i,4]^2))
      dScoreBetai[1,2] <- dScoreBetai[2,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]-Pi^2*Z[i,2]))
      dScoreBetai[1,3] <- dScoreBetai[3,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,3]-Pi^2*Z[i,3]))
      dScoreBetai[1,4] <- dScoreBetai[4,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,4]-Pi^2*Z[i,4]))
      dScoreBetai[2,3] <- dScoreBetai[3,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]*Z[i,3]-Pi^2*Z[i,2]*Z[i,3]))
      dScoreBetai[2,4] <- dScoreBetai[4,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]*Z[i,4]-Pi^2*Z[i,2]*Z[i,4]))
      dScoreBetai[3,4] <- dScoreBetai[4,3] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,3]*Z[i,4]-Pi^2*Z[i,3]*Z[i,4]))
      
      for(j in 5:10) dScoreBetai[1,j] <- dScoreBetai[j,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))-Pi^2*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))))
      for(j in 5:10) dScoreBetai[2,j] <- dScoreBetai[j,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))-Pi^2*Z[i,2]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))))
      for(j in 5:10) dScoreBetai[3,j] <- dScoreBetai[j,3] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,3]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))-Pi^2*Z[i,3]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))))
      for(j in 5:10) dScoreBetai[4,j] <- dScoreBetai[j,4] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,4]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))-Pi^2*Z[i,4]*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))))
      for(j in 5:10) dScoreBetai[j,j] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(v/sqrt(vari)*sigmahat[j-4,j-4]-v*(vari)^(-3/2)*(betaBhat%*%sigmahat[,j-4])^2)+Pi*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))^2-Pi^2*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))^2))
      for(j in 5:9) for(k in (j+1):10) dScoreBetai[j,k] <- dScoreBetai[k,j] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(v/sqrt(vari)*sigmahat[j-4,k-4]-v*(vari)^(-3/2)*(betaBhat%*%sigmahat[,(j-4)])*(betaBhat%*%sigmahat[,(k-4)]))+Pi*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))*(bhat[k-4]+v*betaBhat%*%sigmahat[,k-4]/(sqrt(vari)))-Pi^2*(bhat[j-4]+v*betaBhat%*%sigmahat[,j-4]/(sqrt(vari)))*(bhat[k-4]+v*betaBhat%*%sigmahat[,k-4]/(sqrt(vari)))))
      dScoreBetas <- c(dScoreBetas, list(dScoreBetai))
      
      for(j in 1:Ni1) Scoretau2[i,1] <- Scoretau2[i,1] - Xi1[j,4]/(2*tauhat2[1]) + Xi1[j,4]/(2*tauhat2[1]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^2)-(1-Xi1[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pmax(pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1])),1e-200)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/(2*tauhat2[1]^(3/2)))
      for(j in 1:Ni2) Scoretau2[i,2] <- Scoretau2[i,2] - Xi2[j,4]/(2*tauhat2[2]) + Xi2[j,4]/(2*tauhat2[2]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^2)-(1-Xi2[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pmax(pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2])),1e-200)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/(2*tauhat2[2]^(3/2)))
      for(j in 1:Ni3) Scoretau2[i,3] <- Scoretau2[i,3] - Xi3[j,4]/(2*tauhat2[3]) + Xi3[j,4]/(2*tauhat2[3]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^2)-(1-Xi3[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pmax(pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3])),1e-200)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/(2*tauhat2[3]^(3/2)))
      
      for(j in 1:Ni1) dScoretau2i[1,1] <- dScoretau2i[1,1] + Xi1[j,4]/(2*tauhat2[1]^2) - Xi1[j,4]/(tauhat2[1]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^2)-(1-Xi1[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pmax(pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1])),1e-200)*3/4*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/tauhat2[1]^(5/2)+dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pmax(pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1])),1e-200)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^3/(4*tauhat2[1]^(7/2))+(dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pmax(pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1])),1e-200)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/(2*tauhat2[1]^(3/2)))^2))
      for(j in 1:Ni2) dScoretau2i[2,2] <- dScoretau2i[2,2] + Xi2[j,4]/(2*tauhat2[2]^2) - Xi2[j,4]/(tauhat2[2]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^2)-(1-Xi2[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pmax(pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2])),1e-200)*3/4*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/tauhat2[2]^(5/2)+dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pmax(pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2])),1e-200)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^3/(4*tauhat2[2]^(7/2))+(dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pmax(pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2])),1e-200)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/(2*tauhat2[2]^(3/2)))^2))
      for(j in 1:Ni3) dScoretau2i[3,3] <- dScoretau2i[3,3] + Xi3[j,4]/(2*tauhat2[3]^2) - Xi3[j,4]/(tauhat2[3]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^2)-(1-Xi3[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pmax(pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3])),1e-200)*3/4*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/tauhat2[3]^(5/2)+dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pmax(pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3])),1e-200)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^3/(4*tauhat2[3]^(7/2))+(dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pmax(pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3])),1e-200)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/(2*tauhat2[3]^(3/2)))^2))
      dScoretau2s <- c(dScoretau2s,list(dScoretau2i))
    }
    
    if(sum(X1[(X1[,1]==i),4],X2[(X2[,1]==i),4],X3[(X3[,1]==i),4])/(Ni1+Ni2+Ni3)<0.2){
      printitem <- capture.output({Sample <- MCMCmetrop1R(Sampler, theta.init=c( 1.77783525, -0.07435516, 3.8382696, -0.4193622, 2.1533845, -0.5824151), burnin = 50, mcmc = n.MC,seed=sample(1:10000000,1),Y=Y[i],X12=Xi1[,2],X13=Xi1[,3],X14=Xi1[,4],X22=Xi2[,2],X23=Xi2[,3],X24=Xi2[,4],X32=Xi3[,2],X33=Xi3[,3],X34=Xi3[,4],Z=Z[i,],betaAhat=betaAhat,betaBhat=betaBhat,tauhat2=tauhat2,Bhat=Bhat,Dhat=Dhat)})
      Pi <- (1+exp(as.numeric(-Z[i,]%*%betaAhat)-Sample%*%betaBhat))^(-1)
      
      Bhati[i,] <- colMeans(Sample)
      for(j in 1:n.MC) Dhati <- Dhati + t(Sample[j,]-t(Bhat))%*%(Sample[j,]-t(Bhat))/n.MC
      Dhats <- c(Dhats,list(Dhati))
      
      ScoreBetaA[i,1] <- Y[i]-mean(Pi)
      ScoreBetaA[i,2] <- Y[i]*Z[i,2]-mean(Pi*Z[i,2])
      ScoreBetaA[i,3] <- Y[i]*Z[i,3]-mean(Pi*Z[i,3])
      ScoreBetaA[i,4] <- Y[i]*Z[i,4]-mean(Pi*Z[i,4])
      for(j in 1:6) ScoreBetaB[i,j] <- Y[i]*Bhati[i,j]-mean(Pi*Sample[,j])
      
      dScoreBetai[1,1] <- -mean(Pi*(1-Pi))
      dScoreBetai[2,2] <- -mean(Z[i,2]^2*Pi*(1-Pi))
      dScoreBetai[3,3] <- -mean(Z[i,3]^2*Pi*(1-Pi))
      dScoreBetai[4,4] <- -mean(Z[i,4]^2*Pi*(1-Pi))
      dScoreBetai[1,2] <- dScoreBetai[2,1] <- -mean(Z[i,2]*Pi*(1-Pi))
      dScoreBetai[1,3] <- dScoreBetai[3,1] <- -mean(Z[i,3]*Pi*(1-Pi))
      dScoreBetai[1,4] <- dScoreBetai[4,1] <- -mean(Z[i,4]*Pi*(1-Pi))
      dScoreBetai[2,3] <- dScoreBetai[3,2] <- -mean(Z[i,2]*Z[i,3]*Pi*(1-Pi))
      dScoreBetai[2,4] <- dScoreBetai[4,2] <- -mean(Z[i,2]*Z[i,4]*Pi*(1-Pi))
      dScoreBetai[3,4] <- dScoreBetai[4,3] <- -mean(Z[i,3]*Z[i,4]*Pi*(1-Pi))
      
      for(j in 5:10) dScoreBetai[1,j] <- dScoreBetai[j,1] <- -mean(Pi*(1-Pi)*Sample[,j-4])
      for(j in 5:10) dScoreBetai[2,j] <- dScoreBetai[j,2] <- -mean(Pi*(1-Pi)*Sample[,j-4]*Z[i,2])
      for(j in 5:10) dScoreBetai[3,j] <- dScoreBetai[j,3] <- -mean(Pi*(1-Pi)*Sample[,j-4]*Z[i,3])
      for(j in 5:10) dScoreBetai[4,j] <- dScoreBetai[j,4] <- -mean(Pi*(1-Pi)*Sample[,j-4]*Z[i,4])
      for(j in 5:10) dScoreBetai[j,j] <- -mean(Pi*(1-Pi)*Sample[,j-4]^2)
      for(j in 5:9) for(k in (j+1):10) dScoreBetai[j,k] <- dScoreBetai[k,j] <- -mean(Pi*(1-Pi)*Sample[,j-4]*Sample[,k-4])
      dScoreBetas <- c(dScoreBetas, list(dScoreBetai))
      
      for(j in 1:Ni1) Scoretau2[i,1] <- Scoretau2[i,1] + mean(-Xi1[j,4]/(2*tauhat2[1]) + Xi1[j,4]/(2*tauhat2[1]^2)*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^2-(1-Xi1[j,4])*dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/(2*tauhat2[1]^(3/2)))
      for(j in 1:Ni2) Scoretau2[i,2] <- Scoretau2[i,2] +  mean(-Xi2[j,4]/(2*tauhat2[2]) + Xi2[j,4]/(2*tauhat2[2]^2)*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^2-(1-Xi2[j,4])*dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/(2*tauhat2[2]^(3/2)))
      for(j in 1:Ni3) Scoretau2[i,3] <- Scoretau2[i,3] +  mean(-Xi3[j,4]/(2*tauhat2[3]) + Xi3[j,4]/(2*tauhat2[3]^2)*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^2-(1-Xi3[j,4])*dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/(2*tauhat2[3]^(3/2)))
      
      for(j in 1:Ni1) dScoretau2i[1,1] <- dScoretau2i[1,1] + mean(Xi1[j,4]/(2*tauhat2[1]^2) - Xi1[j,4]/(tauhat2[1]^3)*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^2-(1-Xi1[j,4])*(-dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*3/4*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/tauhat2[1]^(5/2)+dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^3/(4*tauhat2[1]^(7/2))+(dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/(2*tauhat2[1]^(3/2)))^2))
      for(j in 1:Ni2) dScoretau2i[2,2] <- dScoretau2i[2,2] + mean(Xi2[j,4]/(2*tauhat2[2]^2) - Xi2[j,4]/(tauhat2[2]^3)*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^2-(1-Xi2[j,4])*(-dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*3/4*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/tauhat2[2]^(5/2)+dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^3/(4*tauhat2[2]^(7/2))+(dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/(2*tauhat2[2]^(3/2)))^2))
      for(j in 1:Ni3) dScoretau2i[3,3] <- dScoretau2i[3,3] + mean(Xi3[j,4]/(2*tauhat2[3]^2) - Xi3[j,4]/(tauhat2[3]^3)*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^2-(1-Xi3[j,4])*(-dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*3/4*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/tauhat2[3]^(5/2)+dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^3/(4*tauhat2[3]^(7/2))+(dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/(2*tauhat2[3]^(3/2)))^2))
      dScoretau2s <- c(dScoretau2s,list(dScoretau2i))
      
    }
  }
  
  BhatUpdate <- colMeans(Bhati)		#Update for B
  DhatUpdate <- Reduce('+',Dhats)/n	#Update for D
  
  dScoreBeta <- Reduce('+',dScoreBetas)
  dScoreTau2 <- Reduce('+',dScoretau2s)
  
  BetaUpdate <- c(c(betaAhat,betaBhat)+solve(-(dScoreBeta))%*%apply(cbind(ScoreBetaA,ScoreBetaB),2,sum))
  Tau2Update <- c(tauhat2-solve(dScoreTau2)%*%apply(Scoretau2,2,sum))
  
  return(list(BetaUpdate[1:4],BetaUpdate[5:10],tau2Update,BhatUpdate,DhatUpdate))
}


