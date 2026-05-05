#####################################################################################
#NOTE: PLEASE DO NOT USE CODE FOR ANY PUBLICATION PURPOSES WITHOUT FIRST CONTACTING #
#      PAUL BERNHARDT AT PAUL.BERNHARDT@VILLANOVA.EDU                               #
#####################################################################################

##################Simulation Code for Three Covariates Subject to DLs####################
#	Important Functions: 					                 	#
#	1. EM - calculates one-step approximate EM updates for the parameters; 		#
#	   requires InitPar to be run first to get initial values		        #
#	2. InitPar - finds initial values for algorithm as described in paper	        #
#	3. Maxb - used by nlm/optim to find mode and variance of f(b|data); called      #
#	   within EM							                #
#	4. Sampler - used by MCMCmetrop1R (in MCMC package) to obtain random samples    #
#	   from f(b|data)						                #
#	5. Datagen - generate longitudinal covariate data, baseline covariate data,     #
#	   and binary response							        #
#	6. VarEM - calculate variance estimates for parameter estimates		        #
#########################################################################################



###IMPORTANT NOTE: Occasionally, with strict convergence restrictions (small epsilon), due to the approximate nature of
###the algorithm, the intercept and random effects estimates continue to slowly change for many iterations;
###the estimates of main interest - those associated with the longitudinal covariates and possibly the baseline 
###covariates in the logistic model - perform well. Thus, if issues arise, convergence can be defined based on these 
###parameters only (if inference is only important for these parameters) or by using less strict convergence distances 
###for other parameters. 

####################Necessary Packages######################
library(mnormt) 	#generate/evaluate (multivariate) normal distributions
library(corpcor) 	#check positive definiteness of matrices
library(lme4)  		#fit linear mixed models
library(MCMCpack) 	#generate samples for calculating expectations using MCMC
library(lmec) 		#fit linear mixed models with censored responses


##################Simulation Parameters######################
N <- 1000		#number of simulated data sets
n <- 500		#number of individuals in each simulated data set
numtimes <- 8		#maximum number of longitudinal observations per individual
n.MC <- 10000		#number of random draws used to calculate expecation based on MCMC
epsilon <- 0.01		#convergence distance for maximum absolute relative distance
epsilon2 <- 0.005		#convergence distance for maximum absolute distance

#Censoring Levels
d  <- rep(-100,3)	#0% Censoring
#d <- c(0.05,-0.2,-1.9) #25% Censoring
#d<- c(1.5,1,0.52) 	#50% Censoring

#1-Missingness levels
miss <- c(0.9,0.95,0.95,0.95) #10% of all long. covariates at a given obs. time are missing, plus 5% of additional obs. for each long. cov.

#Integration Points (Gaussian-Hermite Quadrature with 9 evaluation points)
u <-c(1.4806867,1.1904181,1.0813192,1.0328036,1.0185664,1.0328036,1.0813192,1.1904181,1.4806867)
v <- c(4.5127459, 3.205429, 2.076848, 1.0232557,0,-1.023256, -2.076848,-3.205429,-4.512746)


#####################Model Parameters#########################
BetaA <- c(-1,0.01)							#logistic regression coefficients for intercept and Z
BetaB <- c(0.2,-0.1,0.4,-0.3,0.25,-0.3)					#logistic regression coefficients for random effects b0-b5
tau2 <- c(1,1.5,1.25)							#variance parameters in linear mixed models for X1-X3
B <- c(1, 0.2,3,-0.1,-2,0.3)						#mean(b0-b5)
D <- matrix(c(1,0.05,0.1,0.1,0.1,0.1,0.05,0.3,0.1,0.1,0.1,0.1,		#var(b0-b5)
              0.1,0.1,1,0.05,0.1,0.1,0.1,0.1,0.05,0.3,0.1,0.1,
              0.1,0.1,0.1,0.1,1,0.05,0.1,0.1,0.1,0.1,0.05,0.3),6,6)


#########################Simulation###########################
Data <- Datagen(N,n,BetaA,BetaB,Tau2,B,D,numtimes,d,miss)  #Generates data for all N datasets 
Xgen1 <- Data[[1]]		#matrix of all X1 data
Xgen2 <- Data[[2]]		# " X2 "
Xgen3 <- Data[[3]]		# " X2 "
Zgen <- Data[[4]]		# " Z "
Ygen <- Data[[5]]		# " Y "
Nni1 <- Data[[6]]		#number of total X1 observations
Nni2 <- Data[[7]]		# " X2 "
Nni3 <- Data[[8]]		# " X3 "


#Initializing vecotors/matrices
Final <- matrix(0,N,53)	#Matrix of final parameter estimates
Iter <- rep(0,N)	#Vector to keep track of current iteration in EM algorithm
SE <- matrix(0,N,8)	#Stores standard errors for Beta parameters of primary interest
SE2 <- matrix(0,N,8)	#Stores sandwich standard errors for Beta parameters of primary interest
Cov <- matrix(0,N,8)	#Stores coverage probabilities for Beta estimates
Cov2 <- matrix(0,N,8)

#Loop through N datasets, finding parameter estimates with EM algorithm
for(k in 1:N){
  
  #Selecting Y (response), X1,X2,X3 (longitudinal covariates), and Z (baseline covariate) for kth of N simulated datasets
  Y <- Ygen[,k]
  X1 <- Xgen1[1:Nni1[k],(4*k-3):(4*k)]
  X2 <- Xgen2[1:Nni2[k],(4*k-3):(4*k)]
  X3 <- Xgen3[1:Nni3[k],(4*k-3):(4*k)]
  
  Z <- Zgen[,(2*k-1):(2*k)]
  
  #Obtaining initial values for EM algorithm
  Pars <- InitPar(X1,X2,X3,Z,Y)
  
  #EM Algorithm: 1. define ParsMatrix (each row is the update after an EM iteration), starting with a "dummy" row consisting 
  #all 1's and then a row of initial estimates (remaining rows are updates, i.e. 3rd row is first update, 4th row is second, etc.)
  #2. define ways for EM algorithm to stop: a.) in while loop, stops if max absolute relative difference (MARD) among parameters is <0.01
  #b.) if at at least 10th iteration, if the MARD average of the last five iterations has not changed from the previous five more than 0.01
  #c.) the max absolute difference (or mean difference of two sets of 5 consecutive iterations) is <0.005
  #d.) if we have completed 50 iterations (this is just to reduce programming time, should be upped in an actual example)  
  #Note: if the initial parameters fail to lead to convergence, this program uses the true parameter values as starting values
  #In practice, alternative starting values could be tried manually until convergence is achieved or more complex starting values
  #can be derived using more computational expensive techniques, perhaps on a subset of the data
  ParsMatrix <- rbind(rep(1,53),c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])))
  STOP <- FALSE
  Iter[k] <- 1 	#defining current iteration 
  while(max(abs((ParsMatrix[(Iter[k]+1),]-ParsMatrix[(Iter[k]),])/ParsMatrix[(Iter[k]),]))>epsilon & max(abs((ParsMatrix[(Iter[k]+1),]-ParsMatrix[(Iter[k]),])))>epsilon2 & Iter[k] < 50 & STOP==FALSE){
    tryCatch(Pars <- EM(Y,X1,X2,X3,Z,Pars[[1]],Pars[[2]],Pars[[3]],Pars[[4]],Pars[[5]],n.MC), error=function(...) Pars <<- EM(Y,X1,X2,X3,Z,BetaA,BetaB,tau2,B,D,n.MC))
    ParsMatrix <- rbind(ParsMatrix,c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])))
    Iter[k] <- Iter[k] + 1
    if(Iter[k]>=10) if(max(abs((colMeans(ParsMatrix[(Iter[k]-3):(Iter[k]+1),])-colMeans(ParsMatrix[(Iter[k]-8):(Iter[k]-4),]))/colMeans(ParsMatrix[(Iter[k]-8):(Iter[k]-4),])))<epsilon | max(abs((colMeans(ParsMatrix[(Iter[k]-3):(Iter[k]+1),])-colMeans(ParsMatrix[(Iter[k]-8):(Iter[k]-4),]))))<epsilon2) STOP <- TRUE
  }
  
  #defines final parameter estimate as last update unless the iteration was stopped because the last five iterations are similar to the previous five
  if(STOP==FALSE) Final[k,] <- c(as.vector(Pars[[1]]),as.vector(Pars[[2]]),as.vector(Pars[[3]]),as.vector(Pars[[4]]),as.vector(Pars[[5]])) else Final[k,] <- colMeans(ParsMatrix[(Iter[k]-3):(Iter[k]+1),])		
  
  #Calculates variance of estimate using Louis (1982) method with MCMC samples
  Variances <- VarEM(Y,X1,X2,X3,Z,as.vector(Final[k,1:2]),as.vector(Final[k,3:8]),as.vector(Final[k,9:11]),as.vector(Final[k,12:17]),matrix(Final[k,18:53],6,6),2000,3000)
  SE[k,] <- sqrt(diag(Variances[[1]]))	#Variance estimator based on information matrix
  SE2[k,] <- sqrt(diag(Variances[[2]]))	#Sandwich variance estimator
  Cov[k,] <- (Final[k,1:8]-qnorm(.975)*SE[k,]< c(BetaA,BetaB) & Final[k,1:8]+qnorm(.975)*SE[k,]> c(BetaA,BetaB))
  Cov2[k,] <- (Final[k,1:8]-qnorm(.975)*SE2[k,]< c(BetaA,BetaB) & Final[k,1:8]+qnorm(.975)*SE2[k,]> c(BetaA,BetaB))
  
}	


#######################Functions##########################


#Function called in an nlm or optim (maximization function built in R) to find posterior mode and variance of f(b|data;thetahat)
Maxb <- function(Y,X1,X2,X3,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat){
  like <- function(bv){
    llike <- -log(dbinom(Y,1,(exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat)/(1+exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat))))*dmnorm(c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6]),Bhat,Dhat))-sum(log((X1[,4]*dnorm(X1[,3],(bv[1] + bv[2]*X1[,2]),sqrt(tauhat2[1]))+(1-X1[,4])*pnorm(X1[,3],(bv[1] + bv[2]*X1[,2]),sqrt(tauhat2[1])))))-sum(log((X2[,4]*dnorm(X2[,3],(bv[3] + bv[4]*X2[,2]),sqrt(tauhat2[2]))+(1-X2[,4])*pnorm(X2[,3],(bv[3] + bv[4]*X2[,2]),sqrt(tauhat2[2])))))-sum(log((X3[,4]*dnorm(X3[,3],(bv[5] + bv[6]*X3[,2]),sqrt(tauhat2[3]))+(1-X3[,4])*pnorm(X3[,3],(bv[5] + bv[6]*X3[,2]),sqrt(tauhat2[3])))))	
    return(llike)}
  like
}

#Function (whose argument is proportional to the likelihood of f(b|data;thetahat)) called in MCMCmetrop1R for obtaining samples from f(b|data;thetahat)
Sampler <- function(bv,Y,X1,X2,X3,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat){
  log(dbinom(Y,1,(exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat)/(1+exp(Z%*%betaAhat+c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6])%*%betaBhat))))*dmnorm(c(bv[1],bv[2],bv[3],bv[4],bv[5],bv[6]),Bhat,Dhat))+sum(log((X1[,4]*dnorm(X1[,3],(bv[1] + bv[2]*X1[,2]),sqrt(tauhat2[1]))+(1-X1[,4])*pnorm(X1[,3],(bv[1] + bv[2]*X1[,2]),sqrt(tauhat2[1])))))+sum(log((X2[,4]*dnorm(X2[,3],(bv[3] + bv[4]*X2[,2]),sqrt(tauhat2[2]))+(1-X2[,4])*pnorm(X2[,3],(bv[3] + bv[4]*X2[,2]),sqrt(tauhat2[2])))))+sum(log((X3[,4]*dnorm(X3[,3],(bv[5] + bv[6]*X3[,2]),sqrt(tauhat2[3]))+(1-X3[,4])*pnorm(X3[,3],(bv[5] + bv[6]*X3[,2]),sqrt(tauhat2[3])))))	
}


#EM algorithm update function: takes data and current parameter estimates and updates the parameter estimates using one EM step
EM <- function(Y,X1,X2,X3,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat,n.MC){
  
  #Printing current parameter values
  print(c(betaAhat,betaBhat,tauhat2,Bhat,Dhat))
  
  #Initializing vectors and matrices
  Dhats <- dScoreBetas <- dScoretau2s <- NULL
  Bhati <-  matrix(0,n,6)
  ScoreBeta <- matrix(0,n,8)
  Scoretau2 <- matrix(0,n,3)
  
  #Loop through each individual to calculate Score and dScores needed for parameter updates
  for(i in 1:n){
    
    #Defining X matrix for ith individual
    Xi1 <- X1[(X1[,1]==i),]
    Xi2 <- X2[(X2[,1]==i),]
    Xi3 <- X3[(X3[,1]==i),]
    Ni1 <- length(Xi1)/4
    Ni2 <- length(Xi2)/4
    Ni3 <- length(Xi3)/4
    
    #Ensuring that R reads Xi1-Xi3 as matrices
    if(Ni1==1) Xi1 <- matrix(Xi1,1,4)
    if(Ni2==1) Xi2 <- matrix(Xi2,1,4)
    if(Ni3==1) Xi3 <- matrix(Xi3,1,4)
    
    #Initializing individual specific matrices
    Dhati <-  matrix(0,6,6)
    dScoretau2i <- matrix(0,3,3)
    dScoreBetai <- matrix(0,8,8)
    
    #First, use adaptive Gaussian-Hermite Quadrature to estimate scores and information matrices (if not too many censored observations)
    if(sum(X1[(X1[,1]==i),4],X2[(X2[,1]==i),4],X3[(X3[,1]==i),4])/(Ni1+Ni2+Ni3)>=0.2){
      #Finding mode and variance of random effects
      tryCatch(f <-nlm(Maxb(Y[i],Xi1,Xi2,Xi3,Z[i,],betaAhat,betaBhat,tauhat2,Bhat,Dhat),Bhat,hessian=TRUE), error=function(...) f <<- optim(Maxb(Y[i],Xi1,Xi2,Xi3,Z[i,],betaAhat,betaBhat,tauhat2,Bhat,Dhat),inits,method="BFGS",hessian=TRUE))
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
      
      #calculating constributions to B- and D-hat to ith individual
      Bhati[i,] <- bhat
      Dhati <- sigmahat+(bhat-Bhat)%*%t(bhat-Bhat)
      Dhats <- c(Dhats,list(Dhati))
      
      #Calculating scores for Beta parameters
      ScoreBeta[i,1] <- Y[i]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi)
      ScoreBeta[i,2] <- Y[i]*Z[i,2]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*Z[i,2])
      for(j in 3:8) ScoreBeta[i,j] <- Y[i]*Bhati[i,j-2]-sum(u/sqrt(2*pi)*exp(-v^2/2)*Pi*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari))))
      
      #Calculating derivative of Q' for Beta parameters
      dScoreBetai[1,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi-Pi^2))
      dScoreBetai[2,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]^2-Pi^2*Z[i,2]^2))
      dScoreBetai[1,2] <- dScoreBetai[2,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]-Pi^2*Z[i,2]))
      for(j in 3:8) dScoreBetai[1,j] <- dScoreBetai[j,1] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))-Pi^2*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))))
      for(j in 3:8) dScoreBetai[2,j] <- dScoreBetai[j,2] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*Z[i,2]*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))-Pi^2*Z[i,2]*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))))
      for(j in 3:8) dScoreBetai[j,j] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(v/sqrt(vari)*sigmahat[j-2,j-2]-v*(vari)^(-3/2)*(betaBhat%*%sigmahat[,j-2])^2)+Pi*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))^2-Pi^2*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))^2))
      for(j in 3:7) for(k in (j+1):8) dScoreBetai[j,k] <- dScoreBetai[k,j] <- -sum(u/sqrt(2*pi)*exp(-v^2/2)*(Pi*(v/sqrt(vari)*sigmahat[j-2,k-2]-v*(vari)^(-3/2)*(betaBhat%*%sigmahat[,(j-2)])*(betaBhat%*%sigmahat[,(k-2)]))+Pi*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))*(bhat[k-2]+v*betaBhat%*%sigmahat[,k-2]/(sqrt(vari)))-Pi^2*(bhat[j-2]+v*betaBhat%*%sigmahat[,j-2]/(sqrt(vari)))*(bhat[k-2]+v*betaBhat%*%sigmahat[,k-2]/(sqrt(vari)))))
      dScoreBetas <- c(dScoreBetas, list(dScoreBetai))
      
      #Calculating scores for Tau2 parameters
      for(j in 1:Ni1) Scoretau2[i,1] <- Scoretau2[i,1] - Xi1[j,4]/(2*tauhat2[1]) + Xi1[j,4]/(2*tauhat2[1]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^2)-(1-Xi1[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/(2*tauhat2[1]^(3/2)))
      for(j in 1:Ni2) Scoretau2[i,2] <- Scoretau2[i,2] - Xi2[j,4]/(2*tauhat2[2]) + Xi2[j,4]/(2*tauhat2[2]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^2)-(1-Xi2[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/(2*tauhat2[2]^(3/2)))
      for(j in 1:Ni3) Scoretau2[i,3] <- Scoretau2[i,3] - Xi3[j,4]/(2*tauhat2[3]) + Xi3[j,4]/(2*tauhat2[3]^2)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^2)-(1-Xi3[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/(2*tauhat2[3]^(3/2)))
      
      #Calculating derivatives of Q' for Tau2 paramters
      for(j in 1:Ni1) dScoretau2i[1,1] <- dScoretau2i[1,1] + Xi1[j,4]/(2*tauhat2[1]^2) - Xi1[j,4]/(tauhat2[1]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^2)-(1-Xi1[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))*3/4*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/tauhat2[1]^(5/2)+dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)^3/(4*tauhat2[1]^(7/2))+(dnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/sqrt(tauhat2[1]))*(Xi1[j,3]-mu21[j]-sqrt(vari21[j])*v)/(2*tauhat2[1]^(3/2)))^2))
      for(j in 1:Ni2) dScoretau2i[2,2] <- dScoretau2i[2,2] + Xi2[j,4]/(2*tauhat2[2]^2) - Xi2[j,4]/(tauhat2[2]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^2)-(1-Xi2[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))*3/4*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/tauhat2[2]^(5/2)+dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)^3/(4*tauhat2[2]^(7/2))+(dnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/sqrt(tauhat2[2]))*(Xi2[j,3]-mu22[j]-sqrt(vari22[j])*v)/(2*tauhat2[2]^(3/2)))^2))
      for(j in 1:Ni3) dScoretau2i[3,3] <- dScoretau2i[3,3] + Xi3[j,4]/(2*tauhat2[3]^2) - Xi3[j,4]/(tauhat2[3]^3)*sum(u/sqrt(2*pi)*exp(-v^2/2)*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^2)-(1-Xi3[j,4])*sum(u/sqrt(2*pi)*exp(-v^2/2)*(-dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))*3/4*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/tauhat2[3]^(5/2)+dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)^3/(4*tauhat2[3]^(7/2))+(dnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/sqrt(tauhat2[3]))*(Xi3[j,3]-mu23[j]-sqrt(vari23[j])*v)/(2*tauhat2[3]^(3/2)))^2))
      dScoretau2s <- c(dScoretau2s,list(dScoretau2i))
    }
    
    #If >80% censoring on long. covariates, use MCMC methods to estimate scores and information matrices (slower)
    if(sum(X1[(X1[,1]==i),4],X2[(X2[,1]==i),4],X3[(X3[,1]==i),4])/(Ni1+Ni2+Ni3)<0.2){
      
      #obtaining MCMC samples 
      inits <- c(0,0,0,0,0,0)
      printitem <- capture.output({Sample <- MCMCmetrop1R(Sampler, theta.init=inits, burnin = 50, mcmc = n.MC,seed=sample(1:10000000,1),Y=Y[i],X1=Xi1,X2=Xi2,X3=Xi3,Z=Z[i,],betaAhat=betaAhat,betaBhat=betaBhat,tauhat2=tauhat2,Bhat=Bhat,Dhat=Dhat)})
      
      #Defining 'Pi' since it appears in formulas frequently
      Pi <- (1+exp(as.numeric(-Z[i,]%*%betaAhat)-Sample%*%betaBhat))^(-1)
      
      #calculating constributions to B- and D-hat to ith individual (Dhat being cumbersome)
      Bhati[i,] <- colMeans(Sample)
      for(j in 1:n.MC) Dhati <- Dhati + t(Sample[j,]-t(Bhat))%*%(Sample[j,]-t(Bhat))/n.MC
      Dhats <- c(Dhats,list(Dhati))
      
      #Calculating scores for Beta parameters
      ScoreBeta[i,1] <- Y[i]-mean(Pi)
      ScoreBeta[i,2] <- Y[i]*Z[i,2]-mean(Pi*Z[i,2])
      for(j in 3:8) ScoreBeta[i,j] <- Y[i]*Bhati[i,j-2]-mean(Pi*Sample[,j-2])
      
      #Calculating derivative of Q' for Beta parameters
      dScoreBetai[1,1] <- -mean(Pi*(1-Pi))
      dScoreBetai[2,2] <- -mean(Z[i,2]^2*Pi*(1-Pi))
      dScoreBetai[1,2] <- dScoreBetai[2,1] <- -mean(Z[i,2]*Pi*(1-Pi))
      for(j in 3:8) dScoreBetai[1,j] <- dScoreBetai[j,1] <- -mean(Pi*(1-Pi)*Sample[,j-2])
      for(j in 3:8) dScoreBetai[2,j] <- dScoreBetai[j,2] <- -mean(Pi*(1-Pi)*Sample[,j-2]*Z[i,2])
      for(j in 3:8) dScoreBetai[j,j] <- -mean(Pi*(1-Pi)*Sample[,j-2]^2)
      for(j in 3:7) for(k in (j+1):8) dScoreBetai[j,k] <- dScoreBetai[k,j] <- -mean(Pi*(1-Pi)*Sample[,j-2]*Sample[,k-2])
      dScoreBetas <- c(dScoreBetas, list(dScoreBetai))
      
      #Calculating scores for Tau2 parameters
      for(j in 1:Ni1) Scoretau2[i,1] <- Scoretau2[i,1] + mean(-Xi1[j,4]/(2*tauhat2[1]) + Xi1[j,4]/(2*tauhat2[1]^2)*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^2-(1-Xi1[j,4])*dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/(2*tauhat2[1]^(3/2)))
      for(j in 1:Ni2) Scoretau2[i,2] <- Scoretau2[i,2] +  mean(-Xi2[j,4]/(2*tauhat2[2]) + Xi2[j,4]/(2*tauhat2[2]^2)*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^2-(1-Xi2[j,4])*dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/(2*tauhat2[2]^(3/2)))
      for(j in 1:Ni3) Scoretau2[i,3] <- Scoretau2[i,3] +  mean(-Xi3[j,4]/(2*tauhat2[3]) + Xi3[j,4]/(2*tauhat2[3]^2)*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^2-(1-Xi3[j,4])*dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/(2*tauhat2[3]^(3/2)))
      
      #Calculating derivatives of Q' for Tau2 paramters
      for(j in 1:Ni1) dScoretau2i[1,1] <- dScoretau2i[1,1] + mean(Xi1[j,4]/(2*tauhat2[1]^2) - Xi1[j,4]/(tauhat2[1]^3)*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^2-(1-Xi1[j,4])*(-dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*3/4*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/tauhat2[1]^(5/2)+dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))^3/(4*tauhat2[1]^(7/2))+(dnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))/pnorm((Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/sqrt(tauhat2[1]))*(Xi1[j,3]-Sample[,1:2]%*%c(1,Xi1[j,2]))/(2*tauhat2[1]^(3/2)))^2))
      for(j in 1:Ni2) dScoretau2i[2,2] <- dScoretau2i[2,2] + mean(Xi2[j,4]/(2*tauhat2[2]^2) - Xi2[j,4]/(tauhat2[2]^3)*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^2-(1-Xi2[j,4])*(-dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*3/4*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/tauhat2[2]^(5/2)+dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))^3/(4*tauhat2[2]^(7/2))+(dnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))/pnorm((Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/sqrt(tauhat2[2]))*(Xi2[j,3]-Sample[,3:4]%*%c(1,Xi2[j,2]))/(2*tauhat2[2]^(3/2)))^2))
      for(j in 1:Ni3) dScoretau2i[3,3] <- dScoretau2i[3,3] + mean(Xi3[j,4]/(2*tauhat2[3]^2) - Xi3[j,4]/(tauhat2[3]^3)*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^2-(1-Xi3[j,4])*(-dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*3/4*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/tauhat2[3]^(5/2)+dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))^3/(4*tauhat2[3]^(7/2))+(dnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))/pnorm((Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/sqrt(tauhat2[3]))*(Xi3[j,3]-Sample[,5:6]%*%c(1,Xi3[j,2]))/(2*tauhat2[3]^(3/2)))^2))
      dScoretau2s <- c(dScoretau2s,list(dScoretau2i))
      
    }
  }
  
  #Calculate the parameter updates:
  BhatUpdate <- colMeans(Bhati)		#Update for B
  DhatUpdate <- Reduce('+',Dhats)/n	#Update for D
  
  #Update for Beta
  dScoreBeta <- Reduce('+',dScoreBetas)
  betaUpdate <- c(c(betaAhat,betaBhat)-solve(dScoreBeta)%*%apply(ScoreBeta,2,sum))
  
  #Update for Tau2
  dScoreTau2 <- Reduce('+',dScoretau2s)
  tau2Update <- c(tauhat2-solve(dScoreTau2)%*%apply(Scoretau2,2,sum))
  
  return(list(betaUpdate[1:2],betaUpdate[3:8],tau2Update,BhatUpdate,DhatUpdate))
}


#Function to generate the data (also outputs some coefficients 
Datagen <- function(N,n,BetaA,BetaB,Tau2,B,D,numtimes,d,miss){
  #initializing matrices 
  Xgen1 <- matrix(0,n*numtimes,4*N) #for all longitudinal covariate X1 data (all individuals)
  Xgen2 <- matrix(0,n*numtimes,4*N) #for all longitudinal covariate X2 data (all individuals)
  Xgen3 <- matrix(0,n*numtimes,4*N) #for all longitudinal covariate X3 data (all individuals)
  Ygen <- matrix(0,n,N)		    #for all binary response Y data (all individuals)
  Zgen <- matrix(0,n,2*N)		    #for all baseline covariate Z data (all individuals)
  Nni1 <- rep(0,N)		          #vector of number of X1 observations for all individuals
  Nni2 <- rep(0,N)		 	    #vector of number of X2 observations for all individuals
  Nni3 <- rep(0,N)			    #vector of number of X3 observations for all individuals
  
  for(k in 1:N){
    X1 <- X2 <- X3 <- NULL  #matrices to hold all longitudinal covariate data
    
    #Step 1: generate number of observations from 1 to 8 for each individual
    totaltimes <- sample(1:numtimes, n, replace=TRUE)
    
    #Step 2: generate random effects for all individuals
    b <- rmnorm(n, B, D)
    
    #Step 3: generate a fully observed, baseline "Z" variable representing age from the GenIMS dataset
    Z <- rbeta(n,3,2)*84+18  
    Z <- cbind(1,Z)
    
    #Step 4: generate X1, X2, X3, the longitudinal covariates subject to DLs
    for(i in 1:n){
      
      #Step 4a: generate a vector of missingness indicators for ith individual corresponding 
      #	    to missingness accross all three long. covariates at a given obs. time  (0 = missing)
      p <- rbinom(totaltimes[i],1,miss[1]) 
      
      
      #Step 4b: generate X1
      #Step 4b1: generate a vector of missingness indicators for ith individual specific to covariate X1 (0 = missing)
      p1 <- rbinom(totaltimes[i],1,miss[2])
      
      #Step 4b2: check to make sure that at least one observation is made for X1 for ith individual if no observation
      #          is made, let missingness vectors p and p1 be all 1's (this will make missingness % slightly < indicated)
      if(sum(p1*p*rep(1,totaltimes[i]))<1) {p <- rep(1,totaltimes[i])
      p1 <- rep(1, totaltimes[i])}
      
      #Step 4b3: generate latent X1 covariate values and multiply by p and p1 so that we observe "0" when value is missing
      Xvals <- rnorm(totaltimes[i], (b[i,1]+b[i,2]*c(1:totaltimes[i])),sqrt(tau2[1]))*p*p1
      
      #Step 4b4: replace latent X1 by DL when it is less than DL (and not missing)
      Xvals[Xvals < d[1] & Xvals!=0] <- d[1]
      
      #Step 4b5: censoring indicator vector ( = 1 implies value above DL)
      Obs <- rep(1,totaltimes[i])
      Obs[Xvals==d[1]] <- 0
      
      #Step 4b6: create matrix of individual id number in column 1, observation time in column 2, observed X1 values 
      #	     (=DL when latent X1<DL and =0 if missing) in column 3,  and censoring indicator vector in column 4
      X1 <- rbind(X1, cbind(rep(i,totaltimes[i]),1:totaltimes[i],Xvals,Obs))
      
      
      #Step 4c: generate X2 (see 4b1-4b6 for line-by-line explanation)
      p2 <- rbinom(totaltimes[i],1,miss[3])
      Obs <- rep(1,totaltimes[i])
      if(sum(p2*p*rep(1,totaltimes[i]))<1) {p <- rep(1,totaltimes[i])
      p2 <- rep(1, totaltimes[i])}
      Xvals <- rnorm(totaltimes[i], (b[i,3]+b[i,4]*c(1:totaltimes[i])),sqrt(tau2[2]))*p*p2
      Xvals[Xvals < d[2] & Xvals!=0] <- d[2]
      Obs[Xvals==d[2]] <- 0
      X2 <- rbind(X2, cbind(rep(i,totaltimes[i]),1:totaltimes[i],Xvals,Obs))
      
      
      #Step 4d: generate X3 (see 4b1-4b6 for line-by-line explanation)
      p3 <- rbinom(totaltimes[i],1,miss[4])
      Obs <- rep(1,totaltimes[i])
      if(sum(p3*p*rep(1,totaltimes[i]))<1) {p <- rep(1,totaltimes[i])
      p3 <- rep(1, totaltimes[i])}
      Xvals <- rnorm(totaltimes[i], (b[i,5]+b[i,6]*c(1:totaltimes[i])),sqrt(tau2[3]))*p*p3
      Xvals[Xvals < d[3] & Xvals!=0] <- d[3]
      Obs[Xvals==d[3]] <- 0
      X3 <- rbind(X3, cbind(rep(i,totaltimes[i]),1:totaltimes[i],Xvals,Obs))
    }
    
    #Step 4e: delete observations that were not observated (we assume missing at random in analysis)
    X1 <- X1[-which(X1[,3]==0),]
    X2 <- X2[-which(X2[,3]==0),]
    X3 <- X3[-which(X3[,3]==0),]
    
    #Step 5: generate Y
    Y <- rbinom(n,1,exp(Z%*%BetaA+b%*%BetaB)/(1+exp(Z%*%BetaA+b%*%BetaB))) #based on inverse logit
    
    #Step 6: combin X1-X3, Y, and Z for all individuals into matrices containing all individuals
    Xgen1[1:length(X1[,1]),(k*4-3):(k*4)] <- X1
    Xgen2[1:length(X2[,1]),(k*4-3):(k*4)] <- X2
    Xgen3[1:length(X3[,1]),(k*4-3):(k*4)] <- X3
    Ygen[,k] <- Y
    Zgen[,(2*k-1):(2*k)] <- Z
    Nni1[k] <- length(X1[,1])
    Nni2[k] <- length(X2[,1])
    Nni3[k] <- length(X3[,1])
  }
  
  return(list(Xgen1,Xgen2,Xgen3,Zgen,Ygen,Nni1,Nni2,Nni3))
}

#Function used to obtain initial values for the EM algorithm; inputs data and outputs initial values
InitPar <- function(X1,X2,X3,Z,Y){
  #Initializing matrices
  Dhat <- matrix(0,6,6)
  Bhat <- rep(0,6)
  tauhat2 <- rep(0,3)
  bEst <- NULL
  
  #Adding intercept term in a matrix containing the covariate data (necessary for syntax of lmer and lmec)
  X21 <- cbind(1,X1[,2])
  X22 <- cbind(1,X2[,2])
  X23 <- cbind(1,X3[,2])
  
  #Running LMM for X1 covariate if no censored observations (quicker to use lmer than lmec if applicable); 
  #initial values for B, tau2, and D found as described in paper
  if(sum(1-X1[,4])==0) {vals1 <- lmer(X1[,3] ~ X1[,2] + (X1[,2] | X1[,1])) 
  Bhat[1:2] <- fixef(vals1)[1:2]
  tauhat2[1] <- attr(VarCorr(vals1), "sc")^2
  Dhat[1:2,1:2] <- matrix(c(as.vector(VarCorr(vals1)$X[,1]),as.vector(VarCorr(vals1)$X[,2])),2,2)
  bEst <- cbind(bEst,as.matrix(ranef(vals1)$X1))}
  
  #Running LMM for X1 if some covariate values are censored
  if(sum(1-X1[,4])>0) {vals1 <- lmec(yL=X1[,3],cens=c(X1[,3]==d[1]),X=X21,Z=X21,cluster=X1[,1],epsstop=0.1)
  Bhat[1:2] <- vals1$beta
  tauhat2[1] <- vals1$sigma^2
  Dhat[1:2,1:2] <- vals1$Psi
  bEst <- cbind(bEst,t(vals1$bi))}
  
  #Running LMM for X2 covariate if no censored observations
  if(sum(1-X2[,4])==0) {vals2 <- lmer(X2[,3] ~ X2[,2] + (X2[,2] | X2[,1])) 
  Bhat[3:4] <- fixef(vals2)[1:2]
  tauhat2[2] <- attr(VarCorr(vals2), "sc")^2
  Dhat[3:4,3:4] <- matrix(c(as.vector(VarCorr(vals2)$X[,1]),as.vector(VarCorr(vals2)$X[,2])),2,2)
  bEst <- cbind(bEst,as.matrix(ranef(vals2)$X2))}
  
  #Running LMM for X2 if some covariate values are censored
  if(sum(1-X2[,4])>0)  {vals2 <- lmec(yL=X2[,3],cens=c(X2[,3]==d[2]),X=X22,Z=X22,cluster=X2[,1],epsstop=0.1)
  Bhat[3:4] <- vals2$beta
  tauhat2[2] <- vals2$sigma^2
  Dhat[3:4,3:4] <- vals2$Psi
  bEst <- cbind(bEst,t(vals2$bi))}
  
  #Running LMM for X3 covariate if no censored observations
  if(sum(1-X3[,4])==0) {vals3 <- lmer(X3[,3] ~ X3[,2] + (X3[,2] | X3[,1]),verbose=TRUE)
  Bhat[5:6] <- fixef(vals3)[1:2]
  tauhat2[3] <- attr(VarCorr(vals3), "sc")^2
  Dhat[5:6,5:6] <- matrix(c(as.vector(VarCorr(vals3)$X[,1]),as.vector(VarCorr(vals3)$X[,2])),2,2)
  bEst <- cbind(bEst,ranef(vals3)$X3)}
  
  #Running LMM for X3 if some covariate values are censored
  if(sum(1-X3[,4])>0)  {vals3 <- lmec(yL=X3[,3],cens=c(X3[,3]==d[3]),X=X23,Z=X23,cluster=X3[,1],epsstop=0.1)
  Bhat[5:6] <- vals3$beta
  tauhat2[3] <- vals3$sigma^2
  Dhat[5:6,5:6] <- vals3$Psi	
  bEst <- cbind(bEst,t(vals3$bi))}	
  
  coefs <- coef(glm(Y ~ Z[,2] + bEst[,1] + bEst[,2]+ bEst[,3] + bEst[,4] + bEst[,5] + bEst[,6],family="binomial"))
  
  Pars <- list(coefs[1:2],coefs[3:8],tauhat2,Bhat,Dhat)
  
  return(Pars)
}

#Function to obtain variances estimates for parameters; here, for simplicity, this function only calculates for 
#Beta parameters, ignoring Tau2, B, and D parameters (which will affect variance estimation slightly)
VarEM <- function(Y,X1,X2,X3,Z,betaAhat,betaBhat,tauhat2,Bhat,Dhat,n.samples,n.samples2){
  #initializing matrices
  dScoreBeta <- ScoreBetaVar <- ScoreBetaVar2 <-  matrix(0,8,8)
  
  #Looping through each individiual to calculates Q'' and var(Q')
  for(i in 1:n){
    dScoreBetai  <-  matrix(0,8,8)
    ScoreBeta <- rep(0,8)
    
    
    #Defining X matrix for ith individual
    Xi1 <- matrix(X1[(X1[,1]==i),],length(X1[(X1[,1]==i),])/4,4)
    Xi2 <- matrix(X2[(X2[,1]==i),],length(X2[(X2[,1]==i),])/4,4)
    Xi3 <- matrix(X3[(X3[,1]==i),],length(X3[(X3[,1]==i),])/4,4)
    Ni1 <- dim(Xi1)[1]
    Ni2 <- dim(Xi2)[1]
    Ni3 <- dim(Xi3)[1]
    
    #Generating an MCMC sample from an approximated normal if at least 20% of longitudinal covariates are observed above DL
    if((sum(Xi1[,4])+sum(Xi2[,4])+sum(Xi3[,4]))/(Ni1+Ni2+Ni3)>=0.2){
      FinalScore <- matrix(0,n.samples,8) 
      
      #Finding mode and variance of random effects
      f <-nlm(Maxb(Y[i],Xi1,Xi2,Xi3,Z[i,],betaAhat,betaBhat,tauhat2,Bhat,Dhat),Bhat,hessian=TRUE)
      Est <- f$estimate
      Vari <- solve(f$hessian)
      Sample <- rmnorm(n.samples,Est,Vari)
    }
    
    #Generating an MCMC sample from f(b|data) if less than 20% of longitudinal covariates are observed above DL
    if((sum(Xi1[,4])+sum(Xi2[,4])+sum(Xi3[,4]))/(Ni1+Ni2+Ni3)<0.2){
      FinalScore <- matrix(0,n.samples2,8)
      inits <- c(0,0,0,0,0,0)			
      printitem <- capture.output({Sample <- MCMCmetrop1R(Sampler, inits, burnin = 50, mcmc = n.samples2,seed=sample(1:10000000,1),Y=Y[i],X1=Xi1,X2=Xi2,X3=Xi3,Z=Z[i,],betaAhat=betaAhat,betaBhat=betaBhat,tauhat2=tauhat2,Bhat=Bhat,Dhat=Dhat)})
    }
    
    #Quantity used in calculations below frequently
    Pi <- (1+exp(as.numeric(-Z[i,]%*%betaAhat)-Sample%*%betaBhat))^(-1)
    
    #Calculting Score functions for ith inidividual
    ScoreBeta[1] <- Y[i]-mean(Pi)
    ScoreBeta[2] <- Y[i]*Z[i,2]-mean(Pi*Z[i,2])
    ScoreBeta[3] <- Y[i]*mean(Sample[,1])-mean(Pi*Sample[,1])
    ScoreBeta[4] <- Y[i]*mean(Sample[,2])-mean(Pi*Sample[,2])
    ScoreBeta[5] <- Y[i]*mean(Sample[,3])-mean(Pi*Sample[,3])
    ScoreBeta[6] <- Y[i]*mean(Sample[,4])-mean(Pi*Sample[,4])
    ScoreBeta[7] <- Y[i]*mean(Sample[,5])-mean(Pi*Sample[,5])
    ScoreBeta[8] <- Y[i]*mean(Sample[,6])-mean(Pi*Sample[,6])
    
    #Calulating centered scores for ith individual for each sample point
    FinalScore[,1] <- (Y[i]-Pi -ScoreBeta[1])
    FinalScore[,2] <- (Y[i]*Z[i,2]-Pi*Z[i,2]-ScoreBeta[2])
    FinalScore[,3] <- (Y[i]*Sample[,1]-Pi*Sample[,1]-ScoreBeta[3])
    FinalScore[,4] <- (Y[i]*Sample[,2]-Pi*Sample[,2]-ScoreBeta[4])
    FinalScore[,5] <- (Y[i]*Sample[,3]-Pi*Sample[,3]-ScoreBeta[5])
    FinalScore[,6] <- (Y[i]*Sample[,4]-Pi*Sample[,4]-ScoreBeta[6])
    FinalScore[,7] <- (Y[i]*Sample[,5]-Pi*Sample[,5]-ScoreBeta[7])
    FinalScore[,8] <- (Y[i]*Sample[,6]-Pi*Sample[,6]-ScoreBeta[8])
    
    #Calculating Q''
    dScoreBetai[1,1] <- -mean(Pi*(1-Pi))
    dScoreBetai[2,2] <- -mean(Z[i,2]^2*Pi*(1-Pi))
    dScoreBetai[1,2] <- dScoreBetai[2,1] <- -mean(Z[i,2]*Pi*(1-Pi))
    for(j in 3:8) dScoreBetai[1,j] <- dScoreBetai[j,1] <- -mean(Pi*(1-Pi)*Sample[,j-2])
    for(j in 3:8) dScoreBetai[2,j] <- dScoreBetai[j,2] <- -mean(Pi*(1-Pi)*Sample[,j-2]*Z[i,2])
    for(j in 3:8) dScoreBetai[j,j] <- -mean(Pi*(1-Pi)*Sample[,j-2]^2)
    for(j in 3:7) for(k in (j+1):8) dScoreBetai[j,k] <- dScoreBetai[k,j] <- -mean(Pi*(1-Pi)*Sample[,j-2]*Sample[,k-2])
    dScoreBeta <- dScoreBeta + dScoreBetai	
    
    #Adding across all individuals
    ScoreBetaVar <- ScoreBetaVar + t(FinalScore)%*%FinalScore/length(FinalScore[,1])				
    ScoreBetaVar2 <- ScoreBetaVar2 + ScoreBeta%*%t(ScoreBeta)
  }
  
  #Calculating info and then returning both regular and sandwich variance matrices
  Info <- -dScoreBeta - ScoreBetaVar
  return(list(solve(Info),solve(Info)%*%ScoreBetaVar2%*%solve(Info)))
  
}

