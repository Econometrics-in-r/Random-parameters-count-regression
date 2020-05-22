rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)

# reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)
#Sample Size
N <- nrow(dataset) 
#Dependent variable: crash counts
DVar <- dataset$Crash 

# Halton Draws 
preparedraws=function()
{
  d=1
  while(d<(length(normaldraws)+1))
  {
    draws1[,normaldraws[d]]<<- qnorm(draws1[,normaldraws[d]])
    d=d+1
  }
}

Ndraws=500      # set number of draws 
dimensions=2    # define number of random parameters in the model

# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# assign names to individual sets of draws - need one entry per dimension
colnames(draws1)=c("HRbeta1","HRbeta2")
# define whether any draws should be transformed to Normals, which is also needed for e.g. lognormals (leave empty if not)
normaldraws=c("HRbeta1","HRbeta2")

# preparing draws for estimation - this may take a while
preparedraws()

# fixing parameters across grouped observations i.e. grouped random parameters
# Do not use if there is no panel
block = length(unique(dataset[,'ID']))
ngroup = length(unique(dataset[,'Group']))
for (i in 1:Ndraws){
  tempInd = ((i-1)*block*ngroup) + (1:block)
  for (ii in 2:ngroup){
    draws1[tempInd+(ii-1)*block,] = draws1[tempInd,]
  }
}

## data preparation
# separating the variables with fixed parameters 
dataF =  as.matrix(data.frame(1,log(dataset$Length)))
# separating the variables with random parameters 
dataR = as.matrix(data.frame(log(dataset$AADT),dataset$LWIDTH))

dataR2=NULL
Dvar2 = NULL
for(i in 1:Ndraws){
  dataR2=rbind(dataR2,dataR)
  Dvar2 = c(Dvar2,DVar)
}

draws1 = draws1[,1:dimensions]

# Likelihood function
LL <- function(params){  
  disp <- params[1] # Dispersion Parameter
  Fbeta <- params[2:3] # Fixed parameters in the mean Function
  MRbeta <- params[4:5]  # Mean of Random parameters in the mean function
  SDRbeta <- params[6:7]  # Std of Random parameters in the mean function
  
  # vector of indipendent variables with fixed parameters
  offset = rep.int(dataF%*%as.matrix(Fbeta,ncol=1),Ndraws)
  # simulating random parameters from their means and standard deviation
  beta = t( t(draws1)*SDRbeta + MRbeta )
  # constructing the mean function
  mu <- exp(offset+rowSums(dataR2*beta))
  # simulated maximum loglikelihood for negative binomial distribution
  loglik <-  sum(log(rowMeans(matrix(dnbinom(Dvar2,size=disp,mu=mu,log = F), ncol = Ndraws))))
  
  return(loglik)
}

# initial values for optimization
init <- c(2,#dispersion parameter
          -5.5,0.66,#fixed parameters
          0.17,-0.14,#mean of random parameters
          0.05,0.08)#standard deviation of random parameters

# optimization (maximization of likelihood function)
fit1 <- maxLik(LL,start=init,method="BFGS")

summary(fit1)

# Predictions, Residuals and Measures of Fit 

params <- fit1$estimate
Fbeta <- params[2:3] # Fixed parameters in Mu Function
MRbeta <- params[4:5]  # Mean of Random parameters in Mu function

offset = as.vector(dataF%*%as.matrix(Fbeta,ncol=1))
mu <- exp(offset+log(dataset$AADT)*MRbeta[1]+dataset$LWIDTH*MRbeta[2])

# fitted values
muFit <- mu
# residuals
RESIDS <- DVar-muFit
# absolute residuals
ABSRES <- abs(RESIDS)
# squared residuals
SQRES <- (RESIDS)^2
# number of estimated parameters in the model
P <- NROW(params)
# mean absolute deviance
MAD <- sum(ABSRES)/N
# mean squared predictive error
MSPE <- sum(SQRES)/N
# goodness of fit
GOF <- data.frame(MAD=MAD,MSPE=MSPE)
  
