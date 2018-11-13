
#Brief intro: estimated Glomerular filtration rate (eGFR) is a measure of kidney function. Lower value, poorer kidney function.
#Given a patient, his eGFR measurements can be extracted from hospital electronic health record (EHR) database.
#Aim: to model kidney function trajectory for each patient based on his longitudinal eGFR measurements by using Penalized Thin-plate Spline Model under Bayesian framework.

#clear working directory
rm(list=ls())


#load required R packages
library(SemiPar)
library(mgcv)
library(rjags)
library(invgamma)
library(pracma)
library(matlib)
library(ggplot2)

#load defined function for math transformation that will be used for modeling thin-plate spline under bayesian framework
source("MathTransform.r")
source("ConvergenceDiag.r")
source("PosteriorPlot.r")
source("DerivativeCurve.r")



#import patient data: 
##time: hospital visit time
##interval_year: interval from his first visit (by year)
##eGFR: a measure of kidney function
patient.df=read.csv("sample_data.csv", header=T, stringsAsFactors=F)

#show no. of eGFR measures for this patient
n = dim(patient.df)[1]

attach(patient.df)

# specify knots that allow us model patient's eGFR trajectory as a bending curve
knots = seq(1, floor(max(interval_year)))

# no. of knots
K = length(knots)

# do math transformation
# X: a 
X = MathTransform(interval_year, knots)[[1]]
Y = MathTransform(interval_year, knots)[[2]]
msqrt_inv = MathTransform(interval_year, knots)[[3]]

#----------------------Begin: Running row-rank thin-plate spline regression under Frequentist framework-----------------------------
fit.spm=spm(eGFR~f(interval_year, knots=knots))
summary(fit.spm)

#display a plot for the modeled thin-plate spline  
## prepare range of follow-up time in x-axis  	  
min_x=min(0, min(interval_year))
max_x=ceiling(max(interval_year))

## prepare range of eGFR in y-axis	  
min_y = min(min(eGFR), max(eGFR)-50)
max_y=max(eGFR)
	  
#Take a glance at fitted curve under Frequentist framework
#plot(fit.spm, ylim=c(min_y, max_y), xlim=c(min_x, max_x), xlab="follow-up year", ylab="eGFR", main="Modeling trajectory under Frequentist ")
#points(interval_year, eGFR)
#----------------------End: Running row-rank thin-plate spline regression under Frequentist framework-----------------------------


#----------------------Begin: Running penalized thin-plate spline regression under Bayesian framework (JAGs)-----------------------------------------------------
# Brief intro: eGFR trajectory is modeled by a low-rank thin-plate spline regression. 
# To avoid overfitting (i.e., having too complex model), we penalize the number of knots in the regression.
# Li et al. 2012 and Crainiceanu et al. 2005 found that this penalized thin-plate spline regression mathematically
# is equivalent to a linear mixed-model: eGFR = X * a + Z * b + eps, where
## X: a 2-by-n matrix (derived from "interval_year")
## unknown parameter a=[a1, a2]: a 1-by-2 parameter for fixed-effects following 2-dimensional Normal distribution
## Norm(mu_a, a.var), where "mu_a" is mean vector and "a.var" is covariance matrix.
## Z: a K-by-K matrix (obtained from MathTransform outputs)
## unknown parameter b: a 1-by-K vector following a K-dimensional Normal distribution Norm(0, Variance_b),
### variance_b = sigma_b * Lambda
#### sigma_b is a unknown value that need be estimated
#### Lambda is a K-by-K matrix (derived from data), and 
## eps: independent residual errors following a n-dimensional Normal distribution Norm(0, sigma_eps*I_n). 
#### In summary: we need estimate parameters
#### mu_a: 2-by-1 vector
#### variance_a: 2-by-1 vector
#### sigma_b: 1-by-1 
#### sigma_eps: 1-by-1

## Step 1. prepare priors for parameters a, sigma_b and sigma_eps
### for parameter a=[a1,a2], its prior is assumed to be normal distribution with
mu_a= fit.spm$fit$coef$fixed
variance_a = c(fit.spm$fit$varFix[1,1]*1e+4, fit.spm$fit$varFix[2,2]*1e+4) 
precision_a = 1/variance_a

#### for sigma_b, prior is a normal distribution with
mean_b = 0                        
variance_b = 1e+10                
precision_b = 1/variance_b

#### for sigma_eps, its prior is a normal distribution with
mean_eps = 0
variance_eps =1e+10
precision_eps = 1/variance_eps

## Step 2. Prepare JAGs dataset
jags.data = list(n = n, K= K, X = X, Y = Y, eGFR = eGFR, msqrt_inv= msqrt_inv, mu_a = mu_a, precision_a = precision_a, precision_eps = precision_eps, precision_b = precision_b)
 
 
## Step 3. Prepare JAGs model
source("PrepareJAGsModel.r") 


## Step 4: Prepare initial values for parameters a, sigma_b, sigma_eps
### We can have multiple sets of initial values to run multiple Markov Chain Monte Carlo (MCMC) chains.
### Here we will run 3 MCMC chains.
jags.inits.1 <- list( a = mu_a, sigma_b = 1, sigma_eps = 2) 
jags.inits.2 <- list( a = mu_a, sigma_b = 2, sigma_eps = 1)
jags.inits.3 <- list( a = mu_a, sigma_b = 3, sigma_eps = 2)

jags.inits =list(jags.inits.1, jags.inits.2, jags.inits.3)

## Step 5: Prepare setting for running JAGs. 
adaptSteps = 100              # Number of steps to "tune" the samplers.
burnInSteps = 2e+4            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
thinSteps= 50                 # Number of steps to "thin" 
nIter = 1e+6                  # Number of iteration for each chain

## Step 6: Now we are ready to run JAGs :)
parameters = c("a", "u", "b", "sigma_eps", "sigma_b")     # The parameter(s) to be monitored.

### randomly selected a seed for replication in future
my.seed = sample(1:1e+6,1)          
print(my.seed)
set.seed(my.seed)

### Create, initialize, and adapt the model:
jagsModel = jags.model(file="modelfile_prior_normal.txt", data=jags.data, inits=jags.inits, n.chains=nChains, n.adapt=adaptSteps)

### Burnin and update model 
print("burning in the MCMC chain...")
update(jagsModel, n.iter=burnInSteps )

### The saved MCMC chain:
print("sampling final MCMC chain...")
bayes.mod.fit.mcmc = coda.samples(jagsModel, variable.names=parameters, n.iter=nIter , thin=thinSteps )

## Step 7: summarize results in screen and save MCMC results to a local folder
### 7.1 Take a look at summary results
summary(bayes.mod.fit.mcmc)

### 7.2 Create a local folder called "Results" for storing output results
out.dir=paste(getwd(), "/Results", sep="")
dir.create(file.path(out.dir), showWarnings = FALSE)

#### 7.3 Save summary of posterior for each parameter monitored to a csv file
summary.stats.file = "./Results/summary_bayes_statistics_.csv"
write.csv(summary(bayes.mod.fit.mcmc)$statistics, summary.stats.file, row.names=T)
	 
#### 7.4 Save quantiles of posterior for each parameter monitored to a csv file
summary.quantiles.file = "./Results/summary_bayes_quantiles.csv"
write.csv(summary(bayes.mod.fit.mcmc)$quantiles,summary.quantiles.file, row.names=T)


### 7.5 Combine 3 chains into one and save it as a matrix.
### Each chain has been stored nIter/thinSteps=2e+4 simulated values for parameters monitored.
mcmc.combi = bayes.mod.fit.mcmc[[1]]
for(c in 2:nChains){
	  mcmc.combi =rbind(mcmc.combi, bayes.mod.fit.mcmc[[c]])
}


## Step 8: Convergence Diagnosis
ConvergenceDiag(bayes.mod.fit.mcmc, out.dir) 

## Step 9: gelman.diagnosis table
gelman.diagnosis.file=paste(out.dir, "/gelman.diagnosis.csv", sep="")
gelman.diagnosis=gelman.diag(bayes.mod.fit.mcmc, multivariate = FALSE)$psrf
write.csv(gelman.diagnosis, gelman.diagnosis.file, row.names=T)

### if(max(gelman.diagnosis[1:8,1]) > 1.1) print("Error! Bayesian MCMC failed for this patient!")


## Step 10: Check posterior distribution with prior distribution for parameters 
PosteriorPlot(mcmc.combi, out.dir)


#----------------------End: Running row-rank thin-plate spline regression under Bayesian framework (JAGs)-----------------------------------------------------

