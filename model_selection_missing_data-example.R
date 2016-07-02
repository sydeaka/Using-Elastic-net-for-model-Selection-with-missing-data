###############################################################################
#### Model Selection with Missing Data and Correlated Covariates
####
#### In this script, I demonstrate the use of elastic net to facilitate
#### model selection given a large number of highly correlated candidate covariates
#### with some covariate data missing.
#### 
#### Details about the method are provided in published manuscript, 
####  "Survival in pulmonary arterial hypertension patients awaiting lung transplantation.",
####   http://www.ncbi.nlm.nih.gov/pubmed/24074527 
#### 
#### Written for R version 3.2.1, but will likely work for other versions.
#### 
#### Note: Given a high performance computing platform, it is trivial to 
####  run the 12 computations (4 imputations x 3 values of alpha) in parallel.
####  The example is written so that it can be easily run on a single computer with
####  4 processors, even on a Windows machine.
#### 
#### Author: Sydeaka Watson, PhD
#### Institution: The University of Chicago, Biostatistics Laboratory
#### Contact: sydeakawatson@gmail.com
###############################################################################

## Set the random seed, load packages
set.seed(102)
library(mvtnorm)
library(mice)
library(glmnet)
library(boot)
library(doSNOW)
library(foreach) 

npreds = 20 # number of candidate predictors
n = 250 # sample size

## Pairwise correlations, mostly low correlations, a few high correlations 
seqcor = seq(0,.7, .05); nseq = length(seqcor)
seqcor = seqcor*sample(c(-1,1), nseq, replace=T)
prob=c(rep(.9/5, 5), rep((1-.9)/(nseq-5),nseq-5))
corrs = sample(seqcor, npreds*(npreds+1)/2, replace=T, prob=prob)

## Fill in correlation matrix
R = matrix(NA, nrow=npreds, ncol=npreds)
diag(R) = 1
q = 1
for (i in 1:npreds) {
  for (j in (i+1):npreds) {
    if (i!=j & j <= npreds) {
      R[i,j] = R[j,i] = corrs[q]
      q = q + 1
    }
  } # end j loop
} # end i loop


## Find the nearest positive definite matrix 
##  (so that it is a proper correlation matrix)
Q <- nearPD(R, posd.tol=1.e-04, corr=T)$mat

## Generate covariate data using multivariate normal distribution
xnames = paste('x', 1:npreds, sep='')
x0 = rmvnorm(n, mean=rep(0,npreds))
x = as.matrix(x0 %*%  chol(Q))
colnames(x) = xnames
x = as.data.frame(x)


## Induce relationship between first four x's and outcome y
y = with(x, 1.74*x1 - 0.68*x2 - 1.5*x3 + 0.93*x4 + rnorm(n,0,2))

## Function to randomly delete some covariate data
rm.data = function(vec, prop) {
  rm.ind = rbinom(n, 1, prop)
  vec[rm.ind==1] = NA
  return(vec)
}

## Randomly delete some covariate data
(prop.missing = sample(seq(0,.4,.01), npreds, replace=T))
xmiss = sapply(1:npreds, function(u) rm.data(x[,u], prop.missing[u]))
colnames(xmiss) = xnames


## Generate 4 imputations using the MICE algorithm 
## (multiple imputations by chained equations)
imp = mice(data.frame(y=y, xmiss), m=4)
dat = complete(imp, 'long')


## Function to perform elastic net on each bootstrap replicate
enet_fun = function(data, indices) {
  yboot = as.matrix(data[indices, 'y'])
  xboot = as.matrix(data[indices, xnames])
  cv.fit = cv.glmnet(xboot, yboot, alpha=alpha)
  fit = glmnet(xboot, yboot, alpha=alpha)
  return(ifelse(coef(fit,s=cv.fit$lambda.1se)[-1] != 0, 1, 0))
} 






## Selected mixing paramers for the penalty term in the elastic net framework
alphas = c(.2, .5, .65, .75, .85)

## Make a cluster using 4 processors
cl<-makeCluster(4)
registerDoSNOW(cl)
clusterExport(cl,ls())

## 200 bootstrap fits of elastic net model for each imputation, for each value of alpha
out = vector('list', length(alphas))
names(out) = paste('alpha =', alphas)
for (a in 1:length(alphas)) {
  alpha = alphas[a]
  cat('alpha = ', alpha, '...')
(run.time = system.time({out[[a]] = 
  foreach(i=1:4) %dopar% {
    library(glmnet)
    library(boot)
    boot(data=subset(dat, .imp==as.character(i)), statistic=enet_fun, R=200)
  }
}))
  cat('done.\n')
  print(run.time)
} # end a loop


## Stop the cluster 
stopCluster(cl)


## Function to compute proportion of times each variable was selected across imputations
summ = function(result) {
  all.results = NULL
  for (i in 1:4) all.results = rbind(all.results, result[[i]]$t)
  selection.props = apply(result[[1]]$t,2,function(u) round(100*mean(u)))
  return(selection.props)
} # end summ function

## Proportion of times each variable was selected across imputations
##  for each value of alpha
selection.props = t(sapply(out, summ))
colnames(selection.props) = xnames
selection.props

## Variables selected for each value of alpha using 85% selection threshold
##  i.e., choose variables that were selected at least 85% of the time.
apply(selection.props,1,function(u) xnames[u>=85])
