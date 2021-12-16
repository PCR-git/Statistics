## STATS 201B
## Homework 4
## Peter Racioppo

# Datasets used in this homework are available at
# https://web.stanford.edu/~hastie/CASI/data.html

## -----------------------------------
## Problem 1
## (1.1) Fit a linear model to the supernova data and recreate Figure 12.1.

# Loading the data:
library(readr)
supernova <- read_table2("~/Desktop/supernova.txt")

# Processing the data:
y = supernova$Magnitude
X = as.matrix(supernova)
n = dim(X)[1]
m = dim(X)[2]
X<-X[,-1]
X<-X[,-(m-1)]
X = matrix(as.numeric(X),nrow=n)

library(glmnet)
# This function computes a LASSO regression,
# using cross validation to choose lambda.
f_Lasso <- function(X,y,nf=13){
  mod_cv <- cv.glmnet(X, y, nfolds=nf)
  # lambd = mod_cv$lambda.min
  lambd = mod_cv$lambda.1se
  coef = coef(mod_cv, mod_cv$lambda.min)
  y_hat = predict(mod_cv, newx=X, s=lambd)
  err = (1/n)*sum((y-y_hat)**2) 
  return(list("y_hat" = y_hat, "err" = err, "mod_cv" = mod_cv, "lambd" = lambd, "coef"= coef))
}

# This function plots y_i vs y_i_hat and
# prints apparent mean squared error.
f_Plot <- function(y,y_hat,err){
  plot(y_hat,y,pch=21,col='blue',bg = "blue",xlab="Predicted Magnitude y_i_hat",ylab="Absolute Magnitude y_i",cex.lab=1.5,cex.axis=1.5)
  legend(1.5, -3, legend=c("Apparent Mean Squared Error =", round(err,3)))
  grid()  
}

Lasso = f_Lasso(X,y) # Compute LASSO
y_hat = Lasso$y_hat # Get y_hat
err = Lasso$err # Get error
coef = Lasso$coef # Get coefficients
mod_cv = Lasso$mod_cv # Get model
lambd1 = Lasso$lambd # Get optimal lambda
f_Plot(y,y_hat,err) # Plot

## -----------------------------------
## (1.2) Remove the five predictors with the smallest (in absolute value)
## regression coefficients, and refit the supernova data using just the
## five remaining predictors.

# This function removes smallest 5 predictors:
f_Remove <- function(X,coef,num = 5){
  coef2 = as.matrix(coef)
  d = dim(coef2)[1]
  ind = order(coef2)[(d-(num-1)):d]
  # coef2[coef2==0] <- NA
  # ind = which(is.na(coef2))
  # Remove 5 smallest coefficients:
  # coef3 = coef2[-ind]
  # Remove corresponding columns of X
  X2 = X[,-ind]
  return(X2)
}

X2 = f_Remove(X,coef) # Remove 5 smallest predictors
Lasso2 = f_Lasso(X2,y) # Compute LASSO
y_hat2 = Lasso2$y_hat # Get y_hat
err2 = Lasso2$err # Get error
mod_cv2 = Lasso2$mod_cv # Get model
lambd2 = Lasso2$lambd # Get optimal lambda
f_Plot(y,y_hat2,err2) # Plot

## -----------------------------------
## (1.3) Compare the two fits in terms of the apparent mean squared error.

# Ratio of MSEs
rato = err2/err
print(err)
print(err2)
print(rato)

# The second model has a smaller MSE.

## -----------------------------------
## (1.4) Compute the leave-one-out cross-validated error (12.21) for
## the two models.

# This function runs f_LASSO using an input
# model and leaving out the ith rows of X & y.
f_LASSO_LOO <- function(modl,X,y,lambd,i){
  Xm = X[-i,]
  ym = y[-i]
  ym_hat = predict(modl, newx=Xm, s=lambd)
  err = mean((ym-ym_hat)**2)
  return(err)
}

# This function computes the leave-one-out
# cross-validation estimate of prediction error.
f_CVE <- function(n,modl,X,y,lambd){
  err_cv = 0
  for(i in 1:n){
    err_cv = err_cv + f_LASSO_LOO(modl,X,y,lambd,i)
  }
  return(err_cv/n)
}

Err_cv_1 = f_CVE(n,mod_cv,X,y,lambd1)
Err_cv_2 = f_CVE(n,mod_cv2,X2,y,lambd2)

print(Err_cv_1)
print(Err_cv_2)

## -----------------------------------
## (1.5) Recompute the cross-validated error following the "remove five smallest"
## rule at each cross-validation step. How does the cross-validated error change?

## ????????????????????????

## -----------------------------------
## -----------------------------------
## Problem 2
## Use the bootstrap to compute the degrees of freedom as in Figure 12.4 and
## (12.65), but now for lowess(x, y, 1/6). Still use sigma^2 = 3.28.

## SURE:

# Loading the data:
library(readr)
kidney <- read_table2("~/Desktop/kidney.txt")

# Processing the data:
x = as.matrix(kidney)[,1]
y = as.matrix(kidney)[,2]

fac = 1/3
# Compute lowess fit:
lowessv <- lowess(x, y, fac)
mu = lowessv$y

# Plot lowess fit:
plot(x,y,pch=21,col='blue',bg = "blue",main="Lowess Fit",xlab="x",ylab="y",cex.lab=1.5,cex.axis=1.5)
lines(x,mu)
grid()

# # Plot fit vs actual values:
# plot(mu,y)
# lines(mu,y)
# grid()

# Compute finite differences:
n = dim(array(y))
dy = y[2:n] - y[1:(n-1)]
dmu = mu[2:n] - mu[1:(n-1)]
dmu_dy = dmu/dy
dmu_dy= c(dmu_dy,(array(dmu_dy))[n-1])

# # Replace outliers:
# m = 10*abs(mean(dmu_dy))
# ind = which(dmu_dy > m)
# dmu_dy[ind] = m
# ind = which(dmu_dy < -m)
# dmu_dy[ind] = -m

# Compute the SURE formula (Stein's unbiased risk estimator):
SURE = sum(dmu_dy)

## ----------------
## Bootstrap:

sd = sqrt(3.28) # Standard dev.
B = 1000 # Number of bootstrap replications

# Compute random normal samples of y:
Y = matrix((1:(n*B))*0,nrow=B)
for(j in 1:n){
  Y[,j] = rnorm(B, mean = mu[j], sd = sd)
}

# Compute lowess fit for each y:
Mu = matrix((1:(n*B))*0,nrow=B)
for(i in 1:B){
  # Compute lowess fit:
  lowessr <- lowess(x, Y[i,], fac)
  Mu[i,] = lowessr$y
}

# Compute means across bootstrap samples:
Y. = array(colSums(Y))/B
Mu. = array(colSums(Mu))/B
Y.M = t(matrix(rep(Y.,B),nrow=dim(Y.)))
Mu.M = t(matrix(rep(Mu.,B),nrow=dim(Mu.)))

# Compute degrees of freedom (eqs. 12.64 and 12.65):
df_v = colSums((Y-Y.M)*(Mu-Mu.M))/B/(sd**2)
df = sum(df_v)

# Plot age vs df estimate, with SURE and bootstrap estimates:
plot(x,dmu_dy,ylim=c(-0.2,0.4),col="blue",xlab="age",ylab="df estimate",cex.lab=1.5,cex.axis=1.5)
lines(x,df_v)
legend(35, 0.4, legend=c("Degrees of Freedom", "SURE =", round(SURE,2), "Bootstrap =", round(df,2)))
grid()

# # Non-parameteric Bootstrap:
# # (sample x and y, with replacement)
# L = 1:n
# B = 5
# Xr = matrix((1:(n*B))*0,nrow=B)
# Yr = matrix((1:(n*B))*0,nrow=B)
# Mu = matrix((1:(n*B))*0,nrow=B)
# for(i in 1:B){
#   Lr = sort(replicate(1,sample(L,n,replace=TRUE)))
#   Xr[i,] = array(x[Lr])
#   Yr[i,] = array(y[Lr])
#   
#   # Compute lowess fit:
#   lowessr <- lowess(Xr[i,], Yr[i,], fac)
#   Mu[i,] = lowessr$y
# }

## -----------------------------------
## -----------------------------------
## Problem 4
## Carry out 1000 Gibbs sampling steps (13.70)-(13.71), as in the blue
## histogram of Figure 13.5.

x = array(c(85,88,88,90,90,93,104,108,110,111,115,120,126,126,128,136,143,151,154,157))
n = length(x)
xb = mean(x)
T1 = var(x)
v = (n-1)/2

B = 1000 # Number of samples
# Prior specifications:
n0 = 1
np = n0 + n
k1 = 1
mu0 = xb
mup = (n0*mu0 + n*xb)/np

# Distributions:
# mu ~ N(mup,tau/np)
# tau ~ Q/G(v+k1+2)

library(Rlab)
library(MCMCpack)

## Computes G():
f_G <- function(v){
  return(rinvgamma(n=1,shape=v,scale=1))
}

## Computes mu:
f_mu <- function(tau,mup,np){
  ## mu ~ N(mup,tau/np)
  return(rnorm(1,mean=mup,sd=(tau/np)))
}

## Computes tau:
f_tau <- function(mu,v,k1,T1,np,xb){
  ## tau ~ Q/G(v+k1+2)
  Q = (v + k1)*T1 + np*(mu-xb)**2
  G = f_G(v+k1+2)
  return(Q*G)
}

b = 10
mu_v = rep(0,b)
tau_v = rep(0,b)
tau = T1
tau_v[1] = tau
for(i in 1:b){
  mu = f_mu(tau,mup,np)
  tau = f_tau(mu,v,k1,T1,np,xb)
  mu_v[i] = mu
  tau_v[i+1] = tau
}
print(mu_v)
print(tau_v)

## -----------------------------------
## -----------------------------------
## Problem 6
## The histogram in Figure 15.9 uses 49 bins, equally spaced between -4 and 4.4.

## Part (1). Compute the histogram counts yi, i = 1, 2, . . . , 49.

library(readr)
# DTI <- read_table2("~/Desktop/DTI.txt")
DTI <- read_csv("~/Desktop/DTI.txt",col_names = FALSE)
DTI = as.matrix(DTI)

# Remove values less than -4:
while(min(DTI) < -4){
  ind = which.min(DTI)
  DTI = DTI[-ind]  
}

# Histogram:
library(OneR)
z <- seq(-4,4.4,l=50)
hist1 <- hist(x=DTI,breaks=z,cex.lab=1.5,cex.axis=1.5,xlab="z-score")
grid()

# -----------------
## Part (2). Fit a Poisson regression glm(y ~ poly(x, 6), Poisson),
## with x the bin center.

# x = (z + (z[2] - z[1])/2)[1:(length(z)-1)] # The bin centers
x = hist1$mids # Mid points of bins on x-axis
y = hist1$counts # Histogram counts
poi.mod <- glm(y ~ poly(x,6), family = "poisson") # Fit Poisson model

# Plot Poisson fit:
y_hat = poi.mod$fitted.values
lines(x,y_hat,lwd=2)

# -----------------
## Part (3). Compute the Poisson deviance residuals (8.41).
## Do you think the fit is satisfactory?

# Deviance:
dev_poi = deviance(poi.mod)
print(dev_poi)

# Compare a Gaussian model:
normal.mod <- fitdistr(DTI,"normal")
para = normal.mod$estimate
scal = sum(hist1$counts)/sum(hist1$density)
curve(scal*dnorm(x, para[1], para[2]),col=2,add=TRUE,lwd=2,lty=2)
y_hat2 = scal*dnorm(x, para[1], para[2])

# Manual Deviance Computations:

# Poisson Deviance:
f_dev_poi <- function(mu1,mu2){
  return(2*mu1*((mu2/mu1-1)-log(mu2/mu1)))
}
S_poi = sum(na.omit(f_dev_poi(y,y_hat)))

# Normal Deviance:
f_dev_normal <- function(mu1,mu2,sig){
  return(((mu1-mu2)/sig)**2)
}
S_norm = sum(f_dev_normal(y,y_hat2,sqrt(para[2])))

# Compare deviance of Poisson and Normal fits:
print(S_poi/S_norm)

# Compare MSEs:
MSE_poi = sum((poi.mod$residuals)**2)/length(poi.mod$residuals)
MSE_norm = sum((y-estim)**2)/length(resid)
print(MSE_poi/MSE_norm)

# df:
print(poi.mod$df.residual)

# -----------------
## Part (4). Plot the equivalent of Figure 15.6.

plot(x=x,y=log(y),xlab="z-value",ylab="log density",cex.lab=1.5,cex.axis=1.5)
grid()

cent = floor(length(y)/2)
intv = 5
polyfit <- lm(log(y)[2:length(y)] ~ poly(x[2:length(x)],4,raw=TRUE))
lines(x[2:length(x)],polyfit$fitted.values,col="black",lwd=2,lty=1)
coef = quadfit$coefficients
quadfit <- lm(log(y)[(cent-intv):(cent+intv)] ~ poly(x[(cent-intv):(cent+intv)],2,raw=TRUE))
quad_vals = -(1/2)*z**2 + coef[1] 
lines(z,quad_vals,col="red",lwd=2,lty=2)
legend(x=-3,y=1,c("4th degree log polynomial","Log null density"),col=c("black","red"),lty=1:2,cex=1.3)

# -----------------
## Part (5). Apply locfdr and comment.

library(locfdr)
# "Compute local false discovery rates, following the definitions
#  and description in references listed below."
fdr = locfdr(DTI)
grid()

# locfdr appears to fit a Poisson model and a normal model,
# and the results are very similar to mine. The purple shade
# section on the right appears to denote the difference in
# the predictions of the two models.

## -----------------------------------
## -----------------------------------
## Problem 7

## Run a simulation to compare the df of best-subset regression and lasso.
## Use p = 30 variables and n = 200 observations to build an X matrix,
## generated from a multivariate Gaussian distribution with non-trivial
## covariance (of your choice). Now pose a response model y = X*beta + eps
## and specify beta in advance. In your simulations hold X and beta fixed, and
## generate new eps at each run. Make a plot similar to the right plot in
## Figure 16.8.

# Generate the data:
n = 200
p = 30
sig = diag(p)*0.8
rand = matrix(runif(p**2,0,0.2),nrow=p)
temp = sig+rand
sig2 = (temp + t(temp))/2

library(mvtnorm)
X = rmvnorm(n, rep(0,length=p), sig2)
beta = runif(30,0,1)

sz = 0.2
eps = runif(200,-sz,sz)
y = X%*%beta + eps

# --------

# Best-subset Regression
library(leaps)
models <- regsubsets(X,y,nvmax=p)
res.sum <- summary(models)
bsr_resid = res.sum$rss
bsr_step = 1:p
bsr_df = 1:p

# LASSO Regression:
library(glmnet)
num = 30
lasso_step = 1:num
# lambd <- 10^seq(10,-2,length=num)
lambd = seq(3,0,length=num)
cv.out <- cv.glmnet(X,y,alpha=0)
lamp <- cv.out$lambda.min
lasso.mod <- glmnet(X,y,alpha=1,lambda=lambd)
lasso_df = lasso.mod$df
# print(lasso_df)

# Compute LASSO residuals for each step
lasso_resid = rep(0,length(lambd))
for(i in 1:length(lambd)){
  lasso.pred <- predict(lasso.mod,s=lambd[i],newx=X)
  lasso_resid[i] = sum((lasso.pred-y)^2) 
}

# Plot Sum of Square Residuals vs Degrees of Freedom
plot(bsr_df,bsr_resid,ylim=c(0,5500),xlab="df",ylab="Sum Squared Residuals",cex.lab=1.5,cex.axis=1.5)
lines(bsr_df,bsr_resid,col="red")
lines(lasso_df,lasso_resid,type="o",col="blue",lty=2)
legend(x=20,y=5250,c("LASSO","Best-subset"),col=c("red","blue"),lty=1:2,cex=1.3)
grid()

# Plot Degrees of Freedom vs Steps
plot(lasso_step,lasso_df,xlab="Step",ylab="df",cex.lab=1.5,cex.axis=1.5)
lines(lasso_step,lasso_df,col="red")
lines(bsr_step,bsr_df,type="o",col="blue",lty=2)
legend(x=2.5,y=30,c("LASSO","Best-subset"),col=c("red","blue"),lty=1:2,cex=1.3)
grid()

## -----------------------------------
## -----------------------------------

