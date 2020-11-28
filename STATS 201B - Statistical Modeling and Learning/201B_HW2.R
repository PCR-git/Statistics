## STATS 201B
## Homework 2
## Peter Racioppo

## -----------------------------------
## Question 2

# Data:
x = c(0:7)
yx = c(7840, 1317, 239, 42, 14, 4, 4, 1)

# Compute E_hat:
E_hat = (1:(length(yx)-1))*0
for (i in 1:(length(yx)-1)){
  E_hat[i] = (i)*yx[i+1]/yx[i] 
}
print(E_hat)

## Compute initial guesses for Gamma distrb.
## parameters (nu and sigma):
## (didn't work)
# mean = mean(yx)
# var = var(yx)
# shape1 = mean**2/var
# scale1 = var/mean
# factor = scale1/shape1
# library(MASS)
# m <- fitdistr(yx/factor, "gamma")
# rate = as.matrix(1/coef(m)["rate"]*factor)[1]
# shape =  as.matrix(1/coef(m)["shape"]*factor)[1]
# scale = 1/rate
# nu = shape
# sigma = scale

LL <- function(sigma,nu){
  x = c(0:7)
  yx = c(7840, 1317, 239, 42, 14, 4, 4, 1)
  gammer = sigma/(1+sigma)
  f = (gammer**(nu+x))*gamma(nu+x)/((sigma**nu)*gamma(nu)*factorial(x))
  LL = -sum(yx*log(f))
  return(LL)
}

library(stats4)
MLE = mle(LL, start = list(sigma=1, nu=1))
sigma = 0.3055760
nu = 0.7014887

f_nu_sigma <- function(nu,sigma,x){
  gammer = sigma/(1+sigma)
  return(((gammer**(nu+x))*gamma(nu+x))/((sigma**nu)*gamma(nu)*factorial(x))) 
}

E_gamma = (x+1)*f_nu_sigma(nu,sigma,x+1)/f_nu_sigma(nu,sigma,x)
print(E_gamma)

## -----------------------------------
## Question 3

x_12 = c(1:12)
y_12 = c(118,74,44,24,29,22,20,19,20,15,12,14)
x_24 = c(1:24)
y_24 = c(118,74,44,24,29,22,20,19,20,15,12,14,6,12,6,9,9,6,10,10,11,5,3,3)

# Calculates expected number of new butterflies
E_Butt <- function(x,y,t){
  E_hat = 0
  for (i in 1:length(x)){
    E_hat = E_hat - y[i]*((-1)**i)*(t**i)
  }
  return(E_hat)
}

sd_Butt <- function(x,y,t){
  i = 1:length(y)
  sd = 0*(1:length(t))
  for (j in 1:length(t)){
    sd[j] = sqrt(sum(y*t[j]**(2*i)))
  }
  return(sd)
}

t = c(0:10)/10
E_hat_24 = E_Butt(x_24,y_24,t)
E_hat_12 = E_Butt(x_12,y_12,t)
print(E_hat_24)
print(E_hat_12)

test = sd_Butt_j(x_24,y_24,0.5)

sd_24 = sd_Butt(x_24,y_24,t)
sd_12 = sd_Butt(x_12,y_12,t)
print(sd_24)
print(sd_12)

# plot(t,E_hat_24,type='l',col='blue')
# lines(t,E_hat_12,col='red',lty=2)
# grid()

## -----------------------------------
## Question 5
# Load and process the data:
baseball_csv <- read.csv(url("https://web.stanford.edu/~hastie/CASI_files/DATA/baseball.txt"))
tbl = data.table(baseball_csv)
N = dim(tbl)[1]
mat = as.numeric(unlist(strsplit(as.matrix(tbl), " ")))
dim2 = length(mat)/N
baseball = matrix(mat, nrow = N, byrow = TRUE)
player = baseball[,1]
p_i = baseball[,2]
truth = baseball[,3]

p_hat_JS <- function(p_i,n,N){
  p_bar = mean(p_i)
  sigma_0 = (p_bar*(1-p_bar)/n)**0.5
  term1 = ((N-3)*sigma_0**2)/sum((p_i-p_bar)**2)
  p_hat_js = p_bar + (1-term1)*(p_i-p_bar)
  return(p_hat_js)
}

p_JS_90 = p_hat_JS(p_i,90,N)
p_JS_180 = p_hat_JS(p_i,180,N)

plot(player,p_i,col='black',type='l',main="Comparison of MLE and Jules Stein Estimates",xlab="Player",ylab="Estimated Batting Average")
lines(player,p_JS_90,col='blue',lty=2)
lines(player,p_JS_180,col='red',lty=4)
grid(); box()
legend(x=12,y=0.325,legend=list('MLE','JS 90','JS 180'),lty=c(1,2,4))

MLE_error = sum((p_i-truth)**2)
JS_90_error = sum((p_JS_90-truth)**2)
JS_180_error = sum((p_JS_180-truth)**2)

## -----------------------------------
## Question 6

csv = read.csv('/Users/peterracioppo/Desktop/diabetes.csv', header = TRUE)
library(data.table)
tbl = data.table(csv)
y = as.matrix(tbl$prog)
y = y - mean(y)
tbl[,prog:=NULL]
x = as.matrix(tbl)
# Standardize x:
for (i in 1:dim(x)[2]){
  x[,i] = (x[,i] - mean(x[,i]))/sd(x[,i])
}

library(glmnet)
lambda = 0
ridge_coeff <- function(lambda,x,y,alpha){
  ridge.mod <- glmnet(x, y, alpha = 0, lambda = lambda, intercept=F)
  coeff = predict(ridge.mod, s = 0, type = 'coefficients', se.fit=T)
  print(coeff)
}

ridge_coeff(0,x,y,alpha)
ridge_coeff(0.1,x,y,alpha)
ridge_coeff(0.2,x,y,alpha)

ols<-lm(y~as.matrix(x))
anova(ols)
mse = 2933
se<-sqrt(mse)
sigmaatrix<-se*solve(t(as.matrix(x))%*%as.matrix(x)+lambda*diag(10))%*%(t(as.matrix(x))%*%as.matrix(x))%*%solve(t(as.matrix(x))%*%as.matrix(x)+lambda*diag(10))

#standard error for betas
for(i in 1:10){
  print(sigmaatrix[i,i])
}

## -----------------------------------
## Question 7

n = 10
y_i = c(0, 0, 0, 3, 6, 6, 5, 9, 9, 10, 10)
x_i = 1:length(y_i)
pi_i = y_i/n
lambda_i = log(pi_i/(1-pi_i))
fit <- glm(pi_i~x_i,family=binomial())
alpha0 = fit$coefficients[1]
alpha1 = fit$coefficients[2]

y_i2 = c(0.1, 0.1, 0.1, 3, 6, 6, 5, 9, 9, 9.9, 9.9)
pi_i2 = y_i2/n
lambda_i2 = log(pi_i2/(1-pi_i2))
lin_fit = lm(lambda_i2~x_i)
alpha0_2 = coefficients(lin_fit)[1]
alpha1_2 = coefficients(lin_fit)[2]

f_pi_hat <- function(alpha0,alpha1,x){
  pi_hat = 1/(1+exp(-(alpha0+alpha1*x)))
  return(pi_hat)
}

x_v = seq(0,11,0.1)
pi_hat = f_pi_hat(alpha0,alpha1,x_v)
pi_hat_2 = f_pi_hat(alpha0_2,alpha1_2,x_v)

plot(x_i,pi_i,xlab="Dose",ylab="Proportion of Deaths")
lines(x_v,pi_hat,type='l',col='red')
# lines(x_v,pi_hat_2,type='l',col='blue',lty=2)
grid(); box()

## -----------------------------------
## Question 8
# Load and process the data:
galaxy_data <- read.csv(url("https://web.stanford.edu/~hastie/CASI_files/DATA/galaxy.txt"))
tbl = data.table(galaxy_data)
dim1 = dim(tbl)[1]
mat = as.matrix(tbl)
mat2 = strsplit(mat, " ")
mat3 = as.numeric(unlist(mat2))
dim2 = length(mat3)/dim1
mat4 = matrix(mat3, nrow = dim1, byrow = TRUE)
counts = mat4[,2:dim2]
counts_flip <- counts[c(dim1:1),]
dim2n = dim(counts_flip)[2]

# Assemble X matrix:
r1 = 1:dim2n
m1 = 1:dim1
R = matrix(rep(r1,dim1),ncol=dim2n,byrow=T)
r = c(R)
m = matrix(rep(m1,dim2n),ncol=1,byrow=T)
r2 = r**2
rm = r*m
m2 = m**2
l = list(r,m,r2,rm,m2)
X = do.call(cbind, l)
# dim(X)
mu = matrix(counts_flip,ncol=1,byrow=T)

# Fit Poisson Distribution:
poisson_fit = glm(mu ~ X, family="poisson")
coeff = as.matrix(poisson_fit$coefficients)

# Calculate predictions:
ones = rep(1, dim(X)[1])
X_ = cbind(ones,X)
lambda_hat = X_%*%coeff
mu_hat = exp(lambda_hat)
counts_hat = matrix(mu_hat,ncol=dim2n,byrow=FALSE)

# Plot:
x_v = 1:dim1
y_v = 1:dim2n
counts_flip <- counts[c(dim1:1),]
persp(x_v, y_v, counts_flip,theta = -45, phi = 25, col='green',main = "Galaxy Data",xlab="Dimmer",ylab="Farther",zlab="Counts")
persp(x_v, y_v, counts_hat,theta = -45, phi = 25, col='green',main = "Poisson GLM Density Estimate",xlab="Dimmer",ylab="Farther",zlab="Counts")

S = deviance(poisson_fit)

summary_fit <- summary(poisson_fit)
resid <-summary_fit$deviance.resid
S1 = sum(resid^2) #S=230.32

f_Poisson_dev <- function(mu_1,mu_2){
  return(2*mu_1*((mu_2/mu_1 - 1) - log(mu_2/mu_1)))
}

D = f_Poisson_dev(mu,mu_hat)
D2 <- na.omit(D)
S2 = sum(abs(D2))

plot(1:length(resid),resid,type='l',xlab="Index",ylab="Residual")
grid(); box()

# par(mfrow=c(1,2))
resid_mat = matrix(resid, nrow=18, byrow = FALSE)
persp(x_v, y_v, resid_mat,theta = -45, phi = 25, col='lightblue',main = "Residuals",xlab="Dimmer",ylab="Farther",zlab="Residual",ticktype="detailed")
persp(x_v, y_v, abs(resid_mat),theta = -45, phi = 25, col='lightblue',main = "Absolute Value of Residuals", xlab="Dimmer",ylab="Farther",zlab="Residual",ticktype="detailed")

## The fit seems to become poorer as dimmness increases and distance
## decreases.

## ------------------------------
## Part (d)

# To improve the fit where it was worst, near the boundaries of the table,
# I tried adding columns to the data matrix X, r*m^2, r^2*m^2, r^3, m^3.
# It didn't seem to help 

rm2 = r*m2
r2m2 = r2*m2
r3 = r**3
m3 = m**3
l2 = list(r,m,r2,rm,m2,rm2,r2m2,r3,m3)
X2 = do.call(cbind, l2)

# Fit Poisson Distribution:
poisson_fit2 = glm(mu ~ X2, family="poisson")
coeff2 = as.matrix(poisson_fit2$coefficients)

# Calculate predictions:
ones = rep(1, dim(X2)[1])
X2_ = cbind(ones,X2)
lambda_hat2 = X2_%*%coeff2
mu_hat2 = exp(lambda_hat2)
counts_hat2 = matrix(mu_hat2,ncol=dim2n,byrow=FALSE)

summary_fit2 <- summary(poisson_fit2)
resid2 <-summary_fit2$deviance.resid
resid_mat2 = matrix(resid2, nrow=18, byrow = FALSE)

# Plot:
persp(x_v, y_v, counts_hat2,theta = -45, phi = 25, col='red',main = "Poisson GLM Density Estimate",xlab="Dimmer",ylab="Farther",zlab="Counts",ticktype="detailed")

persp(x_v, y_v, resid_mat2,theta = -45, phi = 25, col='orange',main = "Residuals",xlab="Dimmer",ylab="Farther",zlab="Residual",ticktype="detailed")
persp(x_v, y_v, abs(resid_mat2),theta = -45, phi = 25, col='orange',main = "Absolute Value of Residuals", xlab="Dimmer",ylab="Farther",zlab="Residual",ticktype="detailed")

