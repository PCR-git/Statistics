## STATS 202B
## Homework 4
## Peter Racioppo

## ===============================
## Question 6
## ===============================
## Part (a)

# Load the data
library(R.matlab)
mnist <- R.matlab::readMat("~/Desktop/mnist.mat")
# mnist_or <- R.matlab::readMat("~/Desktop/mnist_original.mat")
y = mnist$y
X = mnist$X
ytest = mnist$ytest
Xtest = mnist$Xtest
N = 10 # Number of unique labels

# Compute Log-likelihood:
f_LL <- function(X,B,Y){
  num = diag(exp(X%*%B)%*%t(Y))
  P = num/sum(num)
  LL = sum(-log(P))
  return(LL)
}

# Compute Z:
f_Z <- function(m,X,B){
  Z = diag(m)
  for(i in 1:m){
    Z[i,i] = 1/sum(exp(X[i,]%*%B)) 
  }
  return(Z)
}

# Compute one-hot encoding:
f_One_Hot <- function(y){
  nr = dim(y)[2]
  Y = matrix(0,nrow=nr,ncol=10)
  for(k in 1:10){
    yk = y
    yk[yk!=(k-1)]<-11
    yk[yk==(k-1)]<-1
    yk[yk==11]<-0
    Y[,k] = yk
  }
  return(Y) 
}

# Compute the gradient
f_Grad <- function(X,Z,B,Y){
  return(t(X)%*%(Z%*%exp(X%*%B) - Y))
}

# Backtracking line search:
f_Backtrack <- function(a,b,t,X,B,Y,grad){
  lhs = f_LL(X,B-t*grad,Y)
  # print(lhs)
  rhs = f_LL(X,B,Y) - a*t*(norm(grad,"F")**2)
  # print(rhs)
  i = 0
  while(lhs >= rhs){
    t = b*t
    i = i + 1
    if(i>10){break}
  }
  return(t)
}

# Gradient descent:
f_Grad_Descent <- function(t,X,B,Y,it,a,b){
  f = (1:it)*0
  f[1] = f_LL(X,B,Y)
  for(i in 1:it){
    Z = f_Z(1000,X,B)
    grad = f_Grad(X,Z,B,Y)
    t = f_Backtrack(a,b,t,X,B,Y,grad)
    f[i+1] = f_LL(X,B,Y)
    B = B - t*grad
    # t = 0.99*t
    if(i%%100==0){print(i)} # Print iteration #
  }
  return(list("B" = B, "f" = f))
}

# Starting beta
# beta = (1:10)*0
# B = matrix(runif(10000, min = 0, max = 0.5),ncol = 10)
B = matrix(rnorm(10000,0,0.5),ncol=10)

Z = f_Z(1000,X,B)
Y = f_One_Hot(y)
grad = f_Grad(X,Z,B,Y)

# t = 0.01
t = abs(mean(B/grad))/5
a = 0.5 # Alpha
b = 0.5 # Beta
it = 1000 # Iterations

# f = f_LL(X,B,Y)
# f_Backtrack(a,b,t,X,B,Y,grad)
opt = f_Grad_Descent(t,X,B,Y,it,a,b)
B_star = opt$B
f = array(opt$f)
f_0 = f[1]
f_star = f[it]
print("Initial Log-Likelihood:")
print(f_0)
print("Final Log-Likelihood:")
print(f_star)

plot((1:(it+1)),f,type='l',xlab="Iteration",ylab="Minus Log-Likelihood",main="Backstepping Gradient Descent")
grid()

## ===============================
## Part (b)

# Compute column-wise 1-norm of matrix:
f_Mat1Norm <- function(M){
  m = dim(M)[2]
  Mnorm = 0
  for(i in 1:m){
    Mnorm = Mnorm + norm(matrix(M[,i]),"1")
  }
  return(Mnorm)
}

# Regularized Log-likelihood:
f_LL_reg <- function(X,B,Y,lambd){
  num = diag(exp(X%*%B)%*%t(Y))
  P = num/sum(num)
  Bnorm = f_Mat1Norm(B)
  LL = sum(-log(P)) + lambd*Bnorm
  return(LL)
}

# Regularized gradient
f_Grad_reg <- function(X,Z,B,Y,lambd){
  grad_smooth = t(X)%*%(Z%*%exp(X%*%B) - Y)
  Bnorm = f_Mat1Norm(B)
  grad_reg = matrix(0,dim(B)[1],dim(B)[2])
  if(abs(Bnorm) > 1E-7){
    grad_reg = lambd*B/Bnorm
  }
  return(grad_smooth + grad_reg)
}

# Gradient descent with regularized loss:
f_Grad_Descent_reg <- function(t,X,B,Y,it,a,b,lambd){
  f = (1:it)*0
  f[1] = f_LL_reg(X,B,Y,lambd)
  for(i in 1:it){
    Z = f_Z(1000,X,B)
    grad = f_Grad_reg(X,Z,B,Y,lambd)
    t = f_Backtrack(a,b,t,X,B,Y,grad)
    f[i+1] = f_LL_reg(X,B,Y,lambd)
    B = B - t*grad
    # t = 0.99*t
    if(i%%100==0){print(i)} # Print iteration #
  }
  return(list("B" = B, "f" = f))
}

lambd = 0.5 # Test with lambda = 0.5
# f_r = f_LL_reg(X,B,Y,lambd)
# f_Backtrack(a,b,t,X,B,Y,grad)

opt_r = f_Grad_Descent_reg(t,X,B,Y,it,a,b,lambd)
B_r_star = opt_r$B
f_r = array(opt_r$f)
f_r_0 = f_r[1]
f_r_star = f_r[it]
print("Initial Log-Likelihood:")
print(f_r_0)
print("Final Log-Likelihood:")
print(f_r_star)

plot((1:(it+1)),f_r,type='l',xlab="Iteration",ylab="Minus Log-Likelihood",main="Backstepping Gradient Descent (with L1 Column Norm))")
grid()

it = 500
N = 10
library(pracma)
lambd_v = logspace(log10(1E1), log10(1E-2), n = N)

# Run gradient descent with line-search on lambda
f_Lambda <- function(t,X,B,y,Xtest,ytest,it,a,b,lambd_v,N){
  Y = f_One_Hot(y)
  Ytest = f_One_Hot(ytest)
  i = 1
  f_test = (1:N)*0
  for(lambd in lambd_v){
    opt_r = f_Grad_Descent_reg(t,X,B,Y,it,a,b,lambd)
    B_r_star = opt_r$B
    f_test[i] = f_LL_reg(Xtest,B_r_star,Ytest,lambd)
    print(i)
    i = i + 1
  }
  return(f_test)
}

f_test = f_Lambda(t,X,B,y,Xtest,ytest,it,a,b,lambd_v,N)

plot(lambd_v,f_test,log="x",xlab="Lambda",ylab="Testing Error",main="Testing Error vs Lambda (500 Iterations)")
lines(lambd_v,f_test,log="x")
grid()
