## STATS 202B
## Homework 3
## Peter Racioppo

## -------------------
## Question 2

# Define matrix:
A <- matrix(c(10,7,8,7,7,6,6,5,8,6,10,9,7,5,9,10),nrow=4)

## Vector norm:
norm_vec <- function(x){sqrt(sum(x^2))}

## Power Method (one iteration):
f_Power_j <- function(A,tol){
  z = runif(4, min=-1, max=1)
  Az = z
  xj_1 = 0
  for (j in 1:10){
    Az = A%*%Az
    xj = Az/norm_vec(Az)
    if (norm_vec(xj-xj_1)<tol){
      break}
    xj_1 = xj
    xjList <- list("xj" = xj, "iter" = j)
  }
  return(xjList)
}

## Matrix Deflation:
f_Deflate_Matrix <- function(A,x){
  A = A - x%*%t(x)%*%A%*%x%*%t(x)
  return(A)
}

## Power Method:
f_Power_Method  <- function(A,tol){
  n = dim(A)[1]
  Aj = A
  X <- matrix((1:16)*0,nrow=n)
  iter <- c(1:n)*0
  for (i in 1:n){
    xjList = f_Power_j(Aj,tol)
    xj = xjList$xj
    Aj = f_Deflate_Matrix(Aj,xj)
    X[,i] = xj
    iter[i] = xjList$iter
  }
  eigenvecs = X
  eigenvals = colSums((A%*%X)/X)/n
  eigenList <- list("vals" = eigenvals, "vecs" = eigenvecs, "iter" = iter)
  return(eigenList)
}

# Run Power Method:
tol = 0.001
eigenList = f_Power_Method(A,tol)
print("Eigenvalues:")
print(eigenList$vals)
print("Eigenvectors:")
print(eigenList$vecs)
print("Iterations:")
print(eigenList$iter)

# Use inbuilt functions to check our function:
eigen = eigen(A)
# print(eigen$values)
# print(eigen$vectors)

Diff = abs((eigenList$vals - eigen$values)/eigen$values)
print("Percent Errors:")
print(Diff)

## -------------------
## Question 3:

table  = read.table("http://www.stat.ucla.edu/~handcock/202B/datasets/SwissBankNotes.txt", skip=22,header=TRUE)
# X = as.matrix(table)
X = table
X1 = X[1:100,]
X2 = X[101:200,]
dims = dim(X)
# Shift sample mean to zero:
means = array(colSums(X)/dims[1])
X = sweep(X,2,means,FUN="-")
pca_1 = prcomp(X1, scale = TRUE)
pca_2 = prcomp(X2, scale = TRUE)
pca_tot = prcomp(X, scale = TRUE)
# pca_1 = princomp(X1, scale = TRUE)
# pca_2 = princomp(X2, scale = TRUE)
# pca_tot = princomp(X, scale = TRUE)

summary(pca_1,loadings=TRUE)
summary(pca_2,loadings=TRUE)
summary(pca_tot,loadings=TRUE)

plot(c(1:6),pca_1$sdev)
plot(c(1:6),pca_2$sdev)
plot(c(1:6),pca_tot$sdev)

# vecs1 = pca_1$rotation
# vecs2 = pca_2$rotation
# vecs_tot = pca_tot$rotation
# 
# print("Real PCA Vectors:")
# print(vecs1)
# print("Counterfeit PCA Vectors:")
# print(vecs2)
# print("Combined PCA Vectors:")
# print(vecs_tot)

## -------------------
## Question 4

library(MMST)
library(tls)
data(tobacco)
help(tobacco)
head(tobacco)

table <- as.data.frame(tobacco[,c(1,4:9)])
tls_model <- tls(Y1.BurnRate~X1.PercentNitrogen+X2.PercentChlorine+
                 X3.PercentPotassium+X4.PercentPhosphorus+
                 X5.PercentCalcium+X6.PercentMagnesium-1,data=table)

lm_model <- lm(Y1.BurnRate~X1.PercentNitrogen+X2.PercentChlorine+
                   X3.PercentPotassium+X4.PercentPhosphorus+
                   X5.PercentCalcium+X6.PercentMagnesium,data=table)

tls_model
summary(lm_model)

tls_model$coefficient
lm_model$coefficients

## -------------------
## Question 6

library(MASS)
data(pet)
help(pet)
head(pet)

y = as.matrix(pet$y)
x = as.matrix(pet[,1:268])

Ridge = lm.ridge(y ~ x, lambda=seq(0,1,0.1))
# Ridge = lm.ridge(y ~ x, lambda=seq(0,1,0.1),auto=T)
Rp <- Ridge$coef[1:60,]
Ridge$coef = Rp
plot(Ridge)

## -------------------
## Question 7

# library(matrixStats)
#
# # Shift sample mean to zero:
# means = array(colSums(x)/28)
# x = sweep(x,2,means,FUN="-")
# means2 = array(colSums(x)/28)
# print(means2)
# 
# # Divide columns by standard devs:
# for (i in 1:268){
#   x[,i] = x[,i]/sd(x[,i])
# }
# 
# rand_int <- sample(1:28, 28,replace=F)
# print(rand_int)
# 
# L = 28
# d = 4
# n = L/d
# x_shuffle = x*0
# y_shuffle = y*0
# for (i in 1:L){
#   x_shuffle[i,] = x[rand_int[i],]
#   y_shuffle[i,] = y[rand_int[i],]
# }
# 
# x_train = x_shuffle[1:(d-1)*n,]
# x_test = x_shuffle[((d-1)*n+1):L,]
# y_train = y_shuffle[1:(d-1)*n,]
# y_test = y_shuffle[((d-1)*n+1):L,]
# 
# Ridge = lm.ridge(y_train ~ x_train, lambda=0.1)

library(glmnet)

cv.out = cv.glmnet(x,y, lambda = seq(0,1,0.1), nfolds = 7)
cv.out$glmnet.fit
cv.out$lambda
cv.out$lambda.min
cv.out$lambda.1se
plot(cv.out)

## -------------------
