############################################################# 
## Stat 202A - Homework 3
## Author: Peter Racioppo
## Date : 10/20/2019
## Description: This script implements linear regression 
## using the sweep operator
#############################################################

#############################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. You can add examples at the
## end of the script (in the "Optional examples" section) to 
## double-check your work, but MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not use the function "setwd" anywhere
## in your code. If you do, I will be unable to grade your 
## work since R will attempt to change my working directory
## to one that does not exist.
##
## Do not use the following functions for this assignment,
## except when debugging or in the optional examples section:
## 1) lm()
## 2) solve()
#############################################################

 
################################
## Function 1: Sweep operator ##
################################

mySweep <- function(A, k){
  
  # Perform a SWEEP operation on A with the pivot element A[k,k].
  # 
  # A: a square matrix.
  # m: the pivot element is A[k, k].
  # Returns a swept matrix B (which is k by k).
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
    
    # ## For-loop implementation:
    # n <- dim(A)[1] # Number of rows of A
    # B = A*0 # Initialize B
    # m = A[k,k]
    # for(i in 1:n){
    #   for(j in 1:n){
    #     if(i!=k){
    #       if(j!=k){
    #         B[i,j] = A[i,j] - (A[i,k]*A[k,j])/m
    #       }
    #       if(j==k){
    #         B[i,j] = -A[i,j]/m
    #       }
    #     }
    #     if(i==k){
    #       if(j!=k){
    #         B[i,j] = A[i,j]/m
    #       }
    #     }
    #   }
    # }
    # B[k,k] = 1/m
  
  ## Matrix implementation (much faster):
  # n <- dim(A)[1] # Number of rows of A
  m = A[k,k] # Pivot element
  pv = 1/m # m inverse
  # Define B as A minus outerproduct of kth row and column divided by pivot:
  B = A - A[,k]%o%A[k,]*pv
  B[k,] = A[k,]*pv # Overwrite kth row
  B[,k] = -A[,k]*pv # Overwrite kth column
  B[k,k] = pv # Overwrite kkth element
  
  ## The function outputs the matrix B
  return(B)
  
}


############################################################
## Function 2: Linear regression using the sweep operator ##
############################################################

myLinearRegression <- function(X, Y){
  
  # Find the regression coefficient estimates beta_hat
  # corresponding to the model Y = X * beta + epsilon
  # Your code must use the sweep operator you coded above.
  # Note: we do not know what beta is. We are only 
  # given a matrix X and a vector Y and we must come 
  # up with an estimate beta_hat.
  # 
  # X: an 'n row' by 'p column' matrix of input variables.
  # Y: an n-dimensional vector of responses

  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
    n <- nrow(X)
    p <- ncol(X)
    
    # Append column of 1s to X
    Ones = matrix(1, n, 1)
    X = cbind(Ones,X)  # Concatenate horizontally
    
    XpX = t(X)%*%X # X_transpose * X
    XpY = t(X)%*%Y # X_transpose * Y
    YpX = t(Y)%*%X # Y_transpose * X
    YpY = t(Y)%*%Y # Y_transpose * Y
    
    # Construct M matrix:
    M1 = cbind(XpX,XpY) # Concatenate horizontally
    M2 = cbind(YpX,YpY) # Concatenate horizontally
    M = rbind(M1,M2)  # Concatenate vertically
    
    # Sweep M along the rows of XpX
    for(i in 1:(p+1)){
      M = mySweep(M, i)}
    
    # beta_hat
    beta_hat = M[1:(p+1),p+2] # Beta_hat is upper-right block of M
    
    # # Check:
    # beta_hat2 = -t(M[p+2, 1:(p+1)]) # Get beta_hat from lower left of M
    # print(beta_hat - beta_hat2) # Should be vector of zeros
    
  ## Function returns the (p+1)-dimensional vector 
  ## beta_hat of regression coefficient estimates
  return(beta_hat)
  
}

########################################################
## Optional examples (comment out before submitting!) ##
########################################################

testing_Linear_Regression <- function(){

  ## This function is not graded; you can use it to
  ## test out the 'myLinearRegression' function

  # Define parameters
  n    <- 100
  p    <- 3

  ## Simulate data from our assumed model.
  ## We can assume that the true intercept is 0
  X    <- matrix(rnorm(n * p), nrow = n)
  beta <- matrix(1:p, nrow = p)
  Y    <- X %*% beta + rnorm(n)

  ## Save R's linear regression coefficients
  R_coef  <- coef(lm(Y ~ X))
  print(R_coef)

  ## Save our linear regression coefficients
  my_coef <- myLinearRegression(X, Y)
  print(my_coef)

  ## Are these two vectors different?
  sum_square_diff <- sum((R_coef - my_coef)^2)
  if(sum_square_diff <= 0.001){
    return('Both results are identical')
  }else{
    return('There seems to be a problem...')
  }

}

## TESTS
## ------------------------------------
# ## Test that Sweep returns an identity matrix
# ## X is a positive definite matrix
# X = matrix(c(2,-1,0,-1,2,-1,0,-1,2),nrow=3,ncol=3)
# Y = X%*%c(2,1,-1)
# 
# XpX = t(X)%*%X # Positive definite, symmetrix matrix
# XpY = t(X)%*%Y
# YpX = t(Y)%*%X
# YpY = t(Y)%*%Y
# 
# A1 = cbind(XpX,XpY)
# A2 = cbind(YpX,YpY)
# A = rbind(A1,A2)
# 
# B1 = mySweep(A, 1)
# B2 = mySweep(B1, 2)
# B3 = mySweep(B2, 3)
# 
# eye = XpX%*%(B3[1:3,1:3])
# print(eye) # Identity

## ------------------------------------
# ## Test linear regression function
# X = matrix(c(2,-1,0,-1,2,-1,0,-1,2),nrow=3,ncol=3)
# Y = X%*%c(2,1,-1)
# # Y = X%*%c(2,1,-1) + matrix(rnorm(3,0,0.5),nrow=3,ncol=1)
# 
# beta_hat = myLinearRegression(X, Y)
# print(beta_hat)
## ------------------------------------

testing_Linear_Regression()
