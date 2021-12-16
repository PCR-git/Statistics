#########################################################
## Stat 202A - Homework 4
## Author: Peter Racioppo
## Date: 10/27/2019
## Description: This script implements QR decomposition,
## linear regression, and eigen decomposition / PCA 
## based on QR.
#########################################################

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
#############################################################

##################################
## Function 1: QR decomposition ##
##################################

myQR <- function(A){
  
  ## Perform QR decomposition on the matrix A
  ## Input: 
  ## A, an n x m matrix
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  n = dim(A)[1]
  m = dim(A)[2]

  ## Choose Decomposition Method:
  ## 1 for Givens Rotations
  ## 2 for Householder Reflections
  method = 2
  
  if(method == 1){
    ## Using Givens Rotations:
    
    ## Computes the hypotenuse robustly:
    hypot <- function(a,b){
      x = max(abs(a), abs(b))
      y = min(abs(a), abs(b))
      if(x == 0){return(0)}
      t = y/x
      r = x*sqrt(1+t**2)
      return(r)
    }
    
    # Computes Givens Rotation Matrix:
    givens <- function(a,b){
      r = hypot(a,b)
      c = a/r
      s = -b/r
      Qi = matrix(c(c,-s,s,c), nrow = 2, ncol= 2)
      return(Qi)
    }
    
    # Performs QR Decomposition:
    Q = diag(n)
    R = A
    I = diag(n)
    for(j in 1:m){
      for(i in n:-1:(j+1)){
        Gi_t = givens(R[(i-1),j],R[i,j]) # Givens submatrix
        Gi = I # Initialize full Givens matrix
        Gi[(i-1):i,(i-1):i] = Gi_t # Insert 2x2 Givens submatrix
        Ri = matrix(c(R[(i-1),j],R[i,j])) # 2x1 vector to rotate
        R = sign(Ri[1])*round(t(Gi)%*%R,4) # Update R
        Q = round(Q%*%Gi,4) # Update Q
      }
    }
  }
    
  if(method == 2){
    ## Using Householder Reflections:
    R = A # Initialize R
    Q = diag(n) # Initialize Q
    Z = diag(n)*0 # Zero matrix
    
    for(i in 1:min(n,m)){
      xi = matrix(R[i:n,i]) # Lower triangular part of ith column of A
      Ii = diag(length(xi)) # Identity matrix of dimension n+1-i
      xi_n = norm(xi, type="2") # 2-norm of xi
      e1 = matrix(Ii[,1]) # Unit vector with 1 in first entry
      vi_t = xi + sign(xi[1])*xi_n*e1 # ith Householder vector
      c = 2/c(t(vi_t)%*%vi_t) # Constant
      Hi_t = Ii - c*vi_t%*%t(vi_t) # Householder reflection
      Hi = Z # Initialize Hi
      Hi[i:n,i:n] = Hi_t # Nest Hi_tilda in lower right of H_i
      if(i>1){
        Hi[1:i-1,1:i-1] = diag(i-1) # Set upper left of Hi to identity
      }
      Q = round(Q%*%Hi,4) # Update Q
      R = round(Hi%*%R,4) # Update R
    }
  }
  
  ## Function should output a list with Q.transpose and R
  ## Q is an orthogonal n x n matrix
  ## R is an upper triangular n x m matrix
  ## Q and R satisfy the equation: A = Q %*% R
  return(list("Q" = t(Q), "R" = R))
  
}

###############################################
## Function 2: Linear regression based on QR ##
###############################################

myLinearRegression <- function(X, Y){
  
  ## Perform the linear regression of Y on X
  ## Input: 
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Do NOT simulate data in this function. n and p
  ## should be determined by X.
  ## Use myQR inside of this function
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  n = dim(X)[1] # Number of rows of data matrix
  
  # Append column of 1s to X:
  ones = matrix(1, n, 1)
  X = cbind(ones,X) 
  
  p = dim(X)[2] # Number of columns of data matrix
  
  X_list = myQR(X) # Compute QR decomposition of X
  Q_t = X_list$Q # Get Q matrix
  R = X_list$R   # Get R matrix
  
  # We must solve R*b = Q_t*Y:
  RHS = Q_t%*%Y # Right-hand-side
  
  beta_hat = c(1:p)*0 # Initialize beta_hat
  beta_hat[p] = RHS[p]/R[p,p] # Last element of beta_hat
  # Compute elements of beta_hat from bottom up:
  for(i in (p-1):-1:1){
    beta_hat[i] = (RHS[i] - R[i,(i+1):p]%*%beta_hat[(i+1):p])/R[i,i]
  }

  error = (X%*%beta_hat - Y)**2 # Mean squared error
  
  ## Function returns the 1 x (p + 1) vector beta_ls, 
  ## the least squares solution vector
  return(list(beta_hat=beta_hat, error=error))
  
}

# testing <- function(){
#   
#   ## This function is not graded; you can use it to 
#   ## test out the 'myLinearRegression' function 
# 
#   ## Define parameters
#   n    <- 100
#   p    <- 3
# 
#   ## Simulate data from our assumed model.
#   ## We can assume that the true intercept is 0
#   X    <- matrix(rnorm(n * p), nrow = n)
#   beta <- matrix(1:p, nrow = p)
#   Y    <- X %*% beta + rnorm(n)
# 
#   ## Save R's linear regression coefficients
#   R_coef  <- coef(lm(Y ~ X))
#   print(R_coef)
# 
#   ## Save our linear regression coefficients
#   list <- myLinearRegression(X, Y)
#   my_coef = list$beta_hat
#   print(my_coef)
# 
#   ## Are these two vectors different?
#   sum_square_diff <- sum((R_coef - my_coef)^2)
#   if(sum_square_diff <= 0.001){
#     return('Both results are identical')
#   }else{
#     return('There seems to be a problem...')
#   }
#   
# }

## ----------------------------
## TESTS:

# A = t(matrix(c(0.8147, 0.0975, 0.1576, 0.9058, 0.2785,
#                0.9706, 0.1270, 0.5469, 0.9572, 0.9134,
#                0.9575, 0.4854, 0.6324, 0.9649, 0.8003),3))
# list = myQR(A)
# Q = list$Q
# R = list$R
# A2 = t(Q)%*%R
# print(round(A2-A),4)

# testing()
  