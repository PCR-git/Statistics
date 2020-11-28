#########################################################
## Stat 202A - Homework 6
## Author: Peter Racioppo
## Date: 11/10/2019
## Description: See CCLE
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
      Q = Q%*%Hi # Update Q
      R = Hi%*%R # Update R
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

myLM <- function(X, Y){
  
  ## Perform the linear regression of Y on X
  ## Input: 
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Use myQR inside of this function
  
  ########################
  ## FILL IN CODE BELOW ##
  ########################  
  
  n = dim(X)[1] # Number of rows of data matrix
  p = dim(X)[2] # Number of columns of data matrix
  
  X_list = myQR(X) # Compute QR decomposition of X
  Q_t = X_list$Q # Get Q matrix
  R = X_list$R   # Get R matrix
  
  # We must solve R*b = Q_t*Y:
  RHS = Q_t%*%Y # Right-hand-side
  
  beta_ls = c(1:p)*0 # Initialize beta_ls
  beta_ls[p] = RHS[p]/R[p,p] # Last element of beta_ls
  # Compute elements of beta_ls from bottom up:
  for(i in (p-1):-1:1){
    beta_ls[i] = (RHS[i] - R[i,(i+1):p]%*%beta_ls[(i+1):p])/R[i,i]
  }
  
  ## Function returns the 1 x p vector beta_ls, notice this version do not add intercept.
  ## the least squares solution vector
  return(beta_ls)
  
}

######################################
## Function 3: Logistic regression  ##
######################################

## Expit/sigmoid function
expit <- function(x){
  1 / (1 + exp(-x))
}

myLogisticSolution <- function(X, Y){

  ########################
  ## FILL IN CODE BELOW ##
  ########################
  
  num_iter = 20 # Number of iterations
  m = dim(X)[1] # Number of rows of X
  ones = matrix(1, m, 1) # One vector
  X = cbind(ones,X) # Append one vector
  n = dim(X)[2] # Number of columns of X
  
  beta = matrix(0, n, 1) # Initialize weight matrix
  # Initialize loss function vector:
  L_vec = matrix(0, 1, num_iter)
  # Run Newton-Raphson:
  # beta_i+1 = beta_i + (t(X)*Vm*X)^-1 * (t(X)*(Y-p))
  for(i in 1:num_iter){
    p = expit(X%*%beta) # Probability
    Vm = diag(c(p*(1-p))) # Variance matrix
    M = t(X)%*%Vm%*%X
    xtyp = t(X)%*%(Y-p)
    
    # Use QR function to invert M:
    # M^-1 = R^-1 %*% Q_t
    M_list = myQR(M) # Compute QR decomposition of M
    Q_t = M_list$Q # Get Q matrix
    R = M_list$R   # Get R matrix
    Rt = t(R) # R transpose
    Rt_inv = Rt*0 # Initialize R transpose inverse
    # Compute R transpose inverse:
    for(k in 1:n){
      Rt_inv[k,k] = 1/Rt[k,k]
      if(k<n){
        for(i in (k+1):n){
          Rt_inv[i,k] = -Rt[i,k:(i-1)]%*%Rt_inv[k:(i-1),k]/Rt[i,i]
        }
      }
    }
    R_inv = t(Rt_inv) # R inverse
    M_inv = R_inv%*%Q_t # M inverse
    beta_m = beta # Weight matrix before update
    beta = beta + M_inv%*%xtyp # Update weight matrix
    L = t(log(p))%*%Y + t(log(1-p))%*%(1-Y) # Loss at ith step
    L_vec[i] = L # Append loss at ith step to vector
  }
  
  # beta = beta[-1] # Evidently, we don't want the intercept
  return(beta)
}

###################################################
## Function 4: Adaboost  ##
###################################################

myAdaboost <- function(x1, x2, y){

  ## Expit/sigmoid function
  expit <- function(x){
    1 / (1 + exp(-x))
  }
  
  ## Calculate Gradient
  grad <- function(X,Y,W,b){
    m = dim(X)[1]
    ones = matrix(1, 1, m) # One vector
    p = expit(X%*%W + b) # Probability
    db = as.numeric((1/m)*ones%*%(p-Y)) # Gradient
    return(db)
  }
  
  ## Gradient Descent
  gradient_descent <- function(X,Y,W,b,l_rate,num_iter){
    for(i in 1:num_iter){
      db = grad(X,Y,W,b) # Gradient at ith step
      bm = b # Save b value
      b = b - l_rate*db # Update b value
      # If change in b falls below tolerance, break
      if(abs(b-bm)<0.001){
        # print(i)
        # print("break")
        break
      }
    }
    return(b)
  }
  
  ## Create a random ensemble of weak classifiers:
  ensemble <- function(X,Y,range,n_ensmb){
    bt_vec = runif(n_ensmb)*range # Generate random y-intercept values
    # Generate random slopes:
    Wt_vec = matrix(c(runif(n_ensmb*2)),nrow=n_ensmb,ncol=2)
    for(i in 1:n_ensmb){
      # Run gradient descent on the y-intercept value:
      bt = bt_vec[i]
      Wt = Wt_vec[i,]
      bt = gradient_descent(X,Y,Wt,bt,l_rate,num_iter)
      Wt_vec[i,] = Wt
      bt_vec[i] = bt
    }
    alpha = matrix(1, 1, n_ensmb)/n_ensmb
    Ct = list("W_vec"=Wt_vec,"b_vec"=bt_vec,"alpha"=alpha)
    return(Ct)
  }
  
  ## Predict class
  predict <- function(X,Y,C){
    W_vec = C$W_vec
    b_vec = C$b_vec
    alpha = C$alpha
    n_m = length(alpha)
    p = 0
    for(i in 1:n_m){
      W = W_vec[i,]
      b = b_vec[i]
      p = p + expit(X%*%W + b)*alpha[i] # Probability
    }
    p = p/sum(alpha)
    Yp <- as.integer(p >= 0.5) # Choose majority probability
    Yp = 2*Yp - 1 # Change Yp to have 1s or -1s
    return(Yp)
  }
  
  ## Exponential Loss Function
  L_exp <- function(X,Y,C,At,s=1){
    alpha = C$alpha
    n_m = length(alpha)+1
    Yp = predict(X,Y,C)
    arg = -Yp*Y
    if(n_m==1){
      w = 1
    }
    if(n_m>1){
      w_vec = exp(arg)
      if(s==1){
        Yt = predict(X,Y,At)
        Yt = (Yt+1)/2
        Y = (Y+1)/2
        mask = abs(Y-Yt) # Zeros where they're the same
        Y = 2*Y - 1
        w_vec = w_vec*mask
      }
      w = sum(w_vec)
    }
    return(w)
  }
  
  ## Choose the weak classifier that minimizes the total weighted error
  Choose <- function(X,Y,C,Ct){
    Wt_vec = Ct$W_vec
    bt_vec = Ct$b_vec
    alpha_t = Ct$alpha
    n_r = length(alpha_t)
    w_vec = matrix(0,1,n_r)
    for(i in 1:n_r){
      Wt = matrix(Wt_vec[i,],nrow=1,ncol=2)
      bt = bt_vec[i]
      At = list("W_vec"=Wt,"b_vec"=bt,"alpha"=1)
      w = L_exp(X,Y,C,At)
      w_vec[i] = w
    }
    argmin = which.min(w_vec)
    Wt = matrix(Wt_vec[argmin,],nrow=1,ncol=2)
    bt = bt_vec[argmin]
    At = list("W_vec"=Wt,"b_vec"=bt,"alpha"=1)
    return(At)
  }
  
  ## Calculate the new weight
  New_Weight <- function(X,Y,C,At){
    num = L_exp(X,Y,C,At)
    denom = L_exp(X,Y,C,At,s=0)
    eps_m = num/denom
    alpha_m = (1/2)*log((1-eps_m)/eps_m)
    print(alpha_m)
    return(alpha_m)
  }
  
  ## Update Classifier
  Update_Classifier <- function(X,Y,C,At){
    alpha_m = as.numeric(New_Weight(X,Y,C,At))
    if(alpha_m>0){
      Wt_vec = At$W_vec
      bt_vec = At$b_vec
      W_vec = C$W_vec
      b_vec = C$b_vec
      alpha = C$alpha
      W_vec = rbind(W_vec,Wt_vec)
      b_vec = matrix(c(b_vec,as.numeric(bt_vec)))
      alpha = cbind(alpha,alpha_m)
      alpha = alpha/sum(alpha)
      C = list("W_vec"=W_vec,"b_vec"=b_vec,"alpha"=alpha)
    }
    return(C)
  }
  
  ## Preprocessing:
  X = cbind(matrix(x1),matrix(x2)) # Data matrix
  Y = matrix(y) # Class matrix
  Y = 2*Y - 1 # Change Y to have 1s or -1s
  
  ## Hyperparameters:
  l_rate = 0.25 # Learning rate
  num_iter = 200 # Number of iterations in gradient descent
  range = 5 # Range for random values
  n_ensmb = 15 # Number in the ensemble
  
  C = ensemble(X,Y,range,1)
  Yp = predict(X,Y,C)
  acc = (length(which((Y-Yp)==0)))/dim(Y)[1] # Accuracy
  print(acc)
  
  iter = 20
  for(i in 1:iter){
    Ct = ensemble(X,Y,range,n_ensmb)
    At = Choose(X,Y,C,Ct)
    C = Update_Classifier(X,Y,C,At) 
  }
  
  Yp = predict(X,Y,C)
  acc = (length(which((Y-Yp)==0)))/dim(Y)[1] # Accuracy
  print(acc)
  
}

###################################################
## Function 5: XGBoost  ##
###################################################

myXGBoost <- function(x1, x2, y) {

  ## Expit/sigmoid function
  expit <- function(x){
    1 / (1 + exp(-x))
  }
  
  Loss <- function(y,yp){
    Y = 2*y-1
    Yp = 2*yp-1
    L = sum(1-y*yp)
    return(L)
  }
  
  ChooseQuad <- function(x1,x2,x1t,x2t,y){
    mask1 = y
    mask2 = -(y-1)
    
    dx1p = (as.numeric((x1-x1t)>=0))
    dx2p = (as.numeric((x2-x2t)>=0))
    dx1m = (as.numeric((x1-x1t)<0))
    dx2m = (as.numeric((x2-x2t)<0))
    p1 = sum(y)/length(y)
    p2 = 1-p1
    
    q11 = dx1p*dx2p
    q12 = dx1p*dx2m
    q21 = dx1m*dx2p
    q22 = dx1m*dx2m
    
    S11 = sum(q11*mask1)/p1
    S12 = sum(q12*mask1)/p1
    S21 = sum(q21*mask1)/p1
    S22 = sum(q22*mask1)/p1
    P11 = sum(q11*mask2)/p2
    P12 = sum(q12*mask2)/p2
    P21 = sum(q21*mask2)/p2
    P22 = sum(q22*mask2)/p2
    
    v1 = c(S11,P11)
    v2 = c(S12,P12)
    v3 = c(S21,P21)
    v4 = c(S22,P22)
    
    m1 = which(v1==max(v1))
    m2 = which(v2==max(v2))
    m3 = which(v3==max(v3))
    m4 = which(v4==max(v4))
    
    M = c(m1,m2,m3,m4)
    M = -(M-2)
    # print(M)
    return(M)
  }
  
  DecisionTree <- function(x1,x2,y,n_ensmb){
    x1t = runif(1)
    x2t = runif(1)
    M = ChooseQuad(x1,x2,x1t,x2t,y)
    C = list("M"=M,"X"=c(x1t,x2t))
  }
  
  predict <- function(x1,x2,C){
    M = C$M
    X = C$X
    x1t = X[1]
    x2t = X[2]
    
    dx1p = (as.numeric((x1-x1t)>=0))
    dx2p = (as.numeric((x2-x2t)>=0))
    dx1m = (as.numeric((x1-x1t)<0))
    dx2m = (as.numeric((x2-x2t)<0))
    q11 = dx1p*dx2p*M[1]
    q12 = dx1p*dx2m*M[2]
    q21 = dx1m*dx2p*M[3]
    q22 = dx1m*dx2m*M[4]

    yp = q11+q12+q21+q22
    return(yp)
  }
  
  LineSearch <- function(x1,x2,y,yp,rp){
    gamma = 1
    lr = 0.1
    iter = 10
    for(i in 1:iter){
      L1 = Loss(y,yp)
      yp_t = yp
      yp_t = (yp_t + gamma*rp)/(1+gamma)
      L2 = Loss(y,yp)
      if(L2>=L1){gamma=gamma-lr}
      if(L2<L1){gamma=gamma+lr}
    }
    return(gamma)
  }
  
  C = DecisionTree(x1,x2,y)
  yp = predict(x1,x2,C)
  r_index = abs(y-yp)
  index = matrix(which(r_index==1))
  x1r = x1[index]
  x2r = x2[index]
  yr = y[index] # y's that we got wrong
  acc = (length(which((r_index)==0)))/length(y) # Accuracy
  print(acc)
  
  iter = 3
  for(i in 1:iter){
    Cr = DecisionTree(x1r,x2r,yr)
    rp = predict(x1r,x2r,Cr)
    # print(rp)
    yr2 = yp*0
    yr2[index] = rp
    gamma = LineSearch(x1r,x2r,y,yp,yr2)
    yp = (yp + gamma*yr2)/(1+gamma)
    r_index = abs(y-yp) # Indices of y we predicted incorrectly
    index = which(r_index==1)
    x1r = x1[index]
    x2r = x2[index]
    yr = y[index] # y's that we got wrong
    acc = (length(which((r_index)==0)))/length(y) # Accuracy
    print(acc)
    }
}

## Simulation

test <- function() {

  # # Test (1)
  # n <- 5000
  # p <- 4
  # 
  # X    <- matrix(rnorm(n * p), nrow = n)
  # beta <- c(12,-2,-3, 4)
  # Y    <- 1 * (runif(n) < expit(X %*% beta))
  
  ## Our solution
  # logistic_beta <- myLogisticSolution(X, Y)
  # print(logistic_beta)
  
  ## R's solution
  # coef(glm(Y ~ X + 0, family = binomial(link = 'logit')))

  # Test (2, 3)
  num_sample <- 10000

  # Draw from uniform distribution
  # (default is 0 to 1):
  x1 <- runif(num_sample)
  x2 <- runif(num_sample)
  y <- as.integer((x1^2+x2^2 < 1))
  
  # myAdaboost(x1, x2, y)
  myXGBoost(x1, x2, y)
}

###################################################

test()

