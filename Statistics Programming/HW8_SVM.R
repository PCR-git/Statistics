## ----------------------------------
## Support Vector Machine
## Peter Racioppo
## ----------------------------------
## Kernels:

## Linear
Linear <- function(u,v){
  out = t(u)%*%v
  return(out)
}

## Polynomial
Poly <- function(gamma,u,v,coef0){
  out = (gamma*t(u)%*%v+coef0)**2
  return(out)
}

## Radial Basis Function
RBF <- function(gamma,u,v){
  out = exp(-gamma*t(u-v)%*%(u-v))
  return(out)
}

## Sigmoid
Sigmoid <- function(gamma,u,v,coef0){
  out = tanh(gamma*t(u)%*%v+coef0)
  return(out)
}

## Kernel
Kernel <- function(u,v,gamma=0.1,coef0=0,type="linear"){
  if(type == "linear"){
    out = Linear(u,v)
  }
  if(type == "poly"){
    out = Poly(gamma,u,v,coef0)
  }
  if(type == "rbf"){
    out = RBF(gamma,u,v)
  }
  if(type == "sigmoid"){
    out = Sigmoid(gamma,u,v,coef0)
  }
  
  return(out)
}

dK <- function(w,xi,gamma,coef0,type="linear"){
  if(type == "linear"){
    out = xi
  }
  if(type == "poly"){
    out = 2*(t(w)%*%xi)*xi
  }
  if(type == "rbf"){
    norm = t(w-x)%*%(w-x)
    out = -2*gamma*(w-x)*exp(-gamma*norm)
  }
  if(type == "sigmoid"){
    out = (1-tanh(t(w)%*%xi)**2)*xi
  }
  return(out)
}

## ----------------------------------

## Hinge Loss
HingeLoss <- function(y,w,x,lambda){
  L = 0
  n = length(x)
  for(i in 1:n){
    Li = max(0,1 - y[i]*(Linear(w,x[i])))
    L = L + Li
  L = L/n + lambda*t(w)%*%w
  }
  return(L)
}

## Compute subgradient
SubGrad <- function(yi,w,xi,lambda){
  # y_hat = t(w)%*%xi
  y_hat = Kernel(w,xi,type="linear")
  if(yi*y_hat < 1){
    # out = -yi*xi
    out = -yi*dK(w,xi,type="linear") + lambda*w
  } else {out = xi*0}
  return(out)
}

## Sub-Gradient Descent
SGD <- function(y,w,x,lambda,scale,iter){
  n = dim(x)[1] # Length of x
  s = length(w) # Length of w
  for(c in 1:iter){
    j_rand = ceiling(runif(1, min=0, max=s)) # Random element index
    dL = 0
    for(i in 1:n){
      # Compute subgradient:
      dL_i = SubGrad(y[i],w,x[i,],lambda)[j_rand]
      dL = dL + dL_i
    }
    dL = dL/n + lambda*w
    w[j_rand] = w[j_rand] - scale*dL
  }
  return(w)
}

## Coordinate Descent (Using Lagrange Multipliers)
CoordDesc <- function(y,x,lambda,scale,iter){
  n = dim(x)[1]
  s = length(w)
  
  alpha = matrix(1,n)
  # Initialize Lagrange multipliers:
  nu = matrix(0,n)
  rho = matrix(0,n)
  chi = 0
  zeta = 0
  
  for(k in 1:iter){
    for(i in 1:n){
      sum = 0
      for(j in 1:n){
        sum = sum + y[i]*y[j]*(x[i,]%*%x[j,])*alpha[j]
      }
      dalpha = sum - 1 - nu[i] + rho[i] + chi*y[i] - zeta*y[i]
      dnu = -alpha[i] + 1
      drho = alpha[i] - 1/(2*n*lambda) + 1
      dalpha = sum - 1
      
      alpha[i] = alpha[i] - scale*dalpha
      nu[i] = nu[i] - scale*dnu
      rho[i] = rho[i] - scale*drho
    }
    dchi = 0
    dzeta = 0
    for(j in 1:n){
      dchi = dchi + y[j]*alpha[j] + 1
      dzeta = dzeta - y[j]*alpha[j] + 1
    }
    chi = chi - scale*dchi
    zeta = zeta - scale*dzeta
  }
  
  w = matrix(0,s)
  for(i in 1:s){
    w[i] = sum(alpha*y*x[,i])
    print(w[i])
  }
  
  return(w)
}

## Sequential Minimal Optimization
SMO <- function(y,x,lambda,scale,iter){
  n = dim(x)[1] # Length of x
  s = length(w) # Length of w
  C = 1/(2*n*lambda)
  alpha = rbind(matrix(C/n,n/2),matrix(-C/n,n/2))
  for(c in 1:iter){
    j_r = ceiling(runif(1, min=0, max=n))
    samp2 = 1:(j_r-1)
    samp3 = (j_r+1):n
    i_r = sample(c(samp2,samp3),1)
    if(y[i_r]!=y[j_r]){
      L = max(alpha[j_r] - alpha[i_r])
      H = min(C,C+alpha[j_r]-alpha[i_r])
    }
    if(y[i_r]==y[j_r]){
      L = max(alpha[i_r] + alpha[j_r] - C)
      H = min(C,alpha[i_r]+alpha[j_r])
    }
    Ei = t(w)%*%x[i_r,] - y[i_r]
    Ej = t(w)%*%x[j_r,] - y[j_r]
    eta = 2*x[i_r,]%*%x[j_r,] - x[i_r,]%*%x[i_r,] - x[j_r,]%*%x[j_r,]
    alpha_j = alpha[j_r] - y[j_r]*(Ei-Ej)/eta
    if(alpha_j>H){alpha_j=H}
    if(alpha_j<L){alpha_j=L}
    alpha_i = alpha[i_r] + y[i_r]*y[j_r]*(alpha[j_r]-alpha_j)
    
    alpha[j_r] = alpha_j
    alpha[i_r] = alpha_i
  }
  
  w = matrix(0,s)
  for(i in 1:s){
    w[i] = sum(alpha*y*x[,i])
    print(w[i])
  }
  
  return(w)
}

## ----------------------------------
## Generate test data
Num = 100
xc = matrix(c(1:Num)/10)
rand3 = sample(1:6, 1)/2
rand1 = rnorm(Num, mean = 0, sd = 0.3*rand3)
rand2 = rnorm(Num, mean = 0, sd = 0.3*rand3)
yc1 = rand3*xc + rand1
yc2 = rand3*xc + rand2 + 3
x1 = cbind(xc,yc1,-1) # coords
x2 = cbind(xc,yc2,-1)
y1 = matrix(1,Num) # labels
y2 = matrix(-1,Num)
x = rbind(x1,x2)
y = rbind(y1,y2)
n = dim(x)[1]

## ----------------------------------
## Run SGD

lambda = 1
scale = 0.01
iter = 500
w = matrix(c(1,0.1,0))
w = SGD(y,w,x,lambda,scale,iter)
print(w)

## Predict:
y_hat_grad = x%*%w
for(i in 1:n){
  # y_hat_grad[i] = t(w)%*%x[i,]
  y_hat_grad[i] = Kernel(w,x[i,],type="linear")
}
y_pred_grad = -(w[1]/w[2])*xc + w[3]/w[2]
accuracy_grad = sum(as.numeric(y*y_hat_grad>0))/n

# ## ----------------------------------
# ## Run Coordinate Descent
# # scale = 0.01
# iter = 100
# w = CoordDesc(y,x,lambda,scale,iter)
# y_pred_cd = -(w[1]/w[2])*xc + w[3]/w[2]
# 
# ## ----------------------------------
# ## Run Sequential Minimal Optimization
# # scale = 0.01
# iter = 500
# w = SMO(y,x,lambda,scale,iter)
# y_pred_smo = -(w[1]/w[2])*xc + w[3]/w[2]
# # y_hat_smo = x%*%w
# # accuracy_smo = sum(as.numeric(y*y_hat_smo>0))/(2*Num)

# ------------------------
## Plot
plot(xc,yc1,col="blue")
points(xc,yc2,col="red")
lines(xc,y_pred_grad,lty=1,col="black")
# lines(xc,y_pred_smo,lty=2,col="red")
# lines(xc,y_pred_cd,lty=3)
grid()

print(accuracy_grad)

