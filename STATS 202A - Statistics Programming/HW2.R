############################################################# 
## Stat 202A 2019 Fall - Homework 02
## Author: Peter Racioppo
## Date : 10/13/2019
#############################################################

#############################################################
## INSTRUCTIONS: Please fill in the corresponding function. Do not change function names, 
## function inputs or outputs. Do not write anything outside the function. 
## See detailed requirement in python files.
#############################################################

###############################################
## Function 1A                  			 ##
###############################################

sample_uniform <- function(size=1000, low=0, high=1){
    # Function 1A: Uniform (Same as HW1)
    # Detail: This function returns a np.ndarray with given size.
  
    # Constants
    a = 7**5
    b = 0
    M = (2**31)-1
    
    # Set a seed as the last 4 digits of the current time in milliseconds
    # time = as.numeric(as.numeric(Sys.time())*10**3 + as.numeric(Sys.time())*10**6, digits=15)
    time = as.numeric(as.numeric(Sys.time())*10**3, digits=15)
    Y0 = (time*10**3)%%1000

    Y = integer(size+5) # Initialize array, accounting for burn in of 5 points
    Y[1] = Y0 # Set first value of array to seed
    # Generate pseudorandom samples:
    for(i in 2:(size+5)){Y[i] = (a*Y[i-1] + b)%%M}
    U = Y[6:length(Y)]/M # Uniform over [0,1] (include burn in)
    X = U*(high-low) + low  # Uniform over [low,high]
    return(X)
}

###############################################
## Function 1B                  			 ##
###############################################

sample_normal_chain <- function(x0, c, chain_length=100, mean=0, var=1){
    # Function 1B: Normal by Metropolis
    # Detail: This function returns multiple chains by Metropolis sampling from N(mean, var).
    # For every train, the proposal distribution at x is y ~ Uniform[x-c, x+c].
    # The input x0 is a 1 dimension np.ndarray. 
    # Return a np.ndarray X with size x0.shape[0] * chain_length. In X, each row X[i]
    # should be a chain and t-th column X[:, t] corresponding to all different X_t.
  
  X = array(integer(length(x0)*chain_length), dim=c(length(x0),chain_length)) # Initialize X
  for(i in 1:length(x0)){
    x_vec = c(x0[i]) # Starting value of MCMC
    # Sample from a uniform distribution:
    U = sample_uniform(chain_length,0,1)
    # Sample from a different uniform distribution:
    samp = sample_uniform(chain_length,0,1)
    for(j in 1:(chain_length)){
      x = x_vec[length(x_vec)] # X_t
      y = samp[j]*2*c + (x-c) # X_t+1 (sample from proposal distrb.)
      r = exp(((x-mean)**2-(y-mean)**2)/(2*var)) # p(X_t+1)/p(X_t)
      if(min(1, r) > U[j]){
        x = y}
      x_vec = c(x_vec,x)
      X[i,j] = x}}
      
  return(X)}

###############################################
## Function 1C                  			 ##
###############################################
    
metropolis_simulation <- function(num_chain=1000, chain_length=100, mean=0, var=1){
    # Function 1C: Simulate metropolis with different setting.
    # Detail: Try different settings and output movies of histograms.
    
    list_a = c(0,-0.1, -1, -10) # Lower limit for starting sample
    list_b = c(1,0.1, 1, 10)# Upper limit for starting sample
    list_c = c(1, 2) # Sampling parameter
  
    # list_a = c(-1)
    # list_b = c(1)
    # list_c = c(1)

    for(a in list_a){
      for(b in list_b){
        for(c in list_c){
  
          # Sample the starting point from a uniform distribution:
          x0 = sample_uniform(num_chain,a,b)
          
          # Run the Metropolis Algorithm:
          X_normal = sample_normal_chain(x0,c,chain_length,mean,var)
      
          burn_in = ceiling(0.5*chain_length) # Burn in
          X_f = X_normal[1:num_chain,burn_in:chain_length] # Remove values for burn in
          X_f = array(X_f, dim=c(1,(num_chain*(chain_length-burn_in)))) # Flatten array
          
          # Histogram
          sigma = sqrt(var)
          num_bins = 80
          hist(X_f, breaks=num_bins-1,xlim=c(-4,4),panel.first=grid(),
               main="Histogram of Metropolis Samples of Normal Distrb.",
               xlab='X',ylab="Density",col=rgb(0,0,1,0.5))
          
          # Overlaying the target PDF
          x = seq(-4, 4, length=200)
          y = dnorm(x, mean=0, sd=1)
          par(new=TRUE)
          plot(x,y,type="l",lwd=2,xlim=c(-4,4),add=TRUE,lty="dashed",
               main="Histogram of Metropolis Samples of Normal Distrb.",
               xlab = 'X', ylab="Density")
          
        }
      }
    }
    return(X_f)
}

###############################################
## Function 2A                  			 ##
###############################################

gibbs_sample <- function(x0, y0, rho, num_chain=1000, chain_length=100, mean=0, var=1){
  # Function 2A: Bivariate normal with correlation rho
  # Detail :    This function return multiple chains by Gibbs sampling
  # The input x0, y0, rho, num_chain is a number. This time, we use same starting point. 
  # Return a np.ndarray X with size num_chain * chain_length * 2. In X, each row X[i]
  # should be a chain and t-th column X[:, t] corresponding to all different pair (X_t, Y_t).
  
  sample_normal <- function(mean=0,var=1){
    # This function samples from a normal distribution.
    
    # Draw from 2 uniform distributions
    U1 = sample_uniform(2,0,1)[1]
    U2 = sample_uniform(2,0,1)[2]
    
    # Polar coordinates
    theta = 2*pi*U1
    R = sqrt(-2*log(1-U2))
    
    # Cartesian coordinates
    z = R*cos(theta)
    xx = z*sqrt(var) + mean
    
    return(xx)}
  
  X = array(integer(num_chain*chain_length*2), dim=c(num_chain,chain_length,2)) # Initialize X
  for(i in 1:num_chain){
    x_vec = x0
    y_vec = y0
    for(j in 1:chain_length){
      x = x_vec[length(x_vec)]
      y = y_vec[length(y_vec)]
      x = sample_normal(rho*y,1-rho**2)
      y = sample_normal(rho*x,1-rho**2)
      x_vec = c(x_vec,x)
      y_vec = c(y_vec,y)
      X[i,j,1] = x
      X[i,j,2] = y}}
  return(X)}

###############################################
## Function 2B                  			 ##
###############################################

gibbs_simulation <- function(){
  # Function 2B: Simulate Gibbs with different rho and plot
  # Detail : Try different setting and output movie of histgrams. 
  # Discard first 50 steps and output 50~100 steps only.
  
  # Parameters
  mean = 0
  var = 1
  num_chain=1000
  chain_length=100
  x0 = 0
  y0 = 0
  
  list_rho = list(0, 0.2, -1, 1)
  # list_rho = list(0)
  for(rho in list_rho){
    
    # Run Gibbs Sampling
    X = gibbs_sample(x0,y0,rho,num_chain,chain_length,mean,var)
    
    # Discard first 50 steps
    burn_in = 50 # Burn in
    X_f = X[,burn_in:chain_length,] # Remove values
    X_data = X_f[,,1]
    Y_data = X_f[,,2]
    # Flatten:
    X_data = array(X_data, dim=c(1,(num_chain*(chain_length-burn_in))))
    Y_data = array(Y_data, dim=c(1,(num_chain*(chain_length-burn_in))))
    
    # Plot info
    sigma = sqrt(var)
    lim = 4
    num_bins = 80
    
    # Histogram
    sigma = sqrt(var)
    num_bins = 80
    hist(X_data, breaks=num_bins-1,xlim=c(-4,4),panel.first=grid(),
         main="", xlab="",ylab="", col=rgb(0,0,1,0.5))
    par(new=TRUE)
    hist(Y_data, breaks=num_bins-1,xlim=c(-4,4),panel.first=grid(),
         main="", xlab="",ylab="", col=rgb(1,0,0,0.5))
    
    # Overlaying the target PDF
    x = seq(-4, 4, length=200)
    y = dnorm(x, mean=0, sd=1)
    par(new=TRUE)
    plot(x,y,type="l",lwd=2,xlim=c(-4,4),add=TRUE,lty="dashed",
         main="Histogram of Gibbs Samples of Normal Distrb.",
         xlab = 'X, Y', ylab="Density", col = "black")
    
    # Scatter plot of a single chain
    X1 = X_data[1:(chain_length-burn_in)]
    Y1 = Y_data[1:(chain_length-burn_in)]

    plot(X1,Y1, main="Scatterplot of Single Chain, Gibbs Normal Samples",
         xlab="X", ylab="Y", pch=19,panel.first=grid())
  }
}

########################################################
## Optional examples (comment out before submitting!) ##
########################################################

## test()

# X = sample_uniform(1000,0,1)
# X_f = metropolis_simulation(1000,100,0,1)
# gibbs_simulation()

