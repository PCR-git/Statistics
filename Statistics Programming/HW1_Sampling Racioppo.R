############################################################# 
## Stat 202A 2019 Fall - Homework 01
## Author: Peter Racioppo
## Date : Oct 4, 2019
#############################################################

#############################################################
## INSTRUCTIONS: Please fill in the corresponding function. Do not change function names, 
## function inputs or outputs. Do not write anything outside the function. 
## Do not use any of Python's built in functions for matrix inversion or for linear modeling 
## (except for debugging or in the optional examples section).
#############################################################

###############################################
## Function 1: (2a) Uniform 				 ##
###############################################

sample_uniform <- function(low=0, high=1){
  # Constants
  a = 7**5
  b = 0
  M = (2**31)-1
  
  # Number of samples
  num_samp = 100000
  
  # This function is a uniform random number generator.
  f_rand_unif = function(a,b,M,X0,num_samp){
    X = integer(num_samp) # Initialize array
    # Generate uniform random samples
    X[1] = X0 # Set first value of array to seed
    for(i in 2:num_samp){
      X[i] = (a*X[i-1] + b)%%M}
    return(X)}
  
  # Set a seed as the last 4 digits of the current time in milliseconds
  time = as.numeric(as.numeric(Sys.time())*10**3, digits=15)
  # Hmm, not quite sure what precision time I'm taking here, but it works.
  X0 = (time*10**3)%%1000
  # print(round(X0))
  X = f_rand_unif(a,b,M,X0,num_samp) # Generate uniform random numbers
  U = X/M # Uniform over [0,1]
  
  # Histogram of X_i
  num_bins = 100
  hist(X, breaks=seq(0,max(X),l=num_bins+1),col='white',panel.first=grid())
  
  # Scatter Plot of (X_i, X_i+1)
  x = X[1:num_samp-1]
  y = X[2:num_samp]
  plot(x, y, main="Scatter Plot of (X_i, X_i+1)",panel.first=grid(),
       xlab="X_i ", ylab="X_i+1", pch=".",asp=1)
}

###############################################
## Function 2: (2b) Exponential				 ##
###############################################

sample_exponential <- function(k=1){
  X = (-1/k)*log(U)
  
  # Histogram of Exponential
  num_bins = 200
  hist(X, breaks=seq(0,max(X),l=num_bins+1),col='white', xlim = c(0,8),panel.first=grid())
}

###############################################
## Function 3: (2c) Normal	 				 ##
###############################################

sample_normal <- function(mean=0, var=1){
  # Draw from 2 uniform distributions
  # Use different seeds
  X1 = f_rand_unif(a,b,M,X0,num_samp)
  X2 = f_rand_unif(a,b,M,1917,num_samp)
  U1 = X1/M # Uniform over [0,1]
  U2 = X2/M # Uniform over [0,1]
  
  # Polar coordinates
  theta = 2*pi*U1
  R = mean + sqrt(-2*(var**2)*log(1-U2/var))
  
  # Cartesian coordinates
  xx = R*cos(theta)
  yy = R*sin(theta)
  
  # Scatterplot of Normal(0,1)
  plot(xx, yy, main="Scatter Plot of Normal(0,1)",
       xlab="x",ylab="y",pch=".",panel.first=grid(),asp=1)
  
  # Histogram of T = R^2/2
  T = (R**2)/2
  num_bins = 100
  hist(T, breaks=seq(0,max(T),l=num_bins+1),col='white',xlim = c(0,8),panel.first=grid())
}

###############################################
## Function 4: (3) Monte Carlo 				 ##
###############################################

monte_carlo <- function(d=2){
  # Draw from 2 uniform distributions, using different seeds
  num_samp = 100000
  X_t = f_rand_unif(a,b,M,X0,num_samp)
  Y_t = f_rand_unif(a,b,M,1917,num_samp)
  Ux_t = X_t/M # Uniform over [0,1]
  Uy_t = Y_t/M # Uniform over [0,1]
  
  r2 = integer(num_samp) # Initialize array of squared radii
  for(i in 1:num_samp){
    r2[i] = Ux_t[i]**2 + Uy_t[i]**2}
  ze = length(which(floor(r2) == 0)) # Number of points inside the circle
  ratio = ze/num_samp # Ratio of points inside the circle to the total number of points
  # Pi equals the area of a 2x2 square times the ratio of
  # the area in the square which belongs to the circle.
  Pi = 4*ratio
  cat("pi =", Pi) # Print pi
  
  
  # This function computes the volume of an n-dimensional
  # unit ball, using the Monte Carlo method.
  f_MC_ball_volume = function(a,b,M,X0,d,num_samp){
    # Draw from a uniform distribution to generate seeds
    seeds = f_rand_unif(a,b,M,X0,d)
    
    # Generate stacked vectors of random uniform samples
    B = f_rand_unif(a,b,M,seeds[1],num_samp)
    for(i in 2:d){
      X_i = f_rand_unif(a,b,M,seeds[i],num_samp)
      B = c(B, X_i)}
    A = array(B/M,c(d,num_samp))
    
    r2 = integer(num_samp) # Initialize array for radii
    # Compute squared radii
    for(i in 1:num_samp){
      for(j in 1:d){
        r2[i] = r2[i] + A[j,i]**2}}
    
    ze = length(which(floor(r2) == 0)) # Number of points inside the circle
    ratio = ze/num_samp # Ratio of points inside the circle to the total number of points
    vol = ratio*(2**d) # Volume of the unit ball
    return(vol)}
  vol = f_MC_ball_volume(a,b,M,X0,d,num_samp)
  print('\n')
  cat('Volume of a', d, 'dimensional ball =', vol)
}

########################################################
## Optional examples (comment out before submitting!) ##
########################################################

# sample_uniform(0,1)
# sample_exponential(1)
# sample_normal(0,1)
# monte_carlo(2)
# monte_carlo(3)
# monte_carlo(5)
# monte_carlo(10)

## test()

