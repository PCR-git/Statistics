# Probability - HW 3
# Peter Racioppo
# 103953689

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import scipy.integrate

sigma = 2 # Stdev
mu = 1 # Mean
N = 1000

# # Numerical Integration
# # Integrand
# normalizer = 1/math.sqrt(2*math.pi*pow(sigma,2))
# integrand = lambda x: normalizer*(1/(1+pow(x,2)))*np.exp(-pow((x-mu),2)/(2*pow(sigma,2)))
# # Integration Bounds
# a = -np.inf
# b = np.inf
# integral = scipy.integrate.quad(integrand, a, b) # Integration
# print(integral)
# # integral ~= 0.438

# Computation of numerical integral for 100 values of mu
ints = list(range(101))
mu_vec = [i/10 for i in ints] # List of mu values
integral_vec = [] # Initalize array of integral values
for i in mu_vec:
    mu = i
    normalizer = 1/math.sqrt(2*math.pi*pow(sigma,2))
    integrand = lambda x: normalizer*(1/(1+pow(x,2)))*np.exp(-pow((x-mu),2)/(2*pow(sigma,2)))
    integral = scipy.integrate.quad(integrand, -np.inf, np.inf) # Integration
    integral_vec.append(integral[0])

# # Single Monte Carlo Trial
# x_iid = np.random.normal(loc=mu, scale=sigma, size=N)
# J_Ni = []
# for i in x_iid:
#     J_Ni.append(1/(1+pow(i,2)))
# J_N = sum(J_Ni)/N
# #print(J_N)
# # 0.435

# 100 Monte Carlo trials
J_N_vec = [] # Initalize array of Monte Carlo values
for i in mu_vec:
    J_Ni = [] # Initalize array
    mu = i
    x_iid = np.random.normal(loc=mu, scale=sigma, size=N)
    for i in x_iid:
        J_Ni.append(1/(1+pow(i,2)))
    J_N = sum(J_Ni)/N
    J_N_vec.append(J_N)

# Analytic Computation of E(Y)
# The expected value of a normal distribution is just the mean
# E(Y) = E(1/(1+X^2)) = E(1)/(E(1+X^2)) = 1/(1+(E(X))^2) = 1/(1+mu^2) = 1/2

# Plot
fig = plt.figure() # Create figure
ax = plt.axes() # Create plot axes
axg = fig.gca() # Set GCA
ax.plot(ints,integral_vec,linestyle='-', marker='.', color='b',label='Num Integral') # Plot numeric integral
ax.plot(ints,J_N_vec,linestyle='--', marker='.', color='r',label='Monte Carlo') # Plot Monte Carlo
plt.title("Expectation vs Mean") # Plot title
plt.xlabel("Mean") # x-axis label
plt.ylabel("Expectation"); # y-axis label
plt.grid() # Plot grid
plt.legend(); #Legend
plt.show() # Show the plot
