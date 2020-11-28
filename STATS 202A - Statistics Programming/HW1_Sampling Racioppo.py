# -*- coding: utf-8 -*-
"""

 Stat 202A 2019 Fall - Homework 01
 Author: Peter Racioppo
 Date : Oct 4, 2019

 INSTRUCTIONS: Please fill in the corresponding function. Do not change function names, 
 function inputs or outputs. Do not write anything outside the function. 
 Do not use any of Python's built in functions for matrix inversion or for linear modeling 
 (except for debugging or in the optional examples section).
 
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import rcParams
matplotlib.rcParams.update({'font.size':16})
%matplotlib inline

import random
import math
import time

###############################################
## Function 1: (2a) Uniform 				 ##
###############################################

def sample_uniform(low=0, high=1):

	# Constants
	a = np.power(7,5)
	b = 0
	M = np.power(2,31)-1

	# Number of samples
	num_samp = 100000

	# This function is a uniform random number generator.
	def f_rand_unif(a,b,M,X0,num_samp):
	    X = np.zeros(num_samp) # Initialize array
	    X[0] = X0 # Set first value of array to seed
	    # Generate uniform random samples
	    for i in (np.arange(num_samp-1)+1):
	        X[i] = (a*X[i-1] + b)%M
	    U = X/M # Uniform over [0,1]
	    return X, U

	# Set a seed as the last 4 digits of the current time in microseconds
	X0 = int(time.time()*10e6)%1000
	X, U = f_rand_unif(a,b,M,X0,num_samp) # Generate uniform random numbers

	# Histogram of Xi
	fig, ax = plt.subplots()
	num_bins = 100
	n, bins, patches = plt.hist(X, num_bins, edgecolor = 'black', facecolor='white', alpha=1)
	fig.suptitle('Histogram of Xi, i = 1, ..., n − 1', fontsize=15)
	plt.xlabel('X_i')
	plt.ylabel('Frequency')
	rcParams['figure.figsize'] = 12, 6
	ax.grid()
	plt.show()

	x = X[0:-1]
	y = X[1:]

	# Scatter Plot of (X_i, X_i+1)
	fig, ax = plt.subplots()
	ax.scatter(x, y, s=0.1, marker='.')
	fig.suptitle('Scatterplot of (X_i, X_i+1), i = 1, ..., n − 1', fontsize=15)
	plt.xlabel('X_i')
	plt.ylabel('X_i+1')
	ax.set_aspect(1)
	rcParams['figure.figsize'] = 12, 6
	ax.grid()
	plt.show()

	pass
  
###############################################
## Function 2: (2b) Exponential				 ##
###############################################

def sample_exponential(k=1):
	X = (-1/k)*np.log(U)

	# Histogram of Exponential
	fig, ax = plt.subplots()
	num_bins = 200
	n, bins, patches = plt.hist(X, num_bins, edgecolor = 'black', facecolor='white', alpha=1)
	fig.suptitle('Histogram of Xi, i = 1, ..., n − 1', fontsize=15)
	plt.xlabel('X_i')
	plt.ylabel('Frequency')
	axes = plt.gca()
	axes.set_xlim([0,8])
	rcParams['figure.figsize'] = 12, 6
	ax.grid()
	plt.show()

	pass
  
###############################################
## Function 3: (2c) Normal	 				 ##
###############################################

def sample_normal(mean=0, var=1):

	# Draw from 2 uniform distributions
	# Use different seeds
	X1, U1 = f_rand_unif(a,b,M,X0,num_samp)
	X2, U2 = f_rand_unif(a,b,M,1917,num_samp)

	# Polar coordinates
	pi = math.pi
	theta = 2*pi*U1
	R = mean + np.sqrt(-2*(var**2)*np.log(1-U2/var))

	# Cartesian coordinates
	xx = R*np.cos(theta)
	yy = R*np.sin(theta)

	# Scatterplot of Normal(0,1)
	fig, ax = plt.subplots()
	ax.scatter(xx, yy, s=1, marker='.')
	fig.suptitle('Scatterplot of Normal(0,1)', fontsize=15)
	plt.xlabel('x')
	plt.ylabel('y')
	ax.set_aspect(1)
	rcParams['figure.figsize'] = 12, 6
	ax.grid()
	plt.show()

	T = np.power(R,2)/2

	# Histogram of T = R^2/2
	fig, ax = plt.subplots()
	num_bins = 200
	n, bins, patches = plt.hist(T, num_bins, edgecolor = 'black', facecolor='white', alpha=1)
	fig.suptitle('Histogram of T = R^2/2', fontsize=15)
	plt.xlabel('T')
	plt.ylabel('Frequency')
	axes = plt.gca()
	axes.set_xlim([0,8])
	rcParams['figure.figsize'] = 12, 6
	ax.grid()
	plt.show()
	pass
  
###############################################
## Function 4: (3) Monte Carlo 				 ##
###############################################

def monte_carlo(d=2):

	# Reset the seed
	X0 = int(time.time()*10e6)%1000

	# Draw from 2 uniform distributions, using different seeds
	num_samp = 100000
	_, X_t = f_rand_unif(a,b,M,X0,num_samp)
	_, Y_t = f_rand_unif(a,b,M,1917,num_samp)

	r2 = np.zeros(num_samp) # Initialize array of squared radii
	for i in np.arange(len(X_t)):
	    r2[i] = np.power(X_t[i],2) + np.power(Y_t[i],2) # Squared radii
	ze = num_samp - np.count_nonzero(np.floor(r2)) # Number of points inside the circle
	ratio = ze/num_samp # Ratio of points inside the circle to the total number of points
	Pi = 4*ratio # Pi equals the area of a 2x2 square times the
	# ratio of the area in the square which belongs to the circle.
	print('pi =', Pi) # Print pi

	# This function computes the volume of an n-dimensional
	# unit ball, using the Monte Carlo method.
	def f_MC_ball_volume(a,b,M,X0,d,num_samp):
	    # Draw from a uniform distribution to generate seeds
	    seeds,_ = f_rand_unif(a,b,M,X0,d)

	    # Generate stacked vectors of random uniform samples
	    _, A = f_rand_unif(a,b,M,seeds[0],num_samp)
	    for i in np.arange(d-1):
	        _, X_i = f_rand_unif(a,b,M,seeds[i+1],num_samp)
	        A = np.vstack([A, X_i])

	    r2 = np.zeros(num_samp) # Initialize array for squared radii
	    # Compute squared radii
	    for i in np.arange(num_samp):
	        for j in np.arange(d):
	            r2[i] += np.power(A[j][i],2)
	    ze = num_samp - np.count_nonzero(np.floor(r2)) # Number of points inside the circle
	    ratio = ze/num_samp # Ratio of points inside the circle to the total number of points
	    vol = ratio*np.power(2,d) # Volume of the unit ball
	    return vol

	vol = f_MC_ball_volume(a,b,M,X0,d,num_samp)
	print('Volume of a', d, 'dimensional ball =', vol)

	pass
  
########################################################
## Optional examples (comment out before submitting!) ##
########################################################

## test()

## sample_uniform(0,1)




