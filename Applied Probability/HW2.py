
# Probability - HW 2
# Peter Racioppo
# 103953689

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random
import math

#import scipy
import scipy.stats
from scipy.stats import poisson

import pickle
with open('spike_times.p', 'rb') as fp:
        spike_times = pickle.load(fp)

#print(spike_times)
tper = 5 # Bin size
tmin = min(spike_times) # Min time at which a spike occurs
tmax = max(spike_times) # Max time at which a spike occurs
#print(tmin) # Check
#print(tmax) # Check

# The number of bins of length 5 in between the min and max values of spike_times
length = math.ceil((tmax-tmin)/5)

hold1 = list(np.ones(length)) # Vector of ones
time_per = [i * tmin for i in hold1] # Vector with each element = tmin
for i in range(length): # make time_per = tmin k*tper
    time_per[i] += i*tper
#print(time_per[len(time_per)-1]) # Check
#print(time_per) # Check

bins1 = time_per # Vector of bins
hist = np.histogram(spike_times, bins = bins1) # Histogram
nspikes = hist[0] # Vector of number of spikes per bin

time_per.pop(0) # Remove first entry of time_per.
# This means that we are counting the number of spikes in
# LESS THAN the amount of time in the same element of time_per

#print(time_per) # Check
#print(len(time_per)) # Check
#print(len(nspikes)) # Check

# Plot Histogram
fig = plt.figure() # Create figure
ax = plt.axes() # Create plot axes
axg = fig.gca()
plt.stem(time_per,nspikes) # Plot histogram
plt.title("Histogram of Spikes per Time Interval") # Plot title
plt.xlabel("Time Interval (ms)") # x-axis label
plt.ylabel("Number of Spikes"); # y-axis label
plt.show() # Show the plot

nspikes_min = min(nspikes) # Min number of spikes in a bin
nspikes_max = max(nspikes) # Max number of spikes in a bin
nspikes_mean = sum(nspikes)/len(nspikes) # Mean of nspikes

hold1.pop() # Remove 1 element from old vectors of ones
mean_vec = [x * nspikes_mean for x in hold1] # Vector of nspikes_mean
nspikes_diff = nspikes-mean_vec # Vector of nspikes - nspikes_mean
nspikes_diff_2 = hold1 # Initalize nspikes_diff_2 array
for i in range(len(nspikes_diff)):
    # Vector of (nspikes - nspikes_mean)^2
    nspikes_diff_2[i] = math.pow(nspikes_diff[i],2)
nspikes_var = sum(nspikes_diff_2)/len(nspikes) # Variance of nspikes
nspikes_std = math.sqrt(nspikes_var) # Standard deviation of number of spikes in the bins

# print(nspikes_min)  # Check
# print(nspikes_max)  # Check
# print(nspikes_mean) # Check
# print(nspikes_std)  # Check

bins2 = [i+0.1 for i in range(31)] # Vector of 0 through 30 + 0.1
bins2.insert(0, 0) # Insert a zero before the first element of bins2
# We are interested in whole numbers, which now fall between each
# element-pair in bins2
# print(bins2) # Check
hist2 = np.histogram(nspikes, bins = bins2) # Histogram
spike_count = hist2[0] # Vector of number of number of spikes
p_samp = [i/sum(spike_count) for i in spike_count]
# print(p_samp) # Check
bins2.remove(0) # Remove 0 element from bins2
# Make bins2 a vector of 0 through 30
bins2 = [round(bins2[i]-0.1) for i in range(31)]
#print(bins2) # Check
#print(len(bins2)) # Check
#print(len(spike_count)) # Check

# Computing a Poisson Distribution
mu = nspikes_mean # Mean
p_poisson = [0 for i in range(31)] # Initialize p_poisson
for k in range(30): # Vector of Poisson distribution
    p_poisson[k] = scipy.stats.poisson.pmf(k,mu)
#print(p_poisson) # Check
# Poisson = np.exp(-mu)*pow(mu,k)*math.factorial(k)

# Computing a Geometric Distribution
# rho/(1-rho) = nspikes_mean
# => rho = nspikes_mean/(1-nspikes_mean)
rho = nspikes_mean/(1+nspikes_mean)
p_geometric = [(1-rho)*pow(rho,k) for k in bins2]

# Plot Sample Spike Counts and Distributions
fig2 = plt.figure() # Create figure
ax2 = plt.axes() # Create plot axes
ax2g = fig2.gca()
ax2.plot(bins2,p_samp,linestyle='-', marker='o', color='b',label='Sample') # Plot p_samp
ax2.plot(bins2,p_poisson,linestyle='--', marker='o', color='r',label='Poisson') # Plot p_poisson
ax2.plot(bins2,p_geometric,linestyle='-.', marker='o', color='g',label='Geometric') # Plot p_geometric
plt.title("Histogram 2") # Plot title
plt.xlabel("Number of Spikes") # x-axis label
plt.ylabel("Frequency of Number of Spikes"); # y-axis label
plt.grid() # Plot grid
plt.legend(); #Legend
plt.show() # Show the plot

# The number of counts appears to be roughly exponentially distributed

# Comparing P(x>=30) for the geometric distribution & the sample distribution
pa= sum(p_geometric)-p_geometric[30] # P(x<30) = sum_i_30(P(x=i)) - P(x=30)
p0 = 1 - pa # P(x>=30) = 1 - P(x<30)
# pmax: ratio of P(x>=30) for geom. and sample distrb.
#pmax = nspikes_max*p0 # Incorrect value of pmax
#pmax = p0/(1-(sum(p_samp)-p_samp[30])) # Relative probs of P(x>=30)
pmax = p0/p_samp[30] # This approximates the above result
print(pmax)
# pmax ~= 0.0009434818459561622 ~= 0.001
