# Probability - HW 5
# Peter Racioppo
# QS 8

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import scipy.linalg
import scipy.integrate

# ----------------------------

Lambda = 1/2
Sigma = 10
N = 1000

y_i = -50
y_f = 250
y_f2 = 120

Y = []
for i in range(y_f-y_i+1):
    Y.append(y_i+i)

Y2 = []
for i in range(y_f2-y_i+1):
    Y2.append(y_i+i)

X = []
for i in Y:
    X.append(i*(Lambda/2))

X_d = np.random.exponential(1/Lambda,N)
W_d = np.random.normal(0, Sigma, N)

Y_d = []
for i in range(N):
    Y_d.append(pow(X_d[i],2)+W_d[i])

X_lin = []
for i in Y:
    X_lin.append((Lambda/2)*i)

X_hat = []
for i in Y:
    X_hat.append(pow(np.maximum(0,i),0.5))

X_mmse = []
for i in Y2:
    y = i
    integrand1 = lambda x: x*np.exp(-pow((y-pow(x,2)),2)/(2*pow(Sigma,2)))*np.exp(-Lambda*x)
    integrand2 = lambda x: np.exp(-pow((y-pow(x,2)),2)/(2*pow(Sigma,2)))*np.exp(-Lambda*x)
    integral1 = scipy.integrate.quad(integrand1, 0, np.inf, epsabs=1e-10,epsrel=1e-10)
    integral2 = scipy.integrate.quad(integrand2, 0, np.inf, epsabs=1e-10,epsrel=1e-10)
    X_mmse.append(integral1[0]/integral2[0])

fig = plt.figure()
ax = plt.axes()
axg = fig.gca()
plt.plot(Y,X_lin)
plt.plot(Y,X_hat)
plt.plot(Y2,X_mmse)
plt.scatter(Y_d,X_d, marker='.', color='b')

plt.title("Estimation") # Plot title
plt.xlabel("Y") # x-axis label
plt.ylabel("X"); # y-axis label
# plt.grid() # Plot grid

plt.show()
