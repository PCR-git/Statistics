# Probability - HW 4B
# Peter Racioppo
# 103953689

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import scipy.linalg

# ----------------------------

# (1)

size = 1000

mean1 = []
mean2 = []
for i in range(size):
    mean1.append(2)
    mean2.append(3)
mean = np.squeeze([mean1,mean2])

cov = [[1, 0.5], [0.5, 2]]


# x, y = np.random.multivariate_normal([2,3], cov, size).T
# print(x)
# print(y)

mu, sigma = 0, 1 # mean and standard deviation
z1 = np.random.normal(mu, sigma, size)
z2 = np.random.normal(mu, sigma, size)
z = [z1,z2]

A = scipy.linalg.sqrtm(cov)
xi = np.dot(A,z)
x = np.squeeze(xi) + mean

x1 = []
x2 = []
for i in range(size):
    x1.append(x[0,i])
    x2.append(x[1,i])

fig = plt.figure() # Create figure
ax = plt.axes() # Create plot axes
axg = fig.gca()
plt.scatter(x1, x2, marker='x', color='b')
plt.axis('equal')
plt.title("Qs 1: Multivariate Gaussian Distribution")
plt.xlabel("x1")
plt.ylabel("x2")
plt.grid(linestyle='-', linewidth=0.5)
plt.show()

# ----------------------------

# (2)

size = 1000
U = np.random.choice(3, size, p=[0.5, 0.3, 0.2])

count = np.bincount(U)
u1 = count[0]
u2 = count[1]
u3 = count[2]

mean1 = [0, 0]
mean2 = [1, 1]
mean3 = [-2, 2]

cov1 = [[0.5, 0], [0, 0.5]]
cov2 = [[2, -0.1], [-0.1, 0.25]]
cov3 = [[1, 0.5], [0.5, 1]]

x1, y1 = np.random.multivariate_normal(mean1, cov1, u1).T
x2, y2 = np.random.multivariate_normal(mean2, cov2, u2).T
x3, y3 = np.random.multivariate_normal(mean3, cov3, u3).T

r1 = []
r2 = []
r3 = []
for i in range(u1):
    r1i = pow(x1[i],2) + pow(y1[i],2)
    r1.append(r1i)
for i in range(u2):
    r2i = pow(x2[i],2) + pow(y2[i],2)
    r2.append(r2i)
for i in range(u3):
    r3i = pow(x3[i],2) + pow(y3[i],2)
    r3.append(r3i)
r1_m = sum(r1)/u1
r2_m = sum(r2)/u2
r3_m = sum(r3)/u3
print(r1_m)
print(r2_m)
print(r3_m)

# Empirical Values of E(||ri||^2)
# E(||r1||^2) = 1.0190
# (Theoretical value: 1)
# E(||r2||^2) = 4.3850
# (Theoretical value: 4.25)
# E(||r3||^2) = 10.2055
# (Theoretical value: 10)

fig = plt.figure()
ax = plt.axes()
axg = fig.gca()
plt.scatter(x1, y1, marker='.', color='b')
plt.scatter(x2, y2, marker='.', color='r')
plt.scatter(x3, y3, marker='.', color='g')
plt.axis('equal')
plt.title("Qs 2: Gaussian Mixture Model")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(linestyle='-', linewidth=0.5)
plt.show()

# ----------------------------
