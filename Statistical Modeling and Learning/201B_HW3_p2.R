## STATS 201B
## Homework 3
## Peter Racioppo

## --------------------------------------------
## Problem 7

## Load the data:
# data = read.delim("~/Desktop/student_score.txt", " ")
# data$mech.vecs.alg.analy.stat
# mat = as.matrix(data)

x = c(7, 51, 43, 17, 22, 44, 69, 53, 53, 53, 49, 41, 61, 49, 64, 59, 70, 68, 62, 56, 
      34, 42, 50, 47, 29, 46, 40, 47, 29, 17, 0, 40, 21, 9, 14, 32, 45, 49, 57, 64,
      49, 57, 47, 39, 26, 52, 64, 60, 63, 54, 44, 61, 52, 62, 46, 36, 59, 51, 45, 51,
      42, 60, 54, 49, 33, 5, 30, 44, 36, 18, 22, 58, 53, 56, 41, 18, 51, 40, 56, 30, 
      41, 63, 49, 46, 34, 48, 38, 41, 44, 33, 31, 42, 48, 54, 68, 42, 69, 61, 55, 45,
      46, 49, 53, 59, 37, 63, 63, 65, 70, 63)

X = t(matrix(x,nrow=5,ncol=22))
n = dim(X)[1]
m = round(n/2)

f_eigen_ratio <- function(X){
  eigen_vals = eigen(cor(X))$values
  return(max(eigen_vals)/sum(eigen_vals))
}

library(bcaboot)
bca_jack = bcajack(X, B = 2000, func = f_eigen_ratio, m = m, verbose = FALSE)
eigen_ratio_stats = bca_jack$ustats
eigen_ratio_lims = bca_jack$lims

print("Average Bootstrap Eigen-ratio:")
print(eigen_ratio_stats[1])
print("Bootstrap Standard Error of Eigen-ratio:")
print(eigen_ratio_stats[2])

print("BCa confidence limits:")
print(eigen_ratio_lims)

