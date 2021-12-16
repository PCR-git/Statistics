/*
####################################################
## Stat 202A - Homework 3
## Author: Peter Racioppo
## Date : 10/20/2019
## Description: This script implements QR and Sweep
####################################################
 
###########################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not change your working directory
## anywhere inside of your code. If you do, I will be unable 
## to grade your work since R will attempt to change my 
## working directory to one that does not exist.
## // [[Rcpp::depends(RcppArmadillo)]]
###########################################################
 
*/ 

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
 Sign function for later use 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
double signC(double d){
  return d<0?-1:d>0? 1:0;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 1: Sweep operator 
   ~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
NumericMatrix mySweepC(const NumericMatrix A, int k){
  
  /*
  Perform a SWEEP operation on A with the pivot element A[k,k].
  
  A: a square matrix (mat).
  m: the pivot element is A[k, k]. 
  Returns a swept matrix B (which is k by k).
  
  Note the "const" in front of mat A; this is so you
  don't accidentally change A inside your code.
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  */
  
  int n = A.nrow(); // Rows of A
  NumericMatrix B(n,n); // Initialize B
  
  // Define new K
  // (C++ indices start from 0 while R indices start from 1)
  int K = k-1;
  int m = A(K,K)
  // Perform sweep
  for (int i = 0; i<n; i++){
    for (int j = 0; j<n; j++){
      if (i==K && j==K) {
        B(i,j) = -1/m;
      } else if (i==K || j==K){
        B(i,j) = A(i,j)/m;
      } else {
        B(i,j) = A(i,j) - A(i,K)*A(K,j)/m;
      }
    }
  }
  
  // Return B
  return(B);
  
}

