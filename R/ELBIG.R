#' @title The Parameter Estimation for Re-parameterized Length-Biased Inverse Gaussian Distribution
#'
#' @description The package ELBIG provides functions for parameter estimation for re-parameterized length-biased inverse Gaussian distribution with two estimation methods: the maximum likelihood method, the method of moments.
#' 
#' @param X vector of data
#' @param n a number of observations
#' @param lambda a value of the parameter lambda
#' @param theta a value of the parameter theta
#'
#' @return ?
#'
#' @examples
#' X <- c(7.7,1.2,4,2,0.5)
#' MLE(X)
#'
###############The Parameter Estimation for Re-parameterized Length-Biased Inverse Gaussian Distribution#####################
MME<-function(X){
  n<-length(X)
  va<-var(X)* (n - 1) / n
  me<-mean(X)
  lambhat=((me^2)-(2*va)+(me*sqrt((me^2)+(4*va))))/(2*va)
  thehat=((-me)+sqrt((me^2)+(4*va)))/2
  return(list(Lamdahat=lambhat,Thetahat=thehat))
}
#MME(X)
MLE=function(X){
  n<-length(X)
  me<-mean(X)
  T <- sum(1/X)
  lambhat=n/(T*me-n)
  thehat=(T*me-n)/T
  return(list(Lamdahat=lambhat,Thetahat=thehat))
}
#MLE(X)
Mill_Vibration<- matrix(c(
  05.00,77.01,
  05.01,77.01,
  05.02,68.72,
  05.03,71.38,
  05.04,73.99,
  05.05,66.87,
  05.06,71.43,
  05.07,81.3,
  05.08,82.33,
  05.09,76.78,
  05.10,70.31,
  05.11,74.88,
  05.12,69.25,
  05.13,66.57,
  05.14,68.97,
  05.15,67.45,
  05.16,72.09,
  05.17,72.09,
  05.18,71.19,
  05.19,63.63,
  05.20,93.72,
  05.21,82.32,
  05.22,71.93,
  05.23,77.75,
  05.24,92.91,
  05.25,76.27,
  05.26,91.95,
  05.27,73.65,
  05.28,71.42,
  05.29,76.77,
  05.30,66.5,
  05.31,75.61,
  05.32,71.14,
  05.33,69.7,
  05.34,72.48,
  05.35,71.53,
  05.36,74.37,
  05.37,77.88,
  05.38,69.67,
  05.39,66.85,
  05.40,78.9,
  05.41,74.88,
  05.42,77.37,
  05.43,85.67,
  05.44,66.05,
  05.45,66.06,
  05.46,72.03,
  05.47,74.17,
  05.48,70.65,
  05.49,75.99,
  05.50,72.54,
  05.51,85.45,
  05.52,74.6,
  05.53,76.37,
  05.54,70.02,
  05.55,77.66,
  05.56,88.71,
  05.57,86.86,
  05.58,76.82,
  05.59,68.56),
  nrow = 60, byrow = TRUE)
