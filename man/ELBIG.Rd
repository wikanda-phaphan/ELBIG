\name{ELBIG}
\alias{ELBIG}
\title{The Parameter Estimation for Re-parameterized Length-Biased Inverse Gaussian Distribution}
\usage{
MME(X)
MLE(X)
Mill_Vibration
}
\arguments{
\item{X}{vector of data.}   
\item{n}{a number of observations.}   
\item{lambda}{a value of the parameter lambda.}   
\item{theta}{a value of the parameter theta.}   
}
\value{
\strong{The Parameter Estimation for Re-parameterized Length-Biased Inverse Gaussian Distribution}
\item{MME}{gives the value of parameter estimates by the method of moments.}
\item{MLE}{gives the value of parameter estimates by the maximum likelihood method.}
\item{Mill_Vibration}{gives the vibration of the vertical roller mill in 60 minutes, contains the time (a.m.) and values of mill vibration (um), collected on 10 February 2019 from Phaphan (2021).}
}
\description{
The package ELBIG provides functions for parameter estimation for re-parameterized length-biased inverse Gaussian distribution with two estimation methods: the maximum likelihood method, the method of moments.
}

\details{
\strong{Re-parameterized Length-Biased Inverse Gaussian (LBIG) Distribution}

Define X as a positive random variable with the length-biased inverse Gaussian distribution based on Ahmed et al. (2008). A density function can be written in this form:  

\code{fLBIG=(1/(theta*sqrt(2*pi)))*((theta/X)^(1/2))*exp(-0.5*(sqrt(X/theta)-(lambda*sqrt(theta/X))))},
where lambda, theta>0.

}

\examples{
#install.packages("devtools")
library(devtools)
#install_github("wikanda-phaphan/ELBIG")
library(ELBIG)

? ELBIG
########
X <- c(7.7,1.2,4,2,0.5)
MLE(X)
$Lamdahat
[1] 0.7767996
$Thetahat
[1] 1.733454

MME(X)
$Lamdahat
[1] 1.090257
$Thetahat
[1] 1.473503
}

\references{
Folks, J. L. & Chhikara, R. S. (1978). The inverse Gaussian distribution and its statistical application - a review. Journal of the Royal Statistical Society.Serie B, 40, 263–289.

Ahmed, S.E., Budsaba, K., Lisawadi, S., & Volodin, A. (2008). Parametric estimation for the Birnbaum-Saunders lifetime distribution based on new parametrization. Thailand Statistician, 6(2), 213-240.

Phaphan, W. (2021). R package for the two-parameters crack distribution. International Journal of Mathematics and Computer Science, 16(4), in press.
}