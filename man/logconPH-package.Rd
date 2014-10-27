\name{logconPH-package}
\alias{logconPH-package}
\alias{logconPH}
\docType{package}
\title{
Computing a Cox PH Model with a Log Concave Baseline
}
\description{
Fits the log concave NPMLE with logconcave(), either as univaraite data or as a Cox-PH model with the a logconcave baseline distribution. Estimated densities, probabilities and 
quantiles can be found for a given fit via dLC, pLC and qLC. 
}
\details{
\tabular{ll}{
Package: \tab logconPH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-16\cr
License: \tab NA  \cr
}
Computes the univariate log concave estimator. Built specifically for interval censored data. The times should follow standard interval censored notation, i.e. an n by 2 matrix where [i, 1] is the left side of the interval for subject i and [i,2] is right side. Left censoring is handled by setting the left side to 0, right censoring is handled by setting the right side to inf. Uncensored observations are handled by setting the left side equal to the right side. Alternatively, this package can also handle purely uncensored data. If all the data is uncensored, you may enter data as a length(n) vector (as apposed to using the current status form of the data, cbind(uncensored_data, uncensrored_data) ). Entering uncensored data in this format leads to significantly faster computation. 

In the case of uncensored data, the motivation for the logconcave shape constraint is to obtain valid density estimates without specificying either a parametric family or a smoothing parameter. For interval censored data, while the density estimates should be consistent, they are fairly unstable. However, survival estimates empirically show a better rate of convergence than the unconstrained NPMLE without having to select a parametric family, making for a much more efficient estimator than the NPMLE, while being more flexible than a parametric esitmator. 
}
\author{
Clifford Anderson-Bergman

Maintainer: Clifford Anderson-Bergman <cianders@uci.edu>
}
\references{
Semi- and Non-parametric Methods for Interval Censored
  Data with Shape Constraints (Anderson-Bergman 2014)}
\keyword{ Log Concave, CoxPH, Shape Constrained }

\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)      				
  # Estimates the median

  simData <- sim_Censored(n = 400)
  # Simulates current status data
  
  fit = logconcave(simData)
  # Fits a log concave estimator to an interval censored sample

  pLC(0.5, fit)
  # Estimates the cdf at t = 0.5
  
  plotLC(fit,  'surv')
  # Plots the estimated survival function.
  # Options for second argument are 'pdf', 'cdf' and 'surv'

  simData <- simPH_Censored()
  # Simulates current status data from a Cox-PH model
  
  fit <- logconcave(times = simData$times, covariates = simData$x)
  # Fits a Cox-PH model with a logconcave baseline distribution
  
  plotLC(fit, covars = c(0,0) )
  # Plots the estimated baseline survival function
  
  linesLC(fit, covars = c(1,1), col = 'red')
  # Plots the estimated survival function with x1 = 1, x2 = 1
}
