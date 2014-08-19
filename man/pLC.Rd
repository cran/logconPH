\name{pLC}
\alias{pLC}
\title{Probability estimates from a log concave fit}
\usage{
  pLC(x, fit, covars)
}
\arguments{
  \item{x}{A vector of numeric values for which the estimated probabilities will be calculated}
  \item{fit}{Log concave fit (from the logconcave() function)}
  \item{covars}{A matrix of covariate values. Number of columns must match number in original fit}
}
\description{

	Returns the the estimated probability from a log concave fit
  
}
\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  pLC(0, fit)  						
  # Estimates the cdf at the true mode

  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  pLC(1, fit, covars = c(0,0))
  # Estimates the baseline probability at t = 1
}
\author{
  Clifford Anderson-Bergman
}