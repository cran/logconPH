\name{qLC}
\alias{qLC}
\title{Quantiles estimates from a log concave fit}
\usage{
  qLC(p, fit, covars)
}
\arguments{
  \item{p}{A vector of numeric values for which the estimated quantiles will be calculated}
  \item{fit}{Log concave fit (from the logconcave() function)}
  \item{covars}{A matrix of covariate values. Number of columns must match number in original fit}
}
\description{

	Returns the the estimated probability from a log concave fit
  
}
\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)  						
  # Estimates the median

  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  qLC(0.5, fit, covars = c(0,0))
  # Estimates the baseline median
}
\author{
  Clifford Anderson-Bergman
}