\name{dLC}
\alias{dLC}
\title{Density estimates from a log concave fit}
\usage{
  dLC(x, fit, covars)
}
\arguments{
  \item{x}{A vector of numeric values for which the estimated density will be calculated}
  \item{fit}{Log concave fit (from the logconcave() function)}
  \item{covars}{A matrix of covariate values. Number of columns must match number in original fit}

}
\description{

	Returns the the estimated density from a log concave fit
  
}
\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  dLC(0, fit)  						
  # Estimates the density at the true mode
  
  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  dLC(1, fit, covars = c(0,0))
  # Estimates the baseline density at t = 1
}
\author{
  Clifford Anderson-Bergman
}