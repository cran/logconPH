\name{linesLC}
\alias{linesLC}
\title{Draws Lines for Logconcave Fit}
\usage{
  linesLC(fit, funtype = 'surv', covars = NA, ...)
}
\arguments{
  \item{fit}{A fit object returned from a call to logconcave}
  \item{funtype}{Type of function plotted. Choices are 'surv', 'pdf' or 'cdf'. If missing, 'surv' in default}
  \item{covars}{CoxPH covariates}
\item{...}{Additional arguments to be passed to lines}
  
}
\description{
  
  Plots the fit of a log concave object
  
}
\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
  
  plotLC(fit)    					
  # Plots the estimated survival distribution
  
  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  plotLC(fit, covars = c(0,0), funtype = 'cdf')
  # Plots the estimated baseline cdf
  
  linesLC(fit, covars = c(1,1), funtype = 'cdf', col = 'red')
  # Plots the estimates cdf with covariates c(1,1)}
\author{
  Clifford Anderson-Bergman
}