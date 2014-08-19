\name{simPH_Censored}
\alias{simPH_Censored}
\title{Simulate current status data from Cox-PH model}
\usage{
  simPH_Censored(n = 100, b1 = 0.5, b2 = -0.5, shape = 2)
}
\arguments{
  \item{n}{Number of samples simulated}
  \item{b1}{Value of first regression coefficient}
  \item{b2}{Value of second regression coefficient}
  \item{shape}{Shape parameter of baseline gamma distribution}
}
\description{
  
  Simulates current status data from a Cox-PH model with a gamma baseline distribution. Used for demonstration of use of logconcave function. 
  
}
\examples{
  simData <- simPH_Censored()
  # Simulates censored data from a Cox-PH model
  
  fit <- logconcave(times = simData$times, covariates = simData$x)
  # Fits a Cox-PH model with a logconcave baseline distribution
}
\author{
  Clifford Anderson-Bergman
}