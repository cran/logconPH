\name{sim_Censored}
\alias{sim_Censored}
\title{Simulate current status data from a beta(2,2) distribution}
\usage{
  sim_Censored(n = 100)
}
\arguments{
  \item{n}{Number of samples simulated}
}
\description{
  
  Simulates current status data from a beta(2,2) distribution. Inspection times follow a uniform(0,1) distribution.
  
}
\examples{
  simData <- sim_Censored()
  # Simulates current status data
  
  fit <- logconcave(simData)
  # Fits a log concave fit
}
\author{
  Clifford Anderson-Bergman
}