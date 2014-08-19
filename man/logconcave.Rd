\name{logconcave}
\alias{logconcave}
\title{Cox PH model with Log Concave Baseline}
\usage{
  logconcave(times, covariates, aug = TRUE)
}
\arguments{
  \item{times}{Time of event for survival data. If all times are uncensored, can be entered a vector. If times are censored, should be a nx2 matrix of interval censored data}
  \item{covariates}{Optional covariates argument to be fit a Cox-PH model with a logconcave baseline}
  \item{aug}{Should the baseline logconcave density be augmented? If aug = FALSE, algorithm may fail in certain situations}
}
\description{

Fits a semi-parametric Cox PH model with a log-concave baseline distribution for interval censored or uncensored data. If no covariates are supplied, logconcave() will fit the log concave NPMLE. 

For interval censoring, the data should follow standard interval censored notation, i.e. times are entered as an n by 2 matrix where [i, 1] is the left side of the interval for subject i and [i,2] is right side. Left censoring is handled by setting the left side to 0, right censoring is handled by setting the right side to inf. Uncensored observations are handled by setting the left side equal to the right side. Alternatively, this package can also handle purely uncensored data. If all the data is uncensored, you may enter data as a length(n) vector (as apposed to using the current status form of the data, cbind(uncensored_data, uncensrored_data) ). Entering uncensored data in this format leads to significantly faster computation. 

In the case of uncensored data, the motivation for the logconcave shape constraint is to obtain valid density estimates without specificying either a parametric family or a smoothing parameter. For interval censored data, while the density estimates should be consistent, they are fairly unstable. However, survival estimates empirically show a better rate of convergence than the unconstrained NPMLE without having to select a parametric family, making for a much more efficient estimator than the NPMLE, while being more flexible than a parametric esitmator. 

}
\details{
Uses a sequential quadratic programming algorithm, paired with a univariate 
  optimization step. For more details, see Semi- and Non- Parametric Methods for Interval Censored Data
  with Shape Constraints (Anderson-Bergman 2014)
}
\examples{
  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)      				
  # Estimates the median

  data(menodata)
  # source a current status dataset. This is a 2423 x 2 matrix
  
  fit = logconcave(menodata)
  # Fits a log concave estimator to an interval censored sample

  pLC(50, fit)
  # Estimates the cdf at t = 50
  
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
\author{
  Clifford Anderson-Bergman <cianders@uci.edu>
}