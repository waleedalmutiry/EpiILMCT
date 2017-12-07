\name{SpatialData}
\docType{data}
\alias{SpatialData}
\title{Simulated epidemic data set from a distance-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using distance-based SIR individual-level model (ILM) with power-law distance kernel. The model has one binary susceptible covariate and the infectivity rate is given by:

\deqn{\lambda_{jt} = (\alpha_{0} + \alpha_{1}z_{j}) \sum_{i \in I_{t}}{d_{ij}^{-\beta}}}{\lambda_{jt} = (\alpha_0 + \alpha_1 z_j) \sum_{i in I_t}{d_{ij}^{-\beta}}}

The infectious period is assumed to follow exponential distribution with rate \eqn{\lambda}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 0.08}{\alpha_0 = 0.08}, \eqn{\alpha_{1} = 2}{\alpha_1 = 2}, \eqn{\beta = 2} and \eqn{\delta = 2}. 

The data set file contains the XY coordinates of individuals, the Euclidean distance matrix, the covariate matrix and the simulated epidemic data set.
}

\usage{SpatialData}
\format{
A list of the simulated data:    
   \describe{
  \item{loc}{
The XY coordinates of individuals.
  }	
  \item{dis}{
The Euclidean distance matrix.
  }	
  \item{cov}{
A matrix of two columns represents the unity intercept and the binary covariate z.
  }	
  \item{epi}{
The simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
}
}
\keyword{datasets}