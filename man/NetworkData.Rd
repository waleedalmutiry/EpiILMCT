\name{NetworkData}
\docType{data}
\alias{NetworkData}
\title{Simulated epidemic data set from a network-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using network-based SIR individual-level model (ILM) with a contact network that was generated using the power-law model with parameters \eqn{\beta=1.8} and \eqn{\alpha=1}. The model has one binary susceptible covariate and the infectivity rate is given by:

\deqn{\lambda_{jt} = (\alpha_{0} + \alpha_{1}z_{j}) \sum_{i \in I_{t}}{c_{ij}}}{\lambda_{jt} = (\alpha_0 + \alpha_1 z_j) \sum_{i in I_t}{c_{ij}}}

The infectious period is assumed to follow gamma distribution \eqn{\Gamma(4,\delta)}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 0.08}{\alpha_0 = 0.08}, \eqn{\alpha_{1} = 0.5}{\alpha_1 = 0.5} and \eqn{\delta=2}.

The data set file contains the XY coordinates of individuals, the contact network matrix, the covariate matrix and the simulated epidemic data.
}

\usage{data(NetworkData)}
\format{
A list of the simulated data:    
   \describe{
  \item{loc}{
The XY coordinates of individuals.
  }	
  \item{cov}{
A matrix of two columns represents the unity intercept and the binary covariate z.
  }	
  \item{net}{
The undirected binary contact network matrix.
  }	
  \item{epi}{
The simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
}
}
\keyword{datasets}