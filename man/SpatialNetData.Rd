\name{SpatialNetData}
\docType{data}
\alias{SpatialNetData}
\title{Simulated epidemic data set from a distance and network-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using a combined distance and network-based SIR individual-level model (ILM) with power-law distance kernel. The model has one binary susceptible covariate and the infectivity rate is given by:

\deqn{\lambda_{jt} = (\alpha_{0} + \alpha_{1}z_{j}) \sum_{i \in I_{t}}{d_{ij}^{-\beta_{1}}+\beta_{2}c_{ij}}}{\lambda_{jt} = (\alpha_0 + \alpha_1 z_j) \sum_{i in I_t}{d_{ij}^{-\beta_1}+\beta_2 c_{ij}}}

The infectious period is assumed to follow an exponential distribution with rate \eqn{\delta}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 0.008}{\alpha_0 = 0.008}, \eqn{\alpha_{1} = 0.2}{\alpha_1 = 0.2}, \eqn{\beta_{1} = 2}{\beta_1 = 2}, \eqn{\beta_{2} = 0.5}{\beta_2 = 0.5} and \eqn{\delta = 2}. 

The data set file is a list of an object of class "datagen" that contains of \code{type}, \code{kerneltype}, \code{epidat}, \code{location} and \code{network}, and the covariate matrix.
}

\usage{SpatialNetData}
\format{
It is a list of an object of class ``datagen'' that contains the following:
\if{html}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``both'' distance and network kernel functions.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals.
  }	
  \item{network:}{
The undirected binary contact network matrix. 
  }	
}
}
\if{text}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``both'' distance and network kernel functions.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals.
  }	
  \item{network:}{
The undirected binary contact network matrix. 
  }	
}
}
\if{latex}{
 \itemize{
 \item[]
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``both'' distance and network kernel functions.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals.
  }	
  \item{network:}{
The undirected binary contact network matrix. 
  }	
}
}
}
and a \eqn{50 \times 2}{50 by 2} matrix of the covariates represents the unity intercept and the binary covariate z.
}

\keyword{datasets}