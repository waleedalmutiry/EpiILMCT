\name{NetworkData}
\docType{data}
\alias{NetworkData}
\title{Simulated epidemic data set from a network-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using a network-based SIR individual-level model (ILM) with a contact network that was generated using the power-law model with parameters \eqn{\beta=1.8} and \eqn{\alpha=1}. The model has one binary susceptible covariate and the infectivity rate is given by:

\deqn{\lambda_{jt} = (\alpha_{0} + \alpha_{1}z_{j}) \sum_{i \in I_{t}}{c_{ij}}}{\lambda_{jt} = (\alpha_0 + \alpha_1 z_j) \sum_{i in I_t}{c_{ij}}}

The infectious period follows a gamma distribution \eqn{\Gamma(4,\delta)}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 0.08}{\alpha_0 = 0.08}, \eqn{\alpha_{1} = 0.5}{\alpha_1 = 0.5} and \eqn{\delta=2}.

The data set file is a list of an object of class "datagen" that contains of \code{type}, \code{kerneltype}, \code{epidat}, \code{location} and \code{network}, and the covariate matrix.
}

\usage{data(NetworkData)}
\format{
It is a list of an object of class ``datagen'' that contains the following:
\if{html}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``network'' kernel function.
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
The ``network'' kernel function.
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
The ``network'' kernel function.
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