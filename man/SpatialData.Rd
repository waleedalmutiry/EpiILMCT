\name{SpatialData}
\docType{data}
\alias{SpatialData}
\title{Simulated epidemic data set from a distance-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using a distance-based SIR individual-level model (ILM) with power-law distance kernel. The model has no susceptible covariates and the infectivity rate is given by:

\deqn{\lambda_{jt} = \alpha_{0} \sum_{i \in I_{t}}{d_{ij}^{-\beta}}}{\lambda_{jt} = \alpha_0 \sum_{i in I_t}{d_{ij}^{-\beta}}}

The infectious period is assumed to follow an exponential distribution with rate \eqn{\lambda}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 2.0}{\alpha_0 = 2.0}, \eqn{\beta = 3} and \eqn{\delta = 2}. 

The data set file is an object of class "datagen" that has a list of \code{type}, \code{kerneltype}, \code{epidat} and \code{location}.

}

\usage{data(SpatialData)}
\format{
An object of class ``datagen'' that contains the following:
\if{html}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``distance'' kernel function.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals. 
  }
  \item{network:}{
A zero matrix of contact network. 
  }
}
}
\if{text}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``distance'' kernel function.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals. 
  }
  \item{network:}{
A zero matrix of contact network. 
  }
}
}
\if{latex}{
   \describe{
  \item{type:}{
The ``SIR'' compartmental framework.
  }	
  \item{kerneltype:}{
The ``distance'' kernel function.
  }	
  \item{epidat:}{
A matrix of the simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
  \item{location:}{
A matrix of the XY coordinates of individuals. 
  }
  \item{network:}{
A zero matrix of contact network. 
  }
}
}
}
\keyword{datasets}