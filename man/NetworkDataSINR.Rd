\name{NetworkDataSINR}
\docType{data}
\alias{NetworkDataSINR}
\title{Simulated epidemic data set from a network-based SINR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using a network-based SINR individual-level model (ILM) with a contact network that was generated using the power-law model with parameters \eqn{\beta=1.8} and \eqn{\alpha=1}. The model has one binary susceptible covariate and the infectivity rate is given by:

\deqn{\lambda_{jt} = (\alpha_{0} + \alpha_{1}z_{j}) [\sum_{i: \in I_{t}}{c_{ij}}+ \gamma \sum_{i: \in N_{t}}{c_{ij}}]}{\lambda_{jt} = (\alpha_0 + \alpha_1 z_j) [\sum_{i in I_t}{c_{ij}}+ \gamma \sum_{i in N_t}{c_{ij}}]}

The incubation period is assumed to follow a gamma distribution \eqn{\Gamma(4,\delta^{(inc)})}, and the delay period is assumed to follow an exponential distribution \eqn{\delta^{(delay)}}{\delta^(delay)}. The epidemic was simulated with the following parameter values: \eqn{\alpha_{0} = 0.08}{\alpha_0 = 0.08}, \eqn{\alpha_{1} = 2}{\alpha_1 = 2}, \eqn{\gamma = 1}{\gamma = 1}, \eqn{\delta^{(inc)}=2}{\delta^(inc)=2} and \eqn{\delta^{(delay)}=2}{\delta^(delay)=2}.

The data set file contains the XY coordinates of individuals, the contact network matrix, the covariate matrix and the simulated epidemic data.
}

\usage{data(NetworkDataSINR)}
\format{
A list of the simulated data:    
\if{html}{
   \describe{
  \item{loc}{
The XY coordinates of individuals.
  }	
  \item{net}{
The undirected binary contact network matrix.
  }	
  \item{cov}{
A matrix of two columns represents the unity intercept and the binary covariate z.
  }	
  \item{epi}{
The simulated epidemic with six columns as: the id numbers of individuals, removal times, delay periods, notification times, incubation periods, and infection times.
  }	
}
}
\if{text}{
   \describe{
  \item{loc}{
The XY coordinates of individuals.
  }	
  \item{net}{
The undirected binary contact network matrix.
  }	
  \item{cov}{
A matrix of two columns represents the unity intercept and the binary covariate z.
  }	
  \item{epi}{
The simulated epidemic with six columns as: the id numbers of individuals, removal times, delay periods, notification times, incubation periods, and infection times.
  }	
}
}

\if{latex}{
   \describe{
  \item{loc:}{
The XY coordinates of individuals.
  }	
  \item{net:}{
The undirected binary contact network matrix.
  }	
  \item{cov:}{
A matrix of two columns represents the unity intercept and the binary covariate z.
  }	
  \item{epi:}{
The simulated epidemic with six columns as: the id numbers of individuals, removal times, delay periods, notification times, incubation periods, and infection times.
  }	
}
}

}
\keyword{datasets}