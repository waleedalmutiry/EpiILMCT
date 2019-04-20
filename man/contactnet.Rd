\name{contactnet}
\alias{contactnet}

\title{
Generate binary adjacency matrix contact networks
}
\description{
This function allows the user to generate contact network with binary adjacency matrices from spatial (power-law or Cauchy) and random network models.
}
\usage{
contactnet(type, num.id = NULL, location = NULL, beta, nu = NULL)
}

\arguments{

  \item{type}{
type of contact network model, with the choice of ``powerlaw'' for the power-law model, ``Cauchy'' for the Cauchy model, and ``random'' for randomly generated contact network from a Bernoulli distribution with parameter beta.
}
  \item{num.id}{
	 the number of individuals (nodes). It has to be specified when \code{type} is set to ``random''.
}
  \item{location}{
	 a matrix of the XY coordinates of individuals. Required when \code{type} is set to ``powerlaw'' or ``Cauchy''. Default is NULL.
}
  \item{beta}{
spatial parameter of the spatial contact network model or probability parameter of the random network model (>0).
}
  \item{nu}{
scale parameter of the power-law contact network model (>0). Default value is 1.
}

}


\details{
The contact networks considered here are undirected, ``unweighted'' and have binary adjacency matrix, so that \eqn{c_{ij}=c_{ji}} for \eqn{i \ne j; i,j = 1,..,N}, and each element of the contact network is defined as \eqn{c_{ij} = 1} if a connection exist between individuals \eqn{i} and \eqn{j}, and 0 otherwise. We can either generate spatial networks which have connections more likely to be present between two individuals closed together than far apart, or random contact networks. 

For spatial contact networks, there are two options of specifying the model for the probability of connections between individuals. We use a generalized version of the power-law contact network model of Bifolchi et al. (2013) for the first option, in which
the probability of a connection between individuals \eqn{i} and \eqn{j} is given by:
\deqn{p(c_{ij}=1) = 1- e^{-\nu(d_{ij}^{-\beta})} , \quad \nu, \beta > 0,}{p(c_{ij}=1) = 1- exp(-\nu(d_{ij}^{-\beta})) , \nu, \beta > 0,}
where \eqn{d_{ij}} is the Euclidean distance between individuals \eqn{i} and \eqn{j}; \eqn{\beta}{\beta} is the spatial parameter; and \eqn{\nu}{\nu} is the scale parameter. 

The second option is a Cauchy model as used in Jewell et al., (2009). The probability of a connection between individuals \eqn{i} and \eqn{j} is given by:

\deqn{p(c_{ij}=1) = 1- e^{-\beta/(d_{ij}^{2} + \beta^{2})} , \quad \beta > 0,}{p(c_{ij}=1) = 1- exp(-\beta/(d_{ij}^{2} + \beta^{2})), \beta > 0,}

where \eqn{d_{ij}} is the Euclidean distance between individual \eqn{i} and \eqn{j}; and \eqn{\beta} is the spatial parameter.

The final option is a random contact network in which the probability of a connection is generated from a Bernoulli distribution with probability equal to \eqn{\beta}.

}


\value{
An object of class ``contactnet'' that has a list of:
\if{html}{
	\describe{
		\item{location}{A matrix of the XY coordinates of individuals.}
		\item{contact.network}{Binary contact network adjacency matrix.}
		\item{type}{The type of contact network model.}
	}
}
\if{text}{
	\describe{
		\item{location}{A matrix of the XY coordinates of individuals.}
		\item{contact.network}{Binary contact network adjacency matrix.}
		\item{type}{The type of contact network model.}
	}
}
\if{latex}{
	 \itemize{
		\item[]
		\describe{
			\item{location}{A matrix of the XY coordinates of individuals.}
			\item{contact.network}{Binary contact network adjacency matrix.}
			\item{type}{The type of contact network model.}
		}
	}
}
}

\references{
  
  Jewell, C. P., Kypraios, T., Neal, P., and Roberts, G. O. (2009). Bayesian analysis for emerging infectious
diseases. \emph{Bayesian Analysis}, 4(3):465-496.
     
  Bifolchi, N., Deardon, R., and Feng, Z. (2013). Spatial approximations of network-based individual level infectious disease models. \emph{Spatial and Spatio-temporal Epidemiology}, 6:59-70.
}

\examples{
set.seed(12345)
loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, nu = 0.5)
net2<- contactnet(type = "Cauchy", location = loc, beta = 0.5)
net3<- contactnet(type = "random", num.id = 50, beta = 0.08)
}
\keyword{Simulation}