\name{contactnet}
\alias{contactnet}

\title{
Generate binary adjacency contact networks
}
\description{
This function allows the user to generate binary contact network adjacency matrices from spatial (power-law and cauchy) and random network models.
}
\usage{
contactnet(type, location=NULL, beta, alpha=NULL, plot=FALSE)
}

\arguments{

  \item{type}{
Type of contact network model, with the choice of ``powerlaw''  for
the power-law model, ``Cauchy'' for the Cauchy model, and ``random''  for
randomly generated contact network from a Bernoulli distribution with parameter beta.
}
  \item{location}{
	 A matrix of XY coordinates of individuals.
}
  \item{beta}{
Spatial parameter of the spatial contact network model or the probability parameter of the random network model (>0).
}
  \item{alpha}{
Scale parameter of the power-law contact network model (>0). Defualt value is 1.
}
  \item{plot}{
If \code{plot} set to TRUE, a graphical representation of the generated contact network is produced. The default is FALSE.
}

}


\details{
The contact networks considered here are undirected and binary, so that \eqn{c_{ij}=c_{ji}} for \eqn{i \ne j; i,j = 1,..,N}, and each element of the contact network is defined as \eqn{c_{ij} = 1} if a connection exist between individuals \eqn{i} and \eqn{j}, and 0 otherwise. Here, we consider spatial networks with connections more likely to be present between two individuals close together than far apart, or random contact networks. 

For spatial contact networks, there are two options of specifying the model for the probability of connections between individuals. We use a generalized version of the power-law contact network model of Bifolchi et al. (2013) for the first option, in which
the probability of a connection between individuals \eqn{i} and \eqn{j} is given by:
\deqn{p(c_{ij}=1) = 1- e^{-\alpha(d_{ij}^{-\beta})} , \quad \alpha, \beta > 0,}{p(c_{ij}=1) = 1- exp(-\alpha(d_{ij}^{-\beta})) , \alpha, \beta > 0,}
where \eqn{d_{ij}} is the Euclidean distance between individuals \eqn{i} and \eqn{j}; \eqn{\beta}{\beta} is the spatial parameter; and \eqn{\alpha}{\alpha} is the scale parameter. 

The second option is a Cauchy model as used in Jewell et al., (2009). The probability of a connection between individuals \eqn{i} and \eqn{j} is given by:

\deqn{p(c_{ij}=1) = 1- e^{-\beta/(d_{ij}^{2} + \beta^{2})} , \quad \beta > 0,}{p(c_{ij}=1) = 1- exp(-\beta/(d_{ij}^{2} + \beta^{2})), \beta > 0,}

where \eqn{d_{ij}} is the Euclidean distance between individual \eqn{i} and \eqn{j}; and \eqn{\beta} is the spatial parameter.

The final option is a random contact network in which the probability of a connection is generated from a Bernoulli distribution with probability equal to \eqn{\beta}.

}


\value{
Binary contact network adjacency matrix
}

\references{
  
  Jewell, C. P., Kypraios, T., Neal, P., and Roberts, G. O. (2009). Bayesian analysis for emerging infectious
diseases. \emph{Bayesian Analysis}, 4(3):465-496.
     
  Bifolchi, N., Deardon, R., and Feng, Z. (2013). Spatial approximations of network-based individual level infectious disease models. \emph{Spatial and Spatio-temporal Epidemiology}, 6:59-70.
}

\examples{
library( EpiILMCT )
loc<- matrix(cbind(runif(50,0,10),runif(50,0,10)),ncol=2,nrow=50)
net1<- contactnet(type="powerlaw",location=loc,beta=1.5,alpha=0.5,
       plot=TRUE)
net2<- contactnet(type="Cauchy",location=loc,beta=0.5,plot=TRUE)
net3<- contactnet(type="random",location=loc,beta=0.08,plot=TRUE)
}
\keyword{Simulation}