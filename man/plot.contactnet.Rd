\name{plot.contactnet}
\alias{plot.contactnet}
\title{ Plot the contact network of \code{contactnet} object}
\description{
  \code{plot.contactnet} is an S3 method that plots the contact network of an S3 object of class \code{contactnet}. 
}
\usage{
  \method{plot}{contactnet}(x, ...)
}
\arguments{
  \item{x}{an S3 object of class \code{contactnet} (i.e. the output of the \code{contactnet} function).}
  \item{...}{additional arguments that are passed to the generic \code{plot} function.}
}
\value{
  Plot of the contact network.
}

\examples{
set.seed(12345)
loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, nu = 0.5)
plot(net1)
net2<- contactnet(type = "Cauchy", location = loc, beta = 0.5)
plot(net2)
net3<- contactnet(type = "random", num.id = 50, beta = 0.08)
plot(net3)
# adding some options to the igraph plot:
plot(net3, vertex.color = "red", vertex.size = 8, edge.color = "black",
vertex.label.cex = 0.7, vertex.label.color = "black") 
}

\seealso{
  \code{\link{contactnet}}.
}
\keyword{plot}