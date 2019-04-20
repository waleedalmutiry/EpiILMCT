\name{plot.epictmcmc}
\alias{plot.epictmcmc}
\title{ Plot the output of \code{epictmcmc} object}
\description{
  \code{plot.epictmcmc} is an S3 method that plots the output of an S3 object of class \code{epictmcmc}. 
}
\usage{
  \method{plot}{epictmcmc}(x, epi = NULL, start = NULL, end = NULL, 
  thin = NULL, ...)
}
\arguments{
  \item{x}{an S3 object of class \code{epictmcmc} (i.e. the output of the \code{epictmcmc} function).}
  \item{epi}{epidemic data of class \code{datagen}.}
  \item{start, end, thin}{options for creating \code{\link{mcmc}} object.}
  \item{...}{additional arguments that are passed to the generic \code{plot} function.}
}
\value{
  Trace and density plots of the posterior distribution of the model parameters. Additionally, in the case of unobserved event times, a plot of the average posterior and 95\% percentile intervals of the inferred event times are produced in new plot window.
}
\seealso{
  \code{\link{epictmcmc}}, \code{\link{print.epictmcmc}}, \code{\link{summary.epictmcmc}}, \code{\link{mcmc}}, \code{\link{plot.mcmc}}.
}
\keyword{plot}