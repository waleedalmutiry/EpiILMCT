\name{summary.epictmcmc}
\alias{summary.epictmcmc}
\alias{print.summary.epictmcmc}
\title{Summary method for \code{epictmcmc} objects}
\description{Summarize a \code{\link{epictmcmc}} object and return an object of class \code{summary.epictmcmc}.
}
\usage{
\method{summary}{epictmcmc}(object, digits = NULL, start = NULL, end = NULL, thin = NULL, \dots)
\method{print}{summary.epictmcmc}(x, digits, start, end, thin, \dots)
}
\arguments{
  \item{x, object}{an S3 object of class \code{epictmcmc} (i.e. the output of the \code{epictmcmc} function).}
  \item{digits}{the number of printed digits of the estimates. Default value is 6. }
  \item{start, end, thin}{options for creating \code{\link{mcmc}} object.}
  \item{\dots}{potential further arguments (require by generic).}
}

\seealso{
  \code{\link{epictmcmc}}, \code{\link{print.epictmcmc}}, \code{\link{plot.epictmcmc}}.
}

