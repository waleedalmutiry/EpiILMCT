\name{print.epictmcmc}
\alias{print.epictmcmc}
\title{Prints the contents of an \code{epictmcmc} object to the console }
\description{
  \code{print.epictmcmc} is an S3 method that prints the content of an S3 object of class \code{epictmcmc} to the console. 
}
\usage{
\method{print}{epictmcmc}(x, digits = 6, ...)
}
\arguments{
  \item{x}{an S3 object of class \code{epictmcmc} (i.e. the output of the \code{epictmcmc} function).}
  \item{digits}{the number of printed digits of the estimates. Default value is 6. }
  \item{...}{additional arguments that are passed to the generic \code{print} function }
}
\value{
  Formatted output of a \code{epictmcmc} object.
}

\seealso{
  \code{\link{epictmcmc}}, \code{\link{summary.epictmcmc}}, \code{\link{plot.epictmcmc}}.
}
