\name{plot.datagen}
\alias{plot.datagen}
\title{S3 method to plot summary of epidemic.}
\description{
\code{plot.datagen} is an S3 method that provides different plots summarizing epidemic of class \code{datagen}.
}

\usage{
\method{plot}{datagen}(x, plottype, time.index = NULL, ...)
}

\arguments{
  \item{x}{an object that has a class of \code{datagen}.}
  \item{plottype}{History or propagation plots of the epidemic.}
  \item{time.index}{used for obtaining propagation plots at specific infection time points rather than at every infection time. A plot at the first time point is included in any chosen set of infection time points. The default is to provide propagation plot at every infection time point. It is not required when \code{plottype} is set to ``history''.}
  \item{...}{additional arguments that are passed to the generic \code{plot} function.}
}

\details{
If \code{plottype} is set to ``history'', the function produces epidemic curves of infection and removal times (SIR model), or infection, notification, and removal times (SINR model). If it is set to ``propagation'', the function produces plots of the propagation of epidemic over time based on the setting option of the argument \code{kerneltype}. With the network kernel, the function plots all the connections between individuals (in gray) and overlays these with the epidemic pathway direction over time. With the distance kernel, the function plots the epidemic dispersion over time. It shows the changes in the individual status that related to the chosen compartmental framework in \code{type}.
}

\value{
plot
}

\seealso{
  \code{\link{contactnet}, \link{datagen}}.
}

\keyword{plot}