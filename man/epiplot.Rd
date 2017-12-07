\name{epiplot}
\alias{epiplot}

\title{
Epidemic summary plots
}
\description{
Provides different plots summarizing epidemic
}
\usage{
epiplot(type, kerneltype, epidat, location, network=NULL, plottype=NULL, time_index=NULL)
}

\arguments{

  \item{type}{
Type of the compartmental framework, with the choice of ``SIR'' for
Susceptible-Infectious-Removed and ``SINR'' for Susceptible-Infectious-Notified-Removed.
}

  \item{kerneltype}{
The kernel type of the ILMs:``distance'', ``network'' or ``both''.
}

  \item{epidat}{
    Epidemic data. Output from \code{\link{datagen}} or input a data frame in the same format.
}

  \item{location}{
	A matrix or data frame of XY coordinates of individuals.
}
  \item{network}{
Contact network matrix. It is required when \code{kerneltype} is set to ``network''.
}
  \item{plottype}{
History or propagation plots of the epidemic.
}
  \item{time_index}{
Used for obtaining propagation plots at specific infection time points rather than at every infection time. A plot at the first time point is included in any chosen set of infection time points. The default is to provide propagation plot at every infection time point. It is not required when \code{plottype} is set to ``history''.
}
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

\examples{
library( EpiILMCT )
## network-based SIR ILMs:
data(NetworkData)
epiplot(type="SIR",kerneltype="network",epidat=NetworkData$epi,
		location=NetworkData$loc,network=NetworkData$net,
		plottype="propagation",time_index=seq(1,33,by=6))
epiplot(type="SIR",kerneltype="network",epidat=NetworkData$epi,
		location=NetworkData$loc,network=NetworkData$net,
		plottype="history")
## network-based SINR ILMs:
data(NetworkDataSINR)
epiplot(type="SINR",kerneltype="network",epidat=NetworkDataSINR$epi,
		location=NetworkDataSINR$loc,network=NetworkDataSINR$net,
		plottype="propagation",time_index=seq(1,34,by=6))
epiplot(type="SINR",kerneltype="network",epidat=NetworkDataSINR$epi,
		location=NetworkDataSINR$loc,network=NetworkDataSINR$net,
		plottype="history")
## distance-based SIR ILMs:
data(SpatialData)
epiplot(type="SIR",kerneltype="distance",epidat=SpatialData$epi,
		location=SpatialData$loc,plottype="propagation",
		time_index=seq(1,33,by=6))
epiplot(type="SIR",kerneltype="distance",epidat=SpatialData$epi,
		location=SpatialData$loc,plottype="history")
## distance and network-based SIR ILMs:
data(SpatialNetData)
epiplot(type="SIR",kerneltype="both",epidat=SpatialNetData$epi,
		location=SpatialNetData$loc,network= SpatialNetData $net,
		plottype="propagation",time_index=seq(1,21,by=6))
epiplot(type="SIR",kerneltype="both",epidat=SpatialNetData$epi,
		location=SpatialNetData$loc,network= SpatialNetData $net,
		plottype="history")
}
\keyword{plot}