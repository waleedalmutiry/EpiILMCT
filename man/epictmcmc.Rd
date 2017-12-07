\name{epictmcmc}
\alias{epictmcmc}

\title{
Markov Chain Monte Carlo-based tool for analyzing (SIR/SINR) distance-based/network-based individual-level disease data.
}
\description{
Runs a Bayesian data augmented MCMC algorithm for fitting specified models (SIR or SINR) with three possible assumptions regarding the data set: 1. a completely observed epidemic; 2. unknown infection times but known removal times; or 3. unknown removal and infection times but known notification times (SINR only). 
}
\usage{

epictmcmc(epidat, type, kerneltype, kernelmatrix, distancekernel = NULL, datatype, 

blockupdate = NULL, nsim, sus.par = NULL, trans.par = NULL, power.sus = NULL,
	
power.trans = NULL, kernel.par = NULL, spark.par = NULL, delta = NULL, gamma.par = NULL, 
	
suscov = NULL, transcov = NULL, periodproposal=NULL )

}

\arguments{
  \item{epidat}{
Epidemic data. Output from \code{\link{datagen}} or input a data frame in the same format.
}
  \item{type}{
 Type of the compartmental framework, with the choice of ``SIR'' for
Susceptible-Infectious-Removed and ``SINR'' for Susceptible-Infectious-Notified-Removed.
}
  \item{kerneltype}{
Type of kernel function with choice of ``distance'' for distance-based ILM, ``network'' for contact network-based ILM and ``both'' for a combination of network- and distance-based ILM.
}
  \item{kernelmatrix}{
A matrix of distances between individuals when \code{kerneltype} is set to ``distance''. A matrix of edge incidence (0 or 1) or edge weights [0,\eqn{\infty}) representing the contact network when \code{kerneltype} is set to ``network''. A list of two matrices (1: distance matrix, 2: contact network matrix) when \code{kerneltype} is set to ``both''.
}
  \item{distancekernel}{
The spatial kernel type when \code{kerneltype} is set to ``distance'' or ``both''. Choices are ``powerlaw'' for a power-law distance kernel or ``Cauchy'' for a Cauchy distance kernel. 
}
  \item{datatype}{
The ILM data set is analyzed based on three possible assumptions. The possible assumptions are ``known epidemic'' for completely observed epidemic; ``known removal'' times for observing only the removal times for SIR model or removal and notified times for SINR model; or ``unknown removal'' times for observing only the notified times which is applicable only to the SINR model.
}
  \item{blockupdate}{
A vector of the number of the initially observed infection times and the size of blocks for updating the infection times/infectious periods (SIR model), or infection times/incubation periods and removal times/delay periods (SINR model). The default is observing only the first infection time and performing the MCMC updating infection time parameters via single parameter updates. 
}
  \item{nsim}{
The number of MCMC iterations.
}
  \item{sus.par}{
A \eqn{n_{s} \times 5}{n_s by 5} matrix of arguments for updating the parameters of the susceptibility function (>0), where \eqn{n_s} is the number of susceptibility parameters. Default = NULL means the model does not include these parameters. See details for description of required arguments.
}
  \item{trans.par}{
A \eqn{n_{t} \times 5}{n_t by 5} matrix of arguments for updating the parameters of the transmissibility function (>0), where \eqn{n_t} is the number of transmissibility parameters. Default = NULL means the model does not include these parameters. See details for description of required arguments.
}
  \item{power.sus}{
A \eqn{n_{s} \times 5}{n_s by 5} matrix of arguments for updating the power parameters of the susceptibility function (>0), where \eqn{n_s} is the number of susceptibility parameters. Default = NULL means the model does not include these parameters. See details for description of required arguments.
}
  \item{power.trans}{
A \eqn{n_{t} \times 5}{n_t by 5} matrix of arguments for updating the power parameters of the transmissibility function (>0), where \eqn{n_t} is the number of transmissibility parameters. Default = NULL means the model does not include these parameters. See details for description of required arguments.
}
  \item{kernel.par}{
A vector of arguments for updating the spatial parameter of the distance-based kernel when \code{kerneltype} is set to ``distance'' (>0). A \eqn{2 \times 5}{2 by 5} matrix of arguments for updating the spatial and network effect parameters when \code{kerneltype} is set to ``both'' (>0). Default = NULL means the model does not include these parameters (Network-based ILMs). See details for description of required arguments.
}
  \item{spark.par}{
A vector of arguments for updating the spark parameter (>=0). It should be defined if \code{kerneltype} is set to ``network'' with \code{datatype}= ``known removal'' or ``unknown removal''. Default = NULL means the model does not include this parameter. See details for description of required arguments.
}
  \item{delta}{
A vector of arguments for updating the infectious periods for SIR model, or a \eqn{2\times 4}{2 by 4} matrix of arguments for updating the incubation and delay periods for SINR model (>0). Default = NULL means these parameters do not need to be updated, especially when \code{datatype} is set to ``known epidemic''. Note that when \code{type}=``SINR'' with \code{datatype}=``known removal'', the user should just input a vector of arguments for updating the incubation period. The arguments for updating each rate are the known shape and initial value of the rate of the (infectious, incubation or delay) period distribution, the shape and rate parameters of the gamma prior distribution for the rate parameter. 
}
  \item{gamma.par}{
A vector of arguments for updating the notification effect parameter (SINR). The default value(s) is 1, and this means no need to update the parameter. See details for description of required arguments.
}
  \item{suscov}{
Covariate matrix of the susceptibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}
  \item{transcov}{
Covariate matrix of the transmissibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}
  \item{periodproposal}{
A vector/matrix of the proposal distribution parameters of the independence sampler for updating the infectious period (SIR model), or the incubation and/or delay periods (SINR model). Here, we set the proposal distribution to gamma distribution. It is required when \code{datatype} is set to ``known removal'' or ``unknown removal''.
}
}


\details{

Bayesian MCMC is performed to estimate the model parameters and latent variables using data augmentation method. When the \code{datatype} is set to ``known removal'' or ``unknown removal'', we assume the infectious periods (SIR model) or the incubation and/or delay periods (SINR model) to follow a gamma distribution with known shape and unknown rate parameters. We assign a gamma prior distribution for the rate parameter with shape \eqn{a} and rate \eqn{b}. This leads the rate parameter(s) to have a gamma conditional distribution.

Under the SIR model:

\deqn{ \delta|\theta,\boldsymbol{I},\boldsymbol{R} \sim \Gamma(m+a_{\delta},M+b_{\delta}), }{\delta|\theta,I,R ~ \Gamma(m + a_\delta, M + b_\delta),}
where \eqn{\delta}{\delta} is the rate of the infectious period distribution; \eqn{m} is the number of infected individuals; \eqn{M=\sum_{i=1}^{m}{(R_{i}-I_{i})}}{M=\sum_{i=1}^{m}(R_i-I_i)}; and \eqn{a_{\delta}}{a_\delta} and \eqn{b_{\delta}}{b_\delta} are the prior parameters for the infectious period rate.

Under the SINR model:
\deqn{ \delta^{(inc)}|\theta,\boldsymbol{I},\boldsymbol{N},\boldsymbol{R} \sim \Gamma(m+a_{\delta^{(inc)}},M_{inc}+b_{\delta^{(inc)}}),}{\delta^(inc)|\theta, I, N, R ~ \Gamma(m + a_{\delta^(inc)},M_inc + b_{\delta^(inc)}),}
where \eqn{\delta^{(inc)}}{\delta^(inc)} is the rate of the incubation period distribution; \eqn{M_{inc}=\sum_{i=1}^{m}{(N_{i}-I_{i})}}{M_inc=\sum_{i=1}^{m}(N_i-I_i)}; and \eqn{a_{\delta^{(inc)}}}{a_{\delta^(inc)}} and \eqn{b_{\delta^{(inc)}}}{b_{\delta^(inc)}} are the prior parameters for the incubation period rate, and

\deqn{ \delta^{(delay)}|\theta,\boldsymbol{I},\boldsymbol{N},\boldsymbol{R} \sim \Gamma(m+a_{\delta^{(delay)}},M_{delay}+b_{\delta^{(delay)}}),}{\delta^(delay)|\theta, I, N, R ~ \Gamma(m+a_{\delta^(delay)},M_delay+b_{\delta^(delay)}),}
where \eqn{\delta^{(delay)}}{\delta^(delay)} is the rate parameter of the delay period distribution; \eqn{M_{delay}=\sum_{i=1}^{m}{(R_{i}-N_{i})}}{M_delay=\sum_{i=1}^{m}(R_i-N_i)}; and \eqn{a_{\delta^{(delay)}}}{a_{\delta^(delay)}} and \eqn{b_{\delta^{(delay)}}}{b_{\delta^(delay)}} are the prior parameters for the delay period rate.

A Gibbs update (i.e., sampling from the conditional posterior distribution) is then used for the infectious period (SIR model) or the incubation and/or delay period rates (SINR model). An independence sampler is used for updating the infection times and infectious periods (SIR model), and the infection and/or removal times, and the incubation and/or delay periods. To do so undr the SIR model, we update each infection time \eqn{I_{i}}{I_i} by proposing infectious period \eqn{D^{*}_{i}}{D^{*}_i} from a proposal distribution with tuning parameter, such that \eqn{D_{i}\sim}{D_i ~} f(.). Then, the new infection time is just the difference between the fixed known removal time and the new infectious period of the \eqn{i^{th}}{ith} individual. Each infection time/infectious period is updated in this way in turn. The same procedure is applied for updating the missing event times and incubation and delay periods for the SINR model, with their corresponding parameters.

Other model parameters (\code{sus.par}, \code{trans.par}, \code{power.sus}, \code{power.trans}, \code{kernel.par}, \code{spark.par} and \code{gamma.par}) are updated in turn using a random-walk Metropolis Hastings algorithm with normal proposal. The normal proposals must be tuned, via the proposal variance, by the user to achieve good mixing properties. To update each parameter, the user must input the following arguments in order: initial value, prior distribution ("gamma", "half normal" (positive half-normal) or "uniform"), a vector of the prior distribution parameters, and the proposal variance. 

Note that, setting the variance of the normal proposal distribution to zero fixes a parameter at its initial value. This option allows the user to fix such a parameter in the MCMC procedure while updating others.

}
\value{

Returns a list of the MCMC samples as follows:
   \describe{
    \item{1}{The MCMC samples of the model parameters.}
    \item{2}{The MCMC samples of the infection times when \code{datatype} set to ``known removal'' or ``unknown removal''.}
    \item{3}{The MCMC samples of the removal times when \code{datatype} set to ``unknown removal''.}
  }

}

\seealso{
  \code{\link{datagen}, \link{loglikelihoodepiILM}}.
}

\examples{
\dontrun{
	
library(EpiILMCT)
## distance-based SIR continuous-time ILMs ##
data(SpatialData)
## performing the MCMC-tool for analyzing the fully observed spatial data 
## under the SIR distance-based continuous ILM:
mcmcres2<-epictmcmc(epidat=SpatialData$epi,type="SIR",kerneltype="distance", 
kernelmatrix=SpatialData$dis,
distancekernel="powerlaw",datatype="known epidemic",nsim=10000,
sus.par=matrix(c(0.08,2,"uniform","uniform",0,0,100,100,0.5,0.5),ncol=5,nrow=2),
kernel.par=c(2,"gamma",1,0.01,0.5),suscov=SpatialData$cov)

## performing the MCMC-tool for analyzing the partially observed spatial data
## (unknown infection times) under the SIR distance-based continuous ILM:
mcmcres22<-epictmcmc(epidat=SpatialData$epi,type="SIR",kerneltype="distance",
kernelmatrix=SpatialData$dis,distancekernel="powerlaw",
datatype="known removal",nsim=50000,blockupdate=c(1,10),
sus.par=matrix(c(0.08,2,"gamma","gamma",1,4,1,2,0.5,0.5),ncol=5,nrow=2),
kernel.par=c(2,"gamma",4,2,0.5),delta=c(1,2,4,2),
suscov=SpatialData$cov,periodproposal=c(1,1))

## distance-based and network-based SIR ILMs ##
data(SpatialNetData)
## performing the MCMC-tool for analyzing the fully observed spatial and network data
## under the SIR distance-based and network-based continuous-time ILM:
mcmcres3<-epictmcmc(epidat=SpatialNetData$epi,type="SIR",kerneltype="both",
kernelmatrix=list(SpatialNetData$dis,SpatialNetData$net),
distancekernel="powerlaw",datatype="known epidemic",nsim=10000,
sus.par=matrix(c(0.08,0.2,"gamma","gamma",1,1,0.01,0.01,0.5,0.5),ncol=5,nrow=2),
kernel.par=matrix(c(5,0.5,"gamma","gamma",1,1,0.01,0.01,0.5,0.5),ncol=5,nrow=2),
suscov=SpatialNetData$cov)

## network-based SIR ILMs ##
data(NetworkData)
## performing the MCMC-tool for analyzing the fully observed network data 
## under the SIR network-based continuous-time ILM:
mcmcres4<-epictmcmc(epidat=NetworkData$epi,type="SIR",kerneltype="network",
kernelmatrix= NetworkData$net,datatype="known epidemic",nsim=10000,
sus.par=matrix(c(0.08,0.5,"gamma","gamma",1,1,1,1,0.5,0.1),ncol=5,nrow=2),
suscov=NetworkData$cov)

## network-based SINR ILMs ##
data(NetworkDataSINR)
## performing the MCMC-tool for analyzing the fully observed network data
## under the SINR network-based continuous-time ILM:
mcmcres5<-epictmcmc(epidat=NetworkDataSINR$epi,type="SINR",kerneltype="network",
kernelmatrix=NetworkDataSINR$net,datatype="known epidemic",nsim=10000,
sus.par=matrix(c(0.08,0.2,"gamma","gamma",1,1,0.01,0.01,0.5,0.5),ncol=5,nrow=2),
suscov= NetworkDataSINR$cov)
}
}
\keyword{Bayesian Inference}