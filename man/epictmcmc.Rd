\name{epictmcmc}
\alias{epictmcmc}

\title{
Markov Chain Monte Carlo-based tool for analyzing (SIR/SINR) distance-based/network-based individual-level disease data.
}
\description{
Runs a Bayesian data augmented MCMC algorithm for fitting specified models (SIR or SINR) with three possible assumptions regarding the data set: 1. a completely observed epidemic; 2. unknown infection times but known removal times; or 3. unknown removal and infection times but known notification times (SINR only).
}
\usage{

epictmcmc(object, distancekernel = NULL, datatype, blockupdate = NULL, nsim,

nchains = NULL, control.sus = NULL, control.trans = NULL, kernel.par = NULL,

spark.par = NULL, delta = NULL, gamma.par = NULL, periodproposal = NULL,

parallel = FALSE)

}

\arguments{
  \item{object}{
an object of class ``datagen'' that can be the output of \code{\link{datagen}} or \code{\link{as.epidat}} functions.
}
  \item{distancekernel}{
the spatial kernel type when \code{kerneltype} is set to ``distance'' or ``both''. Choices are ``powerlaw'' for a power-law distance kernel or ``Cauchy'' for a Cauchy distance kernel.
}
  \item{datatype}{
the ILM data set is analyzed based on three possible assumptions. The possible assumptions are ``known epidemic'' for a completely observed epidemic; ``known removal'' times for observing only the removal times for the SIR model, or removal and notified times for the SINR model; or ``unknown removal'' times for observing only the notified times which is applicable only when fitting the SINR model.
}
  \item{blockupdate}{
a vector of the number of initially observed infection times and the size of blocks for updating the infection times/infectious periods (SIR model), or infection times/incubation periods and removal times/delay periods (SINR model). The default is observing only the first infection time and performing the MCMC updating infection time parameters via single parameter updates.
}
  \item{nsim}{
the number of MCMC iterations to be carried out for a single chain.
}
  \item{nchains}{
the number of MCMC chains to be run. Default is 1 (a single MCMC chain).
}
  \item{control.sus}{
a list of arguments for updating the parameters of the susceptibility function (>0):
   \describe{
    \item{1st:}{the susceptibility parameters (list):
	   \itemize{
		    \item a vector (of length \eqn{n_{s}}{n_s}) or a \eqn{n_{s} \times \code{nchains}}{n_s by nchains} matrix of initial values of the susceptibility parameters,
		    \item a \eqn{n_{s} \times 4}{n_s by 4} matrix of the required arguments for updating the susceptibility parameters,
		  }
		}
	\item{2nd:}{an \eqn{n \times n_{s}}{n by n_s} matrix of the susceptibility covariates.}
	\item{3rd:}{the power parameters of the susceptibility function (list):
	   \itemize{
		    \item a vector (of length \eqn{n_{s}}{n_s}) or a \eqn{n_{s} \times \code{nchains}}{n_s by nchains} matrix of initial values of the power parameters,
		    \item a \eqn{n_{s} \times 4}{n_s by 4} matrix of the required arguments for updating the power parameters,
		  }
		}
  	}

where \eqn{n} and \eqn{n_s} are the number of individuals and number of susceptibility parameters, respectively. Default = NULL means the model does not include these parameters. See details for a description of required arguments.
}

  \item{control.trans}{
a list of arguments with the same structure as \code{control.sus}, but for the transmissibility function (>0).
}
  \item{kernel.par}{
a list of arguments for updating the spatial parameter of the distance-based kernel. The first argument is a vector of initial values of the spatial parameter when \code{kerneltype} is set to ``distance''  (>0) or a \eqn{2 \times \code{nchains}}{2 by nchains} matrix of initial values of the spatial and network effect parameters when \code{kerneltype} is set to ``both'' (>0). The second argument is a vector of arguments for updating the spatial parameter of the distance-based kernel when \code{kerneltype} is set to ``distance'' or a \eqn{2 \times 4}{2 by 4} matrix of arguments for updating the spatial and network effect parameters when \code{kerneltype} is set to ``both''. Default = NULL means the model does not include these parameters (i.e., we are fitting a network-based ILM). See details for a description of required arguments.
}
  \item{spark.par}{
a list of arguments for updating the spark parameter (>=0). The first argument is a vector of the initial values. The second argument is a vector of the required arguments for updating the spark parameter. It should be defined if \code{kerneltype} is set to ``network'' with \code{datatype}= ``known removal'' or ``unknown removal''. Default = NULL means the model does not include this parameter. See details for a description of required arguments.
}
  \item{delta}{
a list of arguments for updating the infectious period rate for an SIR model, or for updating the incubation and delay periods rates for an SINR model (>0). The arguments for updating each rate are a vector of the fixed and known shape parameter(s) for the distribution of the period type, a vector (matrix) of initial values of the rate(s) of the infectious, incubation and/or delay period distributions, and a vector (matrix) of the shape and rate parameters of the gamma prior distribution for the rate parameter(s). Default = NULL means these parameters do not need to be updated, which is the case when \code{datatype} is set to ``known epidemic''.
}
  \item{gamma.par}{
a list of arguments for updating the notification effect parameter (SINR) (>=0). The first argument is a vector of initial values. The second argument is a vector of the required arguments for updating the notification effect parameter. Default = NULL means the parameter is not updated and is assigned a value equal to 1. See details for a description of required arguments.
}
  \item{periodproposal}{
a vector/matrix of the proposal distribution parameters of the independence sampler for updating the infectious period (SIR model), or the incubation and/or delay periods (SINR model). Here, we use a gamma proposal distribution. It is required when \code{datatype} is set to ``known removal'' or ``unknown removal''. Default = NULL means the parameters of the gamma proposal distribution are the known shape and updated rate values based on the argument \code{delta}.
}
  \item{parallel}{
if set to ``TRUE'', multiple MCMC chains are running in parallel. Default is ``FALSE'' means chains (>1) are run sequentially.
}
}


\details{

Bayesian MCMC is performed to obtain posterior estimates of the model parameters and latent variables. When the \code{datatype} is set to ``known removal'' or ``unknown removal'', we assume the infectious periods (SIR model) or the incubation and/or delay periods (SINR model) follow gamma distributions with fixed shape and unknown rate parameters. We assign a gamma prior distribution for the rate parameter with shape \eqn{a} and rate \eqn{b}. This leads the rate parameter(s) to have a gamma conditional distribution.

Under the SIR model, the conditional distribution is:

\deqn{ \delta|\theta,\boldsymbol{I},\boldsymbol{R} \sim \Gamma(m+a_{\delta},M+b_{\delta}), }{\delta|\theta,I,R ~ \Gamma(m + a_\delta, M + b_\delta),}
where \eqn{\delta}{\delta} is the rate of the infectious period distribution; \eqn{m} is the number of infected individuals; \eqn{M=\sum_{i=1}^{m}{(R_{i}-I_{i})}}{M=\sum_{i=1}^{m}(R_i-I_i)}; and \eqn{a_{\delta}}{a_\delta} and \eqn{b_{\delta}}{b_\delta} are the prior parameters for the infectious period rate, respectively.

Under the SINR model, the conditional distribution is:
\deqn{ \delta^{(inc)}|\theta,\boldsymbol{I},\boldsymbol{N},\boldsymbol{R} \sim \Gamma(m+a_{\delta^{(inc)}},M_{inc}+b_{\delta^{(inc)}}),}{\delta^(inc)|\theta, I, N, R ~ \Gamma(m + a_{\delta^(inc)},M_inc +
	b_{\delta^(inc)}),}
where \eqn{\delta^{(inc)}}{\delta^(inc)} is the rate of the incubation period distribution; \eqn{M_{inc}=\sum_{i=1}^{m}{(N_{i}-I_{i})}}{M_inc=\sum_{i=1}^{m}(N_i-I_i)}; and \eqn{a_{\delta^{(inc)}}}{a_{\delta^(inc)}} and \eqn{b_{\delta^{(inc)}}}{b_{\delta^(inc)}} are the prior parameters for the incubation period rate, respectively; and

\deqn{ \delta^{(delay)}|\theta,\boldsymbol{I},\boldsymbol{N},\boldsymbol{R} \sim \Gamma(m+a_{\delta^{(delay)}},M_{delay}+b_{\delta^{(delay)}}),}{\delta^(delay)|\theta, I, N, R ~ \Gamma(m+a_{\delta^(delay)},M_delay+
	b_{\delta^(delay)}),}
where \eqn{\delta^{(delay)}}{\delta^(delay)} is the rate parameter of the delay period distribution; \eqn{M_{delay}=\sum_{i=1}^{m}{(R_{i}-N_{i})}}{M_delay=\sum_{i=1}^{m}(R_i-N_i)}; and \eqn{a_{\delta^{(delay)}}}{a_{\delta^(delay)}} and \eqn{b_{\delta^{(delay)}}}{b_{\delta^(delay)}} are the prior parameters for the delay period rate, respectively.

A Gibbs update (i.e., sampling from the conditional posterior distribution) is then used for the infectious period (SIR model) or the incubation and/or delay period rates (SINR model). The required information for each period distribution are entered via the \code{delta} argument. We assume each period type follows a gamma distribution with fixed shape and unknown rate. Thus, to update the rate parameter of each period, \code{delta} is a list containing a vector of the fixed shape value(s), a vector (matrix) of the initial values of the rate(s), and a vector (matrix) of the parameters of the prior distribution of the rate parameter(s). In the case of incubation and delay periods being estimated, the input of the initial values is a \eqn{2 \times \code{nchains}}{2 by nchains} matrix, and the prior parameters are given by a \eqn{2 \times 2}{2 by 2} matrix in which each row contains the required information for each period rate.

An independence sampler is used for updating the infection times and infectious periods (SIR model), or the infection and/or removal times, and the incubation and/or delay periods (SINR model). Under the SIR model, we update each infection time \eqn{I_{i}}{I_i} by proposing infectious period \eqn{D^{*}_{i}}{D^{*}_i} from a proposal distribution with tuning parameter, such that \eqn{D_{i}\sim}{D_i ~} f(.). Then, the new infection time is just the difference between the fixed known removal time and the new infectious period of the \eqn{i^{th}}{ith} individual. Each infection time/infectious period is updated in this way in turn. The same procedure is applied for updating the missing event times and incubation and delay periods for the SINR model, with their corresponding parameters.

The \code{control.sus} and \code{control.trans} arguments provide all the information needed for the susceptibility and transmissibility functions. These arguments must be defined as a list of three objects. The first object is a list of: 1) a vector or matrix of initial values; and 2) a vector or matrix of the following arguments in order: initial value, prior distribution ("gamma", (positive) "half normal" or "uniform"), a vector of the prior distribution parameters, and the proposal variance. The second object is a matrix of the susceptibility or transmissibility risk factors. The third contains the power parameters of the susceptibility or transmissibility functions and has the same structure as the first object. The other model parameters \code{kernel.par}, \code{spark.par} and \code{gamma.par} also have the same structure as the first object. These model parameters, along with the model parameters of the susceptibility and transmissibility functions (\code{control.sus} and \code{control.trans}), are updated in turn using a random-walk Metropolis Hastings algorithm with normal proposals. The normal proposals must be tuned, via the proposal variance, by the user to achieve good mixing properties.

Note that, setting the variance of the normal proposal distribution to zero fixes a parameter at its initial value. This option allows the user to fix such a parameter in the MCMC procedure while updating others.

For faster run of the \code{epictmcmc} function, an option of running multiple chains on parallel can be performed, controlled via the two options \code{nchains} and \code{parallel}. If \code{parallel} is set to TRUE, the number of required cores is set to the minimum of \code{nchains} and the number of available cores on the user's computer.
}
\value{
Returns an object of class \code{epictmcmc} that contains:

   \describe{
    \item{compart.framework:}{the compartmental framework model used in the analysis.}
    \item{kernel.type:}{the used \code{kernel.type} in the function.}
    \item{data.assumption:}{the \code{datatype}.}
    \item{parameter.samples:}{the MCMC output of the updated model parameters.}
    \item{log.posterior:}{the log posterior densities.}
    \item{acceptance.rate:}{the acceptance rates for MCMC updates of each parameter.}
    \item{number.iteration:}{the number of iterations.}
    \item{number.parameter:}{the number of the unknown model parameters.}
    \item{number.chains:}{the number of the MCMC chains.}
    \item{infection.times.samples:}{the updated infection times when \code{datatype} is set to ``known removal''.}
    \item{Average.infectious.periods:}{the average infectious period when \code{type} is set to ``SIR'' and \code{datatype} is set to ``known removal''.}
    \item{removal.times.samples:}{the updated removal times when \code{datatype} is set to ``unknown removal''.}
    \item{Average.incubation.periods:}{the average incubation period when \code{type} is set to ``SINR'' and \code{datatype} is set to either ``known removal'' or ``unknown removal''.}
    \item{Average.delay.periods:}{the average delay period when \code{type} is set to ``SINR'' and \code{datatype} is set to ``unknown removal''.}

  }
}

\seealso{
  \code{\link{print.epictmcmc}}, \code{\link{summary.epictmcmc}}, \code{\link{plot.epictmcmc}}, \code{\link{datagen}}, \code{\link{loglikelihoodepiILM}}.
}

\examples{
## This is for testing; the number of MCMC iterations needs to be
## mucher higher in practice.

## distance-based SIR continuous-time ILMs ##
data(SpatialData)
## performing the MCMC-tool for analyzing the fully observed spatial data
## under the SIR distance-based continuous ILM:
suspar <- list(NULL)
suspar[[1]] <- list(2, c("gamma", 1, 0.01, 0.5))
suspar[[2]] <- rep(1, length(SpatialData$epidat[,1]))
kernel1 <- list(2, c("gamma", 1, 0.01, 0.5))
mcmcres2 <- epictmcmc(object = SpatialData,
distancekernel = "powerlaw", datatype = "known epidemic", nsim = 50,
control.sus = suspar, kernel.par = kernel1)
summary(mcmcres2)

## performing the MCMC-tool for analyzing the partially observed spatial
## data (unknown infection times) under the SIR distance-based
## continuous ILM:
suspar <- list(NULL)
suspar[[1]]<-list(2,c("gamma", 1, 0.01, 0.8))
suspar[[2]]<- rep(1,length(SpatialData$epidat[,1]))
kernel1 <- list(2, c("gamma", 1, 0.01, 0.5))
mcmcres22 <- epictmcmc(object = SpatialData,
distancekernel = "powerlaw", datatype = "known removal", nsim = 50,
control.sus = suspar, kernel.par = kernel1, delta = list(1,2,c(4,2)))
summary(mcmcres22)

## distance-based and network-based SIR ILMs ##
data(SpatialNetData)
## performing the MCMC-tool for analyzing the fully observed spatial and
## network data under the SIR distance-based and network-based
## continuous-time ILM:
suspar <- list(NULL)
suspar[[1]] <- list(NULL)
suspar[[1]][[1]] <- c(0.08, 0.2)
suspar[[1]][[2]] <- matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.1, 0.5),
ncol = 4, nrow = 2)
suspar[[2]] <- SpatialNetData[[2]]
kernel1 <- list(NULL)
kernel1[[1]] <- c(5, 0.5)
kernel1[[2]] <- matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.5, 1),
ncol = 4, nrow = 2)
mcmcres3 <- epictmcmc(object = SpatialNetData[[1]],
distancekernel = "powerlaw", datatype = "known epidemic", nsim = 50,
control.sus = suspar, kernel.par = kernel1)
summary(mcmcres3)

## network-based SIR ILMs ##
data(NetworkData)
## performing the MCMC-tool for analyzing the fully observed network data
## under the SIR network-based continuous-time ILM:
suspar <- list(NULL)
suspar[[1]] <- list(NULL)
suspar[[1]][[1]] <- c(0.08, 0.5)
suspar[[1]][[2]] <- matrix(c("gamma", "gamma", 1, 1, 1, 1, 0.1, 0.5),
ncol = 4, nrow = 2)
suspar[[2]] <- NetworkData[[2]]
mcmcres4 <- epictmcmc(object = NetworkData[[1]], datatype = "known epidemic",
nsim = 50, control.sus = suspar)
summary(mcmcres4)
}
\keyword{Bayesian Inference}
