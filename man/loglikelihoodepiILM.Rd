\name{loglikelihoodepiILM}
\alias{loglikelihoodepiILM}

\title{
Calculates the log likelihood 
}
\description{
Calculates the log likelihood for the specific compartmental framework of the continuous-time ILMs.
}
\usage{

loglikelihoodepiILM(type, epidat, kerneltype, kernelmatrix, distancekernel = NULL,

		suspar = NULL, powersus = NULL, transpar = NULL, powertrans = NULL, 
		
		kernel.par = NULL, spark = NULL, gamma = NULL, delta = NULL,
			
		suscov = NULL, transcov = NULL)
}

\arguments{

  \item{type}{
 Type of the compartmental framework, with the choice of ``SIR'' for
Susceptible-Infectious-Removed and ``SINR'' for Susceptible-Infectious-Notified-Removed.
}
  \item{epidat}{
Epidemic data. Output from \code{\link{datagen}} or input a data frame in the same format.
}
  \item{kerneltype}{
Type of kernel function with choice of ``distance'' for a distance-based ILM, ``network'' for a contact network-based ILM and ``both'' for a combination of network- and distance-based ILM.
}
  \item{kernelmatrix}{
A matrix of distances between individuals when \code{kerneltype} is set to ``distance''. A matrix of edge incidence (0 or 1) or edge weights [0,\eqn{\infty}) representing contact network when \code{kerneltype} is set to ``network''. A list of two matrices (1: distance matrix, 2: contact network matrix) when \code{kerneltype} is set to ``both''.
}
  \item{distancekernel}{
The spatial kernel type when \code{kerneltype} is set to ``distance'' or ``both''. Choices are ``powerlaw'' for a power-law distance kernel or ``Cauchy'' for a Cauchy distance kernel. 
}

  \item{suspar}{
Parameter(s) for the susceptibility function (>0).
}
  \item{powersus}{
The power parameter(s) of the susceptibility function (>0). Default value(s) is 1.
}
  \item{transpar}{
Parameter(s) for the transmissibility function (>0).
}
  \item{powertrans}{
The power parameter(s) of the transmissibility function (>0). Default value(s) is 1.
}
  \item{kernel.par}{
A scalar spatial parameter for the distance-based kernel (>0), or a vector of the spatial and network effect parameters of the network and distance-based kernel (both). It is not required when the \code{kerneltype} is set to ``network''.
}
  \item{spark}{
Spark parameter (>=0), representing random infections that are unexplained by other parts of the model. Default value is zero.
}
  \item{gamma}{
The notification effect parameter for SINR model. The default value is 1.
}
  \item{delta}{
A vector of the shape and rate parameters of the gamma-distributed infectious period (SIR) or a 2 \eqn{\times}{by} 2 matrix of the shape and rate parameters of the gamma-distributed incubation and delay periods (SINR).
}
  \item{suscov}{
Covariate matrix of the susceptibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}
  \item{transcov}{
Covariate matrix of the transmissibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}

}

\details{

We label the \eqn{m} infected individuals \eqn{i = 1, 2, \dots, m} corresponding to their infection (\eqn{I_{i}}{I_i}) and removal (\eqn{R_{i}}{R_i}) times; whereas the \eqn{N-m} individuals who remain uninfected are labeled \eqn{i=m+1, m+2, \dots, N} with \eqn{I_{i}= R_{i} = \infty}{I_i= R_i = \infty}. We then denote infection and removal time vectors for the population as \eqn{\boldsymbol{I} = \{I_{1}, \dots, I_{m}\}}{I = {I_1, ..., I_m}} and \eqn{\boldsymbol{R} = \{R_{1}, \dots, R_{m}\}}{R = {R_1,..., R_m}}, respectively. We assume that infectious periods follow a gamma distribution with shape and rate \eqn{\delta}. The likelihood of the general SIR continuous-time ILMs is then given as follows: 

\if{html}{\figure{likelihood-sir.jpg}{options: width=600 }} 
\if{latex}{
\tabular{lcl}{
  \eqn{L(I,R|\theta)}{} \tab \eqn{=}{} \tab \eqn{\prod_{j=2}^{m}{\left(\epsilon +\sum_{i:I_{i} < I_{j} \leq R_{i}}{\lambda_{ij}^{-}}\right)}}{} \cr
    \tab \eqn{\times}{} \tab \eqn{exp \left\{ -\sum_{i=1}^{m}{\left(\sum_{j=1}^{N}{( (R_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{-}}\right)} \right\}}{} \cr
   \tab \eqn{\times}{} \tab \eqn{exp \left(- \epsilon \sum_{i=1}^{N}{\left[(t_{obs} \wedge I_{i}) - I_{1}\right]} \right) \prod_{i=1}^{m}{f(D_{i};\delta)}}{} \cr
 }
}

where \eqn{\theta} is the vector of unknown parameters; f(.;\eqn{\delta}) indicates the density of the infectious period distribution; and \eqn{D_{i}}{D_i} is the infectious period of infected individual \eqn{i} defined as \eqn{D_{i}= R_{i}-I_{i}}{D_i= R_i-I_i}. The likelihood of the general SINR continuous-time ILMs is given by: 

\if{html}{\figure{likelihood-sinr.jpg}{options: width=600 }} 
\if{latex}{
\tabular{lcl}{
  \eqn{L(I, N, R|\theta)} \tab \eqn{=} \tab \eqn{\prod_{j=2}^{m}{\left(\epsilon +\sum_{i:I_{i} < I_{j} \leq N_{i}}{\lambda_{ij}^{-}}+\sum_{i:N_{i} < I_{j} \leq R_{i}}{\lambda_{ij}^{+}}\right)}} \cr
    \tab \eqn{\times} \tab \eqn{exp \left\{ -\sum_{i=1}^{n_{I}}{\left(\sum_{j=1}^{N}{( (t_{obs} \wedge N_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{-} } \right)} \right\}} \cr
    \tab \eqn{\times} \tab \eqn{exp \left\{- \sum_{i=1}^{n_{N}}{\left(\sum_{j=1}^{N}{( (t_{obs} \wedge R_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) - ( (t_{obs} \wedge N_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{+}}\right)} \right\}} \cr
   \tab \eqn{\times} \tab \eqn{exp \left(- \epsilon \sum_{i=1}^{N}{\left[(t_{obs} \wedge I_{i}) - I_{1}\right]} \right)} \cr
   \tab \eqn{\times} \tab \eqn{\prod_{i=1}^{m}{f(D^{(inc)}_{i};\delta^{(inc)})} \prod_{i=1}^{m}{f(D^{(delay)}_{i};\delta^{(delay)})}} \cr
 }
}	 

where \eqn{D^{inc}_i} and \eqn{D^{delay}_i} are the incubation and delay periods such that \eqn{D^{inc}_i = N_i - I_i} and \eqn{D^{delay}_i = R_i - N_i}, and \deqn{\lambda_{ij}^{-} = \Omega_{S}(j) \Omega_{T}(i) \kappa(i,j),}{\lambda_{ij}^{-} = \Omega_{S}(j) \Omega_{T}(i) k(i,j),} for \eqn{i \in I(t), j \in S(t)}{i in I(t), j in S(t)}, and
\deqn{\lambda_{ij}^{+} = \gamma (\Omega_{S}(j) \Omega_{T}(i) \kappa(i,j)),}{\lambda_{ij}^{+} = \gamma (\Omega_{S}(j) \Omega_{T}(i) k(i,j)),} for \eqn{i \in N(t), j \in S(t)}{i in N(t), j in S(t)}.

Note, \eqn{\lambda_{ij}^{+}} is used only under the SINR model.	
}

\value{
Returns the log likelihood value. 
}
\seealso{
  \code{\link{contactnet}, \link{datagen}, \link{epictmcmc}}.
}
\keyword{Likelihood}