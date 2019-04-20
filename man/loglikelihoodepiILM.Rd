\name{loglikelihoodepiILM}
\alias{loglikelihoodepiILM}

\title{
Calculates the log likelihood 
}
\description{
Calculates the log likelihood for the specific compartmental framework of the continuous-time ILMs.
}
\usage{
loglikelihoodepiILM(object, distancekernel = NULL, control.sus = NULL,

control.trans = NULL, kernel.par = NULL, spark = NULL, gamma = NULL,

delta = NULL)
}

\arguments{
  \item{object}{
an object of class ``datagen'' that can be the output of \code{\link{datagen}} or \code{\link{as.epidat}} functions.
}
  \item{distancekernel}{
the spatial kernel type when \code{kerneltype} is set to ``distance'' or ``both''. Choices are ``powerlaw'' for a power-law distance kernel or ``Cauchy'' for a Cauchy distance kernel. 
}
  \item{control.sus}{
a list of values of the susceptibility function (>0):
   \describe{
    \item{1st:}{a vector of values of the susceptibility parameters,}
	\item{2nd:}{an \eqn{n \times n_{s}}{n by n_s} matrix of the susceptibility covariates,}
	\item{3rd:}{a vector of values of the power parameters of the susceptibility function,}
		}
where, \eqn{n} and \eqn{n_s} are the number of individuals and number of susceptibility parameters, respectively. Default = NULL means the model does not include these parameters.
}
  \item{control.trans}{
it has the same structure as the \code{control.sus}, but for the transmissibility function (>0).
}
  \item{kernel.par}{
a scalar spatial parameter for the distance-based kernel (>0), or a vector of the spatial and network effect parameters of the network and distance-based kernel (both). It is not required when the \code{kerneltype} is set to ``network''.
}
  \item{spark}{
spark parameter (>=0), representing random infections that are unexplained by other parts of the model. Default value is zero.
}
  \item{gamma}{
the notification effect parameter for SINR model. The default value is 1.
}
  \item{delta}{
a vector of the shape and rate parameters of the gamma-distributed infectious period (SIR) or a 2 \eqn{\times}{by} 2 matrix of the shape and rate parameters of the gamma-distributed incubation and delay periods (SINR).
}

}

\details{

We label the \eqn{m} infected individuals \eqn{i = 1, 2, \dots, m} corresponding to their infection (\eqn{I_{i}}{I_i}) and removal (\eqn{R_{i}}{R_i}) times; whereas the \eqn{N-m} individuals who remain uninfected are labeled \eqn{i=m+1, m+2, \dots, N} with \eqn{I_{i}= R_{i} = \infty}{I_i= R_i = \infty}. We then denote infection and removal time vectors for the population as \eqn{\boldsymbol{I} = \{I_{1}, \dots, I_{m}\}}{I = {I_1, ..., I_m}} and \eqn{\boldsymbol{R} = \{R_{1}, \dots, R_{m}\}}{R = {R_1,..., R_m}}, respectively. We assume that infectious periods follow a gamma distribution with shape and rate \eqn{\delta}. The likelihood of the general SIR continuous-time ILMs is then given as follows: 

\if{html}{\figure{likelihood-sir.jpg}{options: width=500}} 
\if{latex}{
\tabular{lcl}{
  \eqn{L(I,R|\theta)}{} \tab \eqn{=}{} \tab \eqn{\prod_{j=2}^{m}{\left(\epsilon +\sum_{i:I_{i} < I_{j} \leq R_{i}}{\lambda_{ij}^{-}}\right)}}{} \cr
    \tab \eqn{\times}{} \tab \eqn{exp \left\{ -\sum_{i=1}^{m}{\left(\sum_{j=1}^{N}{( (R_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{-}}\right)} \right\}}{} \cr
   \tab \eqn{\times}{} \tab \eqn{exp \left(- \epsilon \sum_{i=1}^{N}{\left[(t_{obs} \wedge I_{i}) - I_{1}\right]} \right) \prod_{i=1}^{m}{f(D_{i};\delta)}}{} \cr
 }
}
\if{text}{(see the PDF or HTML manuals)}

where \eqn{\theta} is the vector of unknown parameters; f(.;\eqn{\delta}) indicates the density of the infectious period distribution; and \eqn{D_{i}}{D_i} is the infectious period of infected individual \eqn{i} defined as \eqn{D_{i}= R_{i}-I_{i}}{D_i= R_i-I_i}. The likelihood of the general SINR continuous-time ILMs is given by: 

\if{html}{\figure{likelihood-sinr.jpg}{options: width=500}} 
\if{latex}{
\tabular{lcl}{
  \eqn{L(I, N, R|\theta)} \tab \eqn{=} \tab \eqn{\prod_{j=2}^{m}{\left(\epsilon +\sum_{i:I_{i} < I_{j} \leq N_{i}}{\lambda_{ij}^{-}}+\sum_{i:N_{i} < I_{j} \leq R_{i}}{\lambda_{ij}^{+}}\right)}} \cr
    \tab \eqn{\times} \tab \eqn{exp \left\{ -\sum_{i=1}^{n_{I}}{\left(\sum_{j=1}^{N}{( (t_{obs} \wedge N_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{-} } \right)} \right\}} \cr
    \tab \eqn{\times} \tab \eqn{exp \left\{- \sum_{i=1}^{n_{N}}{\left(\sum_{j=1}^{N}{( (t_{obs} \wedge R_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) - ( (t_{obs} \wedge N_{i} \wedge I_{j}) - (I_{i} \wedge I_{j})) \lambda_{ij}^{+}}\right)} \right\}} \cr
   \tab \eqn{\times} \tab \eqn{exp \left(- \epsilon \sum_{i=1}^{N}{\left[(t_{obs} \wedge I_{i}) - I_{1}\right]} \right)} \cr
   \tab \eqn{\times} \tab \eqn{\prod_{i=1}^{m}{f(D^{(inc)}_{i};\delta^{(inc)})} \prod_{i=1}^{m}{f(D^{(delay)}_{i};\delta^{(delay)})}} \cr
 }
}
\if{text}{(see the PDF or HTML manuals)}

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