<p align="center">
	<img src="https://user-images.githubusercontent.com/18523406/64757579-2e6a0e80-d53b-11e9-94ec-473247d05c8a.jpg" width="35%"></img> 
</p>

# EpiILMCT: Continuous Time Distance-Based and Network-Based Individual Level Models for Epidemics

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Downloads](https://cranlogs.r-pkg.org/badges/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Rdoc](http://www.rdocumentation.org/badges/version/EpiILMCT)](http://www.rdocumentation.org/packages/EpiILMCT)

The **R** package **EpiILMCT** provides tools for simulating from continuous-time individual level models of disease transmission, and carrying out infectious disease data analyses with the same models. The epidemic models considered are distance-based and contact network-based models within Susceptible-Infectious-Removed (SIR) or Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks.

## Installing
You can install the **EpiILMCT** version from
[CRAN](https://cran.r-project.org/package=EpiILMCT).

```s
install.packages('EpiILMCT', dependencies = TRUE)
```

You can install the **development** version from
[Github](https://github.com/waleedalmutiry/EpiILMCT-package)

```s
# install.packages("devtools")
devtools::install_github("waleedalmutiry/EpiILMCT-package")
```


## Features
### Simulation
#### Contact network
Different types of undirected unweighted contact networks can be generated through the function **_contactnet_**. This function function has a 'type' option to specify the type of network to be used. It has three available options ("powerlaw", "Cauchy", and "random") for the network model, where the first two options simulate spatial contact networks in which the probability of connections between individuals are based on required XY coordinate input. The output of this function is stored as an object of class 'contactnet'. Also, an S3 method plot function is introduced which uses this object as its input to provide a fancy plot of the contact network.  

```s
library(EpiILMCT)
set.seed(22)

# to generate the XY coordinates of 50 individuals:

loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)
```
##### power-law spatial contact network
```s
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, 
	nu = 0.5)
plot(net1)
```
##### Cauchy spatial contact network
```s
net2<- contactnet(type = "Cauchy", location = loc, beta = 0.5)
plot(net2)
```
##### random contact network
```s
net3<- contactnet(type = "random", num.id = 50, beta = 0.08)
plot(net3)  # the default options in igraph package.
plot(net3, vertex.color = "red", vertex.size = 15, edge.color = "black",
vertex.label.cex = 1, vertex.label.color = "black") 
```
<p align="center">
<img src="https://user-images.githubusercontent.com/18523406/64754143-ae3dac00-d52e-11e9-81e6-5735046fd877.jpg"></img> 
</p>

#### Epidemic data:

The function **_datagen_** allows the user to generate epidemics from the continuous time ILMs under the SIR or SINR compartmental frameworks. This function can generate epidemics based on different kernels through the _kerneltype_ argument which takes one of three options: _"distance"_ for **distance-based**, _"network"_ for **network-based**, or _"both"_ for **distance and network-based**. The appropriate kernel matrix must also be provided via the **_kernelmatrix_** argument. If _"distance"_ is chosen as the **_kerneltype_**, the user must choose a spatial kernel (_"powerlaw"_ or _"Cauchy"_) through
the **_distancekernel_** argument. Here is an example of generating an epidemic from the SIR continuous time network-based ILM.

```s
library("EpiILMCT")
set.seed(91938)

# To simulate the XY coordinate of 50 individuals and their corresponding binary covariate values:
loc <- matrix(cbind(runif(50, 0, 10), runif(50, 0, 10)), ncol = 2, nrow = 50)
cov <- cbind(rep(1, 50), rbinom(50, 1, 0.5))

# To simulate the contact network:
net <- contactnet(type = "powerlaw", location = loc, beta = 1.8, nu = 1)

# To simulate the epidemic:
epi <- datagen(type = "SIR", kerneltype = "network", kernelmatrix = net, suspar = c(0.08, 0.5), delta = c(4, 2), 
   suscov = cov, seedval =  498643)
epi
```
The output of the **_datagen_** function is stored as an _datagen_ object which takes a list of:
1. type
2. kerneltype
3. epidat (event times)
4. location (XY coordinates of individuals)
5. network (contact network matrix), in the case of setting the **_kerneltype_** to _distance_, a zero contact network matrix will be created for the network option. 

The package also contains an S3 method **_plot.datagen_** function, which illustrates disease spread through the epidemic timeline. This function can be used for either **_distance-based_** or **_network-based_** ILMs. The object of this function has to be of class _datagen_. The plot S3 function has a **_plottype_** argument that can be set to _"history"_, to produce epidemic curves of infection and removal times, or set to _"propagation"_ to produce plots of the epidemic propagation over time. The following commands are to produce the below two graphs for the generated epidemic above.

```s
plot(epi, plottype = "propagation", time.index = seq_len(6))
plot(epi, plottype = "history")
```
<p align="center">
<img src="https://user-images.githubusercontent.com/18523406/64774791-8c0f5280-d55d-11e9-8375-6e253fb303ac.jpg"></img> 
</p>

<p align="center">
<img src="https://user-images.githubusercontent.com/18523406/64774792-8c0f5280-d55d-11e9-9f51-ddea23907fb4.jpg"></img> 
</p>



### Analyzing

Metropolis-Hastings MCMC is performed to estimate the joint posterior of the model parameters and latent variables (the latter if various event times are assumed unknown). This is achieved using the function **_epictmcmc_**. The below figure illustrates the structure of most arguments of this function. 

<p align="center">
<img src="https://user-images.githubusercontent.com/18523406/64754037-1a6be000-d52e-11e9-80c0-1864b828591f.jpg"></img> 
</p>


The output of this function is an object of class _epictmcmc_. There are S3 methods: **_print.epictmcmc_**, **_summary.epictmcmc_** and **_plot.epictmcmc_** that depend on the **coda** package. The latter function produced the trace plots of the posterior distributions of the model parameters with the same options of the **_plot.mcmc_** function in the **coda** package, such as, _start_, _thin_, and _density_. 

In the case of datatype is set to either _"known removal"_ (unknown infection times only) or _"unknown removal"_ (unknown infection and removal times), plots of the average posterior and 95% CI of the unobserved event times can also be produced.

The class _epictmcmc_ contains the MCMC samples of the model parameters and the missing information (in case datatype is not set to _known epidemic_), and other useful information to be used in other functions such as the above S3 methods. For example, when **_datatype_** = _known epidemic_, the class _epictmcmc_ has a list contained the following:
1. "compart.framework"
2. "kernel.type"
3. "data.assumption"
4. "parameter.samples"
5. "log.likelihood"
6. "acceptance.rate"
7. "number.iteration"
8. "number.parameter"
9. "number.chains"

The following commands are to perform the MCMC for analyzing the above epidemic data set using the **_epictmcmc_** function.

```s
suscov <- list(NULL)
suscov[[1]] <- list(c(0.01, 0.1), matrix(c("gamma", "gamma", 1, 1, 0.1, 0.1, 0.5, 1), ncol = 4, nrow = 2))
suscov[[2]] <- epi[[2]]
mcmc1 <- epictmcmc(object = epi, datatype = "known epidemic", nsim = 150000, control.sus = suscov, 
seedval = 524837))
```
The estimates of the model parameters can be then obtained either through using S3 **_summary_** function of the **_epictmcmc_** for the _mcmc1_, or using the **_summary_** function of **coda** package for _mcmc1$parameter.samples_, for example. The posterior means and 95% credible intervals of the model parameters using the former **_summary_** can be obtained via the following command:

```s
summary(mcmc1, start = 10000)
********************************************************* 
Model: SIR network-based continuous-time ILM 
Method: Markov chain Monte Carlo (MCMC) 
Data assumption: fully observed epidemic 
number.chains : 1 chains 
number.iteration : 140000 iterations 
number.parameter : 2 parameters 
********************************************************* 
 1. Empirical mean and standard deviation for each variable,
plus standard error of the mean:
                Mean        SD    Naive SE Time-series SE
Alpha_s[1] 0.0973068 0.0325623 8.70262e-05    0.000320181
Alpha_s[2] 0.5335298 0.1295704 3.46290e-04    0.000986805
 2. Quantiles for each variable:
                2.5%       25%       50%      75%    97.5%
Alpha_s[1] 0.0445579 0.0735079 0.0936119 0.117225 0.170676
Alpha_s[2] 0.3020837 0.4428192 0.5257345 0.617407 0.806669
 3. Empirical mean, standard deviation, and quantiles for the log likelihood,
          Mean             SD       Naive SE Time-series SE 
  -48.45599501     1.02806188     0.00274760     0.00963328 
    2.5%      25%      50%      75%    97.5% 
-51.2095 -48.8603 -48.1447 -47.7187 -47.4439 
 4. acceptance.rate : 
Alpha_s[1] Alpha_s[2] 
  0.131621   0.219375 
```

## Built With

* [R](https://cran.r-project.org) - The Comprehensive R Archive Network

## Authors

* **Waleed Almutiry** - *Author and Maintainer*
* **Rob Deardon** - *Author*
* **Vineetha Warriyar K. V.** - *Contributer*

## License

This project is licensed under the GNU General Public License,  version 3 - see the (http://www.r-project.org/Licenses/GPL-3) file for details.
