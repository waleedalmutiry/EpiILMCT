# EpiILMCT: Continuous Time Distance-Based and Network-Based Individual Level Models for Epidemics
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Downloads](https://cranlogs.r-pkg.org/badges/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Rdoc](http://www.rdocumentation.org/badges/version/EpiILMCT)](http://www.rdocumentation.org/packages/EpiILMCT)

The R package *EpiILMCT* provides tools for simulating from continuous-time individual level models of disease transmission, and carrying out infectious disease data analyses with the same models. The epidemic models considered are distance-based and contact network-based models within Susceptible-Infectious-Removed (SIR) or Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks.

## Features
### Simulation
#### Contact network
A function (contactnet) is included to generate undirected unweighted contact networks. It can simulate both spatial networks where connections are more likely to occur between individuals closer in space ("spatial contact networks"), as well as random contact networks. The function contactnet has three available options ("powerlaw", "Cauchy", and "random") for the network model, where the first two options simulate spatial contact networks in which the probability of connections between individuals are based on required XY coordinate input.
```s
library(EpiILMCT)
set.seed(22)

# to generate the XY coordinates of 50 individuals:

loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)

# Spatial contact network:
# power-law model:
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, 
	nu = 0.5)
plot(net1)

# Cauchy model:
net2<- contactnet(type = "Cauchy", location = loc, beta = 0.5)
plot(net2)

# random contact network:
net3<- contactnet(type = "random", num.id = 50, beta = 0.08)
plot(net3)  # the default options in igraph package.
plot(net3, vertex.color = "red", vertex.size = 15, edge.color = "black",
vertex.label.cex = 1, vertex.label.color = "black") 

```
#### Epidemic data:
```s
library(EpiILMCT)

set.seed(22)
loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, 
	nu = 0.5)

#################################
#################################
#################################
########       SIR       ########
#################################
#################################
#################################

#################################
########  spatial ILMs:  ########
#################################

# first: with power-law kernel:

cov1 <- cbind(runif(50, 0, 50), runif(50, 0, 5))
cov2 <- cbind(runif(50, 0, 50), runif(50, 0, 5))

# with specifying all the options of susceptibility and transmissibility functions:

epi.dist.po1 <- datagen(type = "SIR", kerneltype = "distance", kernelmatrix = loc,
				initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
				distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				kernel.par = 8, delta = c(6,2), suscov = cov1, transcov = cov2,
				seedval = 299)

epi.dist.po1$epidat

loglikelihoodepiILM(object = epi.dist.po1, distancekernel = "powerlaw", 
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)), 
kernel.par = 8, delta = c(6,2))

epi.dist.po2 <- datagen(type = "SIR", kerneltype = "network", kernelmatrix = net1$contact.network,
				initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
				suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				delta = c(6,2), suscov = cov1, transcov = cov2,
				seedval = 299)

epi.dist.po2$epidat

loglikelihoodepiILM(object = epi.dist.po2,
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)), 
delta = c(6,2))

epi.dist.po3 <- datagen(type = "SIR", kerneltype = "both", kernelmatrix = list(loc,net1$contact.network),
				initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
				distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				kernel.par = c(8,0.3), delta = c(6,2), suscov = cov1, transcov = cov2,
				seedval = 299)

epi.dist.po3$epidat

loglikelihoodepiILM(object = epi.dist.po3, distancekernel = "powerlaw", 
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)), 
kernel.par = c(8,0.3), delta = c(6,2))

data(NetworkDataSINR)
names(NetworkDataSINR)


netSIR<-as.epidat(type = "SIR", kerneltype = "network", inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "powerlaw")

netSIR<-as.epidat(type = "SIR", kerneltype = "distance", inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "Cauchy")

#################################
#################################
#################################
########       SINR      ########
#################################
#################################
#################################

#################################
########  spatial ILMs:  ########
#################################

# first: with power-law kernel:

cov1 <- cbind(runif(50, 0, 50), runif(50, 0, 5))
cov2 <- cbind(runif(50, 0, 50), runif(50, 0, 5))

# with specifying all the options of susceptibility and transmissibility functions:

epi.dist.po1 <- datagen(type = "SINR", kerneltype = "distance", kernelmatrix = loc,
				initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 4,
				distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				kernel.par = 8, delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE), 
				suscov = cov1, transcov = cov2,
				seedval = 299)

epi.dist.po1$epidat

loglikelihoodepiILM(object = epi.dist.po1, distancekernel = "powerlaw", 
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)), 
kernel.par = 8, delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE))


epi.dist.po2 <- datagen(type = "SINR", kerneltype = "network", kernelmatrix = net1,
				initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 2,
				suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE), 
				suscov = cov1, transcov = cov2,
				seedval = 299)
epi.dist.po2$epidat

loglikelihoodepiILM(object = epi.dist.po2, distancekernel = "powerlaw", 
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)), 
delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE))

epi.dist.po3 <- datagen(type = "SINR", kerneltype = "both", kernelmatrix = list(loc,net1$contact.network),
				initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 2,
				distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
				powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
				kernel.par = c(8,0.3), delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE), 
				suscov = cov1, transcov = cov2,
				seedval = 299)

epi.dist.po3$epidat

data(NetworkDataSINR)
names(NetworkDataSINR)


netSINR<-as.epidat(type = "SINR", kerneltype = "network", incub.time = NetworkDataSINR$epi[,4], inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "powerlaw")


netSINR<-as.epidat(type = "SINR", kerneltype = "distance", incub.time = NetworkDataSINR$epi[,4], inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "Cauchy")
```
### Analyzing


### Plotting


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing
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

## Deployment


## Built With

* [R](https://cran.r-project.org) - The Comprehensive R Archive Network

## Contributing


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/waleedalmutiry/EpiILMCT/tags). 

## Authors

* **Waleed Almutiry** - *Maintainer*
* **Rob Deardon** - *Maintainer*

## License

This project is licensed under the GNU General Public License,  version 3 - see the (http://www.r-project.org/Licenses/GPL-3) file for details.

## Acknowledgments
This project was funded by the Ontario Ministry of Agriculture, Food and Rural Affairs (OMAFRA), the Natural Sciences and Engineering Research Council of Canada (NSERC), Qassim University through the Saudi Arabian Cultural Bureau in Canada, and the University of Calgary Eyes High Post Doctoral Scholarship scheme.
