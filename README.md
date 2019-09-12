<p align="center">
	<img src="https://user-images.githubusercontent.com/18523406/64757579-2e6a0e80-d53b-11e9-94ec-473247d05c8a.jpg" width="35%"></img> 
</p>

# EpiILMCT: Continuous Time Distance-Based and Network-Based Individual Level Models for Epidemics

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Downloads](https://cranlogs.r-pkg.org/badges/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Rdoc](http://www.rdocumentation.org/badges/version/EpiILMCT)](http://www.rdocumentation.org/packages/EpiILMCT)

The R package *EpiILMCT* provides tools for simulating from continuous-time individual level models of disease transmission, and carrying out infectious disease data analyses with the same models. The epidemic models considered are distance-based and contact network-based models within Susceptible-Infectious-Removed (SIR) or Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks.

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

We dene an object of class \datagen" to take a list of values needed for the use of other
functions, such as, plot.datagen and epictmcmc. This list contains: type, kerneltype,
epidat (event times), location (XY coordinates of individuals), and network (contact
network matrix). In the case of setting the kerneltype to \distance", a zero contact
network matrix will be created for the network option. The package has also a separate
function as.epidat that generates an object of class \datagen" for a given epidemic data
set (See Appendix B that has a brief example of using this function).
The package also contains an S3 method plot.datagen function, which illustrates disease
spread through the epidemic timeline. This function can be used for either distancebased
or network-based ILMs. The object of this function has to be of class \datagen".
If the plottype argument is set to "history", the function produces epidemic curves of
infection and removal times. Example plots are shown in Figure 3. Conversely, setting
this argument to "propagation" produces plots of the epidemic propagation over time.
With the latter option, exactly which plots are output varies by kernel. With the network
kernel, the function plots all the connections between individuals and overlays these
with the epidemic pathway direction over time. This path direction consists of directed
edges from all infectious individuals connected to a given newly infected individual i with
infection time Ii (one per plot). Thus, this produces directed networks showing possible
pathways of the disease propagation. Figure 4 shows a propagation plot for the rst six
infection events of the simulated network ILM generated epidemic. With the distance
kernel, the function plots the spatial epidemic dispersion over time. It shows the changes
in the individual status that related to the chosen compartmental framework. To avoid
displaying too many plots, the time.index argument allows user to obtain propagation
plots at specific infection time points rather than at every infection time.

### Analyzing

<p align="center">
<img src="https://user-images.githubusercontent.com/18523406/64754037-1a6be000-d52e-11e9-80c0-1864b828591f.jpg"></img> 
</p>


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
