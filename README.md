# EpiILMCT: Continuous Time Distance-Based and Network-Based Individual Level Models for Epidemics
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Downloads](https://cranlogs.r-pkg.org/badges/EpiILMCT)](https://cran.r-project.org/package=EpiILMCT)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Rdoc](http://www.rdocumentation.org/badges/version/EpiILMCT)](http://www.rdocumentation.org/packages/EpiILMCT)

The R package *EpiILMCT* provides tools for simulating from continuous-time individual level models of disease transmission, and carrying out infectious disease data analyses with the same models. The epidemic models considered are distance-based and contact network-based models within Susceptible-Infectious-Removed (SIR) or Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks.

## Features
### Simulation
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
