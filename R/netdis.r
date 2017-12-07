######################################################################
# FUNCTION: contactnet
# AUTHORS:
#         Waleed Almutiry <walmutir@uoguelph.ca>, 
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# DESCRIPTION:
#
#     To generate undirected binary contact networks
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License, 
#     version 3,  as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful, 
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License,  version 3,  for more details.
# 
#     A copy of the GNU General Public License,  version 3,  is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/EpiILMCT package
#
######################################################################

contactnet <- function(type, location=NULL, beta, alpha=NULL, plot=FALSE) {

# Specifying the network model:
 
	if (type == "powerlaw") {
		anum  <-  1
        if (is.null(alpha)) {
			alpha  <-  1
        }

    } else if (type == "Cauchy") {
        anum  <-  2
        alpha  <-  1

    } else if (type == "random") {
        anum  <-  3
        alpha  <-  1
    } else {
		stop("Specify type as \"powerlaw\" \"Cauchy\" or \"random\" ",  call. = FALSE)
    }


# To calculate the Euclidean distance:

	if (!is.null(location)) {
		n  <-  length(location[, 1])
		distance  <-  as.matrix(dist(location,  method = "euclidean"))
	} else {
		stop("Error: the individual locations must be specified", call.=FALSE)
	}

	n  <-  length(location[, 1])
	net <- matrix(0, ncol=n, nrow=n)

# To generate the contact network:

if (type == "powerlaw") {
for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		pp <- 1-exp(-alpha*(distance[i, j]^(-beta)))
		u <- runif (1)
		if (u<=pp) {
			net[i, j] <- 1
			net[j, i] <- 1
		} else {
			net[i, j] <- 0
			net[j, i] <- 0
		}
	}	
} 
} else if (type == "Cauchy") {
for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		pp <- 1-exp(-beta/( (distance[i, j]^(2))+(beta^(2)) ) )
		u <- runif (1)
		if (u<=pp) {
			net[i, j] <- 1
			net[j, i] <- 1
		} else {
			net[i, j] <- 0
			net[j, i] <- 0
		}
	}	
} 

} else if (type == "random") {

for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		net[i, j] <- rbinom(1, 1, beta)
		net[j, i] <- net[i, j]
	}
} 

}

for (i in 1:n) {
	net[i, i] <- 0
}
    
# Option for producing fancy plot of the contact network:

    result2  <-  net

	if (plot == TRUE) {

		plot(location, axes=F, ylab="", xlab="", type="n")
		title(main="", sub=paste("(",type, ")"))
		for (j in 1:n) {
			for (i in j:n) {
				if (result2[j, i]==1) {
					lines(c(location[j, 1], location[i, 1]), c(location[j, 2], location[i, 2]))
				}
			}
		}

	mn  <-  max(apply(result2, 1, sum))/5
	points(location, cex=apply(result2, 1, sum)/mn, col="red", pch=20)
	}

    return(result2)
    
}
