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

contactnet <- function(type, num.id = NULL, location = NULL, beta, nu = NULL) {

    # Specifying the network model:

    if (type == "powerlaw") {
        anum  <-  1
        if (is.null(nu)) {
            nu  <-  1
        }

    } else if (type == "Cauchy") {
        anum  <-  2
        nu  <-  1

    } else if (type == "random") {
        anum  <-  3
        nu  <-  1
    } else {
        stop("Specify type as \"powerlaw\" \"Cauchy\" or \"random\" ",  call. = FALSE)
    }


    # To calculate the Euclidean distance:

    if (!is.null(location)) {
        if (is.null(num.id)){
            n  <-  length(location[, 1])
            distance  <-  as.matrix(dist(location,  method = "euclidean"))
        } else if (!is.null(num.id)) {
            if (num.id != length(location[,1])) {
                stop("The number of individuals in the option \"num.id\" is not equal to the number of XY coordinates of individuals", call. = FALSE)
            }
            n <- num.id
        }
    } else if (all(is.null(location) & type == "random") == TRUE) {
        if (is.null(num.id)) {
            stop("The option \"num.id\" has to be specified", call. = FALSE)
        }
        n <- num.id
        location <- matrix(0, ncol = 2, nrow = n)
    } else {
        stop("Error: the individual locations must be specified", call.=FALSE)
    }

    #    n  <-  length(location[, 1])
    net <- matrix(0, ncol=n, nrow=n)

    # To generate the contact network:

    if (type == "powerlaw") {
        for (i in 1:(n-1)) {
            for (j in (i+1):n) {
                pp <- 1-exp(-nu*(distance[i, j]^(-beta)))
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

    outnet <- list(location = location, contact.network = net, type = type)

    # Naming the class of the object:

    class(outnet) <- "contactnet"

    outnet

}



# to plot the contact network:

plot.contactnet <- function(x, ...) {

    # Option for producing fancy plot of the contact network:
    if (is(x, "contactnet")) {

        if (x$type == "powerlaw" | x$type == "Cauchy"){
            n <- dim(x$contact.network)[1]
            plot(x$location, axes=FALSE, ylab="", xlab="", type="n", ...)
            title(main="", sub=paste("(",x$type, ")"))
            for (j in 1:n) {
                for (i in j:n) {
                    if (x$contact.network[j, i]==1) {
                        lines(c(x$location[j, 1], x$location[i, 1]), c(x$location[j, 2], x$location[i, 2]))
                    }
                }
            }

            mn  <-  max(apply(x$contact.network, 1, sum))/5
            points(x$location, cex=apply(x$contact.network, 1, sum)/mn, col="red", pch=20)
        } else if (x$type == "random"){

            n <- dim(x$contact.network)[1]
            edges1 <- list(NULL)
            for (i in 1:(n-1)) {
                mn <- NULL
                for (j in (i+1):n) {
                    if (x$contact.network[i,j]==1) {
                        mn <- c(mn,c(i,j))
                    }
                }
                edges1[[i]] <- mn
            }

            combine.edges <- NULL
            for (i in 1:length(edges1)) {
                if (length(edges1[[i]]) > 1) {
                    combine.edges <- c(combine.edges, edges1[[i]])
                }
            }
            graph.random <- igraph::graph(edges = combine.edges, directed = FALSE)
            plot(graph.random, ...)
        }

    } else {
        stop("the input does not have the same class", call.=TRUE)
    }
}
