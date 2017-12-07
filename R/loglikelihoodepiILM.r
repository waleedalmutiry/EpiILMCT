loglikelihoodepiILM<-function(type, epidat, kerneltype, kernelmatrix, distancekernel=NULL, suspar=NULL, powersus=NULL, transpar=NULL, powertrans = NULL, kernel.par=NULL, spark=NULL, gamma=NULL, delta=NULL, suscov=NULL, transcov=NULL) {

# check the compartmental framework:
	if (type == "SIR") {

		kernelpar <- vector(mode="double", length=2)		

# check the epidemic data:
		
		if (dim(epidat)[2]!= 4) {
            stop("Error in the format of the epidat. Please check help file", call.=FALSE)
		}		

# check the kernel function and its required parameters:

        if (kerneltype == "distance") {

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num = 2
			} else if (distancekernel=="Cauchy") {
				num = 3
			}

			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par
			kernelpar[2] <- 0

			n   <- length(kernelmatrix[, 1])
			network  <- matrix(0, ncol=n, nrow=n)
            distance <- kernelmatrix

		} else if (kerneltype == "network") {

			if (dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}
			
			num <- 1
				        
   			n   <- length(kernelmatrix[, 1])
			distance  <- matrix(0, ncol=n, nrow=n)
			network   <- kernelmatrix
	        
			kernelpar[1] <- 0
			kernelpar[2] <- 0

		} else if (kerneltype == "both") {
						
            if (!is.list(kernelmatrix) ) {
                stop("Error: kernelmatrix must be a list of two matrices: distance and contact network: kernelmatrix", call. = FALSE)
            }

            if ( length(kernelmatrix)!=2) {
                stop("Error: kernelmatrix must be a list of two matrices: distance and contact network: kernelmatrix", call. = FALSE)
            }


            if ( dim(kernelmatrix[[1]])[1] != dim(kernelmatrix[[1]])[2] ) {
				stop("The distance matrix must be square", call. = FALSE)
			}

			if (dim(kernelmatrix[[2]])[1] != dim(kernelmatrix[[2]])[2] ) {
				stop("The contact network matrix must be square", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of the distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num = 4
			} else if (distancekernel=="Cauchy") {
				num = 5
			}

			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par[1]
			kernelpar[2] <- kernel.par[2]

			n   <- length(kernelmatrix[[1]][, 1])
            distance <- kernelmatrix[[1]]
			network  <- kernelmatrix[[2]]
		
		}

# obtain the number of infected in the epidemic:
		
		ninfected <- sum(epidat[, 2]!=Inf)

# check the susceptibility function and its covariates and parameters:
	   		
		if (is.null(suspar) & is.null(suscov)) {
			suspar <- 1
			suscov <- matrix(rep(1, n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (!is.null(suspar) & is.null(suscov) ) {
			suscov  <- matrix(rep(1, n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (is.null(suspar) & !is.null(suscov) ) {
			stop("Specify the parameters of the susceptibility function: suspar", call. = FALSE)
		} else {
			if (any(suscov<0)) {
                stop("Covariate(s) values of the susceptibility function must be positive: suscov", call.=FALSE)
			}
			nsuspar <- length(suscov[1, ])
		}

# check the transmissibility function and its covariates and parameters:

		if (is.null(transpar) & is.null(transcov)) {
			transpar <- 1
			transcov <- matrix(rep(1, n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (!is.null(transpar) & is.null(transcov) ) {
			transcov <- matrix(rep(1, n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (is.null(transpar) & !is.null(transcov) ) {
            stop("Specify the transmissibility parameters: transpar", call. = FALSE)
		} else {
			if (any(transcov<0)) {
                stop("Covariate(s) values of the transmissibility function must be positive: transcov", call.=FALSE)
			}
			ntranspar <- length(transcov[1, ])
		}

		if (is.null(powersus)) {
			powersus <- rep(1, nsuspar)
		}

		if (is.null(powertrans)) {
			powertrans <- rep(1, ntranspar)
		}

# check the spark term:

		if (is.null(spark)) {
			spark <- 0
		}

# check the parameters of the distribution of the infectious period:
 
        if (!is.matrix(delta) & length(delta)!=2) {
            stop("Error in entering the parameters of the infectious period distribution: delta", call.=FALSE)
        } else if (is.matrix(delta)) {
            if (dim(delta)[1]!=1 & dim(delta)[2]!=2) {
            	stop("Error in entering the parameters of the infectious period distribution: delta", call.=FALSE)
            }
            deltain1 <- delta[1, 1]
            deltain2 <- delta[1, 2]
        } else {
            deltain1 <- delta[1]
			deltain2 <- delta[2]
		}

		datloglik1<-.Fortran("loglikcontilm", 
		n=as.integer(n),  ninfected =as.integer(ninfected), num=as.integer(num), 
		nsuspar=as.integer(nsuspar), ntranspar=as.integer(ntranspar), 
		cc= matrix(as.double(network), n, n), d333=matrix(as.double(distance), n, n), 
		epidat=matrix(as.numeric(epidat), ncol=4, nrow=n), 
		suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n), 
		transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n), 
		suspar=as.numeric(as.double(suspar)), powersus=as.numeric(as.double(powersus)), 
		transpar=as.numeric(as.double(transpar)), powertrans=as.numeric(as.double(powertrans)), 
		kernelpar=as.vector(kernelpar, mode="double"), spark=as.double(spark), 
		deltain1=as.double(deltain1), deltain2=as.double(deltain2), 
		likk=as.double(0), NAOK = TRUE
		)
   
		result1 <- datloglik1$likk
	
		return(result1)

	} else if (type == "SINR") {

		kernelpar <- vector(mode="double", length=2)		
				
# check the epidemic data:

		if (dim(epidat)[2]!= 6) {
			stop("Error in the format of the epidat. Please check help file", call.=FALSE)
		}		
		
# check the kernel function and its required parameters:

        if (kerneltype == "distance") {

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num = 2
			} else if (distancekernel=="Cauchy") {
				num = 3
			}

			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par
			kernelpar[2] <- 0

			n   <- length(kernelmatrix[, 1])
			network  <- matrix(0, ncol=n, nrow=n)
            distance <- kernelmatrix
            
		} else if (kerneltype == "network") {

			if (dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}

			num <- 1
				        
   			n   <- length(kernelmatrix[, 1])
			distance  <- matrix(0, ncol=n, nrow=n)
			network   <- kernelmatrix

			kernelpar[1] <- 0
			kernelpar[2] <- 0
	        
		} else if (kerneltype == "both") {
						
            if (!is.list(kernelmatrix) ) {
                stop("Error: kernelmatrix must be a list of two matrices: distance and contact network: kernelmatrix", call. = FALSE)
            }

            if ( length(kernelmatrix)!=2) {
                stop("Error: kernelmatrix must be a list of two matrices: distance and contact network: kernelmatrix", call. = FALSE)
            }


            if ( dim(kernelmatrix[[1]])[1] != dim(kernelmatrix[[1]])[2] ) {
				stop("The distance matrix must be square", call. = FALSE)
			}

			if (dim(kernelmatrix[[2]])[1] != dim(kernelmatrix[[2]])[2] ) {
				stop("The contact network matrix must be square", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num = 4
			} else if (distancekernel=="Cauchy") {
				num = 5
			}

			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par[1]
			kernelpar[2] <- kernel.par[2]

			n   <- length(kernelmatrix[[1]][, 1])
            distance  <- kernelmatrix[[1]]
            network   <- kernelmatrix[[2]]
		
		}

# obtain the number of infected in the epidemic:

		ninfected <- sum(epidat[, 2]!=Inf)

# check the susceptibility function and its covariates and parameters:

		if (is.null(suspar) & is.null(suscov)) {
			suspar <- 1
			suscov <- matrix(rep(1, n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (!is.null(suspar) & is.null(suscov) ) {
			suscov  <- matrix(rep(1, n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (is.null(suspar) & !is.null(suscov) ) {
            stop("Specify the susceptibility parameters: suspar", call. = FALSE)
		} else {
			if (any(suscov<0)) {
                stop("Covariate(s) values of the susceptibility function must be positive: suscov", call.=FALSE)
			}
			nsuspar <- length(suscov[1, ])
		}

# check the transmissibility function and its covariates and parameters:

		if (is.null(transpar) & is.null(transcov)) {
			transpar <- 1
			transcov <- matrix(rep(1, n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (!is.null(transpar) & is.null(transcov) ) {
			transcov <- matrix(rep(1, n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (is.null(transpar) & !is.null(transcov) ) {
            stop("Specify the transmissibility parameters: transpar", call. = FALSE)
		} else {
			if (any(transcov<0)) {
                stop("Covariate(s) values of the transmissibility function must be positive: transcov", call.=FALSE)
			}
			ntranspar <- length(transcov[1, ])
		}

		if (is.null(powersus)) {
			powersus <- rep(1, nsuspar)
		}

		if (is.null(powertrans)) {
			powertrans <- rep(1, ntranspar)
		}

# check the spark term:

		if (is.null(spark)) {
			spark <- 0
		}
# check the notification effect parameter:

		if (is.null(gamma)) {
			gamma <- 1
		}

# check the parameters of the distributions of the incubation and delay periods:

		if (is.null(delta)) {
			stop("Specify the incubation and delay periods distribution parameters: delta", call. = FALSE)
		} else if (!is.matrix(delta)) {
			stop("The parameters of the incubation and delay periods distributions must be entered as a 2 by 2 matrix,  where each row represent the scale and shape of the gamma distribution: delta", call.=FALSE)
		} else if (dim(delta)[1]!=2 | dim(delta)[2]!=2) {
			stop("The parameters of the incubation and delay periods distributions must be entered as a 2 by 2 matrix,  where each row represent the scale and shape of the gamma distribution: delta", call.=FALSE)
		} else {
			deltain1 <- delta[1, 1]
			deltain2 <- delta[1, 2]
			deltanr1 <- delta[2, 1]
			deltanr2 <- delta[2, 2]
		}

		datloglikk2<-.Fortran("loglikcontilmsinr", 
		n=as.integer(n), ninfected=as.integer(ninfected), num=as.integer(num), 
		nsuspar=as.integer(nsuspar), ntranspar=as.integer(ntranspar), 
		cc= as.matrix(as.double(network), n, n), d3=as.matrix(as.double(distance), n, n), 
		epidat=as.matrix(as.double(epidat), ncol=6, nrow=n), 
		suscov=as.matrix(as.double(suscov), ncol=nsuspar, nrow=n), 
		transcov=as.matrix(as.double(transcov), ncol=ntranspar, nrow=n), 
		suspar=as.numeric(as.double(suspar)), powersus=as.numeric(as.double(powersus)), 
		transpar=as.numeric(as.double(transpar)), powertrans=as.numeric(as.double(powertrans)), 
		kernelpar=as.vector(kernelpar, mode="double"), spark=as.double(spark), 
		gamma=as.double(gamma), deltain1=as.double(deltain1), 
		deltain2=as.double(deltain2), deltanr1=as.double(deltanr1), deltanr2=as.double(deltanr2), likk=as.double(0), NAOK = TRUE
		)
   
		result2 <- datloglikk2$likk
	
		return(result2)
	
	}

}
