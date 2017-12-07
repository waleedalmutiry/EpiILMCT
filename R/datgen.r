datagen <- function(type, kerneltype, kernelmatrix, distancekernel=NULL, initialepi=NULL, tmax=NULL, suspar=NULL, transpar=NULL, powersus=NULL, powertrans=NULL, kernel.par=NULL, spark=NULL, gamma=NULL, delta, suscov=NULL, transcov=NULL) {
    
    if (type == "SIR") {

####### For SIR ILMs:

		kernelpar <- vector(mode="double", length=2)		

	# checking the type of the kernel function and specifying its parameters:
	
		if (kerneltype == "distance") {
			anum <- 1			

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num <- 2
			} else if (distancekernel=="Cauchy") {
				num <- 3
			}

			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par
			kernelpar[2] <- 0
						
			n   <- length(kernelmatrix[,1])
			distance <- kernelmatrix
			network  <- matrix(0, ncol=n, nrow=n)

		} else if (kerneltype == "network") {

			if(dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
                stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}

	    	if (is.null(spark)) {
		        anum <- 2
        	} else {
	        	anum <- 1
	        }
			
			num <- 1
				        
   			n   <- length(kernelmatrix[,1])
			distance  <- matrix(0, ncol=n, nrow=n)
			network   <- kernelmatrix
			
			kernelpar[1] <- 0
			kernelpar[2] <- 0
	        
		} else if (kerneltype == "both") {
			
			anum <- 1
			
			if (!is.list(kernelmatrix) ) {
				stop("Error: kernelmatrix must be a list of two matrices: distance and contact network", call. = FALSE)
			}
			
			if ( length(kernelmatrix)!=2) {
                stop("Error: kernelmatrix must be a list of two matrices: distance and contact network", call. = FALSE)
			}

			if ( dim(kernelmatrix[[1]])[1] != dim(kernelmatrix[[1]])[2] ) {
				stop("The distance matrix must be square", call. = FALSE)
			}

			if (dim(kernelmatrix[[2]])[1] != dim(kernelmatrix[[2]])[2] ) {
				stop("The contact network matrix must be square", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num <- 4
			} else if (distancekernel=="Cauchy") {
				num <- 5
			}

			n   <- length(kernelmatrix[[1]][,1])
			distance <- kernelmatrix[[1]]
			network  <- kernelmatrix[[2]]

			if (is.null(kernel.par)) {
				stop("Specify the value of the kernel parameters: kernel.par", call. = FALSE)
			}

			kernelpar[1] <- kernel.par[1]
			kernelpar[2] <- kernel.par[2]
					
		}

	# checking the format of the entered infectious period parameters:

		if (is.null(delta)) {
			stop("Specify the values of the parameters of the infectious period distribution: delta", call. = FALSE)
		}
		
        if (!is.matrix(delta) & length(delta)!=2) {
            stop("Error in entering the parameters of the infectious period distribution: delta", call.=FALSE)
        } else if (is.matrix(delta)) {
            if (dim(delta)[1]!=1 & dim(delta)[2]!=2) {
                stop("Error in entering the parameters of the infectious period distribution: delta", call.=FALSE)
            }
            deltain1 <- delta[1,1]
            deltain2 <- delta[1,2]
        } else {
            deltain1 <- delta[1]
			deltain2 <- delta[2]
		}

	# checking whether initial infected individuals are specified or not:

		if (is.null(initialepi)) {
			initial      <- sample(seq(1, n), 1)
			observednum  <- 1
			observedepi  <- c(initial,0,0,0)
		} else {
			if (length(initialepi[1,])!= 4) {
				stop("Error: the initial epidemic must be a matrix of 4 columns: id number(s) of individual(s), removal time(s), infectious period(s) and infection time(s) of the initial infected individual(s)", call.=FALSE)
			}
			observedepi  <- initialepi
			observednum  <- length(initialepi[,1])
		}

	# checking the sucseptibility function information:
		
		if (is.null(suspar) & is.null(suscov)) {
			suspar <- 1
			suscov <- matrix(rep(1,n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (!is.null(suspar) & is.null(suscov) ) {
			suscov  <- matrix(rep(1,n), ncol=1, nrow=n)
			nsuspar <- 1
		} else if (is.null(suspar) & !is.null(suscov) ) {
            stop("Specify the susceptibility parameters: suspar", call. = FALSE)
		} else {
			if (any(suscov<0)) {
                stop("Covariate(s) values of the susceptibility function must be positive: suscov", call.=FALSE)
			}
			nsuspar <- length(suscov[1,])
		}

	# checking the transmissibility function information:

		if (is.null(transpar) & is.null(transcov)) {
			transpar <- 1
			transcov <- matrix(rep(1, n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (!is.null(transpar) & is.null(transcov) ) {
			transcov <- matrix(rep(1,n), ncol=1, nrow=n)
			ntranspar <- 1
		} else if (is.null(transpar) & !is.null(transcov) ) {
            stop("Specify the transmissibility parameters: transpar", call. = FALSE)
		} else {
			if (any(transcov<0)) {
                stop("Covariate(s) values of the transmissibility function must be positive: transcov", call.=FALSE)
			}
			ntranspar <- length(transcov[1,])
		}

		if (is.null(powersus)) {
			powersus <- rep(1, nsuspar)
		}

		if (is.null(powertrans)) {
			powertrans <- rep(1, ntranspar)
		}

	# checking the spark term:

		if (is.null(spark)) {
			spark <- 0
		}

	# checking the maximum infection time of the epidemic:

		if (is.null(tmax)) {
			tmax <- 1000.0
		}

	# defining the fortran function:

		datgg<-.Fortran("datasimulation",	
		n=as.integer(n),anum=as.integer(anum),num=as.integer(num),observednum=as.integer(observednum),
		observedepi=as.matrix(as.double(observedepi),ncol=4,nrow=observednum),tmax=as.double(tmax),
		suspar=as.vector(suspar,mode="double"),nsuspar=as.integer(nsuspar),powersus=as.vector(powersus,mode="double"),
		transpar=as.vector(transpar,mode="double"),ntranspar=as.integer(ntranspar),
		powertrans=as.vector(powertrans,mode="double"),
		kernelpar=as.vector(kernelpar,mode="double"),spark=as.double(spark),delta1=as.double(deltain1),
		delta2=as.double(deltain2),suscov=as.matrix(as.double(suscov),ncol=nsuspar,nrow=n),
		transcov=as.matrix(as.double(transcov),ncol=ntranspar,nrow=n),cc=as.matrix(as.double(network),n,n),
		d3=as.matrix(as.double(distance),n,n),epidat=matrix(0,ncol=4,nrow=n)
		)
	  
	# The output of the fortran function:

		result3 <- datgg$epidat
		return(result3)

    } else if (type == "SINR") {

####### For the SINR ILMs:

		kernelpar <- vector(mode="double", length=2)		

	# checking the type of the kernel function and specifying its parameters:

		if (kerneltype == "distance") {
			anum <- 1			

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num <- 2
			} else if (distancekernel=="Cauchy") {
				num <- 3
			}

			n   <- length(kernelmatrix[,1])
			network  <- matrix(0, ncol=n, nrow=n)
			distance <- kernelmatrix
			
			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
			}
			
			kernelpar[1] <- kernel.par
			kernelpar[2] <- 0

		} else if (kerneltype == "network") {

			if (dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}

	    	if (is.null(spark)) {
		        anum <- 2
        	} else {
	        	anum <- 1
	        }
			
			num <- 1
				        
   			n   <- length(kernelmatrix[,1])
			distance  <- matrix(0, ncol=n, nrow=n)
			network   <- kernelmatrix
			
			kernelpar[1] <- 0
			kernelpar[2] <- 0

		} else if (kerneltype == "both") {
			
			anum <- 1
			
			if (!is.list(kernelmatrix) ) {
				stop("Error: kernelmatrix must be a list of two matrices: distance and contact network", call. = FALSE)
			}
			
			if ( length(kernelmatrix)!=2) {
				stop("Error: kernelmatrix must be a list of two matrices: distance and contact network", call. = FALSE)
			}

			if ( dim(kernelmatrix[[1]])[1] != dim(kernelmatrix[[1]])[2] ) {
				stop("The distance matrix must be square", call. = FALSE)
			}

			if (dim(kernelmatrix[[2]])[1] != dim(kernelmatrix[[2]])[2] ) {
				stop("The contact network matrix must be square", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\"", call. = FALSE)
			} else if (distancekernel=="powerlaw") {
				num <- 4
			} else if (distancekernel=="Cauchy") {
				num <- 5
			}

			n   <- length(kernelmatrix[[1]][,1])
			network  <- kernelmatrix[[2]]
			distance <- kernelmatrix[[1]]
			
			if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameters: kernel.par", call. = FALSE)
			}
			kernelpar[1] <- kernel.par[1]
			kernelpar[2] <- kernel.par[2]
			
		
		}

	# checking the format of the entered incubation and delay periods parameters:

        if (!is.matrix(delta)) {
            stop("The parameters of the incubation and delay periods distributions (delta) must be entered as a 2 by 2 matrix", call.=FALSE)
        } else {
            
            if (dim(delta)[1]!=2 & dim(delta)[2]!=2) {
                stop("Error in entering the parameters of the incubation and delay periods distributions: delta", call.=FALSE)
            }
            
            deltain1 <- delta[1,1]
            deltain2 <- delta[1,2]
            deltanr1 <- delta[2,1]
            deltanr2 <- delta[2,2]

        }
        
	# checking whether initial infected individuals are specified or not:

		if (is.null(initialepi)) {
			initial      <- sample(seq(1, n), 1)
			observednum  <- 1
			observedepi  <- c(initial, 0, 0, 0, 0, 0)
		} else {
			if (is.vector(initialepi)) {
				if (length(initialepi)!= 6) {
					stop("Error: the initial epidemic must be a vector of the id number of individual, removal time, delay period, notification time, incubation period, and infection time of the initial infected individual", call.=FALSE)
				}
			} else {			
				if (length(initialepi[1,])!= 6) {
                    stop("Error: the initial epidemic must be a matrix of 6 columns: id number(s) of individual(s), removal time(s), delay period(s), notification time(s), incubation period(s), and infection time(s) of the initial infected individual(s)", call.=FALSE)
				}
			}
			observedepi  <- initialepi
			observednum  <- length(initialepi[,1])
		}

	# checking the sucseptibility function information:
			
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
			nsuspar <- length(suscov[1,])
		}

	# checking the transmissibility function information:

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
			ntranspar <- length(transcov[1,])
		}

		if (is.null(powersus)) {
			powersus <- rep(1, nsuspar)
		}

		if (is.null(powertrans)) {
			powertrans <- rep(1, ntranspar)
		}

	# checking the spark term:

		if (is.null(spark)) {
			spark <- 0
		}

	# checking the control measure parameter:

		if (is.null(gamma)) {
			gamma <- 1
		}

	# checking the maximum infection time of the epidemic:

		if (is.null(tmax)) {
			tmax <- 1000.0
		}

	# defining the fortran function:

		datgg<-.Fortran("datasimulationsinr",	
		n=as.integer(n),anum=as.integer(anum),num=as.integer(num),observednum=as.integer(observednum),
		observedepi=as.matrix(as.double(observedepi),ncol=6,nrow=observednum),tmax=as.double(tmax),
		suspar=as.vector(suspar,mode="double"),nsuspar=as.integer(nsuspar),powersus=as.vector(powersus,mode="double"),
		transpar=as.vector(transpar,mode="double"),ntranspar=as.integer(ntranspar),
		powertrans=as.vector(powertrans,mode="double"),
		kernelpar=as.vector(kernelpar,mode="double"),
		spark=as.double(spark),gamma=as.double(gamma),deltain1=as.double(deltain1),
		deltain2=as.double(deltain2),deltanr1=as.double(deltanr1),deltanr2=as.double(deltanr2),
		suscov=as.matrix(as.double(suscov),ncol=nsuspar,nrow=n),transcov=as.matrix(as.double(transcov),ncol=ntranspar,nrow=n),
		cc=as.matrix(as.double(network),n,n),d3=as.matrix(as.double(distance),n,n),epidat=matrix(0,ncol=6,nrow=n) )
	  
	# The output of the fortran function:

		result3 <- datgg$epidat
		return(result3)

    } else {
        stop("Error in specifying the compartmental framework of the model: type", call. = FALSE)
    }



}
