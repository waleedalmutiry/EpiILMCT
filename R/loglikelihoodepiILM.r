loglikelihoodepiILM<-function(object, distancekernel=NULL, control.sus=NULL, control.trans=NULL,
kernel.par=NULL, spark=NULL, gamma=NULL, delta=NULL) {

    if (class(object) != "datagen") {
        stop("The epidat object must be in a class of \"datagen\" ", call. = FALSE)
    } else {

        # check the compartmental framework:

        if (object$type == "SIR") {

        kernelpar <- vector(mode="double", length=2)

        # check the epidemic data:

        epidat = object$epidat

        # check the kernel function and its required parameters:

        if (object$kerneltype == "distance") {
            
            n   <-  length(object$location[, 1])
            dis  <- as.matrix(dist(object$location,  method = "euclidean"))
            net  <- matrix(0, ncol = n, nrow = n)
            
            if (is.null(distancekernel)) {
                stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
            } else if (distancekernel=="powerlaw") {
                num  <-  as.integer(2)
            } else if (distancekernel=="Cauchy") {
                num  <-  as.integer(3)
            }

            if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
            }
            
            kernelpar[1] <- kernel.par
            kernelpar[2] <- 0

        } else if (object$kerneltype == "network") {
            
            n   <-  dim(object$network)[1]
            dis  <- matrix(0, ncol = n, nrow = n)
            net  <- object$network
            num <-  1
            
            kernelpar[1] <- 0
            kernelpar[2] <- 0

        } else if (object$kerneltype == "both") {
            
            n   <-  length(object$location[, 1])
            dis  <- as.matrix(dist(object$location,  method = "euclidean"))
            net  <- object$network
            
            if (is.null(distancekernel)) {
                stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
            } else if (distancekernel=="powerlaw") {
                num  <-  4
            } else if (distancekernel=="Cauchy") {
                num  <-  5
            }

            if (is.null(kernel.par)) {
                stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
            }
            
            kernelpar[1] <- kernel.par[1]
            kernelpar[2] <- kernel.par[2]
        
        }

        # obtain the number of infected in the epidemic:
        
        ninfected  <-  sum(object$epidat[, 2]!=Inf)

        # check the susceptibility function and its covariates and parameters:
        
        if (is.null(control.sus)){
            
            nsuspar  <-  1
            suspar <- 1
            suscov  <-  matrix(rep(1, n), ncol= nsuspar, nrow = n)
            powersus <- 1
            
        } else if (!is.list(control.sus)) {
            
            stop("The option control.sus must be a list of values of the susceptibility parameters, susceptibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.sus)) {
            
            len <- length(control.sus)
            
            if (len < 2L) {
                stop("The length of the option control.sus must be at least 2", call. = FALSE)
            } else {
                
                suscov  <- control.sus[[2]]
                
                if (is.matrix(suscov)) {
                    nsuspar <- dim(suscov)[2]
                } else {
                    nsuspar <- 1
                }
                
                if (any(suscov < 0)) {
                    stop("Covariate(s) values of the susceptibility function must be positive: control.sus", call.=FALSE)
                }
                
                if (length(control.sus[[1]]) != nsuspar) {
                    stop("The number of susceptibility parameters must be equal to the number of susceptibility covariates", call. = FALSE)
                } else {
                    suspar    <-  control.sus[[1]]
                }
            }
            
            if (len == 3L) {
                if (length(control.sus[[3]]) != nsuspar) {
                    stop("The number of power parameters of the susceptibility function must be equal to the number of susceptibility covariates", call. = FALSE)
                } else {
                    powersus <-  control.sus[[3]]
                    
                }
            } else {
                powersus <-  rep(1, nsuspar)
            }
        }

        # check the transmissibility function and its covariates and parameters:

        if (is.null(control.trans)){
            
            ntranspar  <-  1
            transpar <- 1
            transcov  <-  matrix(rep(1, n), ncol= ntranspar, nrow = n)
            powertrans <- 1
            
        } else if (!is.list(control.trans)) {
            
            stop("The option control.trans must be a list of values of the transmissibility parameters, transmissibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.trans)) {
            
            len <- length(control.trans)
            
            if (len < 2L) {
                stop("The length of the option control.trans must be at least 2", call. = FALSE)
            } else {
                
                transcov  <- control.trans[[2]]
                
                if (is.matrix(transcov)) {
                    ntranspar <- dim(transcov)[2]
                } else {
                    ntranspar <- 1
                }
                
                if (any(transcov < 0)) {
                    stop("Covariate(s) values of the transmissibility function must be positive: control.trans", call.=FALSE)
                }
                
                if (length(control.trans[[1]]) != ntranspar) {
                    stop("The number of transmissibility parameters must be equal to the number of transceptibility covariates", call. = FALSE)
                } else {
                    transpar    <-  control.trans[[1]]
                }
            }
            
            if (len == 3L) {
                if (length(control.trans[[3]]) != ntranspar) {
                    stop("The number of power parameters of the transmissibility function must be equal to the number of transmissibility covariates", call. = FALSE)
                } else {
                    powertrans <-  control.trans[[3]]
                    
                }
            } else {
                powertrans <-  rep(1, ntranspar)
            }
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
		cc= matrix(as.double(net), n, n), d333=matrix(as.double(dis), n, n),
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

        } else if (object$type == "SINR") {

            kernelpar <- vector(mode="double", length=2)

            # check the epidemic data:

            epidat = object$epidat

            # check the kernel function and its required parameters:

            if (object$kerneltype == "distance") {
                
                n   <-  length(object$location[, 1])
                dis  <- as.matrix(dist(object$location,  method = "euclidean"))
                net  <- matrix(0, ncol = n, nrow = n)
                
                if (is.null(distancekernel)) {
                    stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
                } else if (distancekernel=="powerlaw") {
                    num  <-  as.integer(2)
                } else if (distancekernel=="Cauchy") {
                    num  <-  as.integer(3)
                }
                
                if (is.null(kernel.par)) {
                    stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
                }
                
                kernelpar[1] <- kernel.par
                kernelpar[2] <- 0
                
            } else if (object$kerneltype == "network") {
                
                n   <-  dim(object$network)[1]
                dis  <- matrix(0, ncol = n, nrow = n)
                net  <- object$network
                num <-  1
                
                kernelpar[1] <- 0
                kernelpar[2] <- 0
                
            } else if (object$kerneltype == "both") {
                
                n   <-  length(object$location[, 1])
                dis  <- as.matrix(dist(object$location,  method = "euclidean"))
                net  <- object$network
                
                if (is.null(distancekernel)) {
                    stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)
                } else if (distancekernel=="powerlaw") {
                    num  <-  4
                } else if (distancekernel=="Cauchy") {
                    num  <-  5
                }
                
                if (is.null(kernel.par)) {
                    stop("Specify the value of the kernel parameter: kernel.par", call. = FALSE)
                }
                
                kernelpar[1] <- kernel.par[1]
                kernelpar[2] <- kernel.par[2]
                
            }

            # obtain the number of infected in the epidemic:

            ninfected  <-  sum(object$epidat[, 2]!=Inf)

            # check the susceptibility function and its covariates and parameters:

            if (is.null(control.sus)){
                
                nsuspar  <-  1
                suspar <- 1
                suscov  <-  matrix(rep(1, n), ncol= nsuspar, nrow = n)
                powersus <- 1
                
            } else if (!is.list(control.sus)) {
                
                stop("The option control.sus must be a list of values of the susceptibility parameters, susceptibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
                
            } else if (is.list(control.sus)) {
                
                len <- length(control.sus)
                
                if (len < 2L) {
                    stop("The length of the option control.sus must be at least 2", call. = FALSE)
                } else {
                    
                    suscov  <- control.sus[[2]]
                    
                    if (is.matrix(suscov)) {
                        nsuspar <- dim(suscov)[2]
                    } else {
                        nsuspar <- 1
                    }
                    
                    if (any(suscov < 0)) {
                        stop("Covariate(s) values of the susceptibility function must be positive: control.sus", call.=FALSE)
                    }
                    
                    if (length(control.sus[[1]]) != nsuspar) {
                        stop("The number of susceptibility parameters must be equal to the number of susceptibility covariates", call. = FALSE)
                    } else {
                        suspar    <-  control.sus[[1]]
                    }
                }
                
                if (len == 3L) {
                    if (length(control.sus[[3]]) != nsuspar) {
                        stop("The number of power parameters of the susceptibility function must be equal to the number of susceptibility covariates", call. = FALSE)
                    } else {
                        powersus <-  control.sus[[3]]
                        
                    }
                } else {
                    powersus <-  rep(1, nsuspar)
                }
            }

            # check the transmissibility function and its covariates and parameters:

            if (is.null(control.trans)){
                
                ntranspar  <-  1
                transpar <- 1
                transcov  <-  matrix(rep(1, n), ncol= ntranspar, nrow = n)
                powertrans <- 1
                
            } else if (!is.list(control.trans)) {
                
                stop("The option control.trans must be a list of values of the transmissibility parameters, transmissibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
                
            } else if (is.list(control.trans)) {
                
                len <- length(control.trans)
                
                if (len < 2L) {
                    stop("The length of the option control.trans must be at least 2", call. = FALSE)
                } else {
                    
                    transcov  <- control.trans[[2]]
                    
                    if (is.matrix(transcov)) {
                        ntranspar <- dim(transcov)[2]
                    } else {
                        ntranspar <- 1
                    }
                    
                    if (any(transcov < 0)) {
                        stop("Covariate(s) values of the transmissibility function must be positive: control.trans", call.=FALSE)
                    }
                    
                    if (length(control.trans[[1]]) != ntranspar) {
                        stop("The number of transmissibility parameters must be equal to the number of transceptibility covariates", call. = FALSE)
                    } else {
                        transpar    <-  control.trans[[1]]
                    }
                }
                
                if (len == 3L) {
                    if (length(control.trans[[3]]) != ntranspar) {
                        stop("The number of power parameters of the transmissibility function must be equal to the number of transmissibility covariates", call. = FALSE)
                    } else {
                        powertrans <-  control.trans[[3]]
                        
                    }
                } else {
                    powertrans <-  rep(1, ntranspar)
                }
            }

            # check the spark term:

            if (is.null(spark)) {
                spark <- 0
            }

            # check the parameters of the distribution of the infectious period:

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
            cc= as.matrix(as.double(net), n, n), d3=as.matrix(as.double(dis), n, n),
            epidat=as.matrix(as.double(epidat), ncol=6, nrow=n),
            suscov=as.matrix(as.double(suscov), ncol=nsuspar, nrow=n),
            transcov=as.matrix(as.double(transcov), ncol=ntranspar, nrow=n),
            suspar=as.numeric(as.double(suspar)), powersus=as.numeric(as.double(powersus)),
            transpar=as.numeric(as.double(transpar)), powertrans=as.numeric(as.double(powertrans)),
            kernelpar=as.vector(kernelpar, mode="double"), spark=as.double(spark),
            gamma=as.double(gamma), deltain1=as.double(deltain1),
            deltain2=as.double(deltain2), deltanr1=as.double(deltanr1), deltanr2=as.double(deltanr2),
            likk=as.double(0), NAOK = TRUE
            )

            result2 <- datloglikk2$likk

            return(result2)

        } else {
            stop("Error in specifying the compartmental framework of the model: type", call. = FALSE)
        }
    }
}
