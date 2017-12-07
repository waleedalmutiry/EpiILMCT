epictmcmc <-function(epidat, type, kerneltype, kernelmatrix, distancekernel=NULL, datatype, blockupdate=NULL, nsim, sus.par=NULL, trans.par=NULL, power.sus=NULL, power.trans=NULL, kernel.par=NULL, spark.par=NULL, delta=NULL, gamma.par=NULL, suscov=NULL, transcov=NULL, periodproposal=NULL) {
  
#Declaring variables:

    if (type == "SIR") {

		infperiodproposal  <- vector(mode="double", length=2)		
		delta2prior  <- vector(mode="double", length=2)		

        if (datatype == "known removal") {

				anum11  <- 1

				if (is.null(delta)) {
					stop("Specify the arguments of the parameters of the infectious period distribution: delta",  call. =FALSE)
				} else {
					if ( length(delta)!=4) {
					stop("Error in entering the arguments of the parameters of the infectious period distribution: delta", call.=FALSE)
					}

					delta1 <- delta[1]
					delta2 <- delta[2]
					delta2prior  <- delta[3:4]

				}
			
				if (is.null(periodproposal) ) {
					stop("Specify the parameters of the gamma proposal distribution for updating the infectious period and infection times: periodproposal", call.=FALSE)
				} else {
					if (!is.matrix(periodproposal)) {
						periodproposal  <- matrix(periodproposal, ncol=2, nrow=1)
						infperiodproposal  <- periodproposal[1, ]
					} else if (dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) {
					stop("Error: the parameters of the gamma proposal distribution for updating the infectious periods and infection times should be entered as a 1 by 2 matrix or as a vector: periodproposal", call.=FALSE)
					} else {
						infperiodproposal  <- periodproposal[1, ]			
					}
				}
			 
				if (is.null(blockupdate) ) {
					blockupdate  <- c(1, 1)
				 }
			 
			} else if (datatype == "known epidemic") {

				blockupdate  <- vector(mode="integer", length=2)
				
				anum11  <- 2
				infperiodproposal  <- c(1, 1)
				delta2prior  <- c(0, 0)
				delta1  <- 0.0
				delta2  <- 0.0
				blockupdate  <- c(1, 1)

			} else {
				stop("Specify datatype as \"known removal\" or \"known epidemic\" ",  call. = FALSE)
			}


		kernelpar <- vector(mode="double", length=2)
		kernelparproposalvar <- vector(mode="double", length=2)
		priordistkernelparpar <- vector(mode="integer", length=2)
		kernelparprior <- matrix(0, ncol=2, nrow=2)
		sparkprior <- vector(mode="double", length=2)
		
        if (kerneltype == "distance") {

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("Error: the distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)			
			} else if (distancekernel=="powerlaw") {
				num  <- 2
			} else if (distancekernel=="Cauchy") {
				num  <- 3
			}
			
			n   <- length(kernelmatrix[, 1])
			dis  <- kernelmatrix
			net  <- matrix(0, ncol=n, nrow=n)

			anum55  <- 1

			if (is.null(kernel.par)) {
				stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			if (length(kernel.par)!=5) {
				stop("Error in entering one or more of the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			kernelpar[1] <- kernel.par[1]
			kernelpar[2] <- 0
			kernelparproposalvar[1] <- kernel.par[5]
			kernelparproposalvar[2] <- 0
			kernelparprior[1, ] <- kernel.par[3:4]
			kernelparprior[2, ] <- c(0, 0)

			if (kernel.par[2] == "gamma") {
				priordistkernelparpar[1]  <- 1
			} else if (kernel.par[2] == "half normal") {
				priordistkernelparpar[1]  <- 2
			} else if (kernel.par[2] == "uniform") {
				priordistkernelparpar[1]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}

				priordistkernelparpar[2]  <- 1
				

			if (is.null(spark.par)) {

				spark  <- 0
				anum44  <- 2
				sparkproposalvar  <- 0
				priordistsparkpar  <- 1
				sparkprior  <- c(1, 1)

			} else {

				anum44 <- 1

				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments of the spark parameter: spark.par", call.=FALSE)
				}

				spark <- spark.par[1]
				sparkprior <- spark.par[3:4]
				sparkproposalvar <- spark.par[5]

				if (spark.par[2] == "gamma") {
					priordistsparkpar  <- 1
				} else if (spark.par[2] == "half normal") {
					priordistsparkpar  <- 2
				} else if (spark.par[2] == "uniform") {
					priordistsparkpar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
				}

			}

		} else if (kerneltype == "network") {

			if (dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}
							        
   			n <- length(kernelmatrix[, 1])
			dis <- matrix(0, ncol=n, nrow=n)
			net <- kernelmatrix
			num <- 1
			priordistkernelparpar  <- c(1, 1)
			kernelparproposalvar  <- c(0, 0)
			kernelpar  <- c(0, 0)
			kernelparprior <- matrix(0, ncol=2, nrow=2)

			anum55 <- 3

			if (anum11 == 1) {
				
				anum44   <- 1

				if (is.null(spark.par)) {
					stop("Specify the arguments for updating the spark parameter: spark.par ",  call. =FALSE)
				}

				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments for updating the spark parameter: spark.par", call.=FALSE)
				}

				spark  <- spark.par[1]
				sparkprior <- spark.par[3:4]
				sparkproposalvar <- spark.par[5]

					if (spark.par[2] == "gamma") {
						priordistsparkpar  <- 1
					} else if (spark.par[2] == "half normal") {
						priordistsparkpar  <- 2
					} else if (spark.par[2] == "uniform") {
						priordistsparkpar  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
					}

	
			} else {


				if (is.null(spark.par)) {
			
					spark  <- 0
					anum44  <- 2
					sparkproposalvar <- 0
					priordistsparkpar <- 1
					sparkprior <- c(1, 1)

				} else {
					anum44 <- 1

					if (length(spark.par)!= 5) {
						stop("Error in entering one or more of the arguments for updating the spark parameter: spark.par", call.=FALSE)
					}

					spark <- spark.par[1]
					sparkprior <- spark.par[3:4]
					sparkproposalvar <- spark.par[5]

					if (spark.par[2] == "gamma") {
						priordistsparkpar  <- 1
					} else if (spark.par[2] == "half normal") {
						priordistsparkpar  <- 2
					} else if (spark.par[2] == "uniform") {
						priordistsparkpar  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
					}

				}
	        }
	        
		} else if (kerneltype == "both") {
						
            if (!is.list(kernelmatrix) ) {
                stop("Error: the kernelmatrix must be a list of two matrices: 1) distance and 2) contact network", call. = FALSE)
            }

            if ( length(kernelmatrix)!=2) {
                stop("Error: the kernelmatrix must be a list of two matrices: 1) distance and 2) contact network", call. = FALSE)
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
				num  <- 4
			} else if (distancekernel=="Cauchy") {
				num  <- 5
			}

			n  <- length(kernelmatrix[[1]][, 1])
			dis <- kernelmatrix[[1]]
			net <- kernelmatrix[[2]]
			anum55  <- 2

			if (is.null(kernel.par)) {
				stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			if ( (dim(kernel.par)[1]!=2) & (dim(kernel.par)[2]!=5)) {
				stop("Error in entering one or more of the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			kernelpar[1]  <- kernel.par[1, 1]
			kernelpar[2]  <- kernel.par[2, 1]
			kernelparproposalvar[1]  <- kernel.par[1, 5]
			kernelparproposalvar[2]  <- kernel.par[2, 5]
			kernelparprior[1, ] <- kernel.par[1, 3:4]
			kernelparprior[2, ] <- kernel.par[2, 3:4]
			
			if (kernel.par[1, 2] == "gamma") {
				priordistkernelparpar[1]  <- 1
			} else if (kernel.par[1, 2] == "half normal") {
				priordistkernelparpar[1]  <- 2
			} else if (kernel.par[1, 2] == "uniform") {
				priordistkernelparpar[1]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}
	
			if (kernel.par[2, 2] == "gamma") {
				priordistkernelparpar[2]  <- 1
			} else if (kernel.par[2, 2] == "half normal") {
				priordistkernelparpar[2]  <- 2
			} else if (kernel.par[2, 2] == "uniform") {
				priordistkernelparpar[2]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}

			if (is.null(spark.par)) {
				spark  <- 0
				anum44  <- 2
				sparkproposalvar <- 0
				priordistsparkpar <- 1
				priordistsparkpar <- 1
				sparkprior  <- c(1, 1)

			} else {
				anum44   <- 1

				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments for updating the spark parameter: spark.par", call.=FALSE)
				}

				spark <- spark.par[1]
				sparkprior <- spark.par[3:4]
				sparkproposalvar  <- spark.par[5]

				if (spark.par[2] == "gamma") {
					priordistsparkpar  <- 1
				} else if (spark.par[2] == "half normal") {
					priordistsparkpar  <- 2
				} else if (spark.par[2] == "uniform") {
					priordistsparkpar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
				}
			}		
		}

		# check the number of infected individuals:

		if (anum11 == 1) {
			ninfected  <- sum(epidat[, 2]!=Inf)
		} else if (anum11 == 2) {
			ninfected  <- sum(epidat[, 2]!=Inf)
		} else {
			stop("Error: the epidemic data must be in the same format as in datagen function where the removal and infection times of uninfected individuals are \"Inf\".", call.=FALSE)
		}

		# check Susceptibility terms:


			if (is.null(sus.par) & is.null(suscov)) {

				nsuspar  <- 1
				anum22  <- 2
				suspar  <- 1
				suscov  <- matrix(rep(1, n), ncol= nsuspar, nrow=n)
				susproposalvar  <- 0
				priordistsuspar  <- 1
				priorpar1sus  <- rep(0, nsuspar)
				priorpar2sus  <- rep(0, nsuspar)
				anum66	  <- 2
				powersus <- 1
				powersusproposalvar <- 0
				priordistpowersus <- 1
				priorpar1powersus <- rep(0, nsuspar)
				priorpar2powersus <- rep(0, nsuspar)
  
			} else if (!is.null(sus.par) & is.null(suscov) ) {

				nsuspar  <- 1
				anum22 <- 1
				suscov <- matrix(rep(1, n), ncol= nsuspar, nrow=n)

				if (length(sus.par) != 5) {
					stop("Error in entering one or more of the arguments of the parameters of the susceptibility function: sus.par", call.=FALSE)
				}
				
				suspar <- sus.par[1]
				priorpar1sus <- sus.par[3]
				priorpar2sus <- sus.par[4]
				susproposalvar <- sus.par[5]
				
				anum66 <- 2
				powersus <- 1
				powersusproposalvar  <- 0
				priordistpowersus  <- 1
				priorpar1powersus  <- rep(0, nsuspar)
				priorpar2powersus  <- rep(0, nsuspar)
		
				if (sus.par[2] == "gamma") {
					priordistsuspar  <- 1
				} else if (sus.par[2] == "half normal") {
					priordistsuspar  <- 2
				} else if (sus.par[2] == "uniform") {
					priordistsuspar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: sus.par", call.=FALSE)
				}


			} else if (is.null(sus.par) & !is.null(suscov) ) {

				stop("Specify the arguments of the susceptibility parameters: sus.par", call. = FALSE)

			} else {

				if (any(suscov<0)) {
					stop("Covariate(s) values of the susceptibility function must be positive: suscov", call.=FALSE)
				}

				anum22  	 <- 1
				nsuspar		 <- length(suscov[1, ])
				
				if ((dim(sus.par)[1] != nsuspar) & (dim(sus.par)[2] != 5)) {
					stop("Error in entering one or more of the arguments of the susceptibility parameters: sus.par", call.=FALSE)
				}

				priordistsuspar		 <- vector(mode="integer", length=nsuspar)
				priordistpowersus	 <- vector(mode="integer", length=nsuspar)
				
				suspar          <- sus.par[, 1]
				priorpar1sus    <- sus.par[, 3]
				priorpar2sus    <- sus.par[, 4]
				susproposalvar  <- sus.par[, 5]
				
				for(i in 1:nsuspar) {
					if (sus.par[i, 2] == "gamma") {
						priordistsuspar[i]  <- 1
					} else if (sus.par[i, 2] == "half normal") {
						priordistsuspar[i]  <- 2
					} else if (sus.par[i, 2] == "uniform") {
						priordistsuspar[i]  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: sus.par", call.=FALSE)
					}
				}

				
				if (is.null(power.sus)) {
					anum66 				 <- 2
					powersus			 <- rep(1, nsuspar)
					powersusproposalvar  <- rep(0, nsuspar)
					priordistpowersus 	 <- rep(1, nsuspar)
					priorpar1powersus 	 <- rep(0, nsuspar)
					priorpar2powersus 	 <- rep(0, nsuspar)

				} else {

					anum66  <- 1

					if ((dim(power.sus)[1] != nsuspar) & (dim(power.sus)[2] != 5)) {
						stop("Error in entering one or more of the arguments of the power parameters of the susceptibility function: power.sus", call.=FALSE)
					}

					 powersus              <- power.sus[, 1]
					 priorpar1powersus     <- power.sus[, 3]
					 priorpar2powersus     <- power.sus[, 4]
					 powersusproposalvar   <- power.sus[, 5]
		
					for(i in 1:nsuspar) {
						if (power.sus[i, 2] == "gamma") {
							priordistpowersus[i]  <- 1
						} else if (power.sus[i, 2] == "half normal") {
							priordistpowersus[i]  <- 2
						} else if (power.sus[i, 2] == "uniform") {
							priordistpowersus[i]  <- 3
						} else {
							stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: power.sus", call.=FALSE)
						}
					}
					
				}

			}

		# check Transmissibility terms:
		
			if (is.null(trans.par) & is.null(transcov)) {

				ntranspar 				 <- 1
				
				anum33    				 <- 2
				transpar  				 <- 1
				transcov  				 <- matrix(rep(1, n), ncol= ntranspar, nrow=n)
				transproposalvar 		 <- 0
				priordisttranspar 		 <- 1
				priorpar1trans 			 <- rep(0, ntranspar)
				priorpar2trans 			 <- rep(0, ntranspar)
				anum77	  				 <- 2
				powertrans				 <- 1
				powertransproposalvar 	 <- 0
				priordistpowertrans 	 <- 1
				priorpar1powertrans 	 <- rep(0, ntranspar)
				priorpar2powertrans 	 <- rep(0, ntranspar)
  
			} else if (!is.null(trans.par) & is.null(transcov) ) {

				anum33     <- 1
				ntranspar  <- 1
				transcov   <- matrix(rep(1, n), ncol= ntranspar, nrow=n)

				if (length(trans.par) != 5) {
					stop("Error in entering one or more of the arguments of the transmissibilty parameters: trans.par", call.=FALSE)
				}

				transpar         		 <- trans.par[1]
				priorpartrans    		 <- trans.par[3:4]
				transproposalvar 		 <- trans.par[5]

				anum77	  				 <- 2
				powertrans				 <- 1
				powertransproposalvar 	 <- 0
				priordistpowertrans 	 <- 1
				priorpar1powertrans 	 <- rep(0, ntranspar)
				priorpar2powertrans 	 <- rep(0, ntranspar)

		
				if (trans.par[2] == "gamma") {
					priordisttranspar  <- 1
				} else if (trans.par[2] == "half normal") {
					priordisttranspar  <- 2
				} else if (trans.par[2] == "uniform") {
					priordisttranspar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: trans.par", call.=FALSE)
				}

				priorpar1trans  <- priorpartrans[1]
				priorpar2trans  <- priorpartrans[2]

			} else if (is.null(trans.par) & !is.null(transcov) ) {

				stop("Specify arguments of the transmissibility parameters: trans.par", call. = FALSE)

			} else {

				if (any(transcov<0)) {
					stop("Covariate(s) values of the transmissibility function must be positive: transcov", call.=FALSE)
				}

				anum33   <- 1
				ntranspar <- length(transcov[1, ])

				if ((dim(trans.par)[1] != ntranspar) & (dim(trans.par)[2] != 5)) {
					stop("Error in entering one or more of the arguments of the transmissibility parameters: trans.par", call.=FALSE)
				}

				priordisttranspar		 <- vector(mode="integer", length=ntranspar)
				priordistpowertrans		 <- vector(mode="integer", length=ntranspar)

				transpar         		 <- trans.par[, 1]
				priorpartrans    		 <- trans.par[, 3:4]
				transproposalvar 		 <- trans.par[, 5]
	
				for(i in 1:ntranspar) {
					if (trans.par[i, 2] == "gamma") {
						priordisttranspar[i]  <- 1
					} else if (trans.par[i, 2] == "half normal") {
						priordisttranspar[i]  <- 2
					} else if (trans.par[i, 2] == "uniform") {
						priordisttranspar[i]  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: trans.par", call.=FALSE)
					}
				}

				priorpar1trans 			 <- priorpartrans[, 1]
				priorpar2trans 			 <- priorpartrans[, 2]

				
				if (is.null(power.trans)) {
				
					anum77 					 <- 2
					powertrans				 <- rep(1, ntranspar)
					powertransproposalvar 	 <- rep(0, ntranspar)
					priordistpowertrans 	 <- rep(1, ntranspar)
					priorpar1powertrans 	 <- rep(0, ntranspar)
					priorpar2powertrans 	 <- rep(0, ntranspar)
				
				} else {
					anum77 					 <- 1

					if ((dim(power.trans)[1] != ntranspar) & (dim(power.trans)[2] != 5)) {
						stop("Error in entering the arguments of the power parameters of the transmissibility function: power.trans", call.=FALSE)
					}

					 powertrans              <- power.trans[, 1]
					 prior_dist_trans_power  <- power.trans[, 2]
					 priorparpowertrans      <- power.trans[, 3:4]
					 powertransproposalvar   <- power.trans[, 5]
		
					for(i in 1:ntranspar) {
						if (power.trans[i, 2] == "gamma") {
							priordistpowertrans[i]  <- 1
						} else if (power.trans[i, 2] == "half normal") {
							priordistpowertrans[i]  <- 2
						} else if (power.trans[i, 2] == "uniform") {
							priordistpowertrans[i]  <- 3
						} else {
							stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: power.trans", call.=FALSE)
						}
					}
					
					priorpar1powertrans  <- priorparpowertrans[, 1]
					priorpar2powertrans  <- priorparpowertrans[, 2]
				}

			}
	

			anum2  <- c(anum11, anum22, anum33, anum44, anum55, anum66, anum77)
	

		datmcmc22 <-.Fortran("mcmcsir", 
		n=as.integer(n), 
        nsim=as.integer(nsim), 
        ni=as.integer(ninfected), 
		num=as.integer(num), 
        anum2=as.vector(anum2, mode="integer"), 
        nsuspar=as.integer(nsuspar), 
		ntranspar=as.integer(ntranspar), 
        net=matrix(as.double(net), ncol=n, nrow=n), 
		dis=matrix(as.double(dis), ncol=n, nrow=n), 
        epidat=matrix(as.double(epidat), ncol=4, nrow=n), 
		blockupdate=as.vector(blockupdate, mode="integer"), 
		priordistsuspar=as.vector(priordistsuspar, mode="integer"), 
		priordisttranspar=as.vector(priordisttranspar, mode="integer"), 
		priordistkernelparpar=as.vector(priordistkernelparpar, mode="integer"), 
		priordistsparkpar=as.integer(priordistsparkpar), 
		priordistpowersus=as.vector(priordistpowersus, mode="integer"), 
		priordistpowertrans=as.vector(priordistpowertrans, mode="integer"), 
		suspar=as.vector(suspar, mode="double"), 
        suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n), 
		powersus=as.vector(powersus, mode="double"), 
        transpar=as.vector(transpar, mode="double"), 
        transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n), 
		powertrans=as.vector(powertrans, mode="double"), 
		kernelpar=as.vector(kernelpar, mode="double"), 
        spark=as.double(spark), 
		delta1=as.double(delta1), 
        delta2=as.double(delta2), 
		kernelparproposalvar=as.vector(kernelparproposalvar, mode="double"), 
        sparkproposalvar=as.double(sparkproposalvar), 
		susproposalvar=as.vector(susproposalvar, mode="double"), 
		powersusproposalvar=as.vector(powersusproposalvar, mode="double"), 
		transproposalvar=as.vector(transproposalvar, mode="double"), 
		powertransproposalvar=as.vector(powertransproposalvar, mode="double"), 
		infperiodproposal=as.vector(infperiodproposal, mode="double"), 
        priorpar1sus=as.vector(priorpar1sus, mode="double"), 
		priorpar2sus=as.vector(priorpar2sus, mode="double"), 
        priorpar1powersus=as.vector(priorpar1powersus, mode="double"), 
		priorpar2powersus=as.vector(priorpar2powersus, mode="double"), 
		priorpar1trans=as.vector(priorpar1trans, mode="double"), 
        priorpar2trans=as.vector(priorpar2trans, mode="double"), 
		priorpar1powertrans=as.vector(priorpar1powertrans, mode="double"), 
		priorpar2powertrans=as.vector(priorpar2powertrans, mode="double"), 
		kernelparprior=matrix(as.double(kernelparprior), ncol=2, nrow=2), 
		sparkprior=as.vector(sparkprior, mode="double"), 
        delta2prior=as.vector(delta2prior, mode="double"), 
		susparop= matrix(as.double(0), ncol=nsuspar, nrow=nsim), 
        powersusparop=matrix(as.double(0), ncol=nsuspar, nrow=nsim), 
		transparop= matrix(as.double(0), ncol=ntranspar, nrow=nsim), 
		powertransparop=matrix(as.double(0), ncol=ntranspar, nrow=nsim), 
		kernelparop=matrix(as.double(0), ncol=2, nrow=nsim), 
        sparkop=matrix(as.double(0), ncol=1, nrow=nsim), 
		delta1op=matrix(as.double(0), ncol=1, nrow=nsim), 
        delta2op=matrix(as.double(0), ncol=1, nrow=nsim), 
		epidatmcper=matrix(as.double(0), ncol=n, nrow=nsim), 
		epidatmctim=matrix(as.double(0), ncol=n, nrow=nsim), 
		epidatmcrem=matrix(as.double(0), ncol=n, nrow=nsim), 
		loglik=matrix(as.double(0), ncol=1, nrow=nsim), NAOK = TRUE
		)

			names <-c("Alpha_s[1]", "Alpha_s[2]", "Alpha_s[3]", "Alpha_s[4]", "Alpha_s[5]")
			namet <-c("Alpha_t[1]", "Alpha_t[2]", "Alpha_t[3]", "Alpha_t[4]", "Alpha_t[5]")
			namepowers <-c("Psi_s[1]", "Psi_s[2]", "Psi_s[3]", "Psi_s[4]", "Psi_s[5]")
			namepowert <-c("Psi_t[1]", "Psi_t[2]", "Psi_t[3]", "Psi_t[4]", "Psi_t[5]")
			result77   <- NULL
			namecols  <- NULL
			if (anum2[2]==1) {
			result77 <-cbind(result77, datmcmc22$susparop	)
			namecols <-c(namecols, names[1:nsuspar])
			}
			if (anum2[6]==1) {
			result77 <-cbind(result77, datmcmc22$powersusparop)
			namecols <-c(namecols, namepowers[1:nsuspar])
			}
			if (anum2[3]==1) {
			result77 <-cbind(result77, datmcmc22$transparop)
			namecols <-c(namecols, namet[1:ntranspar])
			}
			if (anum2[7]==1) {
			result77 <-cbind(result77, datmcmc22$powertransparop)
			namecols <-c(namecols, namepowert[1:ntranspar])
			}
			if (anum2[4]==1) {
			result77 <-cbind(result77, datmcmc22$sparkop[, 1])
			namecols <-c(namecols, "Spark")
			}
			if (anum2[5]==1) {
			result77 <-cbind(result77, datmcmc22$kernelparop[, 1])
			namecols <-c(namecols, "Spatial parameter")
			} else if (anum2[5]==2) {
			result77 <-cbind(result77, datmcmc22$kernelparop)
			namecols <-c(namecols, c("Spatial parameter", "Network parameter"))
			}
			if (anum2[1]==1) {
			result77 <-cbind(result77, datmcmc22$delta2op[, 1])
			namecols <-c(namecols, "Rate of infectious period")
			}

			result77 <-cbind(result77, datmcmc22$loglik)
			namecols <-c(namecols, "Log-likelihood")

			result77  <- data.frame(result77)
			colnames(result77)  <- namecols
	
			if (anum2[1]==2) {
				list(result77)
			} else if (anum2[1]==1) {
				list(result77, datmcmc22$epidatmctim)
			}

	} else if (type == "SINR") {

		if (datatype == "known removal") {

			anum66  <- 1

			if (is.null(delta)) {
				stop("Specify the arguments of the parameters of the incubation period distribution: delta",  call. =FALSE)
			} else {
				if ( length(delta)!=4 ) {
				stop("Error in entering the arguments of the parameters of the incubation period distribution: delta", call.=FALSE)
				}

			deltain1 		 <- delta[1]
			deltain2 		 <- delta[2]
			deltanr1		 <- 0
			deltanr2 		 <- 0
			deltain2prior 	 <- delta[3:4]
			deltanr2prior 	 <- c(0, 0)

			}
	
			if (is.null(periodproposal) ) {
				stop("Specify the parameters for the gamma proposal distribution for updating the incubation periods and infection times: periodproposal", call.=FALSE)
			} else {
				if (!is.matrix(periodproposal)) {
					periodproposal 		 <- matrix(periodproposal, ncol=2, nrow=1)
					infperiodproposalin  <- periodproposal[1, ]
					infperiodproposalnr  <- c(0, 0)

				} else if (dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) {
				stop("The parameters of the proposal distribution should be entered as a 1 by 2 matrix or as a vector: periodproposal", call.=FALSE)
				} else {
					infperiodproposalin  <- periodproposal[1, ]			
					infperiodproposalnr  <- c(0, 0)
				}
			}

			if (is.null(blockupdate) ) {
				blockupdate  <- c(1, 1)
			 }
			 
		} else if (datatype == "unknown removal") {

				anum66  <- 2

				if (is.null(delta)) {
					stop("Specify the arguments of the parameters of the incubation and delay periods distributions: delta",  call. =FALSE)
				} else {
					if ( (dim(delta)[1]!=2) & (dim(delta)[2]!=4)) {
					stop("Error in entering one or more of the arguments for updating the parameters of the incubation and delay periods: delta must be a 2 by 4 matrix", call.=FALSE)
					}
				deltain1 		 <- delta[1, 1]
				deltain2 		 <- delta[1, 2]
				deltain2prior 	 <- delta[1, 3:4]
				deltanr1 		 <- delta[2, 1]
				deltanr2 		 <- delta[2, 2]
				deltanr2prior 	 <- delta[2, 3:4]
				}
			
			
				if (is.null(periodproposal) ) {
					stop("Specify the parameters of the gamma proposal distribution for updating the incubation and delay periods and the infection and removal times: periodproposal", call.=FALSE)
				} else {
					if (!is.matrix(periodproposal)) {
						periodproposal 		 <- matrix(periodproposal, ncol=2, nrow=2)
						infperiodproposalin  <- periodproposal[1, ]
						infperiodproposalnr  <- periodproposal[2, ]
					} else if (dim(periodproposal)[1]!=2 & dim(periodproposal)[2]!=2) {
					stop("Enter the proposal distribution for updating the incubation and delay periods as a 2 by 2 matrix: periodproposal", call.=FALSE)
					} else {
						infperiodproposalin  <- periodproposal[1, ]			
						infperiodproposalnr  <- periodproposal[2, ]			
					}
				}

				if (is.null(blockupdate) ) {
					blockupdate  <- c(1, 1)
				 }

			} else if (datatype == "known epidemic") {

				anum66 					 <- 3
				infperiodproposalin 	 <- c(1, 1)
				infperiodproposalnr 	 <- c(1, 1)
				deltain2prior 			 <- c(0, 0)
				deltanr2prior 			 <- c(0, 0)
				deltain1 				 <- 0.0
				deltain2 				 <- 0.0
				deltanr1 				 <- 0.0
				deltanr2 				 <- 0.0
				blockupdate  <- c(1, 1)

			} else {
				stop("Specify the data type as \"known removal\",  \"unknown removal\" or \"known epidemic\": datatype ",  call. = FALSE)
			}

		kernelpar 				 <- vector(mode="double", length=2)
		kernelparproposalvar 	 <- vector(mode="double", length=2)
		priordistkernelparpar 	 <- vector(mode="integer", length=2)
		kernelparprior			 <- matrix(0, ncol=2, nrow=2)


        if (kerneltype == "distance") {

			if ( dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The distance matrix must be square: kernelmatrix", call. = FALSE)
			}

			if (is.null(distancekernel)) {
				stop("Specify the type of distance kernel as \"powerlaw\" or \"Cauchy\": distancekernel", call. = FALSE)			
			} else if (distancekernel=="powerlaw") {
				num  <- 2
			} else if (distancekernel=="Cauchy") {
				num  <- 3
			}

			n    <- length(kernelmatrix[, 1])
			dis  <- kernelmatrix
			net   <- matrix(0, ncol=n, nrow=n)

			anum44  <- 1
			
			if (is.null(kernel.par)) {
				stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			if (length(kernel.par)!=5) {
				stop("Error in entering one or more of the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
			}

			kernelpar[1] 			 <- kernel.par[1]
			kernelpar[2] 			 <- 0
			kernelparproposalvar[1]  <- kernel.par[5]
			kernelparproposalvar[2]  <- 0
			kernelparprior[1, ]		 <- kernel.par[3:4]
			kernelparprior[2, ]		 <- c(0, 0)
			

			if (kernel.par[2] == "gamma") {
				priordistkernelparpar[1]  <- 1
			} else if (kernel.par[2] == "half normal") {
				priordistkernelparpar[1]  <- 2
			} else if (kernel.par[2] == "uniform") {
				priordistkernelparpar[1]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}

			priordistkernelparpar[2]  <- 1

			if (is.null(spark.par)) {
				spark 				 <- 0
				anum33  			 <- 2
				sparkproposalvar 	 <- 0
				priordistsparkpar 	 <- 1
				priordistsparkpar 	 <- 1
				sparkprior 			 <- c(1, 1)

			} else {
			
				anum33  			 <- 1
				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments of the spark parameter: spark.par", call.=FALSE)
				}

				spark            	 <- spark.par[1]
				sparkprior       	 <- spark.par[3:4]
				sparkproposalvar 	 <- spark.par[5]

				if (spark.par[2] == "gamma") {
					priordistsparkpar  <- 1
				} else if (spark.par[2] == "half normal") {
					priordistsparkpar  <- 2
				} else if (spark.par[2] == "uniform") {
					priordistsparkpar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
				}
			}

		} else if (kerneltype == "network") {

			if (dim(kernelmatrix)[1] != dim(kernelmatrix)[2] ) {
				stop("The contact network matrix must be square: kernelmatrix", call. = FALSE)
			}
							        
   			n   					 <- length(kernelmatrix[, 1])
			dis  					 <- matrix(0, ncol=n, nrow=n)
			net   					 <- kernelmatrix

			num     				 <- 1
			kernelpar 				 <- c(0, 0)
			anum44  				 <- 3
			priordistkernelparpar 	 <- c(1, 1)
			kernelparproposalvar 	 <- c(0, 0)
			kernelpar 				 <- c(0, 0)
			kernelparprior			 <- matrix(0, ncol=2, nrow=2)

			if (anum66 == 1 | anum66 == 2) {

				anum33  		 <- 1

				if (is.null(spark.par)) {
					stop("Specify the arguments for updating the spark parameter: spark.par",  call. =FALSE)
				}

				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments of the spark parameter: spark.par", call.=FALSE)
				}

				spark             <- spark.par[1]
				sparkprior        <- spark.par[3:4]
				sparkproposalvar  <- spark.par[5]

					if (spark.par[2] == "gamma") {
						priordistsparkpar  <- 1
					} else if (spark.par[2] == "half normal") {
						priordistsparkpar  <- 2
					} else if (spark.par[2] == "uniform") {
						priordistsparkpar  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
					}
	
			} else {


				if (is.null(spark.par)) {
					spark 				 <- 0
					anum33  			 <- 2
					sparkproposalvar 	 <- 0
					priordistsparkpar 	 <- 1
					priordistsparkpar 	 <- 1
					sparkprior 			 <- c(1, 1)

				} else {
					anum33   <- 1

					if (length(spark.par)!= 5) {
						stop("Error in entering one or more of the arguments of the spark parameter: spark.par", call.=FALSE)
					}

					spark            		 <- spark.par[1]
					sparkprior       		 <- spark.par[3:4]
					sparkproposalvar 		 <- spark.par[5]

					if (spark.par[2] == "gamma") {
						priordistsparkpar  <- 1
					} else if (spark.par[2] == "half normal") {
						priordistsparkpar  <- 2
					} else if (spark.par[2] == "uniform") {
						priordistsparkpar  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
					}
				}
			}
	        
		} else if (kerneltype == "both") {
						
            if (!is.list(kernelmatrix) ) {
                stop("Error: the kernelmatrix must be a list of two matrices: 1) distance and 2) contact network matrix", call. = FALSE)
            }

            if ( length(kernelmatrix)!=2) {
                stop("The kernelmatrix must be a list of two matrices: 1) distance and 2) contact network matrix", call. = FALSE)
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
				num  <- 4
			} else if (distancekernel=="Cauchy") {
				num  <- 5
			}

			n   	 <- length(kernelmatrix[[1]][, 1])
			dis 	 <- kernelmatrix[[1]]
			net  	 <- kernelmatrix[[2]]

			anum44 	 <- 2

			if (is.null(kernel.par)) {
				stop("Specify the arguments for updating the kernel parameters: kernel.par",  call.= FALSE)
			}

			if ((dim(kernel.par)[1]!=2) & (dim(kernel.par)[2]!=5)) {
				stop("Error in entering one or more of the arguments for updating the kernel parameters: kernel.par must be a 2 by 5 matrix.",  call.= FALSE)
			}

			kernelpar[1] 			 <- kernel.par[1, 1]
			kernelpar[2] 			 <- kernel.par[2, 1]
			kernelparproposalvar[1]  <- kernel.par[1, 5]
			kernelparproposalvar[2]  <- kernel.par[2, 5]
			kernelparprior[1, ]		 <- kernel.par[1, 3:4]
			kernelparprior[2, ]		 <- kernel.par[2, 3:4]

			if (kernel.par[1, 2] == "gamma") {
				priordistkernelparpar[1]  <- 1
			} else if (kernel.par[1, 2] == "half normal") {
				priordistkernelparpar[1]  <- 2
			} else if (kernel.par[1, 2] == "uniform") {
				priordistkernelparpar[1]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}

			if (kernel.par[2, 2] == "gamma") {
				priordistkernelparpar[2]  <- 1
			} else if (kernel.par[2, 2] == "half normal") {
				priordistkernelparpar[2]  <- 2
			} else if (kernel.par[2, 2] == "uniform") {
				priordistkernelparpar[2]  <- 3
			} else {
				stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
			}

			if (is.null(spark.par)) {
				
				spark 				 <- 0
				anum33  			 <- 2
				sparkproposalvar 	 <- 0
				priordistsparkpar 	 <- 1
				priordistsparkpar 	 <- 1
				sparkprior 			 <- c(1, 1)

			} else {
				
				anum33   <- 1
				
				if (length(spark.par)!= 5) {
					stop("Error in entering one or more of the arguments of the spark parameter: spark.par", call.=FALSE)
				}

				spark             <- spark.par[1]
				sparkprior        <- spark.par[3:4]
				sparkproposalvar  <- spark.par[5]

				if (spark.par[2] == "gamma") {
					priordistsparkpar  <- 1
				} else if (spark.par[2] == "half normal") {
					priordistsparkpar  <- 2
				} else if (spark.par[2] == "uniform") {
					priordistsparkpar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par",  call. =FALSE)
				}

			}
		
		}



		# check the number of infected individuals:

		if (anum66 == 1) {
			ninfected 	 <- sum(epidat[, 2]!=Inf)
		} else if (anum66 == 2) {
			ninfected 	 <- sum(epidat[, 4]!=Inf)
		} else if (anum66 == 3) {
			ninfected 	 <- sum(epidat[, 2]!=Inf)
		}

		# check Susceptibility terms:
		
			if (is.null(sus.par) & is.null(suscov)) {

				nsuspar  <- 1
				anum11    			 <- 2
				suspar  			 <- 1
				suscov  			 <- matrix(rep(1, n), ncol= nsuspar, nrow=n)
				susproposalvar 		 <- 0
				priordistsuspar 	 <- 1
				priorpar1sus 		 <- rep(0, nsuspar)
				priorpar2sus 		 <- rep(0, nsuspar)
				anum77	  			 <- 2
				powersus			 <- 1
				powersusproposalvar  <- 0
				priordistpowersus 	 <- 1
				priorpar1powersus 	 <- rep(0, nsuspar)
				priorpar2powersus 	 <- rep(0, nsuspar)
  
			} else if (!is.null(sus.par) & is.null(suscov) ) {

				anum11    	 <- 1
				nsuspar 	 <- 1
				suscov  	 <- matrix(rep(1, n), ncol= nsuspar, nrow=n)

				if (length(sus.par) != 5) {
					stop("Error in entering the arguments of the susceptibility parameters: sus.par", call.=FALSE)
				}

				suspar         		 <- sus.par[1]
				priorpar1sus   		 <- sus.par[3]
				priorpar2sus   		 <- sus.par[4]
				susproposalvar 		 <- sus.par[5]
				

				anum77	  			 <- 2
				powersus			 <- 1
				powersusproposalvar  <- 0
				priordistpowersus 	 <- 1
				priorpar1powersus 	 <- rep(0, nsuspar)
				priorpar2powersus 	 <- rep(0, nsuspar)
		
				if (sus.par[2] == "gamma") {
					priordistsuspar  <- 1
				} else if (sus.par[2] == "half normal") {
					priordistsuspar  <- 2
				} else if (sus.par[2] == "uniform") {
					priordistsuspar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: sus.par", call.=FALSE)
				}


			} else if (is.null(sus.par) & !is.null(suscov) ) {

				stop("Specify the arguments of the susceptibility parameters: sus.par", call. = FALSE)

			} else {

				if (any(suscov<0)) {
					stop("Covariate(s) values of the susceptibility function must be positive: suscov", call.=FALSE)
				}

				anum11  		 <- 1
				nsuspar			 <- length(suscov[1, ])

				if ((dim(sus.par)[1] != nsuspar) & (dim(sus.par)[2] != 5)) {
					stop("Error in entering the arguments of the susceptibility parameters: sus.par", call.=FALSE)
				}

				priordistsuspar 	 <- vector(mode="integer", length=nsuspar)
				priordistpowersus 	 <- vector(mode="integer", length=nsuspar)

				suspar         		 <- sus.par[, 1]
				priorpar1sus   		 <- sus.par[, 3]
				priorpar2sus   		 <- sus.par[, 4]
				susproposalvar 		 <- sus.par[, 5]
		
				for(i in 1:nsuspar) {
					if (sus.par[i, 2] == "gamma") {
						priordistsuspar[i]  <- 1
					} else if (sus.par[i, 2] == "half normal") {
						priordistsuspar[i]  <- 2
					} else if (sus.par[i, 2] == "uniform") {
						priordistsuspar[i]  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: sus.par", call.=FALSE)
					}
				}
				
				if (is.null(power.sus)) {
					
					anum77 				 <- 2
					powersus			 <- rep(1, nsuspar)
					powersusproposalvar  <- rep(0, nsuspar)
					priordistpowersus 	 <- rep(1, nsuspar)
					priorpar1powersus 	 <- rep(0, nsuspar)
					priorpar2powersus 	 <- rep(0, nsuspar)

				} else {

					anum77  <- 1

					if ((dim(power.sus)[1] != nsuspar) & (dim(power.sus)[2] != 5)) {
						stop("Error in entering one or more of the arguments of the power parameters of the susceptibility function: power.sus", call.=FALSE)
					}

					 powersus              <- power.sus[, 1]
					 priorpar1powersus     <- power.sus[, 3]
					 priorpar2powersus     <- power.sus[, 4]
					 powersusproposalvar   <- power.sus[, 5]
		
					for(i in 1:nsuspar) {
						if (power.sus[i, 2] == "gamma") {
							priordistpowersus[i]  <- 1
						} else if (power.sus[i, 2] == "half normal") {
							priordistpowersus[i]  <- 2
						} else if (power.sus[i, 2] == "uniform") {
							priordistpowersus[i]  <- 3
						} else {
							stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: power.sus", call.=FALSE)
						}
					}
					
				}

			}

		# check Transmissibility terms:
		
			if (is.null(trans.par) & is.null(transcov)) {

				ntranspar  <- 1
				anum22    				 <- 2
				transpar  				 <- 1
				transcov  				 <- matrix(rep(1, n), ncol= ntranspar, nrow=n)
				transproposalvar 		 <- 0
				priordisttranspar 		 <- 1
				priorpar1trans 			 <- rep(0, ntranspar)
				priorpar2trans 			 <- rep(0, ntranspar)
				anum88	  				 <- 2
				powertrans				 <- 1
				powertransproposalvar 	 <- 0
				priordistpowertrans 	 <- 1
				priorpar1powertrans 	 <- rep(0, ntranspar)
				priorpar2powertrans 	 <- rep(0, ntranspar)
  
			} else if (!is.null(trans.par) & is.null(transcov) ) {

				anum22     <- 1
				ntranspar  <- 1
				transcov   <- matrix(rep(1, n), ncol= ntranspar, nrow=n)

				if (length(trans.par) != 5) {
					stop("Error in entering one or more of the arguments of the transmissibilty parameters: trans.par", call.=FALSE)
				}

				transpar         		 <- trans.par[1]
				priorpar1trans   		 <- trans.par[3]
				priorpar2trans   		 <- trans.par[4]
				transproposalvar 		 <- trans.par[5]
				
				
				anum88	  				 <- 2
				powertrans				 <- 1
				powertransproposalvar 	 <- 0
				priordistpowertrans 	 <- 1
				priorpar1powertrans 	 <- rep(0, ntranspar)
				priorpar2powertrans 	 <- rep(0, ntranspar)

		
				if (trans.par[2] == "gamma") {
					priordisttranspar[1]  <- 1
				} else if (trans.par[2] == "half normal") {
					priordisttranspar[1]  <- 2
				} else if (trans.par[2] == "uniform") {
					priordisttranspar[1]  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: trans.par", call.=FALSE)
				}


			} else if (is.null(trans.par) & !is.null(transcov) ) {

				stop("Specify the arguments of the transmissibility parameters: trans.par", call. = FALSE)

			} else {

				if (any(transcov<0)) {
					stop("Covariate(s) values of the transmissibility function must be positive: transcov", call.=FALSE)
				}

				anum22   <- 1
				ntranspar <- length(transcov[1, ])

				if ((dim(trans.par)[1] != ntranspar) & (dim(trans.par)[2] != 5)) {
					stop("Error in entering one or more of the arguments of the transmissibility parameters: trans.par", call.=FALSE)
				}

				priordisttranspar 		 <- vector(mode="integer", length=ntranspar)
				priordistpowertrans 	 <- vector(mode="integer", length=ntranspar)

				transpar         		 <- trans.par[, 1]
				priorpar1trans   		 <- trans.par[, 3]
				priorpar2trans   		 <- trans.par[, 4]
				transproposalvar 		 <- trans.par[, 5]
	
				for(i in 1:ntranspar) {
					if (trans.par[i, 2] == "gamma") {
						priordisttranspar[i]  <- 1
					} else if (trans.par[i, 2] == "half normal") {
						priordisttranspar[i]  <- 2
					} else if (trans.par[i, 2] == "uniform") {
						priordisttranspar[i]  <- 3
					} else {
						stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: trans.par", call.=FALSE)
					}
				}

				
				if (is.null(power.trans)) {

					anum88 					 <- 2
					powertrans				 <- rep(1, ntranspar)
					powertransproposalvar 	 <- rep(0, ntranspar)
					priordistpowertrans 	 <- rep(1, ntranspar)
					priorpar1powertrans 	 <- rep(0, ntranspar)
					priorpar2powertrans 	 <- rep(0, ntranspar)

				} else {
				
					anum88  <- 1
					
					if ((dim(power.trans)[1] != ntranspar) & (dim(power.trans)[2] != 5)) {
						stop("Error in entering the arguments of the power parameters of the transmissibility function: power.trans", call.=FALSE)
					}

					powertrans              <- power.trans[, 1]
					priorpar1powertrans     <- power.trans[, 3]
					priorpar2powertrans     <- power.trans[, 4]
					powertransproposalvar   <- power.trans[, 5]
		
					for(i in 1:ntranspar) {
						if (power.trans[i, 2] == "gamma") {
							priordistpowertrans[i]  <- 1
						} else if (power.trans[i, 2] == "half normal") {
							priordistpowertrans[i]  <- 2
						} else if (power.trans[i, 2] == "uniform") {
							priordistpowertrans[i]  <- 3
						} else {
							stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: power.trans", call.=FALSE)
						}
					}
				}

			}
	
			if (is.null(gamma.par)) {
			
				anum55 				 <- 2
				gamma  				 <- 1
				priordistgammapar	 <- 1
				gammaproposalvar	 <- 0
				gammaprior			 <- c(1, 1)
			
			} else {
			
				anum55  <- 1
			
				if (length(gamma.par)!=5) {
					stop("Error in entering one or more of the arguments for updating the gamma parameter: gamma.par", call.=FALSE)
				}

				gamma           	 <- gamma.par[1]
				gammaprior			 <- gamma.par[3:4]
				gammaproposalvar	 <- gamma.par[5]
				
				if (gamma.par[2] == "gamma") {
					priordistgammapar  <- 1
				} else if (gamma.par[2] == "half normal") {
					priordistgammapar  <- 2
				} else if (gamma.par[2] == "uniform") {
					priordistgammapar  <- 3
				} else {
					stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: gamma.par", call.=FALSE)
				}

			}

			anum2  <- c(anum11, anum22, anum33, anum44, anum55, anum66, anum77, anum88)
	

		datmcmc22 <-.Fortran("mcmcsinr", 
		n=as.integer(n), nsim=as.integer(nsim), ni=as.integer(ninfected), 
		num=as.integer(num), anum2=as.vector(anum2, mode="integer"), nsuspar=as.integer(nsuspar), 
		ntranspar=as.integer(ntranspar), net=matrix(as.double(net), ncol=n, nrow=n), 
		dis=matrix(as.double(dis), ncol=n, nrow=n), epidat=matrix(as.double(epidat), ncol=6, nrow=n), 
		blockupdate=as.vector(blockupdate, mode="integer"), #sizeblock=as.integer(sizeblock), 
		priordistsuspar=as.vector(priordistsuspar, mode="integer"), 
		priordisttranspar=as.vector(priordisttranspar, mode="integer"), 
		priordistkernelparpar=as.vector(priordistkernelparpar, mode="integer"), 
		priordistsparkpar=as.integer(priordistsparkpar), priordistgammapar=as.integer(priordistgammapar), 
		priordistpowersus=as.vector(priordistpowersus, mode="integer"), 
		priordistpowertrans=as.vector(priordistpowertrans, mode="integer"), 
		suspar=as.vector(suspar, mode="double"), suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n), 
		powersus=as.vector(powersus, mode="double"), 
		transpar=as.vector(transpar, mode="double"), transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n), 
		powertrans=as.vector(powertrans, mode="double"), 
		kernelpar=as.vector(kernelpar, mode="double"), spark=as.double(spark), gamma=as.double(gamma), 
		deltain1=as.double(deltain1), deltain2=as.double(deltain2), 
		deltanr1=as.double(deltanr1), deltanr2=as.double(deltanr2), 
		kernelparproposalvar=as.vector(kernelparproposalvar, mode="double"), sparkproposalvar=as.double(sparkproposalvar), 
		gammaproposalvar=as.double(gammaproposalvar), susproposalvar=as.vector(susproposalvar, mode="double"), 
		powersusproposalvar=as.vector(powersusproposalvar, mode="double"), 
		transproposalvar=as.vector(transproposalvar, mode="double"), 
		powertransproposalvar=as.vector(powertransproposalvar, mode="double"), 
		infperiodproposalin=as.vector(infperiodproposalin, mode="double"), 
		infperiodproposalnr=as.vector(infperiodproposalnr, mode="double"), 
		priorpar1sus=as.vector(priorpar1sus, mode="double"), priorpar2sus=as.vector(priorpar2sus, mode="double"), 
		priorpar1powersus=as.vector(priorpar1powersus, mode="double"), 
		priorpar2powersus=as.vector(priorpar2powersus, mode="double"), 
		priorpar1trans=as.vector(priorpar1trans, mode="double"), priorpar2trans=as.vector(priorpar2trans, mode="double"), 
		priorpar1powertrans=as.vector(priorpar1powertrans, mode="double"), 
		priorpar2powertrans=as.vector(priorpar2powertrans, mode="double"), 
		kernelparprior=matrix(as.double(kernelparprior), ncol=2, nrow=2), sparkprior=as.vector(sparkprior, mode="double"), 
		gammaprior=as.vector(gammaprior, mode="double"), deltain2prior=as.vector(deltain2prior, mode="double"), 
		deltanr2prior=as.vector(deltanr2prior, mode="double"), 
		susparop=matrix(as.double(0), ncol=nsuspar, nrow=nsim), powersusparop=matrix(as.double(0), ncol=nsuspar, nrow=nsim), 
		transparop= matrix(as.double(0), ncol=ntranspar, nrow=nsim), 
		powertransparop=matrix(as.double(0), ncol=ntranspar, nrow=nsim), 
		kernelparop=matrix(as.double(0), ncol=2, nrow=nsim), sparkop=matrix(as.double(0), ncol=1, nrow=nsim), 
		gammaop=matrix(as.double(0), ncol=1, nrow=nsim), 
		deltain2op=matrix(as.double(0), ncol=1, nrow=nsim), 
		deltanr2op=matrix(as.double(0), ncol=1, nrow=nsim), 
		epidatmctim=matrix(as.double(0), ncol=n, nrow=nsim), 
		epidatmcrem=matrix(as.double(0), ncol=n, nrow=nsim), 
		loglik=matrix(as.double(0), ncol=1, nrow=nsim), NAOK = TRUE
		)

			names <-c("Alpha_s[1]", "Alpha_s[2]", "Alpha_s[3]", "Alpha_s[4]", "Alpha_s[5]")
			namet <-c("Alpha_t[1]", "Alpha_t[2]", "Alpha_t[3]", "Alpha_t[4]", "Alpha_t[5]")
			namepowers <-c("Psi_s[1]", "Psi_s[2]", "Psi_s[3]", "Psi_s[4]", "Psi_s[5]")
			namepowert <-c("Psi_t[1]", "Psi_t[2]", "Psi_t[3]", "Psi_t[4]", "Psi_t[5]")
			result77   <- NULL
			namecols  <- NULL
			if (anum2[1]==1) {
			result77 <-cbind(result77, datmcmc22$susparop	)
			namecols <-c(namecols, names[1:nsuspar])
			}
			if (anum2[7]==1) {
			result77 <-cbind(result77, datmcmc22$powersusparop)
			namecols <-c(namecols, namepowers[1:nsuspar])
			}
			if (anum2[2]==1) {
			result77 <-cbind(result77, datmcmc22$transparop)
			namecols <-c(namecols, namet[1:ntranspar])
			}
			if (anum2[8]==1) {
			result77 <-cbind(result77, datmcmc22$powertransparop)
			namecols <-c(namecols, namepowert[1:ntranspar])
			}
			if (anum2[3]==1) {
			result77 <-cbind(result77, datmcmc22$sparkop[, 1])
			namecols <-c(namecols, "Spark")
			}
			if (anum2[4]==1) {
			result77 <-cbind(result77, datmcmc22$kernelparop[, 1])
			namecols <-c(namecols, "Spatial parameter")
			} else if (anum2[4]==2) {
			result77 <-cbind(result77, datmcmc22$kernelparop)
			namecols <-c(namecols, c("Spatial parameter", "Network parameter"))
			}
			if (anum2[5]==1) {
			result77 <-cbind(result77, datmcmc22$gammaop)
			namecols <-c(namecols, "Gamma")
			}
			if (anum2[6]==1) {
			result77 <-cbind(result77, datmcmc22$deltain2op[, 1])
			namecols <-c(namecols, "Rate of incubation period")
			} else if (anum2[6]==2) {
			result77 <-cbind(result77, datmcmc22$deltain2op[, 1])
			namecols <-c(namecols, "Rate of incubation period")
			result77 <-cbind(result77, datmcmc22$deltanr2op[, 1])
			namecols <-c(namecols, "Rate of delay period")
			}

			result77 <-cbind(result77, datmcmc22$loglik)
			namecols <-c(namecols, "Log-likelihood")

			result77  <- data.frame(result77)
			colnames(result77)  <- namecols
	
			if (anum2[6]==3) {
				list(result77)
			} else if (anum2[6]==1) {
				list(result77, datmcmc22$epidatmctim)
			} else if (anum2[6]==2) {
				list(result77, datmcmc22$epidatmctim, datmcmc22$epidatmcrem)
			}

    } else {
        stop("Error in specifying the compartmental framework of the model: type", call. = FALSE)
    }

}
