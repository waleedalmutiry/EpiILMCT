epictmcmc <- function(object, distancekernel = NULL, datatype, blockupdate = NULL, nsim, nchains = NULL, control.sus = NULL, control.trans = NULL, kernel.par = NULL, spark.par = NULL, delta = NULL, gamma.par = NULL, periodproposal = NULL, parallel = FALSE, seedval = NULL) {
    
 
 if (class(object) != "datagen") {
     stop("The epidat object must be in a class of \"datagen\" ", call. = FALSE)
 } else {

 # Declaring variables:
    
    if (object$type == "SIR") {

        if (is.null(nchains)) {
            nchains <- 1
        } else {
            nchains <- nchains
        }
        
        # Set seed value for Fortran random number generator
        if (is.null(seedval)){
            if (nchains == 1) {
                temp1 <- 0
            } else {
                temp1 <- sample(seq(100001, 999999), size = nchains, replace = FALSE)
            }
        } else {
            if (nchains == length(seedval)) {
                temp1 <- seedval
            } else {
                stop("The number of seeds is not equal to the number of chains", call. = FALSE)
            }
        }
        
        initial <- list(NULL)
        infperiodproposal  <-  vector(mode="double", length = 2)
        delta2prior  <-  vector(mode="double", length = 2)
        delta1 <-  vector(mode="double", length = 2)
        
        if (datatype == "known removal") {
            
            anum11  <-  1
            
            if (is.null(delta)) {
                stop("Specify the arguments of the parameters of the infectious period distribution: delta",  call. =FALSE)
            } else {
            
                if (!is.list(delta)) {
				   stop("The argument \"delta\" must be a list of three:\n1) a fixed shape parameter of the infectious period density.\n2) a vector of initial values of the rate parameter of the infectious period density with size equal to \"nchains\". \n3) a vector of the parameter values of the gamma prior distribution for the rate parameter.",  call.= FALSE)
                }

                if (length(delta) != 3) {
                    stop("Error in entering the arguments of the delta parameters of the infectious period distribution: delta", call.=FALSE)
                }

                if ( length(delta[[1]])>1) {
                    stop("Error in entering the arguments of the fixed shape parameter of the infectious period density: delta", call.=FALSE)
                }
    
                initial[[1]] <- matrix(0, ncol = nchains, nrow = 2)
                initial[[1]][1,] <-  rep(delta[[1]], nchains)
				
				if (is.vector(delta[[2]])) {
					if (length(delta[[2]]) == nchains) {
						initial[[1]][2,] <-  delta[[2]]
					} else if (length(delta[[2]]) == 1) {
						initial[[1]] <-  rep(delta[[2]],nchains)
					} else {
						stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
					}
				} else {
					stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
				}
												
				if (length(delta[[3]]) == 2) {
					delta2prior  <-  delta[[3]]
				} else {
					stop("Error in entering the parameter values of the gamma prior distribution of the delta parameter: delta[[3]]",  call.= FALSE)
				}
                
            }
            
            if (is.null(periodproposal) ) {
                infperiodproposal  <-  c(0,0)
            } else {
                if (!is.matrix(periodproposal)) {
                    periodproposal  <-  matrix(periodproposal, ncol = 2, nrow = 1)
                    infperiodproposal  <-  periodproposal[1, ]
                } else if (dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) {
                    stop("Error: the parameters of the gamma proposal distribution for updating the infectious periods and infection times should be entered as a 1 by 2 matrix or as a vector: periodproposal", call.=FALSE)
                } else {
                    infperiodproposal  <-  periodproposal[1, ]
                }
            }
            
            if (is.null(blockupdate) ) {
                blockupdate  <-  c(1, 1)
            }
            
        } else if (datatype == "known epidemic") {
            
            blockupdate  <-  vector(mode="integer", length = 2)

            if (!is.null(delta)) {
                warning("The infectious period rate is not updated as the option of datatype = \"known epidemic\".", call. = TRUE)
            }

            anum11  <-  2
            infperiodproposal  <-  c(0,0)
            delta2prior  <-  c(0, 0)
            initial[[1]] <-  matrix(rep(0.0,nchains), ncol = nchains, nrow = 2)
            blockupdate  <-  c(1, 1)
            
        } else {
            stop("Specify datatype as \"known removal\" or \"known epidemic\" ",  call. = FALSE)
        }
        
        
        kernelpar <-  vector(mode="double", length = 2)
        kernelparproposalvar <-  vector(mode="double", length = 2)
        priordistkernelparpar <-  vector(mode="integer", length = 2)
        kernelparprior <-  matrix(0, ncol = 2, nrow = 2)
        sparkprior <-  vector(mode="double", length = 2)
        
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
            
            anum55  <-  1
            
            if (is.null(kernel.par)) {
                stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
            }
            
            if (!is.list(kernel.par)) {
                stop("The argument \"kernel.par\" must be a list of two:\n1) a vector of initial values of the kernel parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            }
            
            if (length(kernel.par) != 2L) {
                stop("The argument \"kernel.par\" must be a list of two:\n1) a vector of initial values of the kernel parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            } else {
                if (length(kernel.par[[2]]) != 4L) {
                    stop("Error in entering the second list of the argument \"kernel.par\" for updating the kernel parameter.",  call.= FALSE)
                }
                kernelparproposalvar[1] <-  kernel.par[[2]][4]
                kernelparproposalvar[2] <-  0
                kernelparprior[1, ] <-  kernel.par[[2]][2:3]
                kernelparprior[2, ] <-  c(0, 0)

                if (kernel.par[[2]][1] == "gamma") {
                    priordistkernelparpar[1]  <-  1
                } else if (kernel.par[[2]][1] == "half normal") {
                    priordistkernelparpar[1]  <-  2
                } else if (kernel.par[[2]][1] == "uniform") {
                    priordistkernelparpar[1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
                }

                priordistkernelparpar[2]  <-  1
                
                if (length(kernel.par[[1]]) == nchains){
                    initial[[5]] <-  matrix(c(kernel.par[[1]],rep(0,nchains)), ncol = nchains, byrow = TRUE)
                } else if(length(kernel.par[[1]]) == 1) {
                    initial[[5]] <-  matrix(rep(c(kernel.par[[1]],0),nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("Error in entering the initial values of the kernel parameter: kernel.par[[1]]",  call.= FALSE)
                }
                
            }
            
            if (is.null(spark.par)) {
                
                initial[[4]] <-  rep(0,nchains)
                anum44  <-  2
                sparkproposalvar  <-  0
                priordistsparkpar  <-  1
                sparkprior  <-  c(1, 1)
                
            } else {
                
                anum44 <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[4]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        initial[[4]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                    
                }
                
            }
            
        } else if (object$kerneltype == "network") {
            
            n   <-  dim(object$network)[1]
            dis  <- matrix(0, ncol = n, nrow = n)
            net  <- object$network
            
            num <-  1
            priordistkernelparpar  <-  c(1, 1)
            kernelparproposalvar  <-  c(0, 0)
            kernelpar  <-  c(0, 0)
            kernelparprior <-  matrix(0, ncol = 2, nrow = 2)
            initial[[5]] <-  matrix(rep(c(0,0),nchains), ncol = nchains, nrow = 2)
            
            anum55 <-  3
            
            if (anum11 == 1) {
                
                anum44   <-  1
                
                if (is.null(spark.par)) {
                    stop("Specify the arguments for updating the spark parameter: spark.par ",  call. =FALSE)
                }
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[4]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        initial[[4]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                }
                
            } else {
                
                
                if (is.null(spark.par)) {
                    
                    anum44  <-  2
                    sparkproposalvar <-  0
                    priordistsparkpar <-  1
                    sparkprior <-  c(1, 1)
                    initial[[4]] <-  rep(0,nchains)
                    
                } else {
                    anum44 <-  1
                    
                    if (!is.list(spark.par)) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    }
                    
                    if (length(spark.par) != 2L) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    } else {
                        if (length(spark.par[[2]]) != 4L) {
                            stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                        }
                        sparkproposalvar <-  spark.par[[2]][4]
                        sparkprior <-  spark.par[[2]][2:3]
                        
                        if (spark.par[[2]][1] == "gamma") {
                            priordistsparkpar  <-  1
                        } else if (spark.par[[2]][1] == "half normal") {
                            priordistsparkpar  <-  2
                        } else if (spark.par[[2]][1] == "uniform") {
                            priordistsparkpar  <-  3
                        } else {
                            stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                        }
                        
                        if (length(spark.par[[1]]) == nchains){
                            initial[[4]] <-  spark.par[[1]]
                        } else if(length(spark.par[[1]]) == 1) {
                            initial[[4]] <-  rep(spark.par[[1]],nchains)
                        } else {
                            stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                        }

                    }
            
                }
            }
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
            
            anum55  <-  2
            
            if (is.null(kernel.par)) {
                stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
            }
            
            if (!is.list(kernel.par)) {
                stop("The argument \"kernel.par\" must be a list of two matrices:\n1) a matrix of initial values of the two kernel parameters of number of columns equal to \"nchains\". \n2) a matrix of two rows contains prior distribution, prior parameter values, and proposal variance for each parameters.",  call.= FALSE)
            }
            
            if (length(kernel.par) != 2L) {
                stop("The argument \"kernel.par\" must be a list of two matrices:\n1) a matrix of initial values of the two kernel parameters of number of columns equal to \"nchains\". \n2) a matrix of two rows contains prior distribution, prior parameter values, and proposal variance for each parameters.",  call.= FALSE)
            } else {
                if (!is.matrix(kernel.par[[2]])) {
                    stop("The second list of the argument \"kernel.par\" must be a matrix of 2 by 4.",  call.= FALSE)
                } else if (dim(kernel.par[[2]])[1] != 2 & dim(kernel.par[[2]])[2] != 4) {
                    stop("The second list of the argument \"kernel.par\" must be a matrix of 2 by 4.",  call.= FALSE)
                }
                kernelparproposalvar[1] <-  kernel.par[[2]][1,4]
                kernelparproposalvar[2] <-  kernel.par[[2]][2,4]
                kernelparprior[1, ] <-  kernel.par[[2]][1,2:3]
                kernelparprior[2, ] <-  kernel.par[[2]][2,2:3]
                
                if (kernel.par[[2]][1,1] == "gamma") {
                    priordistkernelparpar[1]  <-  1
                } else if (kernel.par[[2]][1,1] == "half normal") {
                    priordistkernelparpar[1]  <-  2
                } else if (kernel.par[[2]][1,1] == "uniform") {
                    priordistkernelparpar[1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                
                if (kernel.par[[2]][2,1] == "gamma") {
                    priordistkernelparpar[2]  <-  1
                } else if (kernel.par[[2]][2,1] == "half normal") {
                    priordistkernelparpar[2]  <-  2
                } else if (kernel.par[[2]][2,1] == "uniform") {
                    priordistkernelparpar[2]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\", \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                
                if (is.matrix(kernel.par[[1]])) {
                    if (dim(kernel.par[[1]])[1] != 2 & dim(kernel.par[[1]])[2] != nchains) {
                        stop("The matrix of the initial values of the kernel parameters must be a 2 by \"nchains\".", call. = FALSE)
                    }
                    initial[[5]] <-  kernel.par[[1]]
                } else if (is.vector(kernel.par[[1]]) & length(kernel.par[[1]]) == 2) {
                    initial[[5]] <-  matrix(rep(kernel.par[[1]],nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("The initial values of the kernel parameters must be a matrix of 2 by \"nchains\" or a vector of two initial values.", call. = FALSE)
                }
                
            }

            if (is.null(spark.par)) {
                
                anum44  <-  2
                sparkproposalvar <-  0
                priordistsparkpar <-  1
                sparkprior <-  c(1, 1)
                initial[[4]] <-  rep(0,nchains)
                
            } else {
                anum44 <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[4]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        initial[[4]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                    
                }
                
            }
        }
        
        # check the number of infected individuals:
        
        if (anum11 == 1 | anum11 == 2) {
            ni  <-  sum(object$epidat[, 2]!=Inf)
        } else {
            stop("Error: the epidemic data must be in the same format as in datagen function where the removal and infection times of uninfected individuals are \"Inf\".", call.=FALSE)
        }
        
        
        # check Susceptibility terms:
        ### new ones
        
        if (is.null(control.sus)){
            
            nsuspar  <-  1
            anum22  <-  2
            initial[[2]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            suscov  <-  matrix(rep(1, n), ncol= nsuspar, nrow = n)
            susproposalvar  <-  0
            priordistsuspar  <-  1
            priorpar1sus  <-  rep(0, nsuspar)
            priorpar2sus  <-  rep(0, nsuspar)
            anum66      <-  2
            initial[[6]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            powersusproposalvar <-  0
            priordistpowersus <-  1
            priorpar1powersus <-  rep(0, nsuspar)
            priorpar2powersus <-  rep(0, nsuspar)
            
        } else if (!is.list(control.sus)) {
            
            stop("The option control.sus must be a list of values of the susceptibility parameters, susceptibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.sus)) {
            
            len <- length(control.sus)
            
            if (len < 2L) {
                stop("The length of the option control.sus must be at least 2", call. = FALSE)
            } else {
                
                anum22       <-  1
                sus.par <- control.sus[[1]]
                suscov  <- control.sus[[2]]
                
                if (is.matrix(suscov)) {
                    nsuspar <- dim(suscov)[2]
                } else {
                    nsuspar <- 1
                }
                
                priordistsuspar    <-  vector(mode="integer", length = nsuspar)
                priordistpowersus  <-  vector(mode="integer", length = nsuspar)
                
                if (any(suscov < 0)) {
                    stop("Covariate(s) values of the susceptibility function must be positive: control.sus", call.=FALSE)
                }
                
                if (!is.list(sus.par)) {
                    stop("The first list of the argument \"control.sus\" (the susceptibility parameters) needs to be entered as a list of two matrices:\n1) a number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                } else {
                    if (length(sus.par) != 2) {
                        stop("The first list of the argument \"control.sus\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                    } else {
                        if (is.matrix(sus.par[[2]])) {
                            if ((dim(sus.par[[2]])[1] != nsuspar) & (dim(sus.par[[2]])[2] != 4L)) {
                                stop("Error in entering the argument of \"control.sus\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            priorpar1sus    <-  sus.par[[2]][, 2]
                            priorpar2sus    <-  sus.par[[2]][, 3]
                            susproposalvar  <-  sus.par[[2]][, 4]
                            
                            for(i in 1:nsuspar) {
                                if (sus.par[[2]][i, 1] == "gamma") {
                                    priordistsuspar[i]  <-  1
                                } else if (sus.par[[2]][i, 1] == "half normal") {
                                    priordistsuspar[i]  <-  2
                                } else if (sus.par[[2]][i, 1] == "uniform") {
                                    priordistsuspar[i]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                }
                            }
                            
                        } else {
                            if (nsuspar > 1L){
                                stop("Error in entering the argument of control.sus, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            if (length(sus.par[[2]]) != 4L){
                                stop("Error in entering one or more of the arguments of the susceptibility parameters in the option control.sus", call. = FALSE)
                            }
                            
                            priorpar1sus    <-  sus.par[[2]][2]
                            priorpar2sus    <-  sus.par[[2]][3]
                            susproposalvar  <-  sus.par[[2]][4]
                            
                            if (sus.par[[2]][1] == "gamma") {
                                priordistsuspar[1]  <-  1
                            } else if (sus.par[[2]][1] == "half normal") {
                                priordistsuspar[1]  <-  2
                            } else if (sus.par[[2]][1] == "uniform") {
                                priordistsuspar[1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                            }
                            
                        }
                        if (is.matrix(sus.par[[1]])) {
                            
                            if ((dim(sus.par[[1]])[1] != nsuspar) & (dim(sus.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                            }
                            initial[[2]]    <-  sus.par[[1]]
                            
                        } else if (is.vector(sus.par[[1]]) & length(sus.par[[1]]) == nsuspar) {
                            initial[[2]]    <-  matrix(rep(sus.par[[1]],nchains), ncol = nchains, nrow = nsuspar)
                        }
                        
                    }
                }
                
                if (len == 3L) {
                    power.sus <- control.sus[[3]]
                    anum66               <-  1

                    if (!is.list(power.sus)) {
                        stop("The third list of the argument \"control.sus\" (power parameters of susceptibility covariates) needs to be entered as a list of two matrices:\n1) A number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                    } else {
                        if (length(power.sus) != 2) {
                            stop("The third list of the argument \"control.sus\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                        } else {
                            if (is.matrix(power.sus[[2]])) {
                                if ((dim(power.sus[[2]])[1] != nsuspar) & (dim(power.sus[[2]])[2] != 4L)) {
                                    stop("Error in entering the argument of \"control.sus\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                
                                priorpar1powersus    <-  power.sus[[2]][, 2]
                                priorpar2powersus    <-  power.sus[[2]][, 3]
                                powersusproposalvar  <-  power.sus[[2]][, 4]
                                
                                
                                for(i in 1:nsuspar) {
                                    if (power.sus[[2]][i, 1] == "gamma") {
                                        priordistpowersus[i]  <-  1
                                    } else if (power.sus[[2]][i, 1] == "half normal") {
                                        priordistpowersus[i]  <-  2
                                    } else if (power.sus[[2]][i, 1] == "uniform") {
                                        priordistpowersus[i]  <-  3
                                    } else {
                                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                    }
                                }
                                
                            } else {
                                if (nsuspar > 1L){
                                    stop("Error in entering the argument of control.sus, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                if (length(power.sus[[2]]) != 4L){
                                    stop("Error in entering one or more of the arguments of the susceptibility parameters in the option control.sus", call. = FALSE)
                                }
                                
                                priorpar1powersus     <-  power.sus[[2]][2]
                                priorpar2powersus     <-  power.sus[[2]][3]
                                powersusproposalvar   <-  power.sus[[2]][4]
                                
                                
                                if (power.sus[[2]][1] == "gamma") {
                                    priordistpowersus[1]  <-  1
                                } else if (power.sus[[2]][1] == "half normal") {
                                    priordistpowersus[1]  <-  2
                                } else if (power.sus[[2]][1] == "uniform") {
                                    priordistpowersus[1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.sus[[1]])) {
                                
                                if ((dim(power.sus[[1]])[1] != nsuspar) & (dim(power.sus[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                                }
                                
                                #                                powersus        <-  power.sus[[1]]
                                initial[[6]]    <-  power.sus[[1]]
                                
                            } else if (is.vector(power.sus[[1]]) & length(power.sus[[1]]) == nsuspar) {
                                #                                powersus        <-  power.sus[[1]]
                                initial[[6]]    <-  matrix(rep(power.sus[[1]],nchains), ncol = nchains, nrow = nsuspar)
                            }
                            
                        }
                        
                    }
                } else if (len == 2L) {
                    anum66               <-  2
                    initial[[6]]         <-  matrix(rep(rep(1, nsuspar),nchains), ncol = nchains, nrow = nsuspar)
                    powersusproposalvar  <-  rep(0, nsuspar)
                    priordistpowersus    <-  rep(1, nsuspar)
                    priorpar1powersus    <-  rep(0, nsuspar)
                    priorpar2powersus    <-  rep(0, nsuspar)
                }
            }
        }
        
        # check transmissibility terms:
        ### new ones
        
        if (is.null(control.trans)){
            
            ntranspar  <-  1
            anum33  <-  2
            #            transpar  <-  1
            initial[[3]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            transcov  <-  matrix(rep(1, n), ncol= ntranspar, nrow = n)
            transproposalvar  <-  0
            priordisttranspar  <-  1
            priorpar1trans  <-  rep(0, ntranspar)
            priorpar2trans  <-  rep(0, ntranspar)
            anum77      <-  2
            #            powertrans <-  1
            initial[[7]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            powertransproposalvar <-  0
            priordistpowertrans <-  1
            priorpar1powertrans <-  rep(0, ntranspar)
            priorpar2powertrans <-  rep(0, ntranspar)
            
        } else if (!is.list(control.trans)) {
            
            stop("The option control.trans must be a list of values of the transmissibility parameters, transmissibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.trans)) {
            
            len <- length(control.trans)
            
            if (len < 2L) {
                stop("The length of the option control.trans must be at least 2", call. = FALSE)
            } else {
                
                anum33       <-  1
                trans.par <- control.trans[[1]]
                transcov  <- control.trans[[2]]
                
                if (is.matrix(transcov)) {
                    ntranspar <- dim(transcov)[2]
                } else {
                    ntranspar <- 1
                }
                
                priordisttranspar    <-  vector(mode="integer", length = ntranspar)
                priordistpowertrans  <-  vector(mode="integer", length = ntranspar)
                
                if (any(transcov < 0)) {
                    stop("Covariate(s) values of the transmissibility function must be positive: control.trans", call.=FALSE)
                }
                
                if (!is.list(trans.par)) {
                    stop("The first list of the argument \"control.trans\" (the transmissibility parameters) needs to be entered as a list of two matrices:\n1) a number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                } else {
                    if (length(trans.par != 2)) {
                        stop("The first list of the argument \"control.trans\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                    } else {
                        if (is.matrix(trans.par[[2]])) {
                            if ((dim(trans.par[[2]])[1] != ntranspar) & (dim(trans.par[[2]])[2] != 4L)) {
                                stop("Error in entering the argument of \"control.trans\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            priorpar1trans    <-  trans.par[[2]][, 2]
                            priorpar2trans    <-  trans.par[[2]][, 3]
                            transproposalvar  <-  trans.par[[2]][, 4]
                            
                            for(i in 1:ntranspar) {
                                if (trans.par[[2]][i, 1] == "gamma") {
                                    priordisttranspar[i]  <-  1
                                } else if (trans.par[[2]][i, 1] == "half normal") {
                                    priordisttranspar[i]  <-  2
                                } else if (trans.par[[2]][i, 1] == "uniform") {
                                    priordisttranspar[i]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                }
                            }
                            
                        } else {
                            if (ntranspar > 1L){
                                stop("Error in entering the argument of control.trans, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            if (length(trans.par[[2]]) != 4L){
                                stop("Error in entering one or more of the arguments of the transmissibility parameters in the option control.trans", call. = FALSE)
                            }
                            
                            priorpar1trans    <-  trans.par[[2]][2]
                            priorpar2trans    <-  trans.par[[2]][3]
                            transproposalvar  <-  trans.par[[2]][4]
                            
                            if (trans.par[[2]][1] == "gamma") {
                                priordisttranspar[1]  <-  1
                            } else if (trans.par[[2]][1] == "half normal") {
                                priordisttranspar[1]  <-  2
                            } else if (trans.par[[2]][1] == "uniform") {
                                priordisttranspar[1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                            }
                            
                        }
                        if (is.matrix(trans.par[[1]])) {
                            
                            if ((dim(trans.par[[1]])[1] != ntranspar) & (dim(trans.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                            }
                            
                            #                            transpar        <-  trans.par[[1]]
                            initial[[3]]    <-  trans.par[[1]]
                            
                        } else if (is.vector(trans.par[[1]]) & length(trans.par[[1]]) == ntranspar) {
                            #                            transpar        <-  trans.par[[1]]
                            initial[[3]]    <-  matrix(rep(trans.par[[1]],nchains), ncol = nchains, nrow = ntranspar)
                        }
                        
                    }
                }
                
                if (len == 3L) {
                    power.trans <- control.trans[[3]]
                    
                    if (!is.list(power.trans)) {
                        stop("The third list of the argument \"control.trans\" (power parameters of transmissibility covariates) needs to be entered as a list of two matrices:\n1) A number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                    } else {
                        if (length(power.trans != 2)) {
                            stop("The third list of the argument \"control.trans\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                        } else {
                            if (is.matrix(power.trans[[2]])) {
                                if ((dim(power.trans[[2]])[1] != ntranspar) & (dim(power.trans[[2]])[2] != 4L)) {
                                    stop("Error in entering the argument of \"control.trans\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                
                                priorpar1powertrans    <-  power.trans[[2]][, 2]
                                priorpar2powertrans    <-  power.trans[[2]][, 3]
                                powertransproposalvar  <-  power.trans[[2]][, 4]
                                
                                
                                for(i in 1:ntranspar) {
                                    if (power.trans[[2]][i, 1] == "gamma") {
                                        priordistpowertrans[i]  <-  1
                                    } else if (power.trans[[2]][i, 1] == "half normal") {
                                        priordistpowertrans[i]  <-  2
                                    } else if (power.trans[[2]][i, 1] == "uniform") {
                                        priordistpowertrans[i]  <-  3
                                    } else {
                                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                    }
                                }
                                
                            } else {
                                if (ntranspar > 1L){
                                    stop("Error in entering the argument of control.trans, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                if (length(power.trans[[2]]) != 4L){
                                    stop("Error in entering one or more of the arguments of the transmissibility parameters in the option control.trans", call. = FALSE)
                                }
                                
                                priorpar1powertrans     <-  power.trans[[2]][, 2]
                                priorpar2powertrans     <-  power.trans[[2]][, 3]
                                powertransproposalvar   <-  power.trans[[2]][, 4]
                                
                                
                                if (power.trans[[2]][1] == "gamma") {
                                    priordistpowertrans[1]  <-  1
                                } else if (power.trans[[2]][1] == "half normal") {
                                    priordistpowertrans[1]  <-  2
                                } else if (power.trans[[2]][1] == "uniform") {
                                    priordistpowertrans[1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.trans[[1]])) {
                                
                                if ((dim(power.trans[[1]])[1] != ntranspar) & (dim(power.trans[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                                }
                                
                                #                                powertrans      <-  power.trans[[1]]
                                initial[[7]]    <-  power.trans[[1]]
                                
                            } else if (is.vector(power.trans[[1]]) & length(power.trans[[1]]) == ntranspar) {
                                #                                powertrans      <-  power.trans[[1]]
                                initial[[7]]    <-  matrix(rep(power.trans[[1]],nchains), ncol = nchains, nrow = ntranspar)
                            }
                            
                        }
                        
                    }
                } else if (len == 2L) {
                    anum77                 <-  2
                    initial[[7]]           <-  matrix(rep(rep(1,ntranspar),nchains), ncol = nchains, nrow = ntranspar)
                    #                    powertrans             <-  rep(1, ntranspar)
                    powertransproposalvar  <-  rep(0, ntranspar)
                    priordistpowertrans    <-  rep(1, ntranspar)
                    priorpar1powertrans    <-  rep(0, ntranspar)
                    priorpar2powertrans    <-  rep(0, ntranspar)
                }
            }
            
            
        }
        
        initial[[8]] <- temp1
        
        anum2  <-  c(anum11, anum22, anum33, anum44, anum55, anum66, anum77)
        
        cat("************************************************","\n")
        cat("* Start performing MCMC for the ", datatype," SIR ILM for","\n")
        #        cat(nsim, "iterations", "\n")
        cat("************************************************","\n")
        
        
        n=as.integer(n);
        nsim=as.integer(nsim);
        ni=as.integer(ni);
        num=as.integer(num);
        anum2=as.vector(anum2, mode="integer");
        temp = as.integer(temp1);
        nsuspar=as.integer(nsuspar);
        ntranspar=as.integer(ntranspar);
        net=matrix(as.double(net), ncol=n, nrow=n);
        dis=matrix(as.double(dis), ncol=n, nrow=n);
        epidat=matrix(as.double(object$epidat), ncol=4, nrow=n);
        blockupdate=as.vector(blockupdate, mode="integer");
        priordistsuspar=as.vector(priordistsuspar, mode="integer");
        priordisttranspar=as.vector(priordisttranspar, mode="integer");
        priordistkernelparpar=as.vector(priordistkernelparpar, mode="integer");
        priordistsparkpar=as.integer(priordistsparkpar);
        priordistpowersus=as.vector(priordistpowersus, mode="integer");
        priordistpowertrans=as.vector(priordistpowertrans, mode="integer");
        suspar=as.vector(initial[[2]][,1], mode="double");
        suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n);
        powersus=as.vector(initial[[6]][,1], mode="double");
        transpar=as.vector(initial[[3]][,1], mode="double");
        transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n);
        powertrans=as.vector(initial[[7]][,1], mode="double");
        kernelpar=as.vector(initial[[5]][,1], mode="double");
        spark=initial[[4]][1];
        delta1=as.vector(initial[[1]][,1], mode = "double");
        kernelparproposalvar=as.vector(kernelparproposalvar, mode="double");
        sparkproposalvar=as.double(sparkproposalvar);
        susproposalvar=as.vector(susproposalvar, mode="double");
        powersusproposalvar=as.vector(powersusproposalvar, mode="double");
        transproposalvar=as.vector(transproposalvar, mode="double");
        powertransproposalvar=as.vector(powertransproposalvar, mode="double");
        infperiodproposal=as.vector(infperiodproposal, mode="double");
        priorpar1sus=as.vector(priorpar1sus, mode="double");
        priorpar2sus=as.vector(priorpar2sus, mode="double");
        priorpar1powersus=as.vector(priorpar1powersus, mode="double");
        priorpar2powersus=as.vector(priorpar2powersus, mode="double");
        priorpar1trans=as.vector(priorpar1trans, mode="double");
        priorpar2trans=as.vector(priorpar2trans, mode="double");
        priorpar1powertrans=as.vector(priorpar1powertrans, mode="double");
        priorpar2powertrans=as.vector(priorpar2powertrans, mode="double");
        kernelparprior=matrix(as.double(kernelparprior), ncol=2, nrow=2);
        sparkprior=as.vector(sparkprior, mode="double");
        delta2prior=as.vector(delta2prior, mode="double");
        susparop= matrix(as.double(0), ncol=nsuspar, nrow=nsim);
        powersusparop=matrix(as.double(0), ncol=nsuspar, nrow=nsim);
        transparop= matrix(as.double(0), ncol=ntranspar, nrow=nsim);
        powertransparop=matrix(as.double(0), ncol=ntranspar, nrow=nsim);
        kernelparop=matrix(as.double(0), ncol=2, nrow=nsim);
        sparkop=matrix(as.double(0), ncol=1, nrow=nsim);
        delta1op=matrix(as.double(0), ncol=1, nrow=nsim);
        delta2op=matrix(as.double(0), ncol=1, nrow=nsim);
        epidatmcper=matrix(as.double(0), ncol=n, nrow=nsim);
        epidatmctim=matrix(as.double(0), ncol=n, nrow=nsim);
        epidatmcrem=matrix(as.double(0), ncol=n, nrow=nsim);
        loglik=matrix(as.double(0), ncol=1, nrow=nsim)
        
        sirmcmc <- list(n,nsim,ni,num,anum2,temp,nsuspar,ntranspar,net,dis,epidat,blockupdate,priordistsuspar,
        priordisttranspar, priordistkernelparpar, priordistsparkpar, priordistpowersus,
        priordistpowertrans,suspar,suscov,powersus,transpar,transcov,powertrans,kernelpar,spark,delta1,
        kernelparproposalvar, sparkproposalvar, susproposalvar, powersusproposalvar, transproposalvar,
        powertransproposalvar, infperiodproposal, priorpar1sus, priorpar2sus, priorpar1powersus,
        priorpar2powersus, priorpar1trans, priorpar2trans, priorpar1powertrans, priorpar2powertrans,
        kernelparprior, sparkprior, delta2prior, susparop, powersusparop, transparop, powertransparop,
        kernelparop, sparkop, delta1op, delta2op, epidatmcper, epidatmctim, epidatmcrem, loglik)

        parallel.function <- function(i) {
            .Fortran("mcmcsir2",
            n=  sirmcmc[[1]], nsim = sirmcmc[[2]],
            ni= sirmcmc[[3]],
            num= sirmcmc[[4]], anum2= sirmcmc[[5]],
            temp = as.integer(initial[[8]][i]),
            nsuspar= sirmcmc[[7]], ntranspar= sirmcmc[[8]], net= sirmcmc[[9]],
            dis= sirmcmc[[10]], epidat= sirmcmc[[11]], blockupdate= sirmcmc[[12]],
            priordistsuspar= sirmcmc[[13]], priordisttranspar= sirmcmc[[14]],
            priordistkernelparpar= sirmcmc[[15]], priordistsparkpar= sirmcmc[[16]],
            priordistpowersus= sirmcmc[[17]], priordistpowertrans= sirmcmc[[18]],
            suspar= as.vector(initial[[2]][,i], mode="double"), suscov= sirmcmc[[20]], 
            powersus= as.vector(initial[[6]][,i], mode="double"),
            transpar= as.vector(initial[[3]][,i], mode="double"), 
            transcov= sirmcmc[[23]], powertrans= as.vector(initial[[7]][,i], mode="double"),
            kernelpar= as.vector(initial[[5]][,i], mode="double"), spark= initial[[4]][i],
            delta1=as.vector(initial[[1]][,i], mode = "double"),
            kernelparproposalvar= sirmcmc[[28]],
            sparkproposalvar= sirmcmc[[29]],
            susproposalvar= sirmcmc[[30]],
            powersusproposalvar= sirmcmc[[31]],
            transproposalvar= sirmcmc[[32]],
            powertransproposalvar= sirmcmc[[33]],
            infperiodproposal= sirmcmc[[34]],
            priorpar1sus= sirmcmc[[35]],
            priorpar2sus= sirmcmc[[36]],
            priorpar1powersus= sirmcmc[[37]],
            priorpar2powersus= sirmcmc[[38]],
            priorpar1trans= sirmcmc[[39]],
            priorpar2trans= sirmcmc[[40]],
            priorpar1powertrans= sirmcmc[[41]],
            priorpar2powertrans= sirmcmc[[42]],
            kernelparprior= sirmcmc[[43]],
            sparkprior= sirmcmc[[44]],
            delta2prior= sirmcmc[[45]],
            susparop=  sirmcmc[[46]],
            powersusparop=  sirmcmc[[47]],
            transparop=  sirmcmc[[48]],
            powertransparop=  sirmcmc[[49]],
            kernelparop=  sirmcmc[[50]],
            sparkop=  sirmcmc[[51]],
            delta1op=  sirmcmc[[52]],
            delta2op=  sirmcmc[[53]],
            epidatmcper=  sirmcmc[[54]],
            epidatmctim=  sirmcmc[[55]],
            epidatmcrem=  sirmcmc[[56]],
            loglik=  sirmcmc[[57]], NAOK = TRUE
            )
        }

        
        
        
        if (nchains > 1L) {
            
            if (parallel == TRUE) {
                no_cores <- min(nchains,detectCores())
                cl <- makeCluster(no_cores)
                varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
                clusterExport(cl, varlist=varlist, envir=environment())
                #                wd <- getwd()
                #                clusterExport(cl, varlist="wd", envir=environment())
                datmcmc22 <- parLapply(cl, seq_len(no_cores), parallel.function)
                stopCluster(cl)
            } else {
                datmcmc22 <- list(NULL)
                for (i in seq_len(nchains)) {
                    datmcmc22[[i]] <- parallel.function(i)
                }
            }
            
        } else if (nchains == 1L) {
            datmcmc22 <- parallel.function(1)
        }

        names <- c(paste("Alpha_s[",seq_len(nsuspar),"]", sep = ""))
        namet <- c(paste("Alpha_t[",seq_len(ntranspar),"]", sep = ""))
        namepowers <- c(paste("Psi_s[",seq_len(nsuspar),"]", sep = ""))
        namepowert <- c(paste("Psi_t[",seq_len(ntranspar),"]", sep = ""))
        
        if (nchains > 1L) {
            
            result77   <-  list(NULL)
            
            for (i in seq_len(nchains)){
                
                result777 <- NULL
                namecols   <-  NULL
                
                if (anum2[2]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$susparop    )
                    namecols <- c(namecols, names[1:nsuspar])
                }
                
                if (anum2[6]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$powersusparop)
                    namecols <- c(namecols, namepowers[1:nsuspar])
                }
                
                if (anum2[3]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$transparop)
                    namecols <- c(namecols, namet[1:ntranspar])
                }
                
                if (anum2[7]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$powertransparop)
                    namecols <- c(namecols, namepowert[1:ntranspar])
                }
                
                if (anum2[4]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$sparkop[, 1])
                    namecols <- c(namecols, "Spark")
                }
                
                if (anum2[5]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$kernelparop[, 1])
                    namecols <- c(namecols, "Spatial parameter")
                } else if (anum2[5]==2) {
                    result777 <- cbind(result777, datmcmc22[[i]]$kernelparop)
                    namecols <- c(namecols, c("Spatial parameter", "Network parameter"))
                }
                
                if (anum2[1]==1) {
                    result777 <- cbind(result777, datmcmc22[[i]]$delta2op[, 1])
                    namecols <- c(namecols, "Infectious period rate")
                }
                
                result777 <- cbind(result777, datmcmc22[[i]]$loglik)
                namecols <- c(namecols, "Log-likelihood")
                
                result777  <-  data.frame(result777)
                colnames(result777)  <-  namecols
                
                if (anum2[1]==2) {
                    result777 <- list(result777)
                } else if (anum2[1]==1) {
                    result777 <- list(result777, datmcmc22[[i]]$epidatmctim)
                }
                
                result77[[i]] <- result777
            }
            
        } else {
            
            result77 <- NULL
            namecols   <-  NULL
            
            if (anum2[2]==1) {
                result77 <- cbind(result77, datmcmc22$susparop)
                namecols <- c(namecols, names[1:nsuspar])
            }
            
            if (anum2[6]==1) {
                result77 <- cbind(result77, datmcmc22$powersusparop)
                namecols <- c(namecols, namepowers[1:nsuspar])
            }
            
            if (anum2[3]==1) {
                result77 <- cbind(result77, datmcmc22$transparop)
                namecols <- c(namecols, namet[1:ntranspar])
            }
            
            if (anum2[7]==1) {
                result77 <- cbind(result77, datmcmc22$powertransparop)
                namecols <- c(namecols, namepowert[1:ntranspar])
            }
            
            if (anum2[4]==1) {
                result77 <- cbind(result77, datmcmc22$sparkop[, 1])
                namecols <- c(namecols, "Spark")
            }
            
            if (anum2[5]==1) {
                result77 <- cbind(result77, datmcmc22$kernelparop[, 1])
                namecols <- c(namecols, "Spatial parameter")
            } else if (anum2[5]==2) {
                result77 <- cbind(result77, datmcmc22$kernelparop)
                namecols <- c(namecols, c("Spatial parameter", "Network parameter"))
            }
            
            if (anum2[1]==1) {
                result77 <- cbind(result77, datmcmc22$delta2op[, 1])
                namecols <- c(namecols, "Infectious period rate")
            }
            
            result77 <- cbind(result77, datmcmc22$loglik)
            namecols <- c(namecols, "Log-likelihood")
            
            result77  <-  data.frame(result77)
            colnames(result77)  <-  namecols
            
            if (anum2[1]==2) {
                result77 <- list(result77)
            } else if (anum2[1]==1) {
                result77 <- list(result77, datmcmc22$epidatmctim)
            }
            
            
        }
        
        
        
        # Creating an epictmcmc object:
        
        if (nchains == 1){
            
            if (length(result77) == 1) {
                dim.results <- dim(result77[[1]])
                mcmcsamp <- as.mcmc(result77[[1]][,-dim.results[2]])
                accpt <- 1-rejectionRate(as.mcmc(mcmcsamp))
                num.iter <- niter(as.mcmc(mcmcsamp))
                num.par <- nvar(as.mcmc(mcmcsamp))
                log.likelihood <- as.mcmc(result77[[1]][,dim.results[2]])
                
                out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                acceptance.rate = accpt, number.iteration = num.iter,
                number.parameter = num.par, number.chains = nchains)
                
            } else {
                dim.results <- dim(result77[[1]])
                mcmcsamp <- as.mcmc(result77[[1]][,-dim.results[2]])
                accpt <- 1-rejectionRate(as.mcmc(mcmcsamp))
                num.iter <- niter(as.mcmc(mcmcsamp))
                num.par <- nvar(as.mcmc(mcmcsamp))
                log.likelihood <- as.mcmc(result77[[1]][,dim.results[2]])
                num.inf <- sum(object$epidat[,2]!=Inf)
                infection.times.samples <- as.mcmc(result77[[2]][,1:num.inf])
                Average.infectious.periods <- as.mcmc(apply(object$epidat[1:num.inf,2]-infection.times.samples,1,mean))
                
                out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                acceptance.rate = accpt, number.iteration = num.iter,
                number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                Average.infectious.periods = Average.infectious.periods)
                
            }
            
        } else {
            
            if (length(result77[[1]]) == 1) {
                dim.results <- dim(result77[[1]][[1]])
                mcmcsamp <- list(NULL)
                log.likelihood <- list(NULL)
                accpt <- list(NULL)
                
                for(i in seq_len(nchains)){
                    mcmcsamp[[i]] <- as.mcmc(result77[[i]][[1]][,-dim.results[2]])
                    accpt[[i]] <- 1-rejectionRate(as.mcmc(mcmcsamp[[i]]))
                    log.likelihood[[i]] <- as.mcmc(result77[[i]][[1]][,dim.results[2]])
                }
                
                num.iter <- niter(as.mcmc(mcmcsamp[[1]]))
                num.par  <- nvar(as.mcmc(mcmcsamp[[1]]))
                
                out <- list(compart.framework = object$type, kernel.type = object$kerneltype,
                data.assumption = datatype, parameter.samples = mcmcsamp,
                log.likelihood = log.likelihood,
                acceptance.rate = accpt, number.iteration = num.iter,
                number.parameter = num.par, number.chains = nchains)
                
            } else {
                
                dim.results <- dim(result77[[1]][[1]])
                mcmcsamp <- list(NULL)
                log.likelihood <- list(NULL)
                accpt <- list(NULL)
                num.inf <- sum(object$epidat[,2]!=Inf)
                infection.times.samples <- list(NULL)
                Average.infectious.periods <- list(NULL)
                
                for(i in seq_len(nchains)){
                    mcmcsamp[[i]] <- as.mcmc(result77[[i]][[1]][,-dim.results[2]])
                    accpt[[i]] <- 1-rejectionRate(as.mcmc(mcmcsamp[[i]]))
                    log.likelihood[[i]] <- as.mcmc(result77[[i]][[1]][,dim.results[2]])
                    infection.times.samples[[i]] <- as.mcmc(result77[[i]][[2]][,seq_len(num.inf)])
                    Average.infectious.periods[[i]] <- as.mcmc(apply(object$epidat[seq_len(num.inf),2]-infection.times.samples[[i]],1,mean))
                }
                
                num.iter <- niter(as.mcmc(mcmcsamp[[1]]))
                num.par  <- nvar(as.mcmc(mcmcsamp[[1]]))
                
                out <- list(compart.framework = object$type, kernel.type = object$kerneltype,
                data.assumption = datatype, parameter.samples = mcmcsamp,
                log.likelihood = log.likelihood,
                acceptance.rate = accpt, number.iteration = num.iter,
                number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                Average.infectious.periods = Average.infectious.periods)
                
            }
            
        }

    } else if (object$type == "SINR") {
        
        if (is.null(nchains)) {
            nchains <- 1
        } else {
            nchains <- nchains
        }
        
        # Set seed value for Fortran random number generator
        if (is.null(seedval)){
            if (nchains == 1) {
                temp1 <- 0
            } else {
                temp1 <- sample(seq(100001, 999999), size = nchains, replace = FALSE)
            }
        } else {
            if (nchains == length(seedval)) {
                temp1 <- seedval
            } else {
                stop("The number of seeds is not equal to the number of chains", call. = FALSE)
            }
        }
        
        initial <- list(NULL)
        infperiodproposalin  <-  vector(mode="double", length = 2)
        infperiodproposalnr  <-  vector(mode="double", length = 2)
        deltain2prior  <-  vector(mode="double", length = 2)
        deltanr2prior  <-  vector(mode="double", length = 2)
        initial[[6]] <- list(NULL)
        
        if (datatype == "known removal") {
            
            anum66  <-  1
            
            if (is.null(delta)) {
                stop("Specify the arguments of the parameters of the incubation period distribution: delta",  call. =FALSE)
            } else {
            
                if (!is.list(delta)) {
				   stop("The argument \"delta\" must be a list of three:\n1) a scalar value of fixed shape parameter of the incubation period density.\n2) a vector of initial values of the rate parameter of the incubation period density with size equal to \"nchains\". \n3) a vector of the parameter values of the gamma prior distribution for the rate parameter.",  call.= FALSE)
                }

                if (length(delta) != 3) {
                    stop("Error in entering the arguments of the delta parameters of the incubation period distribution: delta", call.=FALSE)
                }

                if ( length(delta[[1]])>1) {
                    stop("Error in entering the arguments of the fixed shape parameter of the incubation period density: delta", call.=FALSE)
                }
    
                initial[[6]][[2]]     <- matrix(0, ncol = nchains, nrow = 2)
	            deltanr2prior         <- c(0, 0)
                initial[[6]][[1]]     <- matrix(0, ncol = nchains, nrow = 2)
                initial[[6]][[1]][1,] <- rep(delta[[1]], nchains)
				
				if (is.vector(delta[[2]])) {
					if (length(delta[[2]]) == nchains) {
						initial[[6]][[1]][2,] <-  delta[[2]]
					} else if (length(delta[[2]]) == 1) {
						initial[[6]][[1]] <-  rep(delta[[2]],nchains)
					} else {
						stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
					}
				} else {
					stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
				}
												
				if (is.vector(delta[[3]]) & length(delta[[3]]) == 2) {
					deltain2prior  <-  delta[[3]]
				} else {
					stop("Error in entering the parameter values of the gamma prior distribution of the incubation rate parameter: delta[[3]]",  call.= FALSE)
				}
            }            

            if (is.null(periodproposal) ) {
                infperiodproposalin  <-  c(0, 0)
                infperiodproposalnr  <-  c(0, 0)
            } else {
                if (!is.matrix(periodproposal)) {
                    periodproposal          <-  matrix(periodproposal, ncol=2, nrow=1)
                    infperiodproposalin  <-  periodproposal[1, ]
                    infperiodproposalnr  <-  c(0, 0)
                    
                } else if (dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) {
                    stop("The parameters of the proposal distribution should be entered as a 1 by 2 matrix or as a vector: periodproposal", call. = FALSE)
                } else {
                    infperiodproposalin  <-  periodproposal[1, ]
                    infperiodproposalnr  <-  c(0, 0)
                }
            }
            
            if (is.null(blockupdate) ) {
                blockupdate  <-  c(1, 1)
            }

        } else if (datatype == "unknown removal") {
            
            anum66  <-  2
            
            if (is.null(delta)) {
                stop("Specify the arguments of the parameters of the incubation and delay period distributions: delta",  call. =FALSE)
            } else {
            
                if (!is.list(delta)) {
				   stop("The argument \"delta\" must be a list of three:\n1) a vector of the fixed shape parameters of the incubation and delay period densities.\n2) a (2 by nchains) matrix of initial values of the incubation and delay rate parameters. \n3) a (2 by 2) matrix of the parameter values of the gamma prior distribution for the incubation and delay rate parameters.",  call.= FALSE)
                }

                if (length(delta) != 3) {
                    stop("Error in entering the arguments of the delta parameters of the incubation and delay period distributions: delta", call.=FALSE)
                }

                if (length(delta[[1]]) != 2) {
                    stop("Error in entering the arguments of the fixed shape parameters of the incubation and delay period densities: delta", call.=FALSE)
                }
    
                initial[[6]][[2]]     <- matrix(0, ncol = nchains, nrow = 2)
                initial[[6]][[2]][1,] <- rep(delta[[1]][2], nchains)
                initial[[6]][[1]]     <- matrix(0, ncol = nchains, nrow = 2)
                initial[[6]][[1]][1,] <- rep(delta[[1]][1], nchains)
				
				if (is.vector(delta[[2]]) & length(delta[[2]])==2) {
					initial[[6]][[1]][2,] <-  rep(delta[[2]][1], nchains)
					initial[[6]][[2]][2,] <-  rep(delta[[2]][2], nchains)
				} else if(is.matrix(delta[[2]])) {
					if ((dim(delta[[2]])[1] ==2) & (dim(delta[[2]])[2] == nchains)) {
							initial[[6]][[1]][2,] <-  delta[[2]][1,]
							initial[[6]][[2]][2,] <-  delta[[2]][2,]
						} else {
							stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
						}
				} else {
					stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
				}
												
				if (is.matrix(delta[[3]])) {
					if ((dim(delta[[3]])[1] ==2) & (dim(delta[[3]])[2] == 2)) {
						deltain2prior  <-  delta[[3]][1,]
						deltanr2prior  <-  delta[[3]][2,]
					} else {
					stop("The prior parameters of the incubation and delay periods must be entered as a 2 by 2 matrix: delta[[3]]",  call.= FALSE)
					}
				} else {
					stop("The prior parameters of the incubation and delay periods must be entered as a 2 by 2 matrix: delta[[3]]",  call.= FALSE)
				}
            }
            
            
            if (is.null(periodproposal) ) {
                infperiodproposalin  <-  c(0, 0)
                infperiodproposalnr  <-  c(0, 0)
            } else {
                if (!is.matrix(periodproposal)) {
                    periodproposal          <-  matrix(periodproposal, ncol=2, nrow=2)
                    infperiodproposalin  <-  periodproposal[1, ]
                    infperiodproposalnr  <-  periodproposal[2, ]
                } else if (dim(periodproposal)[1]!=2 & dim(periodproposal)[2]!=2) {
                    stop("Enter the proposal distribution for updating the incubation and delay periods as a 2 by 2 matrix: periodproposal", call. = FALSE)
                } else {
                    infperiodproposalin  <-  periodproposal[1, ]
                    infperiodproposalnr  <-  periodproposal[2, ]
                }
            }
            
            if (is.null(blockupdate) ) {
                blockupdate  <-  c(1, 1)
            }
            
        } else if (datatype == "known epidemic") {
            
            blockupdate  <-  vector(mode="integer", length = 2)
            
            if (!is.null(delta)) {
                warning("The incubation and infectious period rates are not updated as the option of datatype is \"known epidemic\".", call. = TRUE)
            }
            
            anum66                   <-  3
            infperiodproposalin      <-  c(0, 0)
            infperiodproposalnr      <-  c(0, 0)
            deltain2prior            <-  c(0, 0)
            deltanr2prior            <-  c(0, 0)
            blockupdate              <-  c(1, 1)
            initial[[6]][[1]]  <-  matrix(rep(0.0,nchains), ncol = nchains, nrow = 2)
            initial[[6]][[2]]  <-  matrix(rep(0.0,nchains), ncol = nchains, nrow = 2)
            
        } else {
            stop("Specify the data type as \"known removal\",  \"unknown removal\" or \"known epidemic\": datatype ",  call. = FALSE)
        }
        
        
        kernelpar <-  vector(mode="double", length = 2)
        kernelparproposalvar <-  vector(mode="double", length = 2)
        priordistkernelparpar <-  vector(mode="integer", length = 2)
        kernelparprior <-  matrix(0, ncol = 2, nrow = 2)
        sparkprior <-  vector(mode="double", length = 2)
        
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
            
            anum44  <-  1
            
            if (is.null(kernel.par)) {
                stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
            }
            
            if (!is.list(kernel.par)) {
                stop("The argument \"kernel.par\" must be a list of two:\n1) a vector of initial values of the kernel parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            }
            
            if (length(kernel.par) != 2L) {
                stop("The argument \"kernel.par\" must be a list of two:\n1) a vector of initial values of the kernel parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            } else {
                if (length(kernel.par[[2]]) != 4L) {
                    stop("Error in entering the second list of the argument \"kernel.par\" for updating the kernel parameter.",  call.= FALSE)
                }
                kernelparproposalvar[1] <-  kernel.par[[2]][4]
                kernelparproposalvar[2] <-  0
                kernelparprior[1, ] <-  kernel.par[[2]][2:3]
                kernelparprior[2, ] <-  c(0, 0)
                
                if (kernel.par[[2]][1] == "gamma") {
                    priordistkernelparpar[1]  <-  1
                } else if (kernel.par[[2]][1] == "half normal") {
                    priordistkernelparpar[1]  <-  2
                } else if (kernel.par[[2]][1] == "uniform") {
                    priordistkernelparpar[1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
                }
                
                priordistkernelparpar[2]  <-  1
                
                if (length(kernel.par[[1]]) == nchains){
                    initial[[4]] <-  matrix(c(kernel.par[[1]],rep(0,nchains)), ncol = nchains, byrow = TRUE)
                } else if(length(kernel.par[[1]]) == 1) {
                    initial[[4]] <-  matrix(rep(c(kernel.par[[1]],0),nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("Error in entering the initial values of the kernel parameter: kernel.par[[1]]",  call.= FALSE)
                }
                
            }
            
            if (is.null(spark.par)) {
                
                initial[[3]] <-  rep(0,nchains)
                anum33  <-  2
                sparkproposalvar  <-  0
                priordistsparkpar  <-  1
                sparkprior  <-  c(1, 1)
                
            } else {
                
                anum44 <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[3]] <-  spark.par[[1]]
                    } else if(length(kernel.par[[1]]) == 1) {
                        initial[[3]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                    
                }
                
            }
            
        } else if (object$kerneltype == "network") {
            
            n   <-  dim(object$network)[1]
            dis  <- matrix(0, ncol = n, nrow = n)
            net  <- object$network
            
            num <-  1
            priordistkernelparpar  <-  c(1, 1)
            kernelparproposalvar  <-  c(0, 0)
            kernelpar  <-  c(0, 0)
            kernelparprior <-  matrix(0, ncol = 2, nrow = 2)
            initial[[4]] <-  matrix(rep(c(0,0),nchains), ncol = nchains, nrow = 2)
            
            anum44 <-  3
            
            if (anum66 == 1 | anum66 == 2) {
                
                anum33   <-  1
                
                if (is.null(spark.par)) {
                    stop("Specify the arguments for updating the spark parameter: spark.par ",  call. =FALSE)
                }
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[3]] <-  spark.par[[1]]
                    } else if(length(kernel.par[[1]]) == 1) {
                        initial[[3]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                }
                
            } else {
                
                
                if (is.null(spark.par)) {
                    
                    anum33  <-  2
                    sparkproposalvar <-  0
                    priordistsparkpar <-  1
                    sparkprior <-  c(1, 1)
                    initial[[3]] <-  rep(0,nchains)
                    
                } else {
                    anum33 <-  1
                    
                    if (!is.list(spark.par)) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    }
                    
                    if (length(spark.par) != 2L) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    } else {
                        if (length(spark.par[[2]]) != 4L) {
                            stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                        }
                        sparkproposalvar <-  spark.par[[2]][4]
                        sparkprior <-  spark.par[[2]][2:3]
                        
                        if (spark.par[[2]][1] == "gamma") {
                            priordistsparkpar  <-  1
                        } else if (spark.par[[2]][1] == "half normal") {
                            priordistsparkpar  <-  2
                        } else if (spark.par[[2]][1] == "uniform") {
                            priordistsparkpar  <-  3
                        } else {
                            stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                        }
                        
                        if (length(spark.par[[1]]) == nchains){
                            initial[[3]] <-  spark.par[[1]]
                        } else if(length(kernel.par[[1]]) == 1) {
                            initial[[3]] <-  rep(spark.par[[1]],nchains)
                        } else {
                            stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                        }
                        
                    }
                    
                }
            }
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
            
            anum44  <-  2
            
            if (is.null(kernel.par)) {
                stop("Specify the arguments for updating the kernel parameter: kernel.par",  call.= FALSE)
            }
            
            if (!is.list(kernel.par)) {
                stop("The argument \"kernel.par\" must be a list of two matrices:\n1) a matrix of initial values of the two kernel parameters of number of columns equal to \"nchains\". \n2) a matrix of two rows contains prior distribution, prior parameter values, and proposal variance for each parameters.",  call.= FALSE)
            }
            
            if (length(kernel.par) != 2L) {
                stop("The argument \"kernel.par\" must be a list of two matrices:\n1) a matrix of initial values of the two kernel parameters of number of columns equal to \"nchains\". \n2) a matrix of two rows contains prior distribution, prior parameter values, and proposal variance for each parameters.",  call.= FALSE)
            } else {
                if (!is.matrix(kernel.par[[2]])) {
                    stop("The second list of the argument \"kernel.par\" must be a matrix of 2 by 4.",  call.= FALSE)
                } else if (dim(kernel.par[[2]])[1] != 2 & dim(kernel.par[[2]])[2] != 4) {
                    stop("The second list of the argument \"kernel.par\" must be a matrix of 2 by 4.",  call.= FALSE)
                }
                kernelparproposalvar[1] <-  kernel.par[[2]][1,4]
                kernelparproposalvar[2] <-  kernel.par[[2]][2,4]
                kernelparprior[1, ] <-  kernel.par[[2]][1,2:3]
                kernelparprior[2, ] <-  kernel.par[[2]][2,2:3]
                
                if (kernel.par[[2]][1,1] == "gamma") {
                    priordistkernelparpar[1]  <-  1
                } else if (kernel.par[[2]][1,1] == "half normal") {
                    priordistkernelparpar[1]  <-  2
                } else if (kernel.par[[2]][1,1] == "uniform") {
                    priordistkernelparpar[1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                
                if (kernel.par[[2]][2,1] == "gamma") {
                    priordistkernelparpar[2]  <-  1
                } else if (kernel.par[[2]][2,1] == "half normal") {
                    priordistkernelparpar[2]  <-  2
                } else if (kernel.par[[2]][2,1] == "uniform") {
                    priordistkernelparpar[2]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\", \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                
                if (is.matrix(kernel.par[[1]])) {
                    if (dim(kernel.par[[1]])[1] != 2 & dim(kernel.par[[1]])[2] != nchains) {
                        stop("The matrix of the initial values of the kernel parameters must be a 2 by \"nchains\".", call. = FALSE)
                    }
                    initial[[4]] <-  kernel.par[[1]]
                } else if (is.vector(kernel.par[[1]]) & length(kernel.par[[1]]) == 2) {
                    initial[[4]] <-  matrix(rep(kernel.par[[1]],nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("The initial values of the kernel parameters must be a matrix of 2 by \"nchains\" or a vector of two initial values.", call. = FALSE)
                }
                
            }
            
            if (is.null(spark.par)) {
                
                anum33  <-  2
                sparkproposalvar <-  0
                priordistsparkpar <-  1
                sparkprior <-  c(1, 1)
                initial[[3]] <-  rep(0,nchains)
                
            } else {
                anum33 <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    sparkproposalvar <-  spark.par[[2]][4]
                    sparkprior <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        priordistsparkpar  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        priordistsparkpar  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        priordistsparkpar  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        initial[[3]] <-  spark.par[[1]]
                    } else if(length(kernel.par[[1]]) == 1) {
                        initial[[3]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                    
                }
                
            }
        }
        
        # check the number of infected individuals:
        
        if (anum66 == 1 | anum66 == 3) {
            ni  <-  sum(object$epidat[, 2]!=Inf)
        } else if (anum66 == 2) {
            ni  <-  sum(object$epidat[, 4]!=Inf)
        } else {
            stop("Error: the epidemic data must be in the same format as in datagen function where the removal, notified and infection times of uninfected individuals are \"Inf\".", call. = FALSE)
        }
        
        
        # check Susceptibility terms:
        ### new ones
        
        if (is.null(control.sus)){
            
            nsuspar  <-  1
            anum11  <-  2
            initial[[1]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            suscov  <-  matrix(rep(1, n), ncol= nsuspar, nrow = n)
            susproposalvar  <-  0
            priordistsuspar  <-  1
            priorpar1sus  <-  rep(0, nsuspar)
            priorpar2sus  <-  rep(0, nsuspar)
            anum77      <-  2
            initial[[7]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            powersusproposalvar <-  0
            priordistpowersus <-  1
            priorpar1powersus <-  rep(0, nsuspar)
            priorpar2powersus <-  rep(0, nsuspar)
            
        } else if (!is.list(control.sus)) {
            
            stop("The option control.sus must be a list of values of the susceptibility parameters, susceptibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.sus)) {
            
            len <- length(control.sus)
            
            if (len < 2L) {
                stop("The length of the option control.sus must be at least 2", call. = FALSE)
            } else {
                
                anum11       <-  1
                sus.par <- control.sus[[1]]
                suscov  <- control.sus[[2]]
                
                if (is.matrix(suscov)) {
                    nsuspar <- dim(suscov)[2]
                } else {
                    nsuspar <- 1
                }
                
                priordistsuspar    <-  vector(mode="integer", length = nsuspar)
                priordistpowersus  <-  vector(mode="integer", length = nsuspar)
                
                if (any(suscov < 0)) {
                    stop("Covariate(s) values of the susceptibility function must be positive: control.sus", call.=FALSE)
                }
                
                if (!is.list(sus.par)) {
                    stop("The first list of the argument \"control.sus\" (the susceptibility parameters) needs to be entered as a list of two matrices:\n1) a number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                } else {
                    if (length(sus.par) != 2) {
                        stop("The first list of the argument \"control.sus\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                    } else {
                        if (is.matrix(sus.par[[2]])) {
                            if ((dim(sus.par[[2]])[1] != nsuspar) & (dim(sus.par[[2]])[2] != 4L)) {
                                stop("Error in entering the argument of \"control.sus\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            priorpar1sus    <-  sus.par[[2]][, 2]
                            priorpar2sus    <-  sus.par[[2]][, 3]
                            susproposalvar  <-  sus.par[[2]][, 4]
                            
                            for(i in 1:nsuspar) {
                                if (sus.par[[2]][i, 1] == "gamma") {
                                    priordistsuspar[i]  <-  1
                                } else if (sus.par[[2]][i, 1] == "half normal") {
                                    priordistsuspar[i]  <-  2
                                } else if (sus.par[[2]][i, 1] == "uniform") {
                                    priordistsuspar[i]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                }
                            }
                            
                        } else {
                            if (nsuspar > 1L){
                                stop("Error in entering the argument of control.sus, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            if (length(sus.par[[2]]) != 4L){
                                stop("Error in entering one or more of the arguments of the susceptibility parameters in the option control.sus", call. = FALSE)
                            }
                            
                            priorpar1sus    <-  sus.par[[2]][2]
                            priorpar2sus    <-  sus.par[[2]][3]
                            susproposalvar  <-  sus.par[[2]][4]
                            
                            if (sus.par[[2]][1] == "gamma") {
                                priordistsuspar[1]  <-  1
                            } else if (sus.par[[2]][1] == "half normal") {
                                priordistsuspar[1]  <-  2
                            } else if (sus.par[[2]][1] == "uniform") {
                                priordistsuspar[1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                            }
                            
                        }
                        if (is.matrix(sus.par[[1]])) {
                            
                            if ((dim(sus.par[[1]])[1] != nsuspar) & (dim(sus.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                            }
                            
                            initial[[1]]    <-  sus.par[[1]]

                        } else if (is.vector(sus.par[[1]])) {
                            if (nsuspar == 1 & length(sus.par[[1]]) == nchains) {
                                initial[[1]]    <-  matrix(sus.par[[1]], ncol = nchains, nrow = nsuspar)
                            } else if (nsuspar > 1 & length(sus.par[[1]]) == nsuspar) {
                                initial[[1]]    <-  matrix(rep(sus.par[[1]],nchains), ncol = nchains, nrow = nsuspar)
                            } else {
                                stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                            }
                        }
                    }
                }
                
                if (len == 3L) {
                    power.sus <- control.sus[[3]]
                    anum77  <-  1
                    
                    if (!is.list(power.sus)) {
                        stop("The third list of the argument \"control.sus\" (power parameters of susceptibility covariates) needs to be entered as a list of two matrices:\n1) A number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                    } else {
                        if (length(power.sus) != 2) {
                            stop("The third list of the argument \"control.sus\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                        } else {
                            if (is.matrix(power.sus[[2]])) {
                                if ((dim(power.sus[[2]])[1] != nsuspar) & (dim(power.sus[[2]])[2] != 4L)) {
                                    stop("Error in entering the argument of \"control.sus\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                
                                priorpar1powersus    <-  power.sus[[2]][, 2]
                                priorpar2powersus    <-  power.sus[[2]][, 3]
                                powersusproposalvar  <-  power.sus[[2]][, 4]
                                
                                
                                for(i in 1:nsuspar) {
                                    if (power.sus[[2]][i, 1] == "gamma") {
                                        priordistpowersus[i]  <-  1
                                    } else if (power.sus[[2]][i, 1] == "half normal") {
                                        priordistpowersus[i]  <-  2
                                    } else if (power.sus[[2]][i, 1] == "uniform") {
                                        priordistpowersus[i]  <-  3
                                    } else {
                                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                    }
                                }
                                
                            } else {
                                if (nsuspar > 1L){
                                    stop("Error in entering the argument of control.sus, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                if (length(power.sus[[2]]) != 4L){
                                    stop("Error in entering one or more of the arguments of the susceptibility parameters in the option control.sus", call. = FALSE)
                                }
                                
                                priorpar1powersus     <-  power.sus[[2]][2]
                                priorpar2powersus     <-  power.sus[[2]][3]
                                powersusproposalvar   <-  power.sus[[2]][4]
                                
                                
                                if (power.sus[[2]][1] == "gamma") {
                                    priordistpowersus[1]  <-  1
                                } else if (power.sus[[2]][1] == "half normal") {
                                    priordistpowersus[1]  <-  2
                                } else if (power.sus[[2]][1] == "uniform") {
                                    priordistpowersus[1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.sus[[1]])) {
                                
                                if ((dim(power.sus[[1]])[1] != nsuspar) & (dim(power.sus[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the third list regarding to the power parameters of the susceptibility function in the argument of \"control.sus\".", call. = FALSE)
                                }
                                
                                initial[[7]]    <-  power.sus[[1]]
                                
                            } else if (is.vector(power.sus[[1]])) {
                                if (nsuspar == 1 & length(power.sus[[1]]) == nchains) {
                                    initial[[7]]    <-  matrix(power.sus[[1]], ncol = nchains, nrow = nsuspar)
                                } else if (nsuspar > 1 & length(power.sus[[1]]) == nsuspar) {
                                    initial[[7]]    <-  matrix(rep(power.sus[[1]],nchains), ncol = nchains, nrow = nsuspar)
                                } else {
                                    stop("Error in entering the initial values for the third list regarding to the power parameters of the susceptibility function in the argument of \"control.sus\".", call. = FALSE)
                                }
                                
                            }
                            
                        }
                        
                    }
                } else if (len == 2L) {
                    anum77               <-  2
                    initial[[7]]         <-  matrix(rep(rep(1, nsuspar),nchains), ncol = nchains, nrow = nsuspar)
                    powersusproposalvar  <-  rep(0, nsuspar)
                    priordistpowersus    <-  rep(1, nsuspar)
                    priorpar1powersus    <-  rep(0, nsuspar)
                    priorpar2powersus    <-  rep(0, nsuspar)
                }
            }
        }
        
        # check transmissibility terms:
        ### new ones
        
        if (is.null(control.trans)){
            
            ntranspar  <-  1
            anum22  <-  2
            initial[[2]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            transcov  <-  matrix(rep(1, n), ncol= ntranspar, nrow = n)
            transproposalvar  <-  0
            priordisttranspar  <-  1
            priorpar1trans  <-  rep(0, ntranspar)
            priorpar2trans  <-  rep(0, ntranspar)
            anum88      <-  2
            initial[[8]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            powertransproposalvar <-  0
            priordistpowertrans <-  1
            priorpar1powertrans <-  rep(0, ntranspar)
            priorpar2powertrans <-  rep(0, ntranspar)
            
        } else if (!is.list(control.trans)) {
            
            stop("The option control.trans must be a list of values of the transmissibility parameters, transmissibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.trans)) {
            
            len <- length(control.trans)
            
            if (len < 2L) {
                stop("The length of the option control.trans must be at least 2", call. = FALSE)
            } else {
                
                anum22       <-  1
                trans.par <- control.trans[[1]]
                transcov  <- control.trans[[2]]
                
                if (is.matrix(transcov)) {
                    ntranspar <- dim(transcov)[2]
                } else {
                    ntranspar <- 1
                }
                
                priordisttranspar    <-  vector(mode="integer", length = ntranspar)
                priordistpowertrans  <-  vector(mode="integer", length = ntranspar)
                
                if (any(transcov < 0)) {
                    stop("Covariate(s) values of the transmissibility function must be positive: control.trans", call.=FALSE)
                }
                
                if (!is.list(trans.par)) {
                    stop("The first list of the argument \"control.trans\" (the transmissibility parameters) needs to be entered as a list of two matrices:\n1) a number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                } else {
                    if (length(trans.par != 2)) {
                        stop("The first list of the argument \"control.trans\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                    } else {
                        if (is.matrix(trans.par[[2]])) {
                            if ((dim(trans.par[[2]])[1] != ntranspar) & (dim(trans.par[[2]])[2] != 4L)) {
                                stop("Error in entering the argument of \"control.trans\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            priorpar1trans    <-  trans.par[[2]][, 2]
                            priorpar2trans    <-  trans.par[[2]][, 3]
                            transproposalvar  <-  trans.par[[2]][, 4]
                            
                            for(i in 1:ntranspar) {
                                if (trans.par[[2]][i, 1] == "gamma") {
                                    priordisttranspar[i]  <-  1
                                } else if (trans.par[[2]][i, 1] == "half normal") {
                                    priordisttranspar[i]  <-  2
                                } else if (trans.par[[2]][i, 1] == "uniform") {
                                    priordisttranspar[i]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                }
                            }
                            
                        } else {
                            if (ntranspar > 1L){
                                stop("Error in entering the argument of control.trans, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                            }
                            if (length(trans.par[[2]]) != 4L){
                                stop("Error in entering one or more of the arguments of the transmissibility parameters in the option control.trans", call. = FALSE)
                            }
                            
                            priorpar1trans    <-  trans.par[[2]][2]
                            priorpar2trans    <-  trans.par[[2]][3]
                            transproposalvar  <-  trans.par[[2]][4]
                            
                            if (trans.par[[2]][1] == "gamma") {
                                priordisttranspar[1]  <-  1
                            } else if (trans.par[[2]][1] == "half normal") {
                                priordisttranspar[1]  <-  2
                            } else if (trans.par[[2]][1] == "uniform") {
                                priordisttranspar[1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                            }
                            
                        }

                        if (is.matrix(trans.par[[1]])) {
                            
                            if ((dim(trans.par[[1]])[1] != ntranspar) & (dim(trans.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                            }
                            
                            initial[[2]]    <-  trans.par[[1]]
                            
                        } else if (is.vector(trans.par[[1]])) {
                            if (ntranspar == 1 & length(trans.par[[1]]) == nchains) {
                                initial[[2]]    <-  matrix(trans.par[[1]], ncol = nchains, nrow = ntranspar)
                            } else if (ntranspar > 1 & length(trans.par[[1]]) == ntranspar) {
                                initial[[2]]    <-  matrix(rep(trans.par[[1]],nchains), ncol = nchains, nrow = ntranspar)
                            } else {
                                stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                            }
                            
                        }
                    }
                }
                
                if (len == 3L) {
                    power.trans <- control.trans[[3]]
                    anum88                 <-  1
                    
                    if (!is.list(power.trans)) {
                        stop("The third list of the argument \"control.trans\" (power parameters of transmissibility covariates) needs to be entered as a list of two matrices:\n1) A number of initial values equal to \"nchains\" for each parameter \n2) prior distribution, prior parameters, proposal variance for each parameter", call. = FALSE)
                    } else {
                        if (length(power.trans != 2)) {
                            stop("The third list of the argument \"control.trans\" needs to be entered as a list of two matrices see help(epictmcmc)", call. = FALSE)
                        } else {
                            if (is.matrix(power.trans[[2]])) {
                                if ((dim(power.trans[[2]])[1] != ntranspar) & (dim(power.trans[[2]])[2] != 4L)) {
                                    stop("Error in entering the argument of \"control.trans\", the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                
                                priorpar1powertrans    <-  power.trans[[2]][, 2]
                                priorpar2powertrans    <-  power.trans[[2]][, 3]
                                powertransproposalvar  <-  power.trans[[2]][, 4]
                                
                                
                                for(i in 1:ntranspar) {
                                    if (power.trans[[2]][i, 1] == "gamma") {
                                        priordistpowertrans[i]  <-  1
                                    } else if (power.trans[[2]][i, 1] == "half normal") {
                                        priordistpowertrans[i]  <-  2
                                    } else if (power.trans[[2]][i, 1] == "uniform") {
                                        priordistpowertrans[i]  <-  3
                                    } else {
                                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                    }
                                }
                                
                            } else {
                                if (ntranspar > 1L){
                                    stop("Error in entering the argument of control.trans, the number of parameters must be equal to the number of covariates.", call. = FALSE)
                                }
                                if (length(power.trans[[2]]) != 4L){
                                    stop("Error in entering one or more of the arguments of the transmissibility parameters in the option control.trans", call. = FALSE)
                                }
                                
                                priorpar1powertrans     <-  power.trans[[2]][2]
                                priorpar2powertrans     <-  power.trans[[2]][3]
                                powertransproposalvar   <-  power.trans[[2]][4]
                                
                                
                                if (power.trans[[2]][1] == "gamma") {
                                    priordistpowertrans[1]  <-  1
                                } else if (power.trans[[2]][1] == "half normal") {
                                    priordistpowertrans[1]  <-  2
                                } else if (power.trans[[2]][1] == "uniform") {
                                    priordistpowertrans[1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.trans[[1]])) {
                                
                                if ((dim(power.trans[[1]])[1] != ntranspar) & (dim(power.trans[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the third list regarding to the power parameters of the transmissibility function in the argument of \"control.trans\".", call. = FALSE)
                                }
                                
                                initial[[8]]    <-  power.trans[[1]]
                                
                            } else if (is.vector(power.trans[[1]])) {
                                if (ntranspar == 1 & length(power.trans[[1]]) == nchains) {
                                    initial[[8]]    <-  matrix(power.trans[[1]], ncol = nchains, nrow = ntranspar)
                                } else if (ntranspar > 1 & length(power.trans[[1]]) == ntranspar) {
                                    initial[[8]]    <-  matrix(rep(power.trans[[1]],nchains), ncol = nchains, nrow = ntranspar)
                                } else {
                                    stop("Error in entering the initial values for the third list regarding to the power parameters of the transmissibility function in the argument of \"control.trans\".", call. = FALSE)
                                }
                                
                            }
                        }
                        
                    }
                } else if (len == 2L) {
                    anum88                 <-  2
                    initial[[8]]           <-  matrix(rep(rep(1,ntranspar),nchains), ncol = nchains, nrow = ntranspar)
                    powertransproposalvar  <-  rep(0, ntranspar)
                    priordistpowertrans    <-  rep(1, ntranspar)
                    priorpar1powertrans    <-  rep(0, ntranspar)
                    priorpar2powertrans    <-  rep(0, ntranspar)
                }
            }
            
            
        }
        
        
        if (is.null(gamma.par)) {
            
            anum55              <-  2
            priordistgammapar   <-  1
            gammaproposalvar    <-  0
            gammaprior          <-  c(1, 1)
            initial[[5]]        <-  rep(1,nchains)
            
        } else {
            if (!is.list(gamma.par)) {
                stop("The argument \"gamma.par\" must be a list of two:\n1) a vector of initial values of the gamma parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            }
            
            if (length(gamma.par) != 2L) {
                stop("The argument \"gamma.par\" must be a list of two:\n1) a vector of initial values of the gamma parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
            } else {
                if (length(gamma.par[[2]]) != 4L) {
                    stop("Error in entering one or more of the arguments for updating the gamma parameter: gamma.par", call. = FALSE)
                }
                gammaproposalvar <-  gamma.par[[2]][4]
                gammaprior <-  gamma.par[[2]][2:3]
                
                if (gamma.par[[2]][1] == "gamma") {
                    priordistgammapar  <-  1
                } else if (gamma.par[[2]][1] == "half normal") {
                    priordistgammapar  <-  2
                } else if (gamma.par[[2]][1] == "uniform") {
                    priordistgammapar  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: gamma.par", call.=FALSE)
                }
                
                if (length(gamma.par[[1]]) == nchains){
                    initial[[5]] <-  gamma.par[[1]]
                } else if(length(gamma.par[[1]]) == 1) {
                    initial[[5]] <-  rep(gamma.par[[1]],nchains)
                } else {
                    stop("Error in entering the initial values of the gamma parameter: gamma.par[[1]]",  call.= FALSE)
                }
                
            }
        }
            
            initial[[9]] <- temp1
            
            anum2  <-  c(anum11, anum22, anum33, anum44, anum55, anum66, anum77, anum88)
            
            cat("************************************************","\n")
            cat("* Start performing MCMC for the ", datatype," SINR ILM for","\n")
            cat(nsim, "iterations", "\n")
            cat("************************************************","\n")
            
            
            n=as.integer(n);
            nsim=as.integer(nsim);
            ni=as.integer(ni);
            temp = as.integer(temp1);
            num=as.integer(num);
            anum2=as.vector(anum2, mode="integer");
            nsuspar=as.integer(nsuspar);
            ntranspar=as.integer(ntranspar);
            net=matrix(as.double(net), ncol=n, nrow=n);
            dis=matrix(as.double(dis), ncol=n, nrow=n);
            epidat=matrix(as.double(object$epidat), ncol=6, nrow=n);
            blockupdate=as.vector(blockupdate, mode="integer");
            priordistsuspar=as.vector(priordistsuspar, mode="integer");
            priordisttranspar=as.vector(priordisttranspar, mode="integer");
            priordistkernelparpar=as.vector(priordistkernelparpar, mode="integer");
            priordistsparkpar=as.integer(priordistsparkpar);
            priordistgammapar=as.integer(priordistgammapar);
            priordistpowersus=as.vector(priordistpowersus, mode="integer");
            priordistpowertrans=as.vector(priordistpowertrans, mode="integer");
            suspar=as.vector(initial[[1]][,1], mode="double");
            suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n);
            powersus=as.vector(initial[[7]][,1], mode="double");
            transpar=as.vector(initial[[2]][,1], mode="double");
            transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n);
            powertrans=as.vector(initial[[8]][,1], mode="double");
            kernelpar=as.vector(initial[[4]][,1], mode="double");
            spark=as.double(initial[[3]][1]); gamma=as.double(initial[[5]][1]);
            deltain=as.vector(initial[[6]][[1]][,1], mode = "double");
            deltanr=as.vector(initial[[6]][[2]][,1], mode = "double");
            kernelparproposalvar=as.vector(kernelparproposalvar, mode="double");
            sparkproposalvar=as.double(sparkproposalvar);
            gammaproposalvar=as.double(gammaproposalvar);
            susproposalvar=as.vector(susproposalvar, mode="double");
            powersusproposalvar=as.vector(powersusproposalvar, mode="double");
            transproposalvar=as.vector(transproposalvar, mode="double");
            powertransproposalvar=as.vector(powertransproposalvar, mode="double");
            infperiodproposalin=as.vector(infperiodproposalin, mode="double");
            infperiodproposalnr=as.vector(infperiodproposalnr, mode="double");
            priorpar1sus=as.vector(priorpar1sus, mode="double");
            priorpar2sus=as.vector(priorpar2sus, mode="double");
            priorpar1powersus=as.vector(priorpar1powersus, mode="double");
            priorpar2powersus=as.vector(priorpar2powersus, mode="double");
            priorpar1trans=as.vector(priorpar1trans, mode="double");
            priorpar2trans=as.vector(priorpar2trans, mode="double");
            priorpar1powertrans=as.vector(priorpar1powertrans, mode="double");
            priorpar2powertrans=as.vector(priorpar2powertrans, mode="double");
            kernelparprior=matrix(as.double(kernelparprior), ncol=2, nrow=2);
            sparkprior=as.vector(sparkprior, mode="double");
            gammaprior=as.vector(gammaprior, mode="double");
            deltain2prior=as.vector(deltain2prior, mode="double");
            deltanr2prior=as.vector(deltanr2prior, mode="double");
            susparop=matrix(0, ncol=nsuspar, nrow=nsim);
            powersusparop=matrix(0, ncol=nsuspar, nrow=nsim);
            transparop= matrix(0, ncol=ntranspar, nrow=nsim);
            powertransparop=matrix(0, ncol=ntranspar, nrow=nsim);
            kernelparop=matrix(0, ncol=2, nrow=nsim);
            sparkop=matrix(0,ncol=1, nrow=nsim);
            gammaop=matrix(0,ncol=1, nrow=nsim);
            deltain2op=matrix(0,ncol=1, nrow=nsim);
            deltanr2op=matrix(0,ncol=1, nrow=nsim);
            epidatmctim=matrix(0, ncol=n, nrow=nsim);
            epidatmcrem=matrix(0, ncol=n, nrow=nsim);
            loglik=matrix(0,ncol=1, nrow=nsim)
            
            sinrmcmc<-list(n,nsim,ni,temp,num,anum2,nsuspar,ntranspar,net,dis,
            epidat,blockupdate,priordistsuspar, priordisttranspar, priordistkernelparpar,
            priordistpowersus, priordistpowertrans, priordistsparkpar, priordistgammapar,
            suspar,suscov,powersus,transpar,transcov,powertrans,kernelpar,spark, gamma, deltain, deltanr,
            kernelparproposalvar, sparkproposalvar, gammaproposalvar, susproposalvar,
            powersusproposalvar, transproposalvar, powertransproposalvar, infperiodproposalin,
            infperiodproposalnr, priorpar1sus, priorpar2sus, priorpar1powersus,
            priorpar2powersus, priorpar1trans, priorpar2trans, priorpar1powertrans,
            priorpar2powertrans, kernelparprior, sparkprior, gammaprior, deltain2prior, deltanr2prior,
            susparop, powersusparop, transparop, powertransparop, kernelparop, sparkop, gammaop,
            deltain2op, deltanr2op, epidatmctim, epidatmcrem, loglik)

            parallel.function <- function(i) {
                .Fortran("mcmcsinr2",
                n=sinrmcmc[[1]], nsim=sinrmcmc[[2]], ni=sinrmcmc[[3]],
                temp = as.integer(initial[[9]][i]), num=sinrmcmc[[5]], anum2=sinrmcmc[[6]],
                nsuspar=sinrmcmc[[7]], ntranspar=sinrmcmc[[8]], net=sinrmcmc[[9]],
                dis=sinrmcmc[[10]], epidat=sinrmcmc[[11]], blockupdate=sinrmcmc[[12]],
                priordistsuspar=sinrmcmc[[13]],
                priordisttranspar=sinrmcmc[[14]],
                priordistkernelparpar=sinrmcmc[[15]],
                priordistpowersus=sinrmcmc[[16]],
                priordistpowertrans=sinrmcmc[[17]],
                priordistsparkpar=sinrmcmc[[18]],
                priordistgammapar=sinrmcmc[[19]],
                suspar=as.vector(initial[[1]][,i], mode="double"),
                suscov=sinrmcmc[[21]],
                powersus=as.vector(initial[[7]][,i], mode="double"),
                transpar=as.vector(initial[[2]][,i], mode="double"),
                transcov=sinrmcmc[[24]],
                powertrans=as.vector(initial[[8]][,i], mode="double"),
                kernelpar=as.vector(initial[[4]][,i], mode="double"),
                spark=as.double(initial[[3]][i]), gamma=as.double(initial[[5]][i]),
                deltain=as.vector(initial[[6]][[1]][,i], mode="double"),
                deltanr=as.vector(initial[[6]][[2]][,i], mode="double"),
                kernelparproposalvar=sinrmcmc[[31]],
                sparkproposalvar=sinrmcmc[[32]],
                gammaproposalvar=sinrmcmc[[33]],
                susproposalvar=sinrmcmc[[34]],
                powersusproposalvar=sinrmcmc[[35]],
                transproposalvar=sinrmcmc[[36]],
                powertransproposalvar=sinrmcmc[[37]],
                infperiodproposalin=sinrmcmc[[38]],
                infperiodproposalnr=sinrmcmc[[39]],
                priorpar1sus=sinrmcmc[[40]],
                priorpar2sus=sinrmcmc[[41]],
                priorpar1powersus=sinrmcmc[[42]],
                priorpar2powersus=sinrmcmc[[43]],
                priorpar1trans=sinrmcmc[[44]],
                priorpar2trans=sinrmcmc[[45]],
                priorpar1powertrans=sinrmcmc[[46]],
                priorpar2powertrans=sinrmcmc[[47]],
                kernelparprior=sinrmcmc[[48]],
                sparkprior=sinrmcmc[[49]],
                gammaprior=sinrmcmc[[50]],
                deltain2prior=sinrmcmc[[51]],
                deltanr2prior=sinrmcmc[[52]],
                susparop=sinrmcmc[[53]],
                powersusparop=sinrmcmc[[54]],
                transparop= sinrmcmc[[55]],
                powertransparop=sinrmcmc[[56]],
                kernelparop=sinrmcmc[[57]],
                sparkop=sinrmcmc[[58]],
                gammaop=sinrmcmc[[59]],
                deltain2op=sinrmcmc[[60]],
                deltanr2op=sinrmcmc[[61]],
                epidatmctim=sinrmcmc[[62]],
                epidatmcrem=sinrmcmc[[63]],
                loglik=sinrmcmc[[64]], NAOK = TRUE
                )
            }
            
            
            
            
            if (nchains > 1L) {
                
                if (parallel == TRUE) {
                    no_cores <- min(nchains,detectCores())
                    cl <- makeCluster(no_cores)
                    varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
                    clusterExport(cl, varlist=varlist, envir=environment())
                    wd <- getwd()
                    clusterExport(cl, varlist="wd", envir=environment())
                    datmcmc22 <- parLapply(cl, seq_len(no_cores), parallel.function)
                    stopCluster(cl)
                } else {
                    datmcmc22 <- list(NULL)
                    for (i in seq_len(nchains)) {
                        datmcmc22[[i]] <- parallel.function(i)
                    }
                }
                
            } else if (nchains == 1L) {
                datmcmc22 <- parallel.function(1)
            }
            
            names <- c(paste("Alpha_s[",seq_len(nsuspar),"]", sep = ""))
            namet <- c(paste("Alpha_t[",seq_len(ntranspar),"]", sep = ""))
            namepowers <- c(paste("Psi_s[",seq_len(nsuspar),"]", sep = ""))
            namepowert <- c(paste("Psi_t[",seq_len(ntranspar),"]", sep = ""))
            
            if (nchains > 1L) {
                
                result77   <-  list(NULL)
                
                for (i in seq_len(nchains)){
                    
                    result777 <- NULL
                    namecols   <-  NULL
                    
                    if (anum2[1]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$susparop    )
                        namecols <- c(namecols, names[1:nsuspar])
                    }
                    if (anum2[7]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$powersusparop)
                        namecols <- c(namecols, namepowers[1:nsuspar])
                    }
                    if (anum2[2]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$transparop)
                        namecols <- c(namecols, namet[1:ntranspar])
                    }
                    if (anum2[8]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$powertransparop)
                        namecols <- c(namecols, namepowert[1:ntranspar])
                    }
                    if (anum2[3]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$sparkop)
                        namecols <- c(namecols, "Spark")
                    }
                    if (anum2[4]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$kernelparop[, 1])
                        namecols <- c(namecols, "Spatial parameter")
                    } else if (anum2[4]==2) {
                        result777 <- cbind(result777, datmcmc22[[i]]$kernelparop)
                        namecols <- c(namecols, c("Spatial parameter", "Network parameter"))
                    }
                    if (anum2[5]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$gammaop)
                        namecols <- c(namecols, "Gamma")
                    }
                    if (anum2[6]==1) {
                        result777 <- cbind(result777, datmcmc22[[i]]$deltain2op)
                        namecols <- c(namecols, "Incubation period rate")
                    } else if (anum2[6]==2) {
                        result777 <- cbind(result777, datmcmc22[[i]]$deltain2op)
                        namecols <- c(namecols, "Incubation period rate")
                        result777 <- cbind(result777, datmcmc22[[i]]$deltanr2op)
                        namecols <- c(namecols, "Delay period rate")
                    }
                    
                    result777 <-cbind(result777, datmcmc22[[i]]$loglik)
                    namecols <-c(namecols, "Log-likelihood")
                    
                    result777  <- data.frame(result777)
                    colnames(result777)  <- namecols
                    
                    if (anum2[6]==3) {
                        result777 <- list(result777)
                    } else if (anum2[6]==1) {
                        result777 <- list(result777, datmcmc22[[i]]$epidatmctim)
                    } else if (anum2[6]==2) {
                        result777 <- list(result777, datmcmc22[[i]]$epidatmctim, datmcmc22[[i]]$epidatmcrem)
                    }
                    
                    result77[[i]] <- result777
                }
            } else {
                
                result77   <-  NULL
                namecols  <-  NULL
                if (anum2[1]==1) {
                    result77 <- cbind(result77, datmcmc22$susparop    )
                    namecols <- c(namecols, names[1:nsuspar])
                }
                if (anum2[7]==1) {
                    result77 <- cbind(result77, datmcmc22$powersusparop)
                    namecols <- c(namecols, namepowers[1:nsuspar])
                }
                if (anum2[2]==1) {
                    result77 <- cbind(result77, datmcmc22$transparop)
                    namecols <- c(namecols, namet[1:ntranspar])
                }
                if (anum2[8]==1) {
                    result77 <- cbind(result77, datmcmc22$powertransparop)
                    namecols <- c(namecols, namepowert[1:ntranspar])
                }
                if (anum2[3]==1) {
                    result77 <- cbind(result77, datmcmc22$sparkop)
                    namecols <- c(namecols, "Spark")
                }
                if (anum2[4]==1) {
                    result77 <- cbind(result77, datmcmc22$kernelparop[, 1])
                    namecols <- c(namecols, "Spatial parameter")
                } else if (anum2[4]==2) {
                    result77 <- cbind(result77, datmcmc22$kernelparop)
                    namecols <- c(namecols, c("Spatial parameter", "Network parameter"))
                }
                if (anum2[5]==1) {
                    result77 <- cbind(result77, datmcmc22$gammaop)
                    namecols <- c(namecols, "Gamma")
                }
                if (anum2[6]==1) {
                    result77 <- cbind(result77, datmcmc22$deltain2op)
                    namecols <- c(namecols, "Incubation period rate")
                } else if (anum2[6]==2) {
                    result77 <- cbind(result77, datmcmc22$deltain2op)
                    namecols <- c(namecols, "Incubation period rate")
                    result77 <- cbind(result77, datmcmc22$deltanr2op)
                    namecols <- c(namecols, "Delay period rate")
                }
                
                result77 <-cbind(result77, datmcmc22$loglik)
                namecols <-c(namecols, "Log-likelihood")
                
                result77  <- data.frame(result77)
                colnames(result77)  <- namecols
                
                if (anum2[6]==3) {
                    result77 <- list(result77)
                } else if (anum2[6]==1) {
                    result77 <- list(result77, datmcmc22$epidatmctim)
                } else if (anum2[6]==2) {
                    result77 <- list(result77, datmcmc22$epidatmctim, datmcmc22$epidatmcrem)
                }
            }
            
            # Creating an epictmcmc object:
            
            if (nchains == 1){
                
                if (length(result77) == 1) {
                    dim.results <- dim(result77[[1]])
                    mcmcsamp <- as.mcmc(result77[[1]][,-dim.results[2]])
                    log.likelihood <- as.mcmc(result77[[1]][,dim.results[2]])
                    accpt <- 1-rejectionRate(as.mcmc(mcmcsamp))
                    num.iter <- niter(as.mcmc(mcmcsamp))
                    num.par <- nvar(as.mcmc(mcmcsamp))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains)
                    
                } else if (length(result77) == 2) {
                    
                    dim.results <- dim(result77[[1]])
                    mcmcsamp <- as.mcmc(result77[[1]][,-dim.results[2]])
                    log.likelihood <- as.mcmc(result77[[1]][,dim.results[2]])
                    accpt <- 1-rejectionRate(as.mcmc(mcmcsamp))
                    num.iter <- niter(as.mcmc(mcmcsamp))
                    num.par <- nvar(as.mcmc(mcmcsamp))
                    num.inf <- sum(object$epidat[,4]!=Inf)
                    infection.times.samples <- as.mcmc(result77[[2]][,1:num.inf])
                    Average.incubation.periods <- as.mcmc(apply(object$epidat[1:num.inf,4] - infection.times.samples,1,mean))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                    Average.incubation.periods = Average.incubation.periods)
                    
                } else if (length(result77) == 3) {
                    
                    dim.results <- dim(result77[[1]])
                    mcmcsamp <- as.mcmc(result77[[1]][,-dim.results[2]])
                    log.likelihood <- as.mcmc(result77[[1]][,dim.results[2]])
                    accpt <- 1-rejectionRate(as.mcmc(mcmcsamp))
                    num.iter <- niter(as.mcmc(mcmcsamp))
                    num.par <- nvar(as.mcmc(mcmcsamp))
                    num.inf <- sum(object$epidat[,4]!=Inf)
                    infection.times.samples <- as.mcmc(result77[[2]][,1:num.inf])
                    Average.incubation.periods <- as.mcmc(apply(object$epidat[1:num.inf,4] - infection.times.samples,1,mean))
                    removal.times.samples <- as.mcmc(result77[[3]][,1:num.inf])
                    Average.delay.periods <- as.mcmc(apply(removal.times.samples - object$epidat[1:num.inf,4],1,mean))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                    Average.incubation.periods = Average.incubation.periods,
                    removal.times.samples = removal.times.samples,
                    Average.delay.periods = Average.delay.periods)
                    
                }
                
            } else{
                
                if (length(result77[[1]]) == 1) {
                    dim.results <- dim(result77[[1]][[1]])
                    mcmcsamp <- list(NULL)
                    log.likelihood <- list(NULL)
                    accpt <- list(NULL)
                    
                    for (i in seq_len(nchains)) {
                        mcmcsamp[[i]] <- as.mcmc(result77[[i]][[1]][,-dim.results[2]])
                        accpt[[i]] <- 1-rejectionRate(as.mcmc(mcmcsamp[[i]]))
                        log.likelihood[[i]] <- as.mcmc(result77[[i]][[1]][,dim.results[2]])
                    }
                    
                    num.iter <- niter(as.mcmc(mcmcsamp[[1]]))
                    num.par <- nvar(as.mcmc(mcmcsamp[[1]]))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains)
                    
                } else if (length(result77[[1]]) == 2) {
                    
                    dim.results <- dim(result77[[1]][[1]])
                    mcmcsamp <- list(NULL)
                    log.likelihood <- list(NULL)
                    accpt <- list(NULL)
                    num.inf <- sum(object$epidat[,4]!=Inf)
                    infection.times.samples <- list(NULL)
                    Average.incubation.periods <- list(NULL)
                    
                    for (i in seq_len(nchains)) {
                        mcmcsamp[[i]] <- as.mcmc(result77[[i]][[1]][,-dim.results[2]])
                        accpt[[i]] <- 1-rejectionRate(as.mcmc(mcmcsamp[[i]]))
                        log.likelihood[[i]] <- as.mcmc(result77[[i]][[1]][,dim.results[2]])
                        infection.times.samples[[i]] <- as.mcmc(result77[[i]][[2]][,seq_len(num.inf)])
                        Average.incubation.periods[[i]] <- as.mcmc(apply(object$epidat[seq_len(num.inf),4]-infection.times.samples[[i]],1,mean))
                    }
                    
                    num.iter <- niter(as.mcmc(mcmcsamp[[1]]))
                    num.par <- nvar(as.mcmc(mcmcsamp[[1]]))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                    Average.incubation.periods = Average.incubation.periods)
                    
                } else if (length(result77[[1]]) == 3) {
                    
                    dim.results <- dim(result77[[1]][[1]])
                    mcmcsamp <- list(NULL)
                    log.likelihood <- list(NULL)
                    accpt <- list(NULL)
                    num.inf <- sum(object$epidat[,4]!=Inf)
                    infection.times.samples <- list(NULL)
                    Average.incubation.periods <- list(NULL)
                    removal.times.samples <- list(NULL)
                    Average.delay.periods <- list(NULL)
                    
                    for (i in seq_len(nchains)) {
                        mcmcsamp[[i]] <- as.mcmc(result77[[i]][[1]][,-dim.results[2]])
                        accpt[[i]] <- 1-rejectionRate(as.mcmc(mcmcsamp[[i]]))
                        log.likelihood[[i]] <- as.mcmc(result77[[i]][[1]][,dim.results[2]])
                        infection.times.samples[[i]] <- as.mcmc(result77[[i]][[2]][,seq_len(num.inf)])
                        Average.incubation.periods[[i]] <- as.mcmc(apply(object$epidat[seq_len(num.inf),4] - infection.times.samples[[i]],1,mean))
                        removal.times.samples[[i]] <- as.mcmc(result77[[i]][[3]][,seq_len(num.inf)])
                        Average.delay.periods[[i]] <- as.mcmc(apply(removal.times.samples[[i]] - object$epidat[seq_len(num.inf),4],1,mean))
                    }
                    
                    num.iter <- niter(as.mcmc(mcmcsamp))
                    num.par <- nvar(as.mcmc(mcmcsamp))
                    
                    out <- list(compart.framework = object$type, kernel.type = object$kerneltype, data.assumption = datatype,
                    parameter.samples = mcmcsamp, log.likelihood = log.likelihood,
                    acceptance.rate = accpt, number.iteration = num.iter,
                    number.parameter = num.par, number.chains = nchains, infection.times.samples = infection.times.samples,
                    Average.incubation.periods = Average.incubation.periods,
                    removal.times.samples = removal.times.samples,
                    Average.delay.periods = Average.delay.periods)
                    
                }
                
            }

    } else {
        stop("Error in specifying the compartmental framework of the model: type", call. = FALSE)
    }
 
 }
# Naming the class of the object:

    class(out) <- "epictmcmc"
    out

}


# S3 Method to print the output of the epictmcmc function:
# x is must be an epictmcmc object

print.epictmcmc <- function(x, digits = 6, ...) {

    if (class(x) != "epictmcmc") {
        stop("The object has to be of \"epictmcmc\" class", call. = FALSE)
    }

    if (x$number.chains == 1){
        if (x$compart.framework == "SIR" & x$data.assumption == "known epidemic") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
            cat(paste("Data assumption: fully observed epidemic"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(x$parameter.samples), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(x$log.likelihood), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SIR" & x$data.assumption == "known removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(x$parameter.samples), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(x$log.likelihood), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(x$infection.times.samples), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(x$Average.infectious.periods), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "known epidemic") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SINR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
            cat(paste("Data assumption: fully observed epidemic"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(x$parameter.samples), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(x$log.likelihood), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "known removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SINR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(x$parameter.samples), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(x$log.likelihood), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(x$infection.times.samples), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(x$Average.incubation.periods), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "unknown removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection & removal times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(x$parameter.samples), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(x$log.likelihood), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(x$infection.times.samples), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(x$Average.incubation.periods), digits)
            cat(paste("\n",names(x)[12],": \n"))
            print(head(x$removal.times.samples), digits)
            cat(paste("\n",names(x)[13],": \n"))
            print(head(x$Average.delay.periods), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        }
        
    } else {
        
        if (x$compart.framework == "SIR" & x$data.assumption == "known epidemic") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
            cat(paste("Data assumption: fully observed epidemic"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(mcmc.list(x$parameter.samples)), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(mcmc.list(x$log.likelihood)), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SIR" & x$data.assumption == "known removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(mcmc.list(x$parameter.samples)), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(mcmc.list(x$log.likelihood)), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(mcmc.list(x$infection.times.samples)), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(mcmc.list(x$Average.infectious.periods)), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "known epidemic") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SINR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
            cat(paste("Data assumption: fully observed epidemic"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(mcmc.list(x$parameter.samples)), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(mcmc.list(x$log.likelihood)), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "known removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SINR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(mcmc.list(x$parameter.samples)), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(mcmc.list(x$log.likelihood)), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(mcmc.list(x$infection.times.samples)), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(mcmc.list(x$Average.incubation.periods)), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        } else if (x$compart.framework == "SINR" & x$data.assumption == "unknown removal") {
            cat(paste("\n********************************************************* \n"))
            cat(paste("Model:", "SIR", x$kernel.type, "-based continuous-time ILM \n" ))
            cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
            cat(paste("Data assumption: partially observed epidemic (unknown infection & removal times)"))
            cat(paste("\n********************************************************* \n"))
            cat(paste("\n"))
            cat(paste("Output: \n"))
            cat(paste("\n",names(x)[4],": \n"))
            print(head(mcmc.list(x$parameter.samples)), digits)
            cat(paste("\n",names(x)[5],": \n"))
            print(head(mcmc.list(x$log.likelihood)), digits)
            cat(paste("\n",names(x)[10],": \n"))
            print(head(mcmc.list(x$infection.times.samples)), digits)
            cat(paste("\n",names(x)[11],": \n"))
            print(head(mcmc.list(x$Average.incubation.periods)), digits)
            cat(paste("\n",names(x)[12],": \n"))
            print(head(mcmc.list(x$removal.times.samples)), digits)
            cat(paste("\n",names(x)[13],": \n"))
            print(head(mcmc.list(x$Average.delay.periods)), digits)
            cat("\nAvailable components:\n")
            print(names(x))
            invisible(x)
        }
        
    }
    
}

# S3 Method to print the summary of the epictmcmc outputs:
# x is must be an epictmcmc object

summary.epictmcmc <- function(object, digits = NULL, start = NULL, end = NULL, thin = NULL, ...)
{
    
    if (is.null(start)) {
        start <- 1
    }
    
    if (is.null(end)) {
        end <- object$number.iteration
    }
    
    if (is.null(thin)) {
        thin <- 1
    }
    
    if (is.null(digits)) {
        digits <- 6
    }
    
    print.summary.epictmcmc(object, digits = digits, start = start, end = end, thin = thin, ...)
}

print.summary.epictmcmc <- function(x, digits, start, end, thin, ...) {
    
    if (class(x) != "epictmcmc") {
        stop("The object has to be of \"epictmcmc\" class", call. = FALSE)
    }
    
    if (x$compart.framework == "SIR" & x$data.assumption == "known epidemic") {
        cat(paste("\n********************************************************* \n"))
        cat(paste("Model: SIR ", x$kernel.type, "-based continuous-time ILM \n", sep=""))
        cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
        cat(paste("Data assumption: fully observed epidemic \n"))
        cat(paste(names(x)[9],":",x$number.chains, "chains \n"))
        cat(paste(names(x)[7],":",end-start, "iterations \n"))
        cat(paste(names(x)[8],":",x$number.parameter, "parameters \n"))
        cat(paste("\n********************************************************* \n"))
        cat(paste("\n"))
        cat(paste("\n 1. Empirical mean and standard deviation for each variable,\n"))
        cat(paste("plus standard error of the mean:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 2. Quantiles for each variable:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 3. Empirical mean, standard deviation, and quantiles for the log likelihood,\n"))
        print(summary(window(x$log.likelihood, start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 4.",names(x)[6],": \n"))
        print(x$acceptance.rate, digits = digits)
    } else if (x$compart.framework == "SIR" & x$data.assumption == "known removal") {
        cat(paste("\n********************************************************* \n"))
        cat(paste("Model:", "SIR ", x$kernel.type, "-based continuous-time ILM \n", sep="" ))
        cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
        cat(paste("Data assumption: partially observed epidemic (unknown infection times) \n"))
        cat(paste(names(x)[9],":",x$number.chains, "chains \n"))
        cat(paste(names(x)[7],":",end-start, "iterations \n"))
        cat(paste(names(x)[8],":",x$number.parameter, "parameters"))
        cat(paste("\n********************************************************* \n"))
        cat(paste("\n"))
        cat(paste("\n 1. Empirical mean and standard deviation for each variable,\n"))
        cat(paste("plus standard error of the mean:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 2. Quantiles for each variable:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 3. Empirical mean, standard deviation, and quantiles for the log likelihood,\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 4. Empirical mean, standard deviation, and quantiles for the average infectious periods,\n"))
        print(summary(window(mcmc.list(x$Average.infectious.periods), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$Average.infectious.periods), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 5.",names(x)[6],": \n"))
        print(x$acceptance.rate, digits = digits)
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known epidemic") {
        cat(paste("\n********************************************************* \n"))
        cat(paste("Model:", "SINR ", x$kernel.type, "-based continuous-time ILM \n", sep="" ))
        cat(paste("Method: Markov chain Monte Carlo (MCMC) \n"))
        cat(paste("Data assumption: fully observed epidemic \n"))
        cat(paste(names(x)[9],":",x$number.chains, "chains \n"))
        cat(paste(names(x)[7],":",end-start, "iterations \n"))
        cat(paste(names(x)[8],":",x$number.parameter, "parameters"))
        cat(paste("\n********************************************************* \n"))
        cat(paste("\n"))
        cat(paste("\n 1. Empirical mean and standard deviation for each variable,\n"))
        cat(paste("plus standard error of the mean:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 2. Quantiles for each variable:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 3. Empirical mean, standard deviation, and quantiles for the log likelihood,\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 4.",names(x)[6],": \n"))
        print(x$acceptance.rate, digits = digits)
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known removal") {
        cat(paste("\n********************************************************* \n"))
        cat(paste("Model:", "SINR", x$kernel.type, "-based continuous-time ILM \n" ))
        cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
        cat(paste("Data assumption: partially observed epidemic (unknown infection times) \n"))
        cat(paste(names(x)[9],":",x$number.chains, "chains \n"))
        cat(paste(names(x)[7],":",end-start, "iterations \n"))
        cat(paste(names(x)[8],":",x$number.parameter, "parameters \n"))
        cat(paste("\n********************************************************* \n"))
        cat(paste("\n"))
        cat(paste("\n 1. Empirical mean and standard deviation for each variable,\n"))
        cat(paste("plus standard error of the mean:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 2. Quantiles for each variable:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 3. Empirical mean, standard deviation, and quantiles for the log likelihood,\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 4. Empirical mean, standard deviation, and quantiles for the average incubation periods,\n"))
        print(summary(window(x$Average.incubation.periods, start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$Average.incubation.periods), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 5.",names(x)[6],": \n"))
        print(x$acceptance.rate, digits = digits)
    } else if (x$compart.framework == "SINR" & x$data.assumption == "unknown removal") {
        cat(paste("\n********************************************************* \n"))
        cat(paste("Model:", "SINR ", x$kernel.type, "-based continuous-time ILM \n", sep="" ))
        cat(paste("Method: Data augmented Markov chain Monte Carlo (DA-MCMC) \n"))
        cat(paste("Data assumption: partially observed epidemic (unknown infection & removal times) \n"))
        cat(paste(names(x)[9],":",x$number.chains, "chains \n"))
        cat(paste(names(x)[7],":",end-start, "iterations \n"))
        cat(paste(names(x)[8],":",x$number.parameter, "parameters \n"))
        cat(paste("\n********************************************************* \n"))
        cat(paste("\n"))
        cat(paste("\n 1. Empirical mean and standard deviation for each variable,\n"))
        cat(paste("plus standard error of the mean:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 2. Quantiles for each variable:\n"))
        print(summary(window(mcmc.list(x$parameter.samples), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 3. Empirical mean, standard deviation, and quantiles for the log likelihood,\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$log.likelihood), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 4. Empirical mean, standard deviation, and quantiles for the average incubation periods,\n"))
        print(summary(window(mcmc.list(x$Average.incubation.periods), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$Average.incubation.periods), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 5. Empirical mean, standard deviation, and quantiles for the average delay periods,\n"))
        print(summary(window(mcmc.list(x$Average.delay.periods), start = start, end = end,
        thin = thin))$statistics, digits = digits)
        cat(paste("\n"))
        print(summary(window(mcmc.list(x$Average.delay.periods), start = start, end = end,
        thin = thin))$quantiles, digits = digits)
        cat(paste("\n"))
        cat(paste("\n"))
        cat(paste("\n 6.",names(x)[6],": \n"))
        print(x$acceptance.rate, digits = digits)
    }
    
}

# to plot the traceplots and densities of the model parameters:
# in case of unobserved event times, the below plot function will
# produce plot of the average posterior and 95% of the unobserved event times.

plot.epictmcmc <- function(x, epi = NULL, start = NULL, end = NULL, thin = NULL, ...) {
    
    if (class(x) != "epictmcmc") {
        stop("The object x has to be of \"epictmcmc\" class", call. = FALSE)
    }
    
    if (is.null(start)) {
        start <- 1
    }
    
    if (is.null(end)) {
        end <- x$number.iteration
    }
    
    if (is.null(thin)) {
        thin <- 1
    }
    
    if (x$compart.framework == "SIR" & x$data.assumption == "known epidemic") {
        
        if (x$number.chains == 1) {
            part <- window(x$parameter.samples, start = start, end = end, thin = thin)
            plot(part, ...)
        } else {
            part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
            plot(part, ...)
        }
    } else if (x$compart.framework == "SIR" & x$data.assumption == "known removal") {
        
        if (x$number.chains == 1) {
            
            part <- window(x$parameter.samples, start = start, end = end, thin = thin)
            plot(part, ...)

            if (!is.null(epi) & class(epi) == "datagen") {

                inft <- window(x$infection.times.samples, start = start, end = end, thin = thin)
                k1 <- length(inft[1,])
                minpoint <- min(apply(inft, 2, min))
                maxpoint <- max(epi$epidat[1:k1,2])
                
                plot(epi$epidat[1:k1, 2], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Removal times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)

            }
                
        } else {
            
            part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
            plot(part, ...)
            
            if (!is.null(epi) & class(epi) == "datagen") {

                inft <- window(mcmc.list(x$infection.times.samples), start = start, end = end, thin = thin)
                k1 <- length(inft[[1]][1,])
                minpoint1 <- NULL
                
                for (i in seq_len(x$number.chains)) {
                    minpoint1[i] <- min(apply(inft[[i]], 2, min))
                }
                minpoint <- min(minpoint1)
                maxpoint <- max(epi$epidat[1:k1,2])
                
                plot(epi$epidat[1:k1, 2], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Removal times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
            }
        }
        
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known epidemic") {
        
        if (x$number.chains == 1) {
            part <- window(x$parameter.samples, start = start, end = end, thin = thin)
            plot(part, ...)
        } else {
            part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
            plot(part, ...)
        }
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known removal") {
        
        if (x$number.chains == 1) {
            
            part <- window(x$parameter.samples, start = start, end = end, thin = thin)
            plot(part, ...)
            
            if (!is.null(epi) & class(epi) == "datagen") {
                inft <- window(x$infection.times.samples, start = start, end = end, thin = thin)
                k1 <- length(inft[1,])
                minpoint <- min(apply(inft, 2, min))
                maxpoint <- max(epi$epidat[1:k1,4])
                
                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
            }
            
        } else {
            
            part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
            plot(part, ...)
            
            if (!is.null(epi) & class(epi) == "datagen") {

                inft <- window(mcmc.list(x$infection.times.samples), start = start, end = end, thin = thin)
                k1 <- length(inft[[1]][1,])
                
                minpoint1 <- NULL
                
                for (i in seq_len(x$number.chains)) {
                    minpoint1[i] <- min(apply(inft[[i]], 2, min))
                }
                
                minpoint <- min(minpoint1)
                maxpoint <- max(epi$epidat[1:k1,4])
                
                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
            }
        }
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "unknown removal") {
        
        if (x$number.chains == 1) {
            
            part <- window(x$parameter.samples, start = start, end = end, thin = thin)
            plot(part, ...)

            if (!is.null(epi) & class(epi) == "datagen") {
                inft <- window(x$infection.times.samples, start = start, end = end, thin = thin)
                remt <- window(x$removal.times.samples, start = start, end = end, thin = thin)
                
                k1 <- length(inft[1,])
                minpoint <- min(apply(inft, 2, min))
                maxpoint <- max(epi$epidat[1:k1,4])
                
                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
                
                minpoint <- min(apply(remt, 2, min))
                maxpoint <- max(epi$epidat[1:k1,4])
                
                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Removal times", main = "The average posterior and 95% CI \n of the removal times")
                lines(summary(remt)$statistics[,1], col = "red", type = "o")
                lines(summary(remt)$quantiles[,1], col = "red", lty = 2)
                lines(summary(remt)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of removal times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
            }
        } else {
            
            part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
            plot(part, ...)
            
            if (!is.null(epi) & class(epi) == "datagen") {
                inft <- window(mcmc.list(x$infection.times.samples), start = start, end = end, thin = thin)
                remt <- window(mcmc.list(x$removal.times.samples), start = start, end = end, thin = thin)
                k1 <- length(inft[[1]][1,])

                minpoint1 <- NULL

                for (i in seq_len(x$number.chains)) {
                    minpoint1[i] <- min(apply(inft[[i]], 2, min))
                }

                minpoint <- min(minpoint1)
                maxpoint <- max(epi$epidat[1:k1,4])

                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Infection times", main = "The average posterior and 95% CI \n of the infection times")
                lines(summary(inft)$statistics[,1], col = "red", type = "o")
                lines(summary(inft)$quantiles[,1], col = "red", lty = 2)
                lines(summary(inft)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of infection times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)

                maxpoint1 <- NULL

                for (i in seq_len(x$number.chains)) {
                    maxpoint1[i] <- max(apply(remt[[i]], 2, max))
                }

                minpoint <- min(epi$epidat[1:k1,4])
                maxpoint <- max(maxpoint1)

                plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                ylab = "Removal times", main = "The average posterior and 95% CI \n of the removal times")
                lines(summary(remt)$statistics[,1], col = "red", type = "o")
                lines(summary(remt)$quantiles[,1], col = "red", lty = 2)
                lines(summary(remt)$quantiles[,5], col = "red", lty = 2)
                legend("bottomright", legend=c("Notification times","Average posterior of removal times", "95% CI"),
                col=c("black","red","red"), lty=c(1,1,2), cex = 0.5)
            }
        }
        
    }
    
}




