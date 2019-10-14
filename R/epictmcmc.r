epictmcmc <- function(object, distancekernel = NULL, datatype, blockupdate = NULL, nsim, nchains = NULL, control.sus = NULL, control.trans = NULL, kernel.par = NULL, spark.par = NULL, delta = NULL, gamma.par = NULL, periodproposal = NULL, parallel = FALSE, seedval = NULL) {
    
    
    if (class(object) != "datagen") {
        stop("The epidat object must be in a class of \"datagen\" ", call. = FALSE)
    } else {
        
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

        spark <- list(NULL)
        kernel <- list(NULL)
        kernel[[1]] <-  vector(mode="double", length = 2)
        kernel[[2]] <-  vector(mode="integer", length = 2)
        kernel[[3]] <-  matrix(0, ncol = 2, nrow = 2)
        
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
            
            kernel[[4]]  <-  1
            
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
                kernel[[1]][1] <-  kernel.par[[2]][4]
                kernel[[1]][2] <-  0
                kernel[[3]][1, ] <-  kernel.par[[2]][2:3]
                kernel[[3]][2, ] <-  c(0, 0)
                
                if (kernel.par[[2]][1] == "gamma") {
                    kernel[[2]][1]  <-  1
                } else if (kernel.par[[2]][1] == "half normal") {
                    kernel[[2]][1]  <-  2
                } else if (kernel.par[[2]][1] == "uniform") {
                    kernel[[2]][1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par", call.=FALSE)
                }
                
                kernel[[2]][2]  <-  1
                
                if (length(kernel.par[[1]]) == nchains){
                    kernel[[5]] <-  matrix(c(kernel.par[[1]],rep(0,nchains)), ncol = nchains, byrow = TRUE)
                } else if(length(kernel.par[[1]]) == 1) {
                    kernel[[5]] <-  matrix(rep(c(kernel.par[[1]],0),nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("Error in entering the initial values of the kernel parameter: kernel.par[[1]]",  call.= FALSE)
                }
                
            }
            
            if (is.null(spark.par)) {
                
                spark[[5]] <-  rep(0,nchains)
                spark[[4]]  <-  2
                spark[[1]]  <-  0
                spark[[2]]  <-  1
                spark[[3]]  <-  c(1, 1)
                
            } else {
                
                spark[[4]] <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    spark[[1]] <-  spark.par[[2]][4]
                    spark[[3]] <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        spark[[2]]  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        spark[[2]]  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        spark[[2]]  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        spark[[5]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        spark[[5]] <-  rep(spark.par[[1]], nchains)
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
            kernel[[2]]  <-  c(1, 1)
            kernel[[1]]  <-  c(0, 0)
            kernel[[3]] <-  matrix(0, ncol = 2, nrow = 2)
            kernel[[5]] <-  matrix(rep(c(0,0),nchains), ncol = nchains, nrow = 2)
            
            kernel[[4]] <-  3
            
            if (datatype == "known removal" | datatype == "unknown removal") {
                
                spark[[4]]   <-  1
                
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
                    spark[[1]] <-  spark.par[[2]][4]
                    spark[[3]] <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        spark[[2]]  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        spark[[2]]  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        spark[[2]]  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        spark[[5]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        spark[[5]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                }
                
            } else {
                
                
                if (is.null(spark.par)) {
                    
                    spark[[4]]  <-  2
                    spark[[1]] <-  0
                    spark[[2]] <-  1
                    spark[[3]] <-  c(1, 1)
                    spark[[5]] <-  rep(0,nchains)
                    
                } else {
                    spark[[4]] <-  1
                    
                    if (!is.list(spark.par)) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    }
                    
                    if (length(spark.par) != 2L) {
                        stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                    } else {
                        if (length(spark.par[[2]]) != 4L) {
                            stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                        }
                        spark[[1]] <-  spark.par[[2]][4]
                        spark[[3]] <-  spark.par[[2]][2:3]
                        
                        if (spark.par[[2]][1] == "gamma") {
                            spark[[2]]  <-  1
                        } else if (spark.par[[2]][1] == "half normal") {
                            spark[[2]]  <-  2
                        } else if (spark.par[[2]][1] == "uniform") {
                            spark[[2]]  <-  3
                        } else {
                            stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                        }
                        
                        if (length(spark.par[[1]]) == nchains){
                            spark[[5]] <-  spark.par[[1]]
                        } else if(length(spark.par[[1]]) == 1) {
                            spark[[5]] <-  rep(spark.par[[1]],nchains)
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
            
            kernel[[4]]  <-  2
            
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
                kernel[[1]][1] <-  kernel.par[[2]][1,4]
                kernel[[1]][2] <-  kernel.par[[2]][2,4]
                kernel[[3]][1, ] <-  kernel.par[[2]][1,2:3]
                kernel[[3]][2, ] <-  kernel.par[[2]][2,2:3]
                
                if (kernel.par[[2]][1,1] == "gamma") {
                    kernel[[2]][1]  <-  1
                } else if (kernel.par[[2]][1,1] == "half normal") {
                    kernel[[2]][1]  <-  2
                } else if (kernel.par[[2]][1,1] == "uniform") {
                    kernel[[2]][1]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                
                if (kernel.par[[2]][2,1] == "gamma") {
                    kernel[[2]][2]  <-  1
                } else if (kernel.par[[2]][2,1] == "half normal") {
                    kernel[[2]][2]  <-  2
                } else if (kernel.par[[2]][2,1] == "uniform") {
                    kernel[[2]][2]  <-  3
                } else {
                    stop("The possible choices of the prior distribution are \"gamma\", \"half normal\" or \"uniform\" distributions: kernel.par[[2]]", call. = FALSE)
                }
                
                if (is.matrix(kernel.par[[1]])) {
                    if (dim(kernel.par[[1]])[1] != 2 & dim(kernel.par[[1]])[2] != nchains) {
                        stop("The matrix of the initial values of the kernel parameters must be a 2 by \"nchains\".", call. = FALSE)
                    }
                    kernel[[5]] <-  kernel.par[[1]]
                } else if (is.vector(kernel.par[[1]]) & length(kernel.par[[1]]) == 2) {
                    kernel[[5]] <-  matrix(rep(kernel.par[[1]],nchains), ncol = nchains, nrow = 2)
                } else {
                    stop("The initial values of the kernel parameters must be a matrix of 2 by \"nchains\" or a vector of two initial values.", call. = FALSE)
                }
                
            }
            
            if (is.null(spark.par)) {
                
                spark[[4]]  <-  2
                spark[[1]] <-  0
                spark[[2]] <-  1
                spark[[3]] <-  c(1, 1)
                spark[[5]] <-  rep(0,nchains)
                
            } else {
                spark[[4]] <-  1
                
                if (!is.list(spark.par)) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                }
                
                if (length(spark.par) != 2L) {
                    stop("The argument \"spark.par\" must be a list of two:\n1) a vector of initial values of the spark parameter of size equal to \"nchains\". \n2) a vector of prior distribution, prior parameter values, and proposal variance.",  call.= FALSE)
                } else {
                    if (length(spark.par[[2]]) != 4L) {
                        stop("Error in entering the second list of the argument \"spark.par\" for updating the spark parameter.",  call.= FALSE)
                    }
                    spark[[1]] <-  spark.par[[2]][4]
                    spark[[3]] <-  spark.par[[2]][2:3]
                    
                    if (spark.par[[2]][1] == "gamma") {
                        spark[[2]]  <-  1
                    } else if (spark.par[[2]][1] == "half normal") {
                        spark[[2]]  <-  2
                    } else if (spark.par[[2]][1] == "uniform") {
                        spark[[2]]  <-  3
                    } else {
                        stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: spark.par", call.=FALSE)
                    }
                    
                    if (length(spark.par[[1]]) == nchains){
                        spark[[5]] <-  spark.par[[1]]
                    } else if(length(spark.par[[1]]) == 1) {
                        spark[[5]] <-  rep(spark.par[[1]],nchains)
                    } else {
                        stop("Error in entering the initial values of the spark parameter: spark.par[[1]]",  call.= FALSE)
                    }
                    
                }
                
            }
        }
        
        # check the number of infected individuals:
        
        if (datatype == "known epidemic" | datatype == "known removal" | datatype == "unknown removal") {
            ni  <-  sum(object$epidat[, 2]!=Inf)
        } else {
            stop("Error: the epidemic data must be in the same format as in datagen function where the removal and infection times of uninfected individuals are \"Inf\".", call.=FALSE)
        }
        
        
        # check Susceptibility terms:
        sus <- list(NULL)
        suspower <- list(NULL)
        
        if (is.null(control.sus)){
            
            nsuspar  <-  1
            sus[[4]]  <-  2
            sus[[5]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            sus[[6]]  <-  matrix(rep(1, n), ncol= nsuspar, nrow = n)
            sus[[1]]  <-  0
            sus[[2]]  <-  1
            sus[[3]] <- list(NULL)
            sus[[3]][[1]]  <-  rep(0, nsuspar)
            sus[[3]][[2]]  <-  rep(0, nsuspar)
            suspower[[4]]      <-  2
            suspower[[5]] <- matrix(rep(1,nchains), ncol = nchains, nrow = nsuspar)
            suspower[[1]] <-  0
            suspower[[2]] <-  1
            suspower[[3]] <- list(NULL)
            suspower[[3]][[1]] <-  rep(0, nsuspar)
            suspower[[3]][[2]] <-  rep(0, nsuspar)
            
        } else if (!is.list(control.sus)) {
            
            stop("The option control.sus must be a list of values of the susceptibility parameters, susceptibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.sus)) {
            
            len <- length(control.sus)
            
            if (len < 2L) {
                stop("The length of the option control.sus must be at least 2", call. = FALSE)
            } else {
                
                sus[[4]]       <-  1
                sus.par <- control.sus[[1]]
                sus[[6]]  <- control.sus[[2]]
                
                if (is.matrix(sus[[6]])) {
                    nsuspar <- dim(sus[[6]])[2]
                } else {
                    nsuspar <- 1
                }
                
                sus[[2]]    <-  vector(mode="integer", length = nsuspar)
                suspower[[2]]  <-  vector(mode="integer", length = nsuspar)
                
                if (any(sus[[6]] < 0)) {
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
                            sus[[3]] <- list(NULL)
                            sus[[3]][[1]]    <-  sus.par[[2]][, 2]
                            sus[[3]][[2]]    <-  sus.par[[2]][, 3]
                            sus[[1]]  <-  sus.par[[2]][, 4]
                            
                            for(i in 1:nsuspar) {
                                if (sus.par[[2]][i, 1] == "gamma") {
                                    sus[[2]][i]  <-  1
                                } else if (sus.par[[2]][i, 1] == "half normal") {
                                    sus[[2]][i]  <-  2
                                } else if (sus.par[[2]][i, 1] == "uniform") {
                                    sus[[2]][i]  <-  3
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
                            
                            sus[[3]] <- list(NULL)
                            sus[[3]][[1]]    <-  sus.par[[2]][2]
                            sus[[3]][[2]]    <-  sus.par[[2]][3]
                            sus[[1]]  <-  sus.par[[2]][4]
                            
                            if (sus.par[[2]][1] == "gamma") {
                                sus[[2]][1]  <-  1
                            } else if (sus.par[[2]][1] == "half normal") {
                                sus[[2]][1]  <-  2
                            } else if (sus.par[[2]][1] == "uniform") {
                                sus[[2]][1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                            }
                            
                        }
                        if (is.matrix(sus.par[[1]])) {
                            
                            if ((dim(sus.par[[1]])[1] != nsuspar) & (dim(sus.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                            }
                            sus[[5]]    <-  sus.par[[1]]
                            
                        } else if (is.vector(sus.par[[1]]) & length(sus.par[[1]]) == nsuspar) {
                            sus[[5]]    <-  matrix(rep(sus.par[[1]],nchains), ncol = nchains, nrow = nsuspar)
                        } else if (is.vector(sus.par[[1]]) & length(sus.par[[1]]) == nchains) {
                            sus[[5]]    <-  matrix(sus.par[[1]], ncol = nchains, nrow = nsuspar)
                       }
                        
                    }
                }
                
                if (len == 3L) {
                    power.sus <- control.sus[[3]]
                    suspower[[4]]               <-  1
                    
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
                                
                                suspower[[3]]    <-  list(NULL)
                                suspower[[3]][[1]]    <-  power.sus[[2]][, 2]
                                suspower[[3]][[2]]    <-  power.sus[[2]][, 3]
                                suspower[[1]]  <-  power.sus[[2]][, 4]
                                
                                
                                for(i in 1:nsuspar) {
                                    if (power.sus[[2]][i, 1] == "gamma") {
                                        suspower[[2]][i]  <-  1
                                    } else if (power.sus[[2]][i, 1] == "half normal") {
                                        suspower[[2]][i]  <-  2
                                    } else if (power.sus[[2]][i, 1] == "uniform") {
                                        suspower[[2]][i]  <-  3
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
                                
                                suspower[[3]]    <-  list(NULL)
                                suspower[[3]][[1]]     <-  power.sus[[2]][2]
                                suspower[[3]][[2]]     <-  power.sus[[2]][3]
                                suspower[[1]]   <-  power.sus[[2]][4]
                                
                                
                                if (power.sus[[2]][1] == "gamma") {
                                    suspower[[2]][1]  <-  1
                                } else if (power.sus[[2]][1] == "half normal") {
                                    suspower[[2]][1]  <-  2
                                } else if (power.sus[[2]][1] == "uniform") {
                                    suspower[[2]][1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.sus", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.sus[[1]])) {
                                
                                if ((dim(power.sus[[1]])[1] != nsuspar) & (dim(power.sus[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the susceptibility parameters in the argument of \"control.sus\".", call. = FALSE)
                                }
                                
                                suspower[[5]]    <-  power.sus[[1]]
                                
                            } else if (is.vector(power.sus[[1]]) & length(power.sus[[1]]) == nsuspar) {
                                suspower[[5]]    <-  matrix(rep(power.sus[[1]],nchains), ncol = nchains, nrow = nsuspar)
                            } else if (is.vector(power.sus[[1]]) & length(power.sus[[1]]) == nchains) {
                                suspower[[5]]    <-  matrix(sus.par[[1]], ncol = nchains, nrow = nsuspar)
                            }
                            
                        }
                        
                    }
                } else if (len == 2L) {
                    suspower[[4]]               <-  2
                    suspower[[5]]         <-  matrix(rep(rep(1, nsuspar),nchains), ncol = nchains, nrow = nsuspar)
                    suspower[[1]]  <-  rep(0, nsuspar)
                    suspower[[2]]    <-  rep(1, nsuspar)
                    suspower[[3]]    <-  list(NULL)
                    suspower[[3]][[1]]    <-  rep(0, nsuspar)
                    suspower[[3]][[2]]    <-  rep(0, nsuspar)
                }
            }
        }
        
        # check transmissibility terms:
        trans <- list(NULL)
        transpower <- list(NULL)
        
        if (is.null(control.trans)){
            
            ntranspar  <-  1
            trans[[4]]  <-  2
            trans[[5]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            trans[[6]]  <-  matrix(rep(1, n), ncol= ntranspar, nrow = n)
            trans[[1]]  <-  0
            trans[[2]]  <-  1
            trans[[3]] <- list(NULL)
            trans[[3]][[1]]  <-  rep(0, ntranspar)
            trans[[3]][[2]]  <-  rep(0, ntranspar)
            transpower[[4]]      <-  2
            transpower[[5]] <- matrix(rep(1,nchains), ncol = nchains, nrow = ntranspar)
            transpower[[1]] <-  0
            transpower[[2]] <-  1
            transpower[[3]] <- list(NULL)
            transpower[[3]][[1]] <-  rep(0, ntranspar)
            transpower[[3]][[2]] <-  rep(0, ntranspar)
            
        } else if (!is.list(control.trans)) {
            
            stop("The option control.trans must be a list of values of the transmissibility parameters, transmissibility cavariates, values of the non-linearity (power) parameters of the covariates", call. = FALSE)
            
        } else if (is.list(control.trans)) {
            
            len <- length(control.trans)
            
            if (len < 2L) {
                stop("The length of the option control.trans must be at least 2", call. = FALSE)
            } else {
                
                trans[[4]]       <-  1
                trans.par <- control.trans[[1]]
                trans[[6]]  <- control.trans[[2]]
                
                if (is.matrix(trans[[6]])) {
                    ntranspar <- dim(trans[[6]])[2]
                } else {
                    ntranspar <- 1
                }
                
                trans[[2]]    <-  vector(mode="integer", length = ntranspar)
                transpower[[2]]  <-  vector(mode="integer", length = ntranspar)
                
                if (any(trans[[6]] < 0)) {
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
                            
                            trans[[3]] <- list(NULL)
                            trans[[3]][[1]]    <-  trans.par[[2]][, 2]
                            trans[[3]][[2]]    <-  trans.par[[2]][, 3]
                            trans[[1]]  <-  trans.par[[2]][, 4]
                            
                            for(i in 1:ntranspar) {
                                if (trans.par[[2]][i, 1] == "gamma") {
                                    trans[[2]][i]  <-  1
                                } else if (trans.par[[2]][i, 1] == "half normal") {
                                    trans[[2]][i]  <-  2
                                } else if (trans.par[[2]][i, 1] == "uniform") {
                                    trans[[2]][i]  <-  3
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
                            
                            trans[[3]] <- list(NULL)
                            trans[[3]][[1]]    <-  trans.par[[2]][2]
                            trans[[3]][[2]]    <-  trans.par[[2]][3]
                            trans[[1]]  <-  trans.par[[2]][4]
                            
                            if (trans.par[[2]][1] == "gamma") {
                                trans[[2]][1]  <-  1
                            } else if (trans.par[[2]][1] == "half normal") {
                                trans[[2]][1]  <-  2
                            } else if (trans.par[[2]][1] == "uniform") {
                                trans[[2]][1]  <-  3
                            } else {
                                stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                            }
                            
                        }
                        if (is.matrix(trans.par[[1]])) {
                            
                            if ((dim(trans.par[[1]])[1] != ntranspar) & (dim(trans.par[[1]])[2] != nchains)) {
                                stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                            }
                            
                            trans[[5]]    <-  trans.par[[1]]
                            
                        } else if (is.vector(trans.par[[1]]) & length(trans.par[[1]]) == ntranspar) {
                            trans[[5]]    <-  matrix(rep(trans.par[[1]],nchains), ncol = nchains, nrow = ntranspar)
                        } else if (is.vector(trans.par[[1]]) & length(trans.par[[1]]) == nchains) {
                            trans[[5]]    <-  matrix(trans.par[[1]], ncol = nchains, nrow = ntranspar)
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
                                
                                transpower[[3]] <- list(NULL)
                                transpower[[3]][[1]]    <-  power.trans[[2]][, 2]
                                transpower[[3]][[2]]    <-  power.trans[[2]][, 3]
                                transpower[[1]]  <-  power.trans[[2]][, 4]
                                
                                
                                for(i in 1:ntranspar) {
                                    if (power.trans[[2]][i, 1] == "gamma") {
                                        transpower[[2]][i]  <-  1
                                    } else if (power.trans[[2]][i, 1] == "half normal") {
                                        transpower[[2]][i]  <-  2
                                    } else if (power.trans[[2]][i, 1] == "uniform") {
                                        transpower[[2]][i]  <-  3
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
                                
                                transpower[[3]] <- list(NULL)
                                transpower[[3]][[1]]     <-  power.trans[[2]][, 2]
                                transpower[[3]][[2]]     <-  power.trans[[2]][, 3]
                                transpower[[1]]   <-  power.trans[[2]][, 4]
                                
                                
                                if (power.trans[[2]][1] == "gamma") {
                                    transpower[[2]][1]  <-  1
                                } else if (power.trans[[2]][1] == "half normal") {
                                    transpower[[2]][1]  <-  2
                                } else if (power.trans[[2]][1] == "uniform") {
                                    transpower[[2]][1]  <-  3
                                } else {
                                    stop("The possible choices of the prior distribution are \"gamma\" ,  \"half normal\" or \"uniform\" distributions: control.trans", call. = FALSE)
                                }
                                
                            }
                            
                            if (is.matrix(power.trans[[1]])) {
                                
                                if ((dim(power.trans[[1]])[1] != ntranspar) & (dim(power.trans[[1]])[2] != nchains)) {
                                    stop("Error in entering the initial values for the transmissibility parameters in the argument of \"control.trans\".", call. = FALSE)
                                }
                                
                                transpower[[5]]    <-  power.trans[[1]]
                                
                            } else if (is.vector(power.trans[[1]]) & length(power.trans[[1]]) == ntranspar) {
                                transpower[[5]]    <-  matrix(rep(power.trans[[1]],nchains), ncol = nchains, nrow = ntranspar)
                            } else if (is.vector(power.trans[[1]]) & length(power.trans[[1]]) == nchains) {
                                power.trans[[5]]    <-  matrix(power.trans[[1]], ncol = nchains, nrow = ntranspar)
                            }
                            
                        }
                        
                    }
                } else if (len == 2L) {
                    transpower[[4]]                 <-  2
                    transpower[[5]]           <-  matrix(rep(rep(1,ntranspar),nchains), ncol = nchains, nrow = ntranspar)
                    transpower[[1]]  <-  rep(0, ntranspar)
                    transpower[[2]]    <-  rep(1, ntranspar)
                    transpower[[3]] <- list(NULL)
                    transpower[[3]][[1]]    <-  rep(0, ntranspar)
                    transpower[[3]][[2]]    <-  rep(0, ntranspar)
                }
            }
            
            
        }

        # Running MCMC for the specific compartmental framework:
        
        if (object$type == "SIR") {
            
            out <- epictmcmcsir(object, distancekernel, datatype, blockupdate, nsim, nchains, sus, suspower, trans, transpower, kernel, spark, delta, periodproposal, parallel, temp1, n, ni, net, dis, num, nsuspar, ntranspar)

        } else if (object$type == "SINR") {

            out <- epictmcmcsinr(object, distancekernel, datatype, blockupdate, nsim, nchains, sus, suspower, trans, transpower, kernel, spark, delta, gamma.par, periodproposal, parallel, temp1, n, ni, net, dis, num, nsuspar, ntranspar)

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

# x must be an epictmcmc object

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

plot.epictmcmc <- function(x, epi = NULL, plottype = NULL, start = NULL, end = NULL, thin = NULL, ...) {
    
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
        if (plottype == "parameter") {
            if (x$number.chains == 1) {
                part <- window(x$parameter.samples, start = start, end = end, thin = thin)
                plot(part, ...)
            } else {
                part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
                plot(part, ...)
            }
        } else {
            stop("As datatype = \"known epidemic\", the only option to plot is \"parameters\".", call. = TRUE)
        }
    } else if (x$compart.framework == "SIR" & x$data.assumption == "known removal") {
        if (x$number.chains == 1) {
            if (plottype == "parameter") {
                
                part <- window(x$parameter.samples, start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                    
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\" and \"inf.times\".", call. = TRUE)
            }
            
        } else {
            
            if (plottype == "parameter") {
                part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\" and \"inf.times\".", call. = TRUE)
            }
        }
        
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known epidemic") {
        
        if (plottype == "parameter") {
            
            if (x$number.chains == 1) {
                part <- window(x$parameter.samples, start = start, end = end, thin = thin)
                plot(part, ...)
            } else {
                part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
                plot(part, ...)
            }
            
        } else {
            stop("As datatype = \"known epidemic\", the only option to plot is \"parameters\".", call. = TRUE)
        }
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "known removal") {
        
        if (x$number.chains == 1) {
            
            if (plottype == "parameter") {
                part <- window(x$parameter.samples, start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\" and \"inf.times\".", call. = TRUE)
            }
            
        } else {
            
            if (plottype == "parameter") {
                part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\" and \"inf.times\".", call. = TRUE)
            }
        }
        
    } else if (x$compart.framework == "SINR" & x$data.assumption == "unknown removal") {
        
        if (x$number.chains == 1) {
            if (plottype == "parameter") {
                part <- window(x$parameter.samples, start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else if (plottype == "rem.times") {
                if (!is.null(epi) & class(epi) == "datagen") {
                    remt <- window(x$removal.times.samples, start = start, end = end, thin = thin)
                    k1 <- length(inft[1,])
                    
                    minpoint <- min(apply(remt, 2, min))
                    maxpoint <- max(epi$epidat[1:k1,4])
                    
                    plot(epi$epidat[1:k1, 4], type = "o", ylim = c(minpoint, maxpoint),
                    ylab = "Removal times", main = "The average posterior and 95% CI \n of the removal times")
                    lines(summary(remt)$statistics[,1], col = "red", type = "o")
                    lines(summary(remt)$quantiles[,1], col = "red", lty = 2)
                    lines(summary(remt)$quantiles[,5], col = "red", lty = 2)
                    legend("bottomright", legend=c("Notification times","Average posterior of removal times", "95% CI"),
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\", \"inf.times\", and \"rem.times\".", call. = TRUE)
            }
            
        } else {
            
            if (plottype == "parameter") {
                part <- window(mcmc.list(x$parameter.samples), start = start, end = end, thin = thin)
                plot(part, ...)
            } else if (plottype == "inf.times") {
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else if (plottype == "rem.times") {
                if (!is.null(epi) & class(epi) == "datagen") {
                    remt <- window(mcmc.list(x$removal.times.samples), start = start, end = end, thin = thin)
                    k1 <- length(remt[[1]][1,])
                    
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
                    col=c("black","red","red"), lty=c(1,1,2), cex = 0.8)
                }
            } else {
                stop("As datatype = \"known removal\", the only options to plot are \"parameters\", \"inf.times\", and \"rem.times\".", call. = TRUE)
            }
            
        }
        
    }
    
}





