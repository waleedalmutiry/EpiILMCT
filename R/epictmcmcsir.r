epictmcmcsir <- function(object, distancekernel, datatype, blockupdate, nsim, nchains, sus, suspower, trans, transpower, kernel, spark, delta, periodproposal, parallel, temp1, n, ni, net, dis, num, nsuspar, ntranspar) {

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
                } else if (all(dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) == TRUE) {
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

        anum55  <-  kernel[[4]]
        kernelparproposalvar <-  kernel[[1]]
        kernelparprior <-  kernel[[3]]
        priordistkernelparpar  <-  kernel[[2]]
        initial[[5]] <-  kernel[[5]]

        initial[[4]] <-  spark[[5]]
        anum44  <-  spark[[4]]
        sparkproposalvar  <-  spark[[1]]
        priordistsparkpar  <-  spark[[2]]
        sparkprior  <-  spark[[3]]


        anum22  <-  sus[[4]]
        initial[[2]] <- sus[[5]]
        suscov  <-  sus[[6]]
        susproposalvar  <-  sus[[1]]
        priordistsuspar  <-  sus[[2]]
        priorpar1sus  <-  sus[[3]][[1]]
        priorpar2sus  <-  sus[[3]][[2]]
        anum66      <-  suspower[[4]]
        initial[[6]] <- suspower[[5]]
        powersusproposalvar <-  suspower[[1]]
        priordistpowersus <-  suspower[[2]]
        priorpar1powersus <-  suspower[[3]][[1]]
        priorpar2powersus <-  suspower[[3]][[2]]

        anum33  <-  trans[[4]]
        initial[[3]] <- trans[[5]]
        transcov  <-  trans[[6]]
        transproposalvar  <-  trans[[1]]
        priordisttranspar  <-  trans[[2]]
        priorpar1trans  <-  trans[[3]][[1]]
        priorpar2trans  <-  trans[[3]][[2]]
        anum77      <-  transpower[[4]]
        initial[[7]] <- transpower[[5]]
        powertransproposalvar <-  transpower[[1]]
        priordistpowertrans <-  transpower[[2]]
        priorpar1powertrans <-  transpower[[3]][[1]]
        priorpar2powertrans <-  transpower[[3]][[2]]

        initial[[8]] <- temp1

        anum2  <-  c(anum11, anum22, anum33, anum44, anum55, anum66, anum77)

        cat("************************************************","\n")
        cat("* Start performing MCMC for the ", datatype," SIR ILM for","\n")
        cat(nsim, "iterations", "\n")
        cat("************************************************","\n")

        if (nchains > 1L) {

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
            delta2op=matrix(as.double(0), ncol=1, nrow=nsim);
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
            kernelparop, sparkop, delta2op, epidatmctim, epidatmcrem, loglik)

            parallel.function <- function(i) {
                .Fortran("mcmcsir_f",
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
                delta2op=  sirmcmc[[52]],
                epidatmctim=  sirmcmc[[53]],
                epidatmcrem=  sirmcmc[[54]],
                loglik=  sirmcmc[[55]], NAOK = TRUE)
            }

            if (parallel) {
                no_cores <- min(nchains,getOption("cl.cores", detectCores()))
                cl <- makeCluster(no_cores)
                varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
                clusterExport(cl, varlist=varlist, envir=environment())
                wd <- getwd()
                clusterExport(cl, varlist="wd", envir=environment())
                datmcmc22 <- parLapply(cl=cl, X=1:no_cores, fun=parallel.function)
                stopCluster(cl)
            } else {
                datmcmc22 <- list(NULL)
                for (i in seq_len(nchains)) {
                    datmcmc22[[i]] <- parallel.function(i)
                }
            }

        } else if (nchains == 1L) {
            datmcmc22 <- .Fortran("mcmcsir_f",
            n=as.integer(n),
            nsim=as.integer(nsim),
            ni=as.integer(ni),
            num=as.integer(num),
            anum2=as.vector(anum2, mode="integer"),
            temp = as.integer(initial[[8]][1]),
            nsuspar=as.integer(nsuspar),
            ntranspar=as.integer(ntranspar),
            net=matrix(as.double(net), ncol=n, nrow=n),
            dis=matrix(as.double(dis), ncol=n, nrow=n),
            epidat=matrix(as.double(object$epidat), ncol=4, nrow=n),
            blockupdate=as.vector(blockupdate, mode="integer"),
            priordistsuspar=as.vector(priordistsuspar, mode="integer"),
            priordisttranspar=as.vector(priordisttranspar, mode="integer"),
            priordistkernelparpar=as.vector(priordistkernelparpar, mode="integer"),
            priordistsparkpar=as.integer(priordistsparkpar),
            priordistpowersus=as.vector(priordistpowersus, mode="integer"),
            priordistpowertrans=as.vector(priordistpowertrans, mode="integer"),
            suspar=as.vector(initial[[2]][,1], mode="double"),
            suscov=matrix(as.double(suscov), ncol=nsuspar, nrow=n),
            powersus=as.vector(initial[[6]][,1], mode="double"),
            transpar=as.vector(initial[[3]][,1], mode="double"),
            transcov=matrix(as.double(transcov), ncol=ntranspar, nrow=n),
            powertrans=as.vector(initial[[7]][,1], mode="double"),
            kernelpar=as.vector(initial[[5]][,1], mode="double"),
            spark=initial[[4]][1],
            delta1=as.vector(initial[[1]][,1], mode = "double"),
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
            delta2op=matrix(as.double(0), ncol=1, nrow=nsim),
            epidatmctim=matrix(as.double(0), ncol=n, nrow=nsim),
            epidatmcrem=matrix(as.double(0), ncol=n, nrow=nsim),
            loglik=matrix(as.double(0), ncol=1, nrow=nsim), NAOK = TRUE)

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

}
