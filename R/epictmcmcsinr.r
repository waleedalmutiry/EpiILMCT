epictmcmcsinr <- function(object, distancekernel, datatype, blockupdate, nsim, nchains, sus, suspower, trans, transpower, kernel, spark, delta, gamma.par, periodproposal, parallel, temp1, n, ni, net, dis, num, nsuspar, ntranspar) {

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

            if (is.null(periodproposal)) {
                infperiodproposalin  <-  c(0, 0)
                infperiodproposalnr  <-  c(0, 0)
            } else {
                if (!is.matrix(periodproposal)) {
                    periodproposal          <-  matrix(periodproposal, ncol=2, nrow=1)
                    infperiodproposalin  <-  periodproposal[1, ]
                    infperiodproposalnr  <-  c(0, 0)

                } else if (all(dim(periodproposal)[1]!=1 & dim(periodproposal)[2]!=2) == TRUE) {
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
                    if (all((dim(delta[[2]])[1] ==2) & (dim(delta[[2]])[2] == nchains)) == TRUE) {
                        initial[[6]][[1]][2,] <-  delta[[2]][1,]
                        initial[[6]][[2]][2,] <-  delta[[2]][2,]
                    } else {
                        stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
                    }
                } else {
                    stop("Error in entering the initial values of the delta parameter: delta[[2]]",  call.= FALSE)
                }

                if (is.matrix(delta[[3]])) {
                    if (all((dim(delta[[3]])[1] ==2) & (dim(delta[[3]])[2] == 2)) == TRUE) {
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
                } else if (all(dim(periodproposal)[1]!=2 & dim(periodproposal)[2]!=2) == TRUE) {
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

        anum44  <-  kernel[[4]]
        kernelparproposalvar <-  kernel[[1]]
        kernelparprior <-  kernel[[3]]
        priordistkernelparpar  <-  kernel[[2]]
        initial[[4]] <-  kernel[[5]]

        initial[[3]] <-  spark[[5]]
        anum33  <-  spark[[4]]
        sparkproposalvar  <-  spark[[1]]
        priordistsparkpar  <-  spark[[2]]
        sparkprior  <-  spark[[3]]


        anum11  <-  sus[[4]]
        initial[[1]] <- sus[[5]]
        suscov  <-  sus[[6]]
        susproposalvar  <-  sus[[1]]
        priordistsuspar  <-  sus[[2]]
        priorpar1sus  <-  sus[[3]][[1]]
        priorpar2sus  <-  sus[[3]][[2]]
        anum77      <-  suspower[[4]]
        initial[[7]] <- suspower[[5]]
        powersusproposalvar <-  suspower[[1]]
        priordistpowersus <-  suspower[[2]]
        priorpar1powersus <-  suspower[[3]][[1]]
        priorpar2powersus <-  suspower[[3]][[2]]

        anum22  <-  trans[[4]]
        initial[[2]] <- trans[[5]]
        transcov  <-  trans[[6]]
        transproposalvar  <-  trans[[1]]
        priordisttranspar  <-  trans[[2]]
        priorpar1trans  <-  trans[[3]][[1]]
        priorpar2trans  <-  trans[[3]][[2]]
        anum88      <-  transpower[[4]]
        initial[[8]] <- transpower[[5]]
        powertransproposalvar <-  transpower[[1]]
        priordistpowertrans <-  transpower[[2]]
        priorpar1powertrans <-  transpower[[3]][[1]]
        priorpar2powertrans <-  transpower[[3]][[2]]

        initial[[9]] <- temp1

        anum2  <-  c(anum11, anum22, anum33, anum44, anum55, anum66, anum77, anum88)

        cat("************************************************","\n")
        cat("* Start performing MCMC for the ", datatype," SINR ILM for","\n")
        cat(nsim, "iterations", "\n")
        cat("************************************************","\n")


        n=as.integer(n);
        nsim=as.integer(nsim);
        ni=as.integer(ni);
        temp = as.integer(initial[[9]][1]);
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
            .Fortran("mcmcsinr_f",
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
            loglik=sinrmcmc[[64]], NAOK = TRUE)
        }




        if (nchains > 1L) {

            if (parallel) {

                no_cores <- min(nchains, getOption("cl.cores", detectCores()))
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

}
