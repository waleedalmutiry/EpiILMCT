library(EpiILMCT)
set.seed(22)
## distance-based SIR continuous-time ILMs ##
data(SpatialData)
## performing the MCMC-tool for analyzing the fully observed spatial data
## under the SIR distance-based continuous ILM:
suspar <- list(NULL)
suspar[[1]]<-list(2,c("gamma", 1, 0.01, 0.5))
suspar[[2]]<- rep(1,length(SpatialData$epidat[,1]))
kernel1 <- list(2, c("gamma", 1, 0.01, 0.5))

mcmcres2 <- epictmcmc(object = SpatialData,
distancekernel = "powerlaw", datatype = "known epidemic", nsim = 50,
control.sus = suspar,
kernel.par = kernel1)
#plot(mcmcres2, plottype = "parameter")
print(mcmcres2)
summary(mcmcres2)

suspar <- list(NULL)
suspar[[1]]<-list(2,c("gamma", 1, 0.1, 2.5))
suspar[[2]]<- matrix(rep(1,length(SpatialData$epidat[,1])),ncol=1)
kernel1 <- list(0.2, c("gamma", 1, 0.1, 0.01))

mcmcres22 <- epictmcmc(object = SpatialData,
distancekernel = "Cauchy", datatype = "known epidemic", nsim = 50,
control.sus = suspar,
kernel.par = kernel1)
#plot(mcmcres22, plottype = "parameter")
print(mcmcres22)
summary(mcmcres22)

#plot(mcmcres2$log.likelihood)
#plot(mcmcres22$log.likelihood)

## performing the MCMC-tool for analyzing the partially observed spatial
## data (unknown infection times) under the SIR distance-based
## continuous ILM:

suspar <- list(NULL)
suspar[[1]]<-list(2,c("gamma", 1, 0.01, 0.8))
suspar[[2]]<- matrix(rep(1,length(SpatialData$epidat[,1])),ncol=1)
kernel1 <- list(2, c("gamma", 1, 0.01, 0.5))

mcmcres22 <- epictmcmc(object = SpatialData, distancekernel = "powerlaw",
datatype = "known removal", nsim = 50,
control.sus = suspar, kernel.par = kernel1, delta = list(1, 2, c(4, 2)))

#plot(mcmcres22, plottype = "parameter")
print(mcmcres22)
summary(mcmcres22)

## distance-based and network-based SIR ILMs ##
set.seed(22)
data(SpatialNetData)
## performing the MCMC-tool for analyzing the fully observed spatial and
## network data
## under the SIR distance-based and network-based continuous-time ILM:
suspar <- list(NULL)
suspar[[1]]<-list(c(0.08,0.2),matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.1, 0.5),
ncol = 4, nrow = 2))
suspar[[2]]<- SpatialNetData[[2]]
kernel1 <- list(c(5, 0.5), matrix(c("gamma", "gamma", 1, 1,
0.01, 0.01, 0.5, 1), ncol = 4, nrow = 2))

mcmcres3 <- epictmcmc(object = SpatialNetData[[1]], distancekernel = "powerlaw",
datatype = "known epidemic", nsim = 50,
control.sus = suspar, kernel.par = kernel1)
#plot(mcmcres3, plottype = "parameter")
print(mcmcres3)
summary(mcmcres3)

## network-based SIR ILMs ##
set.seed(22)
data(NetworkData)
## performing the MCMC-tool for analyzing the fully observed network data
## under the SIR network-based continuous-time ILM:

suspar <- list(NULL)
suspar[[1]]<-list(c(0.08,0.5),matrix(c("gamma", "gamma", 1, 1, 1, 1, 0.1, 0.5),
ncol = 4, nrow = 2))
suspar[[2]]<- NetworkData[[2]]

mcmcres4 <- epictmcmc(object = NetworkData[[1]], datatype = "known epidemic",
nsim = 50, control.sus = suspar)
#plot(mcmcres4, plottype = "parameter")
print(mcmcres4)
summary(mcmcres4)

## network-based SINR ILMs ##
set.seed(22)
data(NetworkDataSINR)
names(NetworkDataSINR)

netSINR<-as.epidat(type = "SINR", kerneltype = "network", incub.time = NetworkDataSINR$epi[,4], inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "powerlaw")

## performing the MCMC-tool for analyzing the fully observed network data
## under the SINR network-based continuous-time ILM:
suspar <- list(NULL)
suspar[[1]]<-list(c(0.08,0.2),matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.05, 0.5),
ncol = 4, nrow = 2))
suspar[[2]]<- NetworkDataSINR$cov

mcmcres5 <- epictmcmc(object = netSINR, datatype = "known epidemic",
nsim = 500, control.sus = suspar)

mcmcres5
#plot(mcmcres5, plottype = "parameter")
print(mcmcres5)
summary(mcmcres5)

suspar <- list(NULL)
suspar[[1]]<-list(c(0.08,0.2),matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.05, 0.5),
ncol = 4, nrow = 2))
suspar[[2]]<- NetworkDataSINR$cov
delta1<-list(1,2,c(4,2))
spark<-list(1,matrix(c("gamma", 1, 0.01, 0.05), ncol = 4, nrow = 1))

mcmcres5 <- epictmcmc(object = netSINR, datatype = "known removal",
nsim = 500, control.sus = suspar, spark.par = spark, delta = delta1)

print(mcmcres5)
summary(mcmcres5)
#plot(mcmcres5, plottype = "parameter")
#plot(mcmcres5, plottype = "inf.times")

suspar <- list(NULL)
suspar[[1]]<-list(c(0.08,0.2),matrix(c("gamma", "gamma", 1, 1, 0.01, 0.01, 0.05, 0.5),
ncol = 4, nrow = 2))
suspar[[2]]<- NetworkDataSINR$cov
delta1<-list(NULL)
delta1[[1]]<-c(1,1)
delta1[[2]]<-c(2,2)
delta1[[3]]<-matrix(c(4,4,2,2),ncol=2,nrow=2)
spark<-list(1,matrix(c("gamma", 1, 0.01, 0.05), ncol = 4, nrow = 1))

mcmcres5 <- epictmcmc(object = netSINR, datatype = "unknown removal",
nsim = 500, control.sus = suspar, spark.par = spark, delta = delta1)

print(mcmcres5)
summary(mcmcres5)
#plot(mcmcres5, plottype = "parameter")
#plot(mcmcres5, plottype = "inf.times")
#plot(mcmcres5, plottype = "rem.times")
