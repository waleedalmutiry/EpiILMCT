library(EpiILMCT)
set.seed(22)
loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5,
    nu = 0.5)

#################################
#################################
#################################
########       SIR       ########
#################################
#################################
#################################

#################################
########  spatial ILMs:  ########
#################################

# first: with power-law kernel:

cov1 <- cbind(runif(50, 0, 50), runif(50, 0, 5))
cov2 <- cbind(runif(50, 0, 50), runif(50, 0, 5))

# with specifying all the options of susceptibility and transmissibility functions:

epi.dist.po1 <- datagen(type = "SIR", kerneltype = "distance", kernelmatrix = loc,
                initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
                distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                kernel.par = 8, delta = c(6,2), suscov = cov1, transcov = cov2)

epi.dist.po1$epidat

loglikelihoodepiILM(object = epi.dist.po1, distancekernel = "powerlaw",
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)),
kernel.par = 8, delta = c(6,2))

epi.dist.po2 <- datagen(type = "SIR", kerneltype = "network", kernelmatrix = net1$contact.network,
                initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
                suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                delta = c(6,2), suscov = cov1, transcov = cov2)

epi.dist.po2$epidat

loglikelihoodepiILM(object = epi.dist.po2,
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)),
delta = c(6,2))

epi.dist.po3 <- datagen(type = "SIR", kerneltype = "both", kernelmatrix = list(loc,net1$contact.network),
                initialepi = matrix(c(13, 2, 2, 0), ncol = 4, nrow = 1), tmax = 4,
                distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                kernel.par = c(8,0.3), delta = c(6,2), suscov = cov1, transcov = cov2)

epi.dist.po3$epidat

loglikelihoodepiILM(object = epi.dist.po3, distancekernel = "powerlaw",
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)),
kernel.par = c(8,0.3), delta = c(6,2))

data(NetworkDataSINR)
names(NetworkDataSINR)


netSIR<-as.epidat(type = "SIR", kerneltype = "network", inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "powerlaw")

netSIR<-as.epidat(type = "SIR", kerneltype = "distance", inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "Cauchy")

#################################
#################################
#################################
########       SINR      ########
#################################
#################################
#################################

#################################
########  spatial ILMs:  ########
#################################

# first: with power-law kernel:

set.seed(22)

cov1 <- cbind(runif(50, 0, 50), runif(50, 0, 5))
cov2 <- cbind(runif(50, 0, 50), runif(50, 0, 5))

# with specifying all the options of susceptibility and transmissibility functions:

epi.dist.po1 <- datagen(type = "SINR", kerneltype = "distance", kernelmatrix = loc,
                initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 4,
                distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                kernel.par = 8, delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE),
                suscov = cov1, transcov = cov2)

epi.dist.po1$epidat

loglikelihoodepiILM(object = epi.dist.po1, distancekernel = "powerlaw",
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)),
kernel.par = 8, delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE))


epi.dist.po2 <- datagen(type = "SINR", kerneltype = "network", kernelmatrix = net1,
                initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 2,
                suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE),
                suscov = cov1, transcov = cov2)
epi.dist.po2$epidat

loglikelihoodepiILM(object = epi.dist.po2, distancekernel = "powerlaw",
control.sus = list(c(0.01,2), cov1, c(0.5,0.7)), control.trans = list(c(0.03,0.2), cov2, c(0.7,1.3)),
delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE))

epi.dist.po3 <- datagen(type = "SINR", kerneltype = "both", kernelmatrix = list(loc,net1$contact.network),
                initialepi = matrix(c(13, 2, 1, 1, 1, 0), ncol = 6, nrow = 1), tmax = 2,
                distancekernel = "powerlaw", suspar = c(0.01, 2), transpar = c(0.03,0.2),
                powersus = c(0.5, 0.7), powertrans = c(0.7, 1.3),
                kernel.par = c(8,0.3), delta = matrix(c(1,2,6,2), ncol = 2, byrow = TRUE),
                suscov = cov1, transcov = cov2)

epi.dist.po3$epidat

data(NetworkDataSINR)
names(NetworkDataSINR)


netSINR<-as.epidat(type = "SINR", kerneltype = "network", incub.time = NetworkDataSINR$epi[,4], inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "powerlaw")


netSINR<-as.epidat(type = "SINR", kerneltype = "distance", incub.time = NetworkDataSINR$epi[,4], inf.time = NetworkDataSINR$epi[,6], rem.time = NetworkDataSINR$epi[,2], id.individual = NetworkDataSINR$epi[,1], location  = NetworkDataSINR$loc, network = NetworkDataSINR$net, network.type = "Cauchy")
