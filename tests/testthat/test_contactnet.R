library(EpiILMCT)
set.seed(22)

# to generate the XY coordinates of 50 individuals:

loc<- matrix(cbind(runif(50, 0, 10),runif(50, 0, 10)), ncol = 2, nrow = 50)

# Spatial contact network:
# power-law model:
net1<- contactnet(type = "powerlaw", location = loc, beta = 1.5, 
	nu = 0.5)
plot(net1)

# Cauchy model:
net2<- contactnet(type = "Cauchy", location = loc, beta = 0.5)
plot(net2)

# random contact network:
net3<- contactnet(type = "random", num.id = 50, beta = 0.08)
plot(net3)  # the default options in igraph package.
plot(net3, vertex.color = "red", vertex.size = 15, edge.color = "black",
vertex.label.cex = 1, vertex.label.color = "black") 


