######################################################################
# FUNCTION: epiplot
# AUTHORS:
#         Waleed Almutiry <walmutir@uoguelph.ca>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# DESCRIPTION:
#
#     To provide different plots summarizing epidemic
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Founepidation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/EpiILMCT package
#
######################################################################

epiplot<- function(type, kerneltype, epidat, location, network = NULL, plottype = NULL, time_index = NULL) {

    k1 <- sum(epidat[,2] != Inf)
    n <- length(epidat[,1])

	if (is.null(time_index)) {
		ss <- c(seq(1, k1))
	} else {
		if (1 %in% time_index) { 
			ss <- sort(time_index)
		} else {
			ss <- c(1, sort(time_index))
		}
	}

	if (length(ss) < 2) {
		a <- 1
		b <- 2
	} else if (length(ss) <= 4 & length(ss) >=  2) {
		a <- 2
		b <- 2
	} else if (length(ss) <= 6 & length(ss) > 4) {
		a <- 2
		b <- 3
	} else if (length(ss) <= 9 & length(ss) > 6) {
		a <- 3
		b <- 3
	} else if (length(ss) <= 12 & length(ss) > 9) {
		a <- 4
		b <- 3
	} else {
		a <- 4
		b <- 4
	}				

if (plottype == "propagation") {    
    if (kerneltype == "distance") {
        if (type == "SIR") {

			par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
            
            plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
            ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
            main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
            points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
            
            for (m in 2:length(ss)) {
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))),
				panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				infectious <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
				removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 4])
				points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
				points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
				points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
            }
            
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
            
        } else if (type == "SINR") {

			par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
            
            plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
            ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
            main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
            points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
 
            for (m in 2:length(ss)) {
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))),
				panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				infectious <- which(epidat[(1:(ss[m]-1)), 6] < epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 6])
				notified   <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
				removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 6])
				points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
				points(location[epidat[notified, 1], 1], location[epidat[notified, 1], 2], col = "yellow", pch = 19)
				points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
				points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
            }
            
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
        }
               
    } else if (kerneltype == "network") {

        if (is.null(network)) {
            stop("Error: the network matrix is missing", call. = FALSE)
        } else if (all(network %in% 0:1)) {

			if (type == "SIR") {
			par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))

				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)
			
				for (i in 1:(n-1)) {
					for (j in (i+1):n) {
						if (network[i,j] == 1) {
							segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
						}
					}
				}
			
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "red", pch = 19)
						
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
					main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")
				
					for (i in 1:(n-1)) {
						for (j in (i+1):n) {
							if (network[i,j] == 1) {
								segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
							}
						}
					}
				
					points(location[, 1], location[, 2], col = "gray", pch = 16)

					for (i in 1:(ss[m]-1)) {
						infectious <- which(epidat[1:i, 4]  <=  epidat[i, 4] & epidat[1:i, 2] > epidat[i, 4])
						removed    <- which(epidat[1:k1, 2] < epidat[i, 4])
						points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
						points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
						for (j in 1:length(infectious)) {
							if (network[epidat[infectious[j], 1], epidat[i, 1]] == 1) {
								arrows(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
								location[epidat[i, 1], 1], location[epidat[i, 1], 2], code  =  2, length  =  0.09)
							}# network if-condition
						}# j for-loop
					}# i for-loop
				
					infectious <- which(epidat[1:(ss[m]-1), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					removed    <- which(epidat[1:k1, 2] <= epidat[ss[m], 4])
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
				
					for (j in 1:length(infectious)) {
						if (network[epidat[infectious[j], 1], epidat[ss[m], 1]] == 1) {
							arrows(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
							location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", code = 2, length = 0.09)
						}# network if-condition
					}# j for-loop
				}
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"), 
			pch = 20, col = c("gray", "blue", "red", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
					
			} else if (type == "SINR") {

			par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))

				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)
			
				for (i in 1:(n-1)) {
					for (j in (i+1):n) {
						if (network[i,j] == 1) {
							segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
						}
					}
				}
			
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "red", pch = 19)
						
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
					main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")
				
					for (i in 1:(n-1)) {
						for (j in (i+1):n) {
							if (network[i,j] == 1) {
								segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
							}
						}
					}
				
					points(location[, 1], location[, 2], col = "gray", pch = 16)

					for (i in 1:(ss[m]-1)) {
						infectious1 <- which(epidat[1:i, 6]  <=  epidat[i, 6] & epidat[1:i, 2] > epidat[i, 6])
						infectious  <- which(epidat[1:i, 6]  <=  epidat[i, 6] & epidat[1:i, 4] > epidat[i, 6])
						notified    <- which(epidat[1:i, 4]  <=  epidat[i, 4] & epidat[1:i, 2] > epidat[i, 4])
						removed     <- which(epidat[1:k1, 2] < epidat[i, 6])
						points(location[epidat[infectious1, 1], 1], location[epidat[infectious1, 1], 2], col = "red", pch = 19)	
						points(location[epidat[notified, 1], 1], location[epidat[notified, 1], 2], col = "yellow", pch = 19)	
						points(location[epidat[removed, 1], 1],location[epidat[removed, 1], 2], col = "green", pch = 19)	
						for (j in 1:length(infectious)) {
							if (network[epidat[infectious[j], 1], epidat[i, 1]] == 1) {
								arrows(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
								location[epidat[i, 1], 1], location[epidat[i, 1], 2], code = 2, length = 0.09)
							}# network if-condition
						}# j for-loop
					}# i for-loop
				
					notified    <- which(epidat[1:(ss[m]-1), 4]  < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					infectious1 <- which(epidat[1:(ss[m]-1), 6]  < epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 4] > epidat[ss[m], 6])
					infectious  <- which(epidat[1:(ss[m]-1), 6]  < epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 6])
					removed     <- which(epidat[1:k1, 2]    < epidat[ss[m], 6])
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)	
					points(location[epidat[infectious1, 1], 1], location[epidat[infectious1, 1], 2], col = "red", pch = 19)	
					points(location[epidat[notified, 1], 1], location[epidat[notified, 1], 2], col = "yellow", pch = 19)	
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)	
				
					for (j in 1:length(infectious)) {
						if (network[epidat[infectious[j], 1], epidat[ss[m], 1]] == 1) {
							arrows(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
							location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", code = 2, length = 0.09)
						}# network if-condition
					}# j for-loop
				}
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			}

		} else {

			if (type == "SIR") {

			par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
			
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[,1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
			
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))),
					panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
					infectious <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 4])
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
				}
			
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			
			} else if (type == "SINR") {

				par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
		
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
 
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))),
					panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
					infectious <- which(epidat[(1:(ss[m]-1)), 6] < epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 6])
					notified   <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 6])
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[notified, 1], 1], location[epidat[notified, 1], 2], col = "yellow", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
				}
			
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			}
						
		}
    } else if (kerneltype == "both") {
        if (is.null(network)) {
            stop("Error: the network matrix is missing", call. = FALSE)
        } else if (all(network %in% 0:1)) {
				
			if (type == "SIR") {

				par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))

				plot(location[, 1], location[, 2], xlim = c(floor(min(location[,1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)
			
				for (i in 1:(n-1)) {
					for (j in (i+1):n) {
						if (network[i,j] == 1) {
							segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
						}
					}
				}
			
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "red", pch = 19)
						
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
					main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")
				
					for (i in 1:(n-1)) {
						for (j in (i+1):n) {
							if (network[i,j] == 1) {
								segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
							}
						}
					}
				
					points(location[, 1], location[, 2], col = "gray", pch = 16)

					for (i in 1:(ss[m]-1)) {
						infectious <- which(epidat[1:i, 4]  <=  epidat[i, 4] & epidat[1:i, 2] > epidat[i, 4])
						removed    <- which(epidat[1:k1, 2] < epidat[i, 4])
						points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
						points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
						for (j in 1:length(infectious)) {
							if (network[epidat[infectious[j], 1],epidat[i, 1]] == 1) {
								segments(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
								location[epidat[i, 1], 1], location[epidat[i, 1], 2])
							}# network if-condition
						}# j for-loop
					}# i for-loop
				
					infectious <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2]>epidat[ss[m], 4])
					removed    <- which(epidat[1:k1, 2] <=  epidat[ss[m], 4])
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
				
					for (j in 1:length(infectious)) {
						if (network[epidat[infectious[j], 1], epidat[ss[m], 1]] == 1) {
							segments(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
							location[epidat[ss[m], 1], 1],location[epidat[ss[m], 1], 2], col = "blue")
						}# network if-condition
					}# j for-loop
				}
				
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
					
			} else if (type == "SINR") {

				par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))

				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)
			
				for (i in 1:(n-1)) {
					for (j in (i+1):n) {
						if (network[i, j] == 1) {
							segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
						}
					}
				}
			
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "red", pch = 19)
						
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
					main =  paste("Infection time (",m,")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")
				
					for (i in 1:(n-1)) {
						for (j in (i+1):n) {
							if (network[i, j] == 1) {
								segments(location[i, 1], location[i, 2], location[j, 1], location[j, 2], col = "lightgray")
							}
						}
					}
				
					points(location[, 1], location[, 2], col = "gray", pch = 16)

					for (i in 1:(ss[m]-1)) {
						infectious <- which(epidat[1:i, 6]  <=  epidat[i, 6] & epidat[1:i, 2] > epidat[i, 6])
						removed    <- which(epidat[1:k1, 2] < epidat[i, 6])
						points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)	
						points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)	
						for (j in 1:length(infectious)) {
							if (network[epidat[infectious[j], 1], epidat[i, 1]] == 1) {
								segments(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
								location[epidat[i, 1], 1], location[epidat[i, 1], 2])
							}# network if-condition
						}# j for-loop
					}# i for-loop
				
					infectious <- which(epidat[1:(ss[m]-1), 6] <=  epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 6])
					removed    <- which(epidat[1:k1, 2] < epidat[ss[m], 6])
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)	
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)	
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)	
				
					for (j in 1:length(infectious)) {
						if (network[epidat[infectious[j], 1],epidat[ss[m], 1]] == 1) {
							segments(location[epidat[infectious[j], 1], 1], location[epidat[infectious[j], 1], 2],
							location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue")
						}# network if-condition
					}# j for-loop
				}
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			}

		} else {

			if (type == "SIR") {
							
				par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
			
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
			
				for (m in 1:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))),
					panel.first = grid(), main = paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
					infectious <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 4])
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
				}
			
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			
			} else if (type == "SINR") {

				par(mar = c(3.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
			
				plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
				ylim = c(floor(min(location[, 2])), ceiling(max(location[, 2]))), panel.first = grid(),
				main = paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
				points(location[epidat[1, 1], 1], location[epidat[1, 1], 2], col = "blue", pch = 16)
 
				for (m in 2:length(ss)) {
					plot(location[, 1], location[, 2], xlim = c(floor(min(location[, 1])), ceiling(max(location[, 1]))),
					ylim = c(floor(min(location[, 2])), ceiling(max(location[,2]))),
					panel.first = grid(), main = paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
					infectious <- which(epidat[(1:(ss[m]-1)), 6] < epidat[ss[m], 6] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 6])
					notified   <- which(epidat[(1:(ss[m]-1)), 4] < epidat[ss[m], 4] & epidat[(1:(ss[m]-1)), 2] > epidat[ss[m], 4])
					removed    <- which(epidat[(1:(ss[m]-1)), 2] < epidat[ss[m], 6])
					points(location[epidat[infectious, 1], 1], location[epidat[infectious, 1], 2], col = "red", pch = 19)
					points(location[epidat[notified, 1], 1], location[epidat[notified, 1], 2], col = "yellow", pch = 19)
					points(location[epidat[removed, 1], 1], location[epidat[removed, 1], 2], col = "green", pch = 19)
					points(location[epidat[ss[m], 1], 1], location[epidat[ss[m], 1], 2], col = "blue", pch = 19)
				}
			
			opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
			mar = c(0, 0, 0, 0), new = TRUE)
			on.exit(par(opar))
			plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
			legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"), 
			pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
			horiz = TRUE, bty = 'n', cex = 0.9)
			}
						
		}
    } else {
		stop("Specify the type of the kernel function as \"distance\", \"network\" or \"both\"", call.  =  FALSE)			
    }# type if-condition

} else if (plottype == "history") {

	if (type == "SIR") {
		par(mfrow = c(2, 2))
		k1 <- sum(epidat[, 2] != Inf)
		plot(density(epidat[1:k1, 4], from = min(epidat[1:k1, 4])), main = "Epidemic curve of \n the infection times", xlab = "Infection times")
		plot(density(epidat[1:k1, 2], from = min(epidat[1:k1, 2])), main = "Epidemic curve of \n the removal times", xlab = "Removal times")
		plot(epidat[1:k1, 4], type = "l",main = "The epidemic time-lines history",ylim = c(min(epidat[1:k1, 4]),max(epidat[1:k1, 2])), ylab = "event times", xlab = "Time points")
		lines(epidat[1:k1, 2], col = "black")
		
		polygon(c(seq(1, k1), rev(seq(1, k1))), c(epidat[1:k1, 4], rev(epidat[1:k1, 2])), col = "red", border = NA)
	
	} else if (type == "SINR") {
		par(mfrow = c(2, 2))
		k1 <- sum(epidat[, 2] != Inf)
		plot(density(epidat[1:k1, 6], from = min(epidat[1:k1, 6])), main = "Epidemic curve of \n the infection times", xlab = "Infection times")
		plot(density(epidat[1:k1, 4], from = min(epidat[1:k1, 4])), main = "Epidemic curve of \n the notification times", xlab = "Notification times")
		plot(density(epidat[1:k1, 2], from = min(epidat[1:k1, 2])), main = "Epidemic curve of \n the removal times", xlab = "Removal times")
		
		plot(epidat[1:k1, 6], type = "l", main = "The epidemic time-line history", ylim = c(min(epidat[1:k1, 6]), max(epidat[1:k1, 2])), ylab = "event times", xlab = "Time points")
		lines(epidat[1:k1, 4], col = "black")
		lines(epidat[1:k1, 2], col = "black")
		polygon(c(seq(1, k1), rev(seq(1, k1))), c(epidat[1:k1, 6], rev(epidat[1:k1, 4])), col = "red", border = NA)
		polygon(c(seq(1, k1), rev(seq(1, k1))), c(epidat[1:k1, 4], rev(epidat[1:k1, 2])), col = "blue", border = NA)

}

}
    
}
