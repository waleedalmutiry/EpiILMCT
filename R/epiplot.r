###########################################################################
# FUNCTION: S3 method: plot.datagen
# AUTHORS:
#         Waleed Almutiry <walmutir@uoguelph.ca>,
#         Vineetha Warriyar. K. V. <vineethawarriyar.kod@ucalgary.ca> and
#         Rob Deardon <robert.deardon@ucalgary.ca>
#
# DESCRIPTION:
#
#     To provide different plots summarizing epidemic using the S3 plot
#     function. The object x has to be of <datagen> class.
#     The argument plottype has two options: "history" and "propagation".
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
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
###########################################################################

plot.datagen<- function(x, plottype, time.index = NULL, ...) {

    if (is(x, "datagen")) {

        k1 <- sum(x$epidat[,2] != Inf)
        n <- length(x$epidat[,1])

        if (is.null(time.index)) {
            ss <- c(seq(1, k1))
        } else {
            if (1 %in% time.index) {
                ss <- sort(time.index)
            } else {
                ss <- c(1, sort(time.index))
            }
        }

        if (length(ss) <= 2) {
            a <- 1
            b <- 2
        } else if (length(ss) <= 4 & length(ss) >  2) {
            a <- 2
            b <- 2
        } else if (length(ss) <= 6 & length(ss) > 4) {
            a <- 2
            b <- 3
        } else if (length(ss) > 6) {
            a <- 3
            b <- 3
        }

        if (plottype == "propagation") {
            if (x$kerneltype == "distance") {

                op1 <- par(no.readonly = TRUE)

                if (x$type == "SIR") {

                    op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                    on.exit(par(op))

                    plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                    ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                    main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                    points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                    for (m in 2:length(ss)) {
                        plot(x$location[, 1], x$location[, 2],
                        xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))),
                        panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        infectious <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                        removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 4])
                        points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                        points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                        points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                        if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                            if (a+b == 3) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                pch = 20, col = c("gray", "blue", "red", "green"),
                                horiz = TRUE, bty = 'n', cex = 0.8)
                            } else if (a+b == 4) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                pch = 20, col = c("gray", "blue", "red", "green"),
                                horiz = TRUE, bty = 'n', cex = 1.0)
                            } else if (a+b > 4) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                pch = 20, col = c("gray", "blue", "red", "green"),
                                horiz = TRUE, bty = 'n', cex = 1.3)
                            }
                            op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                            on.exit(par(op), add = TRUE)
                        }

                    }

                } else if (x$type == "SINR") {

                    op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                    on.exit(par(op))

                    plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                    ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                    main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                    points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                    for (m in 2:length(ss)) {
                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))),
                        panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        infectious <- which(x$epidat[(1:(ss[m]-1)), 6] < x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 6])
                        notified   <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                        removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 6])
                        points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                        points(x$location[x$epidat[notified, 1], 1], x$location[x$epidat[notified, 1], 2], col = "yellow", pch = 19)
                        points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                        points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                        if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                            opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
                            on.exit(par(opar), add = TRUE)
                            plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                            if (a+b == 3) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                horiz = TRUE, bty = 'n', cex = 0.8)
                            } else if (a+b == 4) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                horiz = TRUE, bty = 'n', cex = 1.0)
                            } else if (a+b > 4) {
                                legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                horiz = TRUE, bty = 'n', cex = 1.3)
                            }

                            op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                            on.exit(par(op), add = TRUE)
                        }

                    }

                }

                on.exit(par(op1))

            } else if (x$kerneltype == "network") {

                if (is.null(x$network)) {
                    stop("the network matrix is missing", call. = FALSE)
                } else if (all(x$network %in% 0:1)) {

                    op1 <- par(no.readonly = TRUE)

                    if (x$type == "SIR") {
                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)

                        for (i in 1:(n-1)) {
                            for (j in (i+1):n) {
                                if (x$network[i,j] == 1) {
                                    segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                }
                            }
                        }

                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "red", pch = 19)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                            main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")

                            for (i in 1:(n-1)) {
                                for (j in (i+1):n) {
                                    if (x$network[i,j] == 1) {
                                        segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                    }
                                }
                            }

                            points(x$location[, 1], x$location[, 2], col = "gray", pch = 16)

                            for (i in 1:(ss[m]-1)) {
                                infectious <- which(x$epidat[1:i, 4]  <=  x$epidat[i, 4] & x$epidat[1:i, 2] > x$epidat[i, 4])
                                removed    <- which(x$epidat[1:k1, 2] < x$epidat[i, 4])
                                points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                                points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                                for (j in 1:length(infectious)) {
                                    if (x$network[x$epidat[infectious[j], 1], x$epidat[i, 1]] == 1) {
                                        arrows(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                        x$location[x$epidat[i, 1], 1], x$location[x$epidat[i, 1], 2], code  =  2, length  =  0.09)
                                    }# network if-condition
                                }# j for-loop
                            }# i for-loop

                            infectious <- which(x$epidat[1:(ss[m]-1), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[1:k1, 2] <= x$epidat[ss[m], 4])
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)

                            for (j in 1:length(infectious)) {
                                if (x$network[x$epidat[infectious[j], 1], x$epidat[ss[m], 1]] == 1) {
                                    arrows(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                    x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", code = 2, length = 0.09)
                                }# network if-condition
                            }# j for-loop

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious","removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }

                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    } else if (x$type == "SINR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)

                        for (i in 1:(n-1)) {
                            for (j in (i+1):n) {
                                if (x$network[i,j] == 1) {
                                    segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                }
                            }
                        }

                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "red", pch = 19)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                            main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")

                            for (i in 1:(n-1)) {
                                for (j in (i+1):n) {
                                    if (x$network[i,j] == 1) {
                                        segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                    }
                                }
                            }

                            points(x$location[, 1], x$location[, 2], col = "gray", pch = 16)

                            for (i in 1:(ss[m]-1)) {
                                infectious1 <- which(x$epidat[1:i, 6]  <=  x$epidat[i, 6] & x$epidat[1:i, 2] > x$epidat[i, 6])
                                infectious  <- which(x$epidat[1:i, 6]  <=  x$epidat[i, 6] & x$epidat[1:i, 4] > x$epidat[i, 6])
                                notified    <- which(x$epidat[1:i, 4]  <=  x$epidat[i, 4] & x$epidat[1:i, 2] > x$epidat[i, 4])
                                removed     <- which(x$epidat[1:k1, 2] < x$epidat[i, 6])
                                points(x$location[x$epidat[infectious1, 1], 1], x$location[x$epidat[infectious1, 1], 2], col = "red", pch = 19)
                                points(x$location[x$epidat[notified, 1], 1], x$location[x$epidat[notified, 1], 2], col = "yellow", pch = 19)
                                points(x$location[x$epidat[removed, 1], 1],x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                                for (j in 1:length(infectious)) {
                                    if (x$network[x$epidat[infectious[j], 1], x$epidat[i, 1]] == 1) {
                                        arrows(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                        x$location[x$epidat[i, 1], 1], x$location[x$epidat[i, 1], 2], code = 2, length = 0.09)
                                    }# network if-condition
                                }# j for-loop
                            }# i for-loop

                            notified    <- which(x$epidat[1:(ss[m]-1), 4]  < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            infectious1 <- which(x$epidat[1:(ss[m]-1), 6]  < x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 4] > x$epidat[ss[m], 6])
                            infectious  <- which(x$epidat[1:(ss[m]-1), 6]  < x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 6])
                            removed     <- which(x$epidat[1:k1, 2]    < x$epidat[ss[m], 6])
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)
                            points(x$location[x$epidat[infectious1, 1], 1], x$location[x$epidat[infectious1, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[notified, 1], 1], x$location[x$epidat[notified, 1], 2], col = "yellow", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)

                            for (j in 1:length(infectious)) {
                                if (x$network[x$epidat[infectious[j], 1], x$epidat[ss[m], 1]] == 1) {
                                    arrows(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                    x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", code = 2, length = 0.09)
                                }# network if-condition
                            }# j for-loop

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }

                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }
                    }

                    on.exit(par(op1))

                } else {

                    op1 <- par(no.readonly = TRUE)

                    if (x$type == "SIR") {

                    op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                    on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[,1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))),
                            panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                            infectious <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 4])
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    } else if (x$type == "SINR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))),
                            panel.first = grid(), main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                            infectious <- which(x$epidat[(1:(ss[m]-1)), 6] < x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 6])
                            notified   <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 6])
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[notified, 1], 1], x$location[x$epidat[notified, 1], 2], col = "yellow", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    }

                    on.exit(par(op1))

                }
            } else if (x$kerneltype == "both") {
                if (is.null(x$network)) {
                    stop("the network matrix is missing", call. = FALSE)
                } else if (all(x$network %in% 0:1)) {

                    op1 <- par(no.readonly = TRUE)

                    if (x$type == "SIR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[,1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)

                        for (i in 1:(n-1)) {
                            for (j in (i+1):n) {
                                if (x$network[i,j] == 1) {
                                    segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                }
                            }
                        }

                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "red", pch = 19)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                            main =  paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")

                            for (i in 1:(n-1)) {
                                for (j in (i+1):n) {
                                    if (x$network[i,j] == 1) {
                                        segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                    }
                                }
                            }

                            points(x$location[, 1], x$location[, 2], col = "gray", pch = 16)

                            for (i in 1:(ss[m]-1)) {
                                infectious <- which(x$epidat[1:i, 4]  <=  x$epidat[i, 4] & x$epidat[1:i, 2] > x$epidat[i, 4])
                                removed    <- which(x$epidat[1:k1, 2] < x$epidat[i, 4])
                                points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                                points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                                for (j in 1:length(infectious)) {
                                    if (x$network[x$epidat[infectious[j], 1],x$epidat[i, 1]] == 1) {
                                        segments(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                        x$location[x$epidat[i, 1], 1], x$location[x$epidat[i, 1], 2])
                                    }# network if-condition
                                }# j for-loop
                            }# i for-loop

                            infectious <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2]>x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[1:k1, 2] <=  x$epidat[ss[m], 4])
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)

                            for (j in 1:length(infectious)) {
                                if (x$network[x$epidat[infectious[j], 1], x$epidat[ss[m], 1]] == 1) {
                                    segments(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                    x$location[x$epidat[ss[m], 1], 1],x$location[x$epidat[ss[m], 1], 2], col = "blue")
                                }# network if-condition
                            }# j for-loop

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    } else if (x$type == "SINR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "dark gray", pch = 19)

                        for (i in 1:(n-1)) {
                            for (j in (i+1):n) {
                                if (x$network[i, j] == 1) {
                                    segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                }
                            }
                        }

                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "red", pch = 19)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                            main =  paste("Infection time (",m,")"), xlab = "x", ylab = "y", pch = 19, col = "dark gray")

                            for (i in 1:(n-1)) {
                                for (j in (i+1):n) {
                                    if (x$network[i, j] == 1) {
                                        segments(x$location[i, 1], x$location[i, 2], x$location[j, 1], x$location[j, 2], col = "lightgray")
                                    }
                                }
                            }

                            points(x$location[, 1], x$location[, 2], col = "gray", pch = 16)

                            for (i in 1:(ss[m]-1)) {
                                infectious <- which(x$epidat[1:i, 6]  <=  x$epidat[i, 6] & x$epidat[1:i, 2] > x$epidat[i, 6])
                                removed    <- which(x$epidat[1:k1, 2] < x$epidat[i, 6])
                                points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                                points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                                for (j in 1:length(infectious)) {
                                    if (x$network[x$epidat[infectious[j], 1], x$epidat[i, 1]] == 1) {
                                        segments(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                        x$location[x$epidat[i, 1], 1], x$location[x$epidat[i, 1], 2])
                                    }# network if-condition
                                }# j for-loop
                            }# i for-loop

                            infectious <- which(x$epidat[1:(ss[m]-1), 6] <=  x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 6])
                            removed    <- which(x$epidat[1:k1, 2] < x$epidat[ss[m], 6])
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)

                            for (j in 1:length(infectious)) {
                                if (x$network[x$epidat[infectious[j], 1],x$epidat[ss[m], 1]] == 1) {
                                    segments(x$location[x$epidat[infectious[j], 1], 1], x$location[x$epidat[infectious[j], 1], 2],
                                    x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue")
                                }# network if-condition
                            }# j for-loop

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    }

                    on.exit(par(op1))

                } else {

                    op1 <- par(no.readonly = TRUE)

                    if (x$type == "SIR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main =  paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                        for (m in 1:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))),
                            panel.first = grid(), main = paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                            infectious <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 4])
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    } else if (x$type == "SINR") {

                        op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                        on.exit(par(op))

                        plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                        ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[, 2]))), panel.first = grid(),
                        main = paste("Infection time (",1,")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                        points(x$location[x$epidat[1, 1], 1], x$location[x$epidat[1, 1], 2], col = "blue", pch = 16)

                        for (m in 2:length(ss)) {
                            plot(x$location[, 1], x$location[, 2], xlim = c(floor(min(x$location[, 1])), ceiling(max(x$location[, 1]))),
                            ylim = c(floor(min(x$location[, 2])), ceiling(max(x$location[,2]))),
                            panel.first = grid(), main = paste("Infection time (",ss[m],")"), xlab = "x", ylab = "y", col = "light gray", pch = 19)
                            infectious <- which(x$epidat[(1:(ss[m]-1)), 6] < x$epidat[ss[m], 6] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 6])
                            notified   <- which(x$epidat[(1:(ss[m]-1)), 4] < x$epidat[ss[m], 4] & x$epidat[(1:(ss[m]-1)), 2] > x$epidat[ss[m], 4])
                            removed    <- which(x$epidat[(1:(ss[m]-1)), 2] < x$epidat[ss[m], 6])
                            points(x$location[x$epidat[infectious, 1], 1], x$location[x$epidat[infectious, 1], 2], col = "red", pch = 19)
                            points(x$location[x$epidat[notified, 1], 1], x$location[x$epidat[notified, 1], 2], col = "yellow", pch = 19)
                            points(x$location[x$epidat[removed, 1], 1], x$location[x$epidat[removed, 1], 2], col = "green", pch = 19)
                            points(x$location[x$epidat[ss[m], 1], 1], x$location[x$epidat[ss[m], 1], 2], col = "blue", pch = 19)

                            if (any(m == seq.int(9, n, 9)) | (m == length(ss))) {
                                opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                                mar = c(0, 0, 0, 0), new = TRUE)
                                on.exit(par(opar), add = TRUE)
                                plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
                                if (a+b == 3) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 0.8)
                                } else if (a+b == 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.0)
                                } else if (a+b > 4) {
                                    legend("bottom", legend = c("susceptible", "newly-infected", "infectious", "notified", "removed"),
                                    pch = 20, col = c("gray", "blue", "red", "yellow", "green"),
                                    horiz = TRUE, bty = 'n', cex = 1.3)
                                }
                                op <- par(mar = c(5.9, 4.0, 1.5, 0.5), omi = c(0.2, 0.05, 0.15, 0.15), mfrow = c(a, b))
                                on.exit(par(op), add = TRUE)
                            }

                        }

                    }

                    on.exit(par(op1))

                }
            } else {
                stop("Specify the type of the kernel function as \"distance\", \"network\" or \"both\"", call.  =  FALSE)
            }# type if-condition

        } else if (plottype == "history") {

            op1 <- par(no.readonly = TRUE)

            if (x$type == "SIR") {
                op <- par(mar = c(5.1, 4.1, 4.1, 2.1), omi = c(0, 0, 0, 0), mfrow = c(2, 2))
                on.exit(par(op))
                k1 <- sum(x$epidat[, 2] != Inf)
                plot(density(x$epidat[1:k1, 4], from = min(x$epidat[1:k1, 4])), main = "Epidemic curve of \n the infection times", xlab = "Infection times")
                plot(density(x$epidat[1:k1, 2], from = min(x$epidat[1:k1, 2])), main = "Epidemic curve of \n the removal times", xlab = "Removal times")
                plot(x$epidat[1:k1, 4], type = "l",main = "The epidemic time-lines history",ylim = c(min(x$epidat[1:k1, 4]),max(x$epidat[1:k1, 2])), ylab = "event times", xlab = "Time points")
                lines(x$epidat[1:k1, 2], col = "black")

                polygon(c(seq(1, k1), rev(seq(1, k1))), c(x$epidat[1:k1, 4], rev(x$epidat[1:k1, 2])), col = "red", border = NA)

            } else if (x$type == "SINR") {
                op <- par(mar = c(5.1, 4.1, 4.1, 2.1), omi = c(0, 0, 0, 0), mfrow = c(2, 2))
                on.exit(par(op))
                k1 <- sum(x$epidat[, 2] != Inf)
                plot(density(x$epidat[1:k1, 6], from = min(x$epidat[1:k1, 6])), main = "Epidemic curve of \n the infection times", xlab = "Infection times")
                plot(density(x$epidat[1:k1, 4], from = min(x$epidat[1:k1, 4])), main = "Epidemic curve of \n the notification times", xlab = "Notification times")
                plot(density(x$epidat[1:k1, 2], from = min(x$epidat[1:k1, 2])), main = "Epidemic curve of \n the removal times", xlab = "Removal times")

                plot(x$epidat[1:k1, 6], type = "l", main = "The epidemic time-line history", ylim = c(min(x$epidat[1:k1, 6]), max(x$epidat[1:k1, 2])), ylab = "event times", xlab = "Time points")
                lines(x$epidat[1:k1, 4], col = "black")
                lines(x$epidat[1:k1, 2], col = "black")
                polygon(c(seq(1, k1), rev(seq(1, k1))), c(x$epidat[1:k1, 6], rev(x$epidat[1:k1, 4])), col = "red", border = NA)
                polygon(c(seq(1, k1), rev(seq(1, k1))), c(x$epidat[1:k1, 4], rev(x$epidat[1:k1, 2])), col = "blue", border = NA)

            }

            on.exit(par(op1))

        }

    } else {
        stop("The class of x must be a class of datagen", call. = FALSE)
    }
}
