
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
#if (file.exists("../plotting/rEEMSplots/")) {
#  install.packages("rEEMSplots", repos=NULL, type="source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}
#detach("package:rEEMSplots", unload=TRUE)
#install.packages("../plotting/rEEMSplots", repos=NULL, type="source")
library("rEEMSplots")

## mcmcpath is a list of three output directories; the results will be averaged
#plotpath <- 'plots'
#mcmcpath <- 'popres-EEMS2-test-sim'
plotpath <- 'data/3x4_uniform/plot'
mcmcpath <- 'data/3x4_uniform/3x4_uniform-EEMS2-test-sim'
longlat <- TRUE

Shat <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDhatJ.txt")))
Sobs <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDobsJ.txt")))
plot(c(Sobs), c(Shat), xlab = "observed sim", ylab = "expected sim", col = "red", lwd = 3)
abline(a=0, b=1, col = "black")

pilogl <- as.matrix(read.table(paste0(mcmcpath, "/mcmcpilogl.txt")))
plot(1:length(pilogl[,1]),  pilogl[,1], xlab = "Iteration (after burn-in and under thinning)", ylab = "prior probability")
plot(1:length(pilogl[,2]),  pilogl[,2], xlab = "Iteration (after burn-in and under thinning)", ylab = "logl")

mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black", m.colscale = c(-0.4, 0.4), q.colscale = c(-0.6, 0.6),
                                                     plot.height=8, plot.width=14, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))


