
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
#if (file.exists("../plotting/rEEMSplots/")) {
#  install.packages("rEEMSplots", repos=NULL, type="source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}

library("rEEMSplots")

## mcmcpath is a list of three output directories; the results will be averaged
plotpath <- 'plots'
mcmcpath <- 'popres-EEMS2-test-sim'
longlat <- FALSE


eems.plots(mcmcpath,plotpath,longlat,
           add.map=TRUE, lwd.map = 1, col.map = "black",
           add.grid=TRUE, add.demes = TRUE, projection.in = '+proj=longlat +datum=WGS84', plot.height = 5, plot.width = 7)

