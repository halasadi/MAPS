
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
longlat <- FALSE

mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black",
                                                     plot.height=8, plot.width=14, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))

