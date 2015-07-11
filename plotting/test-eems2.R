
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
#if (file.exists("./rEEMSplots/")) {
#  install.packages("rEEMSplots", repos=NULL, type="source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}

library(rEEMSplots)


## mcmcpath is a list of three output directories; the results will be averaged
plotpath <- '../data/run_3/plots'
mcmcpath <- '../data/run_3/l_1.3e8_N_2000_m_0.01_mu_1.25e-8_r_1e-8_ts_20_4D-EEMS2-test-sim'
longlat <- TRUE


eems.plots(mcmcpath,plotpath,longlat,
           add.map=FALSE,
           add.grid=TRUE, add.demes = TRUE,plot.height=8, plot.width=8)
