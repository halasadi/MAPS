
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (file.exists("../plotting/rEEMSplots/")) {
  install.packages("rEEMSplots", repos=NULL, type="source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

#library(rEEMSplots)


## mcmcpath is a list of three output directories; the results will be averaged
#mcmcpath <- '../data/l_1e8_N_1000_m_0.01_mu_1.25e-8_r_1e-8_ts_100_12D-EEMS2-test-simno1'
plotpath <- '../data/plots'
mcmcpath <- '../data/run_2/l_1.3e8_N_2000_m_1e-2_mu_1.25e-8_r_1e-8_ts_100_4D-EEMS2-test-sim'
#plotpath <- '../data/plot/l_1.3e8_N_2000_m_1e-2_mu_1.25e-8_r_1e-8_ts_800_4D'
longlat <- TRUE


eems.plots(mcmcpath,plotpath,longlat,
           add.map=FALSE,
           add.grid=TRUE,
           add.samples=TRUE, plot.height=8, plot.width=8)
