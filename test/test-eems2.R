#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
#setwd("~/eems2/test/")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
#plotpath <- 'data/4x5/4x5_two_locations/4_Inf_plot'
#mcmcpath <- 'data/4x5/4x5_two_locations/eems_4_Inf-EEMS2-test-sim'
#plotpath <- 'data/4x5/uniform/plot-6-Inf'
#mcmcpath <- 'data/4x5/uniform/eems_6_Inf-EEMS2-test-sim' 
#plotpath <- 'data/popres/plot-3_6'
#mcmcpath <- 'data/popres/popressard_3_6-EEMS2-test-sim' 
#plotpath <- 'data/12x8/uniform/plot'
#mcmcpath <- 'data/12x8/uniform/eems_4_Inf-EEMS2-test-sim' 
plotpath <- 'data/4x5/recent_barrier/eems_6_Inf_plot-prevpath'
mcmcpath <- 'data/4x5/recent_barrier/eems_6_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes_unevensampling2/eems_4_Inf_plot'
#mcmcpath <- 'data/4x5/heter_popsizes_unevensampling2/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/past_barrier/eems_2_4_plot'
#mcmcpath <- 'data/4x5/past_barrier/eems_2_4-EEMS2-test-sim' 
#plotpath <- 'data/4x5/high_migration/plot'
#mcmcpath <- 'data/4x5/high_migration/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/plot-2'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim-2' 

#q <- 10^(read.table('data/4x5/heter_popsizes_unevensampling2/eems_4_Inf-EEMS2-test-sim/log10qMeanRates.txt'))
#N <- (1/(2*q))
#N.true <- c(1e4,1e4,1e2,1e2,1e4,1e4,1e4,1e2,1e2,1e4,1e4,1e4,1e2,1e2,1e4,1e4,1e4,1e2,1e2,1e4)
#map <- mean(abs(N-N.true)/N.true) * 100

m <- 10^(read.table('data/4x5/uniform/eems_2_6-EEMS2-test-sim/log10mMeanRates.txt'))
mean(abs(m-0.01)/0.01)
longlat=TRUE
oldcontourpath = NA

mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=FALSE, add.demes = TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=4, col.grid="black",
                                                     lwd.map=2, col.map = "black", 
                                                     add.abline=TRUE, scale.by.demes = FALSE,
                                                     plot.height=8, plot.width=14,
                                                     projection.in = "+proj=longlat +datum=WGS84"))
