#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
#setwd("~/eems2/test/")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
plotpath <- 'data/4x5/recent_barrier/6_Inf_plot'
mcmcpath <- 'data/4x5/recent_barrier/eems_6_Inf-EEMS2-test-sim'
#plotpath <- 'data/4x5/uniform/plot'
#mcmcpath <- 'data/4x5/uniform/eems_8_Inf-EEMS2-test-sim' 
#plotpath <- 'data/popres/plot-3_6'
#mcmcpath <- 'data/popres/popressard_3_6-EEMS2-test-sim' 
#plotpath <- 'data/12x8/uniform/plot'
#mcmcpath <- 'data/12x8/uniform/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/recent_barrier/eems_2_6_plot'
#mcmcpath <- 'data/4x5/recent_barrier/eems_2_6-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/eems_4_Inf_plot'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/past_barrier/eems_6_Inf_plot'
#mcmcpath <- 'data/4x5/past_barrier/eems_6_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/high_migration/plot'
#mcmcpath <- 'data/4x5/high_migration/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/plot-2'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim-2' 


longlat=TRUE

mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=FALSE, add.demes = FALSE,
                                                     add.grid=TRUE, add.outline=FALSE, add.expression = TRUE, lwd.grid=4, col.grid="black",
                                                     lwd.map=2, col.map = "black", 
                                                     add.abline=TRUE,
                                                     plot.height=8, plot.width=14,
                                                     projection.in = "+proj=longlat +datum=WGS84"))
