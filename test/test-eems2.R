#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
setwd("~/eems2/test/")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
#plotpath <- 'data/popres/3.5_6/plots_lowerBnd_3.5_upperBnd_6'
#mcmcpath <- 'data/popres/popres_lowerBnd_3.5_upperBnd_6-EEMS2-test-sim'
#plotpath <- 'data/popres/plots_lowerBnd_6_upperBnd_Inf'
#mcmcpath <- 'data/popres/popres_lowerBnd_6_upperBnd_Inf-EEMS2-test-sim'
mcmcpath <- '../data/POBI/out/POBI_10_INF'
plotpath <- '../data/POBI/out/POBI_10_INF'
#plotpath <- 'data/4x5/recent_barrier/6_Inf_plot'
#mcmcpath <- 'data/4x5/recent_barrier/eems_6_Inf-EEMS2-test-sim'
#plotpath <- 'data/4x5/uniform/plot'
#mcmcpath <- 'data/4x5/uniform/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/past_barrier/eems_6_Inf_plot'
#mcmcpath <- 'data/4x5/past_barrier/eems_6_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/past_barrier/eems_2_6_plot'
#mcmcpath <- 'data/4x5/past_barrier/eems_2_6-EEMS2-test-sim' 
#plotpath <- 'data/4x5/high_migration/plot'
#mcmcpath <- 'data/4x5/high_migration/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/plot'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim' 
# FALSE FOR POPRES
longlat <- TRUE


mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     m.standardize = TRUE,
                                                     q.standardize = TRUE,
                                                     m.log10transform = TRUE,
                                                     q.log10transform = TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black", 
                                                     #m.colscale = c(0.09, 0.11), 
                                                     #m.colscale = c(0.000, 0.02),
                                                     add.abline=TRUE,
                                                     #N.colscale = c(8000, 12000),
                                                     plot.height=8, plot.width=10, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))

#eems.voronoi(mcmcpath, plotpath, longlat, add.seeds=FALSE, plot.height=8, plot.width=10)
