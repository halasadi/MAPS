#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
#setwd("~/eems2/test/")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
#plotpath <- 'data/4x5/recent_barrier/eems_2_6_plot'
#mcmcpath <- 'data/4x5/recent_barrier/eems_2_6-EEMS2-test-sim'
#plotpath <- 'data/4x5/uniform/plot_4_8'
#mcmcpath <- 'data/4x5/uniform/eems_4_8-EEMS2-test-sim' 
#plotpath <- 'data/test-df/test_mindf_200d/df_3/rep0/plot_4_Inf'
#mcmcpath <- 'data/test-df/test_mindf_200d/df_3/rep0/popressard_4_Inf/output' 
mcmcpath <- 'data/popres/popressard_3_6-EEMS2-test-sim'
plotpath <- 'data/popres/porpessard_3_6'
#plotpath <- 'data/12x8/uniform/plot'
#mcmcpath <- 'data/12x8/uniform/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/recent_barrier/eems_8_Inf_plot'
#mcmcpath <- 'data/4x5/recent_barrier/eems_8_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/eems_4_Inf_plot'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/unevensample_heter_popsizes/eems_4_Inf_plot'
#mcmcpath <- 'data/4x5/unevensample_heter_popsizes/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/past_barrier/eems_6_Inf_plot'
#mcmcpath <- 'data/4x5/past_barrier/eems_6_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/high_migration/plot'
#mcmcpath <- 'data/4x5/high_migration/eems_4_Inf-EEMS2-test-sim' 
#plotpath <- 'data/4x5/heter_popsizes/plot-2'
#mcmcpath <- 'data/4x5/heter_popsizes/eems_4_Inf-EEMS2-test-sim-2' 


# FALSE FOR POPRES
longlat <- FALSE

eems.plots(mcmcpath, plotpath, longlat, add.map=TRUE, add.demes=TRUE,
           add.grid=FALSE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
           lwd.map=2, col.map = "black",
           add.abline=TRUE,
           plot.height=8, plot.width=14, 
           projection.in = "+proj=longlat +datum=WGS84", scale.by.demes=TRUE,
           prob.levels = c(0.95, .99),
           N.colscale = c(2.5, 8), m.colscale = c(-4, 0))

#mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
#                                                     add.map=TRUE, add.demes = TRUE,
#                                                     add.grid=FALSE, add.outline=FALSE, lwd.grid=5, col.grid="black",
#                                                     lwd.map=2, col.map = "black", 
#                                                     add.abline=TRUE,
#                                                     plot.height=8, plot.width=14,
#                                                     projection.in = "+proj=longlat +datum=WGS84"))

#eems.voronoi.samples(mcmcpath, plotpath = paste(plotpath,"-voronoi-diagrams",sep=""),
#            longlat = longlat, post.draws = 100)
