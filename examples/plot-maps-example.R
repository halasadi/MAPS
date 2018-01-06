## loading in libraries
library(plotmaps)

mcmcpath <- '/Users/halasadi/eems2/examples/2_6/2_6-MAPS-test-sim/'
outpath  <- '/Users/halasadi/eems2/examples/2_6/'

help(plot_maps)

plot_maps(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
         longlat = TRUE, mcmcpath = mcmcpath, 
         outpath = outpath, width = 10, height = 6)




mcmcpath <- '/Users/halasadi/eems2/examples/6_Inf/6_Inf-MAPS-test-sim/'
outpath  <- '/Users/halasadi/eems2/examples/6_Inf/'

plot_maps(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
         longlat = TRUE, mcmcpath = mcmcpath, 
         outpath = outpath, width = 10, height = 6, plot.difference=TRUE)

