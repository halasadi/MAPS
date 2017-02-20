#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
mcmcpath <- 'include_swiss/chain1/popressard_4_10/output'
plotpath <- 'include_swiss/chain1/popressard_4_10/plot'
#mcmcpath <- 'test-df/df_free/popressard_4_Inf/output'
#plotpath <- 'test-df/df_free/popressard_4_Inf/plot'

# FALSE FOR POPRES
longlat <- FALSE


#lwd.map=2

# read the estimates from the previous run 
# important for plotting the plots with 'abs' in the filename
oldcontourpath = 'include_swiss/chain1/popressard_2_4/output/contours.rda'
#oldcontourpath = NA


eems.plots(mcmcpath, plotpath, longlat, add.map=TRUE, add.demes=TRUE, oldcontourpath = oldcontourpath,
                                                     add.grid=FALSE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
                                                     lwd.map=2, col.map = "black", # m.colscale = c(-4, 0),
                                                               add.abline=TRUE, #N.colscale = c(3.5,7),
                                                     plot.height=8, plot.width=14, 
                                                     projection.in = "+proj=longlat +datum=WGS84", scale.by.demes=FALSE)
                                                      #N.colscale = c(2.5, 8), m.colscale = c(-4, 0))
#eems.voronoi(mcmcpath, plotpath, longlat, add.seeds=FALSE, plot.height=8, plot.width=10)
