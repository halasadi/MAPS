#detach("package:rEEMSplots2", unload=TRUE)
install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
library(maps)
library(rEEMSplots2)



## mcmcpath is a list of three output directories; the results will be averaged
mcmcpath <- 'midway/vary-grid/100d/popressard_2_Inf/output'
plotpath <- 'midway/vary-grid/100d/popressard_2_Inf/plot'
d <- c(100, 150, 200, 250, 300, 250, 400, 450)
mcmcpath <- rep(NA, length(d))

for (i in 1:length(d)){
  mcmcpath[i] <- paste0('midway/vary-grid/', d[i], 'd/popressard_2_Inf/output')
}

plotpath <- rep(NA, length(d))
for (i in 1:length(d)){
  plotpath[i] <- paste0('midway/overall/r', inds[i], '/popressard_2_Inf/plot')
}

#plotpath <- 'backup/r20/popressard_2_Inf/plot'
#mcmcpath <- 'backup/r20/popressard_2_Inf/output'

# FALSE FOR POPRES
longlat <- FALSE

# total area = 7955710
#N.scalingfactor = totalarea/ndemes (330)
#m.scalignfactor = dx^2
#m.scalingfactor = 39842
#N.scalingfactor = 24108

# 200 demes
#m.scalingfactor = 76000
#N.scalingfactor = 40000


#lwd.map=2

# read the estimates from the previous run 
# important for plotting the plots with 'abs' in the filename
oldcontourpath = NA
#oldcontourpath = 'include_swiss/chain3-3/popressard_2_8/output/contours.rda'

# 117 demes
m.scalingfactor = 115943
N.scalingfactor = 68000

eems.plots(mcmcpath, plotpath[1], longlat, add.map=TRUE, m.scalingfactor = m.scalingfactor, N.scalingfactor=N.scalingfactor, add.demes=TRUE, oldcontourpath = oldcontourpath,
                                                     add.grid=FALSE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
                                                     lwd.map=1, col.map = "black", # m.colscale = c(-4, 0),
                                                               add.abline=TRUE, #N.colscale = c(3.5,7),
                                                     plot.height=8, plot.width=14, 
                                                     projection.in = "+proj=longlat +datum=WGS84", scale.by.demes=FALSE)
                                                      #N.colscale = c(2.5, 8), m.colscale = c(-4, 0))
#eems.voronoi(mcmcpath, plotpath, longlat, add.seeds=FALSE, plot.height=8, plot.width=10)


eems.plot.difference(mcmcpath = mcmcpath, contourpath1 = "midway/partition_time/r1/popressard_2_8/output/contours.rds",  
                     contourpath2 = "midway/partition_time/r2/popressard_8_Inf/output/contours.rds", plotpath= "diff", 
                     longlat, add.map=TRUE, m.scalingfactor = m.scalingfactor, N.scalingfactor=N.scalingfactor, add.demes=TRUE,
           add.grid=FALSE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
           lwd.map=1, col.map = "black", # m.colscale = c(-4, 0),
           add.abline=TRUE, #N.colscale = c(3.5,7),
           plot.height=8, plot.width=14, 
           projection.in = "+proj=longlat +datum=WGS84", scale.by.demes=FALSE)

# convert -delay 20 -depth 8 tiles/plot-mv* mvoronoi.gif
# convert mvoronoi.gif -coalesce -scale 700x525 -fuzz 2% +dither -remap mvoronoi.gif[0] -layers Optimize mvoronoi1.gif
# eems.voronoi.samples(mcmcpath, plotpath = "through-time-2/r1/popressard_8_Inf/tiles/plot", longlat, post.draws = 100)
