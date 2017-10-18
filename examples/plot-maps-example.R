## Install plotting ##
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")



## loading in libraries

library(rEEMSplots2)
library(maps)
library(geosphere)


## Here we consider the simulation scenario of a recent barrier



### 2-6 cM ###

plotpath <- '2_6/plot_2_6'
mcmcpath <- '2_6/2_6-MAPS-test-sim/'
eems.plots(mcmcpath, plotpath, longlat = TRUE, add.map=FALSE, scale.by.area = FALSE, add.demes=TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
                                                     lwd.map=1, col.map = "black",
                                                               add.abline=TRUE,
                                                     plot.height=8, plot.width=14, 
                                                     projection.in = "+proj=longlat +datum=WGS84")


### 6-Inf cM ####

plotpath <- '6_Inf/plot_6_Inf'
mcmcpath <- '6_Inf/6_Inf-MAPS-test-sim/'
eems.plots(mcmcpath, plotpath, longlat = TRUE, add.map=FALSE, scale.by.area = FALSE, add.demes=TRUE,
           add.grid=TRUE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
           lwd.map=1, col.map = "black",
           add.abline=TRUE,
           plot.height=8, plot.width=14, 
           projection.in = "+proj=longlat +datum=WGS84")


### plotting the difference between the contour maps ###

eems.plot.difference(mcmcpath = mcmcpath, contourpath1 = "2_6/2_6-MAPS-test-sim/contours.rds",  
                     contourpath2 = "6_Inf/6_Inf-MAPS-test-sim/contours.rds", plotpath= "6_Inf/diff", 
                     longlat, add.map=FALSE, add.demes=TRUE,
           add.grid=FALSE, add.outline=FALSE, lwd.grid=0.5, col.grid="black",
           lwd.map=1, col.map = "black", 
           add.abline=TRUE, 
           plot.height=8, plot.width=14, 
           projection.in = "+proj=longlat +datum=WGS84")

