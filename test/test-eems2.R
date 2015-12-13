
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
#if (file.exists("../plotting/rEEMSplots/")) {
#  install.packages("rEEMSplots", repos=NULL, type="source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}
#detach("package:rEEMSplots2", unload=TRUE)
#install.packages("../plotting/rEEMSplots2", repos=NULL, type="source")
library(rEEMSplots2)

## mcmcpath is a list of three output directories; the results will be averaged
#plotpath <- 'data/popres/3.5_6/plots_lowerBnd_3.5_upperBnd_6'
#mcmcpath <- 'data/popres/popres_lowerBnd_3.5_upperBnd_6-EEMS2-test-sim'
#plotpath <- 'data/popres/6_Inf/plots_lowerBnd_6_upperBnd_Inf'
#plotpath <- 'data/3x4/recent_barrier_10gen/8_Inf_plot'
#mcmcpath <- 'data/3x4/recent_barrier_10gen/eems_8_Inf-EEMS2-test-sim'
#plotpath <- 'data/3x4/past_barrier_40gen/3_5_plot'
#mcmcpath <- 'data/3x4/past_barrier_40gen/eems_3_5-EEMS2-test-sim' 
#plotpath <- 'data/3x4/3x4_uniform/plot'
#mcmcpath <- 'data/3x4/3x4_uniform/3x4_uniform-EEMS2-test-sim' 
#plotpath  <- 'data/3x4/barrier_40gen_heter_popsizes/4_Inf_plot'
#mcmcpath  <- 'data/3x4/barrier_40gen_heter_popsizes/eems_4_Inf-EEMS2-test-sim'
plotpath <- 'data/4x5/uniform/plot'
mcmcpath <- 'data/4x5/uniform/eems_4_Inf-EEMS2-test-sim' 
# FALSE FOR POPRES
longlat <- TRUE

#Shat <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDhatJ.txt")))
#Sobs <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDobsJ.txt")))

#loc <- as.matrix(read.table("/Users/halasadi/eems2/test/data/3x4/barrier_40gen_heter_popsizes/eems_4_Inf.demes"))
#plot(loc[,1], loc[,2], xlab = "lat", ylab = "long", bty='L')
#edges <- as.matrix(read.table("/Users/halasadi/eems2/test/data/popres/8_Inf/popres_lowerBnd_8_upperBnd_Inf.edges"))
#segments(loc[edges[,1],1], loc[edges[,1],2], loc[edges[,2],1], loc[edges[,2],2])

#diag(Sobs) <- NA
#labs = matrix(cut(Sobs, breaks = 5, labels = 1:5), nrow = 12, ncol = 12)
#colfunc <- colorRampPalette(c("white", "black"))
#cols = colfunc(5)
#for (i in 1:nrow(Sobs)){
#  if (i+1 < ncol(Sobs)){
#    
#    for (j in (i+1):(ncol(Sobs))){
#      if (Sobs[i,j] != 1){
#        segments(loc[i,1], loc[i,2], loc[j,1], loc[j,2], col = cols[as.numeric(labs[i,j])], lwd = 1)
#      }
#  }
# 
#  }
#}

#cuts<-levels(cut(Sobs,breaks = 5))
#cuts<-gsub(","," - ",cuts)
#cuts<-gsub("\\(","[",cuts)
#par(xpd=TRUE)
#legend("topright", cuts, col = cols, pch = 16, cex = 0.6)


mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     m.standardize = FALSE,
                                                     q.standardize = FALSE,
                                                     m.log10transform = FALSE,
                                                     q.log10transform = FALSE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black", 
                                                     #m.colscale = c(0.005, 0.015), 
                                                     #m.colscale = c(0.000, 0.02),
                                                     add.abline=TRUE,
                                                     #N.colscale = c(8000, 12000),
                                                     plot.height=8, plot.width=10, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))

#eems.voronoi(mcmcpath, plotpath, longlat, add.seeds=FALSE, plot.height=8, plot.width=10)
