
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
#if (file.exists("../plotting/rEEMSplots/")) {
#  install.packages("rEEMSplots", repos=NULL, type="source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}
#detach("package:rEEMSplots", unload=TRUE)
#install.packages("../plotting/rEEMSplots", repos=NULL, type="source")
library("rEEMSplots")

## mcmcpath is a list of three output directories; the results will be averaged
plotpath <- 'data/popres/6_Inf/plots_lowerBnd_6_upperBnd_Inf'
mcmcpath <- 'data/popres/popres_lowerBnd_6_upperBnd_Inf-EEMS2-test-sim'
#plotpath <- 'data/3x4_recent_barrier_2/plot'
#mcmcpath <- 'data/3x4_recent_barrier_2/3x4_recent_barrier_2-EEMS2-test-sim'
#plotpath <- 'data/3x4_uniform_nsamp_40/plot'
#mcmcpath <- 'data/3x4_uniform_nsamp_40/3x4_uniform-EEMS2-test-sim' 

# FALSE FOR POPRES
longlat <- FALSE

Shat <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDhatJ.txt")))
Sobs <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDobsJ.txt")))

loc <- as.matrix(read.table(paste0(mcmcpath, "/demes.txt")))
plot(loc[,1], loc[,2], xlab = "lat", ylab = "long")
color_palette <- c("white", "blue", "red")
col_index = c(1, 2,3)

Sobs[Sobs > 0.01] = 3
Sobs[Sobs > 0 & Sobs <= 0.005] = 2
Sobs[Sobs == 0] = 1

for (i in 1:nrow(Sobs)){
  if (i+1 < ncol(Sobs)){
    
    for (j in (i+1):(ncol(Sobs))){
      if (Sobs[i,j] != 1){
        segments(loc[i,1], loc[i,2], loc[j,1], loc[j,2], col = color_palette[col_index[Sobs[i,j]]], lwd = 0.5)
      }
  }
 
  }
}


mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black", m.colscale = c(-0.5, 0.5), add.abline=TRUE,
                                                     N.colscale = c(-0.5, 0.5),
                                                     plot.height=8, plot.width=10, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))


