
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
#plotpath <- 'data/sim_output/barrier_L_8e+06_M_0.05_N_10000_nsamp_per_deme_10'
#mcmcpath <- 'data/sim_output/barrier_L_8e+06_M_0.05_N_10000_nsamp_per_deme_10-EEMS2-test-sim'
#plotpath <- 'data/3x4_uniform_nsamp_40/plot'
#mcmcpath <- 'data/3x4_uniform_nsamp_40/3x4_uniform-EEMS2-test-sim' 

# FALSE FOR POPRES
longlat <- FALSE

Shat <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDhatJ.txt")))
Sobs <- as.matrix(read.table(paste0(mcmcpath, "/rdistJtDobsJ.txt")))

loc <- as.matrix(read.table(paste0(mcmcpath, "/demes.txt")))
plot(loc[,1], loc[,2], xlab = "lat", ylab = "long", bty='L')

diag(Sobs) <- NA
labs = matrix(cut(Sobs, breaks = 5, labels = 1:5), nrow = 12, ncol = 12)

colfunc <- colorRampPalette(c("white", "black"))
cols = colfunc(5)

for (i in 1:nrow(Sobs)){
  if (i+1 < ncol(Sobs)){
    
    for (j in (i+1):(ncol(Sobs))){
      if (Sobs[i,j] != 1){
        segments(loc[i,1], loc[i,2], loc[j,1], loc[j,2], col = cols[as.numeric(labs[i,j])], lwd = 1)
      }
  }
 
  }
}

cuts<-levels(cut(Sobs,breaks = 5))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","[",cuts)
par(xpd=TRUE)
legend("topright", cuts, col = cols, pch = 16, cex = 0.6)


mapply(eems.plots, mcmcpath, plotpath, MoreArgs=list(longlat,
                                                     add.map=TRUE, add.demes=TRUE,
                                                     add.grid=TRUE, add.outline=FALSE, lwd.grid=0.3, col.grid="black",
                                                     lwd.map=2, col.map = "black", 
                                                     #m.colscale = c(-1, 1), 
                                                     add.abline=TRUE,
                                                     #N.colscale = c(-1, 1),
                                                     plot.height=8, plot.width=10, 
                                                     projection.in = "+proj=longlat +datum=WGS84"))


