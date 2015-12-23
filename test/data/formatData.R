readIBD <- function(infile, nhaploids, lowerCutOff, upperCutOff)
{
  # example
  # infile = data/run_7/l_1.3e8_N_1000_mu_1.25e-8_r_1e-8_1D
  # nhaploids = 80
  # cutOff = 2e6
  
  
  file = paste(infile, ".out.ibd", sep = "")
  ibd = read.table(file, sep = "\t", header=FALSE)
  file = paste(infile, ".out.hbd", sep = "")
  hbd = read.table(file, sep = "\t", header=FALSE)
  ndiploids = nhaploids/2
  
  toIndex = matrix(nrow=ndiploids, ncol = 2, 0)
  cnt = 1
  for (i in 1:ndiploids){
    for (j in 1:2){
      toIndex[i,j] = cnt
      cnt = cnt + 1
    }
  }
  
  n = nrow(ibd)
  ibdM = matrix(nrow=nhaploids, ncol = nhaploids, 0)
  for (i in 1:n){
    ii = toIndex[ibd$V1[i], ibd$V2[i]]
    jj = toIndex[ibd$V3[i], ibd$V4[i]]
    diff = ibd$V7[i]-ibd$V6[i]
    if (diff > lowerCutOff & diff < upperCutOff){
      ibdM[ii,jj] = ibdM[ii,jj] + 1
      ibdM[jj,ii] = ibdM[ii,jj]
    }
  }
  
  n = nrow(hbd)
  for (i in 1:n){
    ii = toIndex[hbd$V1[i], hbd$V2[i]]
    jj = toIndex[hbd$V3[i], hbd$V4[i]]
    diff = hbd$V7[i]-hbd$V6[i]
    if (diff > lowerCutOff & diff < upperCutOff){
      # number of segments
      ibdM[ii,jj] = ibdM[ii,jj] + 1
      ibdM[jj,ii] = ibdM[ii,jj]
    }
  }
  
  return(ibdM)
}

workingDir <- "/Users/halasadi/eems2/test/data/4x5/recent_barrier/"
ibdM = readIBD(paste0(workingDir, "recent_barrier_mt_300_nsamp_per_deme_20.merged"), 400, 4e6, Inf)
write.table(ibdM, file = paste0(workingDir, "eems_4_Inf.sims"), quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)

#### computing variance of IBD empirically
#m = 4
#L = 264
#N = 1e4
#analytical_m <- 100*(25+m*N)/((50+m*N)^2)
#analytical_v <- (100/(N*L))*log(L/m)
#num = 0
#denom = 0
#for (i in 1:22){
#  num = num + 12*log(12/m)
#  denom = denom + 12
#}
#denom = denom^2
#analytical_v <- (100/N)*(num/denom)

#ibdM = readIBD("3x4_barrier5/3x4_l_1.3e8_N_1e4_nsamp_20_rs_10", 240, 4e6)
#write.table(ibdM, file = "3x4_l_1.3e8_N_1e4_nsamp_20_rs_10.sims", quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)