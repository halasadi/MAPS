readIBD <- function(infile, nhaploids, cutOff)
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
    if (diff > cutOff){
      ibdM[ii,jj] = ibdM[ii,jj] + 1
      ibdM[jj,ii] = ibdM[ii,jj]
    }
  }
  
  n = nrow(hbd)
  for (i in 1:n){
    ii = toIndex[hbd$V1[i], hbd$V2[i]]
    jj = toIndex[hbd$V3[i], hbd$V4[i]]
    diff = hbd$V7[i]-hbd$V6[i]
    if (diff > cutOff){
      ibdM[ii,jj] = ibdM[ii,jj] + 1
      ibdM[jj,ii] = ibdM[ii,jj]
    }
  }
  
  return(ibdM)
}

ibdM = readIBD("3x4_uniform/3x4_l_1.3e8_N_1e4_nsamp_50_rs_10", 600, 4e6)
write.table(ibdM, file = "3x4_uniform.sims", quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)

#ibdM = readIBD("3x4_barrier5/3x4_l_1.3e8_N_1e4_nsamp_20_rs_10", 240, 4e6)
#write.table(ibdM, file = "3x4_l_1.3e8_N_1e4_nsamp_20_rs_10.sims", quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)