ndiploids = 120

toIndex = matrix(nrow=ndiploids, ncol = 2, 0)
cnt = 1
for (i in 1:ndiploids){
  for (j in 1:2){
    toIndex[i,j] = cnt
    cnt = cnt + 1
  }
}

toID <- function(id){
  id1 <- strsplit(id, split = "_")[[1]][1]
  hap <- strsplit(id, split = "_")[[1]][2]
  return(toIndex[as.integer(id1), as.integer(hap)])
}

readIBD <- function(infile, nhaploids, lowerCutOff, upperCutOff)
{
  ibd = read.table(infile, sep = " ", header=T, stringsAsFactors = F)
  
  n = nrow(ibd)
  ibdM = matrix(nrow=nhaploids, ncol = nhaploids, 0)
  for (i in 1:n){
    id1 = toID(ibd$id1[i])
    id2 = toID(ibd$id2[i])
    diff = ibd$end[i]-ibd$start[i]
    if (diff > lowerCutOff & diff < upperCutOff){
      ibdM[id1,id2] = ibdM[id1,id2] + 1
      ibdM[id2,id1] = ibdM[id1,id2]
    }
  }
  return(ibdM)
}

workingDir <- "/Users/halasadi/eems2/test/data/4x5/heter_popsizes_unevensampling2/"
nchr = 20
lowerCutOff = 4e6
upperCutoff = Inf
nhaploids = ndiploids*2
ibdM = readIBD(paste0(workingDir, "heter_popsizes_uneven_sampling_mt_300gen_1.out.qc.ibd"), nhaploids, lowerCutOff, upperCutoff)

for (i in 2:nchr){
  ibdM = ibdM + readIBD(paste0(workingDir, "heter_popsizes_uneven_sampling_mt_300gen_", i, ".out.qc.ibd"), nhaploids, lowerCutOff, upperCutoff)
}
write.table(ibdM, file = paste0(workingDir, "eems_4_Inf.sims"), quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)