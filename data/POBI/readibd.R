sampleIDs <- read.table("POBI_SampleLocations/POBI_SampleLocations.WTCCC2_ID.txt", header = F, stringsAsFactor=F)$V3
sample_info <- read.csv(file = "POBI_SampleLocations/Geocoding for Fabian 151113.csv", as.is=T)
sampleIDs <- na.omit(sample_info$id2[match(sampleIDs, sample_info$id2)])

match(sampleIDs, sample_info$id2)

ndip <- length(sampleIDs)
ibdM <- matrix(nrow = 2*ndip, ncol = 2*ndip, 0)

ids <- c(paste0(sampleIDs, "_1"), paste0(sampleIDs, "_2"))
rownames(ibdM) = ids
colnames(ibdM) = ids

MIN_LEN <- 10
MAX_LEN <- Inf

#latlong <- read.table("latlong.txt", header = T, stringsAsFactors = F)
#region_definitions <- read.table("POBI_SampleLocations/RegionDefinitions.csv", header = T, sep = ",", stringsAsFactors = F)


for (chr in 1:22){
  ibd <- read.table(paste0("WTCC2_POBI_illumina.finalqc.parsed.", chr, ".ibd"), header = T, stringsAsFactors = F)
  ibd <- ibd[which(ibd$id1 %in% ids & ibd$id2 %in% ids),]
  lens <- ibd$end - ibd$start
  ibd_keep <- ibd[which(lens > MIN_LEN & lens < MAX_LEN),]
  for (i in 1:nrow(ibd_keep)){
    ibdM[ibd_keep$id1[i], ibd_keep$id2[i]] = ibdM[ibd_keep$id1[i], ibd_keep$id2[i]] + 1
    ibdM[ibd_keep$id2[i], ibd_keep$id1[i]]  = ibdM[ibd_keep$id1[i], ibd_keep$id2[i]] 
  }
}

write.table(ibdM, file = "POBI_10_INF.sims", quote = F, row.names = F, col.names = F)

lat = sample_info$Lat_MGM_U[match(sampleIDs, sample_info$id2)]
long = sample_info$Lon_MGM_U[match(sampleIDs, sample_info$id2)]
lat = c(lat, lat)
long = c(long, long)

write.table(cbind(long, lat), file = "POBI_10_INF.coord", quote = F, row.names = F, col.names = F)
