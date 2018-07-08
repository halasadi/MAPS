## read in the meta data
meta_info_popres <- read.table("novembre2008.txt", sep = "\t", as.is = TRUE, header = TRUE)
meta_info <- data.frame(fid = meta_info_popres$ID, iid = meta_info_popres$ID, long = meta_info_popres$longitude, 
                        lat = meta_info_popres$latitude, stringsAsFactors = FALSE)

## read in the PSC calls. Here, we read in only chromosome 1
## Here, the file was generated using the Snakemake pipeline found here: https://github.com/halasadi/ibd_data_pipeline
## In this pipeline, both the .ibd and .hbd files are utilized (see the relabelIBD rule in the Snakemake file)
ibd_data <- read.table("POPRES_CHR1.merged.nosparse_cm.finalqc.ibd", header = TRUE, stringsAsFactors = FALSE)


ids <- paste0(meta_info$iid, "_", meta_info$fid)

locs <- paste0(meta_info$lat, " ", meta_info$long)

# BEAGLE outputs haplotype data, we treat
# the two different chromosomes of each individual as independent
ids <- c(paste0(ids, "_1"), paste0(ids, "_2"))
locs <- c(locs, locs)


n <- length(ids) 

# set up the matrix
ibd_summary <- matrix(nrow = n, ncol = n, 0)
rownames(ibd_summary) <- ids
colnames(ibd_summary) <- ids

# compute the lengths of the PSC segments
lengths <- ibd_data$end-ibd_data$start

# highlight PSC segments greater than 2cM
lowerBnd <- 2
upperBnd <- Inf
selected_inds <- which(lengths > lowerBnd & lengths < upperBnd)

for (i in 1:length(selected_inds)){
  id1 <- ibd_data$id1[selected_inds[i]]
  id2 <- ibd_data$id2[selected_inds[i]]
  if (id1 %in% ids & id2 %in% ids){
    ibd_summary[id1, id2] = ibd_summary[id1, id2] + 1
    ibd_summary[id2, id1] = ibd_summary[id1, id2]
  }
}


## write similarity matrix and coordinates to a .sims and .coord file respectively
write.table(ibd_summary, file = paste0("maps_", lowerBnd, "_", upperBnd, ".sims"), quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
write.table(data.frame(locs), file = paste0("maps_", lowerBnd, "_", upperBnd, ".coord"), quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
