rm(list = ls())

library(data.table)
ibd_blocklens = read.table("ibd-blocklens.csv", sep = ",", header=TRUE)
ibd_labels = read.csv("ibd-pop-info.csv", stringsAsFactors = FALSE)


loc_info <- read.csv("countries.txt", header=TRUE, sep = ",", stringsAsFactors = FALSE, check.names=FALSE)
row.names(loc_info) <- loc_info$country

# in centimorgan
lowerBound = 8
upperBound = Inf

ibd_blocklens = subset(ibd_blocklens, ibd_blocklens$maplen > lowerBound & ibd_blocklens$maplen < upperBound)

#countries_to_keep = c("Italy", "United Kingdom", "France", "Ireland", "Portugal", "Spain", "Germany", 
#                      "Switzerland", "Austria", "Belgium", "Netherlands", "Czech Republic", "Slovenia")

countries_to_keep = c("France", "Ireland", "Italy", "Portugal", "Spain", "United Kingdom", "Germany", "Austria", "Belgium")

ids_to_keep <- ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF %in% countries_to_keep)]

ibd_blocklens = subset(ibd_blocklens, ibd_blocklens$id1 %in% ids_to_keep & ibd_blocklens$id2 %in% ids_to_keep)

unqids = unique(c(ibd_blocklens[,1], ibd_blocklens[,2]))
nsamp = length(unqids)
ibd_df <- data.frame(matrix(0, nrow = nsamp, ncol = nsamp))
rownames(ibd_df) = unqids
colnames(ibd_df) = unqids

for (i in 1:nrow(ibd_blocklens)){
  id1 = toString(ibd_blocklens[i,]$id1)
  id2 = toString(ibd_blocklens[i,]$id2)
  ibd_df[id1, id2] = ibd_df[id1, id2] + 1
  ibd_df[id2,id1] = ibd_df[id1, id2]
}

write.table(ibd_df, file = paste0("popres_lowerBnd_", lowerBound, "_upperBnd_", upperBound, ".sims"), row.names=FALSE, col.names=FALSE)


coords = matrix(0, nrow = nsamp, ncol = 2)
for (i in 1:nsamp){
  ind = which(ibd_labels$SUBJID == unqids[i])
  lat = loc_info[ibd_labels$COUNTRY_SELF[ind],]$lat
  long = loc_info[ibd_labels$COUNTRY_SELF[ind],]$long
  coords[i,1] = lat
  coords[i,2] = long
}

write.table(coords, file = paste0("popres_lowerBnd_", lowerBound, "_upperBnd_", upperBound, ".coord"), row.names=FALSE, col.names=FALSE)