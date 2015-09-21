rm(list = ls())

library(data.table)
ibd_blocklens = read.table("ibd-blocklens.csv", sep = ",", header=TRUE)
ibd_labels = read.csv("ibd-pop-info.csv", stringsAsFactors = FALSE)


loc_info <- read.csv("countries.txt", header=TRUE, sep = ",", stringsAsFactors = FALSE, check.names=FALSE)
row.names(loc_info) <- loc_info$country

# in centimorgan
lowerBound = 6
upperBound = Inf

ibd_blocklens = subset(ibd_blocklens, ibd_blocklens$maplen > lowerBound & ibd_blocklens$maplen < upperBound)

countries_to_keep = c("Italy", "United Kingdom", "France", "Ireland", "Portugal", "Spain", "Germany", 
                      "Switzerland", "Austria", "Belgium", "Netherlands", "Czech Republic", "Slovenia")
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

# italy_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Italy")]
# uk_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "United Kingdom")]
# france_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "France")]
# ireland_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Ireland")]
# portugal_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Portugal")]
# spain_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Spain")]
# 
# # this is for assigning geographic locations
# total_ids = c(italy_ids, uk_ids, france_ids, ireland_ids, portugal_ids, spain_ids)
# n = length(total_ids)
# 
# sim = matrix(0, nrow = n, ncol = n)
# for (i in 1:n){
#   print(i)
#   id1 = total_ids[i]
#   for (j in (i+1):n){
#     id2 = total_ids[j]
#     inds = which(ibd_blocklens$id1 == id1 & ibd_blocklens$id2 == id2)
#     if (sum(inds) > 0){
#       blocks = ibd_blocklens$maplen[inds]
#       sim[i,j] = sum(blocks > cutOff)
#       sim[j,i] = sim[i,j]
#     }
#   }
# }
# write.table(sim, file = "popres.sim", row.names=FALSE, col.names=FALSE)
# 
# # in same order as total_ids
# lat = c(10, -1.22, 2, -7.6, -8.4, -3.8)
# long = c(46, 52.2, 46, 53, 39.4, 40.7)
# 
# l_vector = c(length(italy_ids), length(uk_ids), length(france_ids), length(ireland_ids),
#              length(portugal_ids), length(spain_ids))
# lat_w = rep(lat, l_vector)
# long_w = rep(long, l_vector)
# 
# df = data.frame(long_w, lat_w)
# write.table(df, file = "popres.coord", row.names=FALSE, col.names=FALSE)
# 
