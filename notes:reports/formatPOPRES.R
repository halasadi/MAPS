library(data.table)
ibd_blocklens = read.csv("data/coop_popres/ibd-blocklens.csv")
ibd_labels = read.csv("data/coop_popres/ibd-pop-info.csv")

# in centimorgan
cutOff = 4

ibd_blocklens = subset(ibd_blocklens, ibd_blocklens$maplen > cutOff)

italy_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Italy")]
uk_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "United Kingdom")]
france_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "France")]
ireland_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Ireland")]
portugal_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Portugal")]
spain_ids = ibd_labels$SUBJID[which(ibd_labels$COUNTRY_SELF == "Spain")]

# this is for assigning geographic locations
total_ids = c(italy_ids, uk_ids, france_ids, ireland_ids, portugal_ids, spain_ids)
n = length(total_ids)

sim = matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  print(i)
  id1 = total_ids[i]
  for (j in (i+1):n){
    id2 = total_ids[j]
    inds = which(ibd_blocklens$id1 == id1 & ibd_blocklens$id2 == id2)
    if (sum(inds) > 0){
      blocks = ibd_blocklens$maplen[inds]
      sim[i,j] = sum(blocks > cutOff)
      sim[j,i] = sim[i,j]
    }
  }
}
write.table(sim, file = "popres.sim", row.names=FALSE, col.names=FALSE)

# in same order as total_ids
lat = c(10, -4.9, 2, -7.6, 8.4, -3.8)
long = c(46, 50.3, 46, 53, -39.4, 40.7)

l_vector = c(length(italy_ids), length(uk_ids), length(france_ids), length(ireland_ids),
             length(portugal_ids), length(spain_ids))
lat_w = rep(lat, l_vector)
long_w = rep(long, l_vector)

df = data.frame(long_w, lat_w)
write.table(df, file = "popres.coord", row.names=FALSE, col.names=FALSE)

