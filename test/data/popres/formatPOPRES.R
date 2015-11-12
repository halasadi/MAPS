rm(list = ls())

library(data.table)
ibd_blocklens = read.table("ibd-blocklens.csv", sep = ",", header=TRUE)
ibd_labels = read.csv("ibd-pop-info.csv", stringsAsFactors = FALSE)


loc_info <- read.csv("countries.txt", header=TRUE, sep = ",", stringsAsFactors = FALSE, check.names=FALSE)
row.names(loc_info) <- loc_info$country

# in centimorgan
lowerBound = 4
upperBound = Inf

ibd_blocklens = subset(ibd_blocklens, ibd_blocklens$maplen > lowerBound & ibd_blocklens$maplen < upperBound)

#countries_to_keep = c("Italy", "United Kingdom", "France", "Ireland", "Portugal", "Spain", "Germany",
#                      "Switzerland", "Austria", "Belgium", "Netherlands", "Czech Republic", "Hungary", "Bosnia", "Slovenia", "Croatia",
#                      "Poland", "Kosovo", "Albania", "Serbia", "Macedonia")

countries_to_keep = c("Spain")

#countries_to_keep = c("France", "Ireland", "Italy", "Portugal", "Spain", "United Kingdom", "Germany", "Austria", "Belgium")

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
  
  # number of segments
  #ibd_df[id1, id2] = ibd_df[id1, id2] + 1
  # total
  ibd_df[id1, id2] = ibd_df[id1, id2] + ibd_blocklens[i,]$maplen
  
  ibd_df[id2,id1] = ibd_df[id1, id2]
}

# divide by 278 cM (total genome size)
ibd_df = ibd_df/278

#write.table(ibd_df, file = paste0("popres_fraction_lowerBnd_", lowerBound, "_upperBnd_", upperBound, ".sims"), row.names=FALSE, col.names=FALSE)


coords = matrix(0, nrow = nsamp, ncol = 2)
for (i in 1:nsamp){
  ind = which(ibd_labels$SUBJID == unqids[i])
  lat = loc_info[ibd_labels$COUNTRY_SELF[ind],]$lat
  long = loc_info[ibd_labels$COUNTRY_SELF[ind],]$long
  coords[i,1] = lat
  coords[i,2] = long
}

#n = nrow(ibd_df)
#range <- 1:n
#covs <- c()
#for (i in 2:n){
#  # fix i
#  print(i)
#  rem_indvs <- range[-i]
#  combs <- combn(rem_indvs, m = 2)
#  X = c()
#  Y = c()
#  for (j in 1:ncol(combs)){
#    X = c(X,ibd_df[i,combs[,j][1]])
#    Y = c(Y, ibd_df[i,combs[,j][2]])
#  }
#  covs <- c(covs, cov(X,Y))
#}

#write.table(coords, file = paste0("popres_fracion_lowerBnd_", lowerBound, "_upperBnd_", upperBound, ".coord"), row.names=FALSE, col.names=FALSE)

L = 278
m = 4
ibd_m = as.matrix(ibd_df)

indep_ibd <- ibd_m[row(ibd_m) == (col(ibd_m)-1)]
hist(indep_ibd, 100)

empirical_var <- var(indep_ibd)
analytical_var <- log(L/m)*(m/L)*mean(indep_ibd)

print(empirical_var)
print(analytical_var)
