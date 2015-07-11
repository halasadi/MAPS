basePath = "/Users/halasadi/Documents/eems2/data/test_prior"
df = read.csv(paste(basePath, "/mcmcthetas.txt", sep =""), header=FALSE)

# should be uniform
hist(log(df$V1))

mcmcmtiles = read.csv(paste(basePath, "/mcmcmtiles.txt", sep =""), header=FALSE, sep = " ")
hist(mcmcmtiles$V1)

mcmcmhyper = read.csv(paste(basePath, "/mcmcmhyper.txt", sep =""), header=FALSE, sep = " ")
mmeans = mcmcmhyper$V2
mvars = mcmcmhyper$V4
# should be uniform
# these are the mean rates
hist(mmeans)
# should be distributed inv-gamma
# the variances
hist(mvars,100)

# each row i is drawn from N(mcmchyper$V1[i], mcmchyper$V2[i])
# but that doesn't seem to be true! 
mcmcmrates = read.csv(paste(basePath, "/mcmcmrates.txt", sep =""), header=FALSE, sep = " ")
pvals = vector()
cnt = 1
for (i in 1:nrow(mcmcmrates)){
  x = as.numeric(mcmcmrates[i,])
  x <- x[!is.na(x)]
  for (j in 1:length(x)){
    pvals[cnt] = dnorm(x[j], mean = mmeans[i], sd = sqrt(mvars[i]))
    cnt = cnt + 1
  }
}
