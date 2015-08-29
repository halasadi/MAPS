basePath = "popres-EEMS2-test-sim"

mcmcmtiles = read.table(paste(basePath, "/mcmcmtiles.txt", sep =""))
hist(mcmcmtiles$V1)

mcmcqtiles = read.table(paste(basePath, "/mcmcqtiles.txt", sep =""))
hist(mcmcqtiles$V1)

mcmcmhyper = read.table(paste(basePath, "/mcmcmhyper.txt", sep =""))
mmeans = mcmcmhyper$V1
mvars = mcmcmhyper$V2
# should be uniform
# these are the mean rates
hist(mmeans)
# should be distributed inv-gamma
# the variances should be uniform on the log scale 
hist(log(mvars))

library(LearnBayes)
cq = 0.001/2
cw = 0.001/2
dq = 1/2
dm = 1/2
D = rigamma(n = 10000, a = cq, b = dq)
hist(log(D))

mcmcqhyper = read.table(paste(basePath, "/mcmcqhyper.txt", sep =""))
qmeans = mcmcqhyper$V1
qvars = mcmcqhyper$V2
# should be uniform on the log scale
# these are the mean rates
hist(qmeans)
# should be distributed inv-gamma
# the variances should be uniform on the log scale 
hist(log(qvars))

# doesn't seem like prior is working

f = file(paste(basePath, "/mcmcmrates.txt", sep =""))
mcmcmrates <- strsplit(readLines(f), " ")
mcmcmrates <- lapply(mcmcmrates, as.numeric)
close(f)

f = file(paste(basePath, "/mcmcqrates.txt", sep =""))
mcmcqrates <- strsplit(readLines(f), " ")
mcmcqrates <- lapply(mcmcqrates, as.numeric)
close(f)

pqvals = vector()
for (i in 1:length(mcmcqrates)){
  for (j in 1:length(mcmcqrates[[i]])){
    #pqvals = c(pqvals, dnorm(log10(mcmcqrates[[i]][j]), mean = qmeans[i], sd = sqrt(qvars[i])))
    pqvals = c(pqvals, dnorm(mcmcqrates[[i]][j], mean = qmeans[i], sd = sqrt(qvars[i])))
    
  }
}


pmvals = vector()
for (i in 1:length(mcmcmrates)){
  for (j in 1:length(mcmcmrates[[i]])){
    print(log10(mcmcmrates[[i]][j]))
    #pmvals = c(pmvals, dnorm(log10(mcmcmrates[[i]][j]), mean = mmeans[i], sd = sqrt(mvars[i])))
    pmvals = c(pmvals, dnorm(mcmcmrates[[i]][j], mean = mmeans[i], sd = sqrt(mvars[i])))
    
  }
}

f = file(paste(basePath, "/mcmcwcoord.txt", sep =""))
qx <- strsplit(readLines(f), " ")
qx <- unlist(lapply(qx, as.numeric))
close(f)

f = file(paste(basePath, "/mcmczcoord.txt", sep =""))
qy <- strsplit(readLines(f), " ")
qy <- unlist(lapply(qy, as.numeric))
close(f)

f = file(paste(basePath, "/mcmcxcoord.txt", sep =""))
mx <- strsplit(readLines(f), " ")
mx <- unlist(lapply(mx, as.numeric))
close(f)

f = file(paste(basePath, "/mcmcycoord.txt", sep =""))
my <- strsplit(readLines(f), " ")
my <- unlist(lapply(my, as.numeric))
close(f)

plot(mx, my)
plot(qx, qy)
