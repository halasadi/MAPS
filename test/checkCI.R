expectedIBD <- as.matrix(read.table("../src/expectedIBD.txt", as.is=T, header=F))
observedIBD <- as.matrix(read.table("data/4x5/uniform/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))
#observedIBD <- as.matrix(read.table("data/4x5/old/uniform_nsamp_50/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))

totalInside = 0
total = 0
for (i in 1:20){
  for (j in i:20){
    total = total + 1
    if (i == j){
      n = choose(20,2)
    }
    else{
      n = 20*20
    }
    
    lowerBnd = observedIBD[i,j] - 1.64*sqrt(expectedIBD[i,j]/n)
    upperBnd = observedIBD[i,j] + 1.64*sqrt(expectedIBD[i,j]/n)
    if (expectedIBD[i,j] > lowerBnd & expectedIBD[i,j] < upperBnd){
      totalInside = totalInside + 1;
    }
    
  }
}

# should be around 90% for 90% CI
frac = totalInside/total