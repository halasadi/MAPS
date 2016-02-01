expectedIBD <- as.matrix(read.table("../src/expectedIBD_m_0.05.txt", as.is=T, header=F))
#observedIBD <- as.matrix(read.table("data/4x5/uniform/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))
observedIBD <- as.matrix(read.table("data/4x5/old/uniform_nsamp_50/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))

m = 4
L = 260
totalInside = 0
total = 0
nsamp <- 50
for (i in 1:20){
  for (j in i:20){
    total = total + 1
    if (i == j){
      n = choose(nsamp,2)
      varterm <- expectedIBD[i,j]/n + (2*(nsamp-2)*(m/L)*expectedIBD[i,j]*expectedIBD[i,j])/n
      print((varterm-(expectedIBD[i,j]/n))/(expectedIBD[i,j]/n))
    }
    else{
      n = nsamp*nsamp
      varterm <- expectedIBD[i,j]/n + ((nsamp+nsamp-2)*(m/L)*expectedIBD[i,j]*expectedIBD[i,j])/(n)
    }
    

    lowerBnd = observedIBD[i,j] - 1.64*sqrt(varterm)
    upperBnd = observedIBD[i,j] + 1.64*sqrt(varterm)
    if (expectedIBD[i,j] > lowerBnd & expectedIBD[i,j] < upperBnd){
      totalInside = totalInside + 1;
    }
    
  }
}

# should be around 90% for 90% CI
frac = totalInside/total