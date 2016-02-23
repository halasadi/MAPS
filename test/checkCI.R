expectedIBD <- as.matrix(read.table("../src/expectedIBD_m_0.01.txt", as.is=T, header=F))
#expectedIBD <- as.matrix(read.table("../src/expectedIBD_m_0.05.txt", as.is=T, header=F))
observedIBD <- as.matrix(read.table("data/4x5/uniform/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))
#observedIBD <- as.matrix(read.table("data/4x5/old/uniform_nsamp_50/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))
#observedIBD <- as.matrix(read.table("data/4x5/uniform_nsamp_100/eems_4_Inf-EEMS2-test-sim/rdistJtDobsJ.txt", as.is=T, header=F))

m = 4;
L = 260;
total = 0;
nsamp <- 20;

computeCI <- function(nsamp, m, L, expectedIBD, observedIBD, alpha){
  totalInside = 0;
  for (i in 1:20){
    for (j in i:20){
      total = total + 1;
      if (i == j){
        varterm <- expectedIBD[i,j]/choose(nsamp,2) #+ (2*(nsamp-2)*(m/L)*expectedIBD[i,j]*expectedIBD[i,j])/choose(nsamp,2);
        #print((varterm-(expectedIBD[i,j]/n))/(expectedIBD[i,j]/n))
      }
      else{
        varterm <- expectedIBD[i,j]/(nsamp*nsamp) #+ ((nsamp+nsamp-2)*(m/L)*expectedIBD[i,j]*expectedIBD[i,j])/(nsamp*nsamp);
      }
      
      
      z <- qnorm(1 - (alpha/2));
      lowerBnd = observedIBD[i,j] - z*sqrt(varterm);
      upperBnd = observedIBD[i,j] + z*sqrt(varterm);
      if (expectedIBD[i,j] > lowerBnd & expectedIBD[i,j] < upperBnd){
        totalInside = totalInside + 1;
      }
      
    }
  }
  
  # should be around 90% for 90% CI
  frac = totalInside/total;
  return(frac);
}

alphas = seq(0.01, 0.99, by = 0.01);
fracs <- vector();
for (ii in 1:length(alphas)){
  fracs[ii] = computeCI(nsamp, m, L, expectedIBD, observedIBD, alphas[ii]);
}
plot(1-alphas, fracs, ylab = "estimated coverage", xlab = "expected coverage", xlim = c(0,1), ylim = c(0,1));
abline(a=0, b=1, col = "red", lwd=3);
mse <- mean(sum((fracs-(1-alphas))^2));
print(mse)
