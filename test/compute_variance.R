library(expm)
library(Bolstad)
#source("myImagePlot.R")
rm(list = ls())

makeQSym <- function(M, q){
  Q = matrix(nrow=4, ncol=4, 0)
  # STATES:
  # (1,1), (1,2), (2,2), C
  Q[1,2] = 2*M[1,2]
  Q[1,4] = q[1]
  
  Q[2,1] = M[2,1]
  Q[2,3] = M[1,2]
  
  Q[3,2] = 2*M[2,1]
  Q[3,4] = q[2]
  
  diag(Q) = -rowSums(Q)
  return(Q)
}

makeQ <- function(M, q){
  Q = matrix(nrow=5, ncol=5, 0)
  # STATES:
  # (1,1), (1,2), (2,1), (2,2), C
  Q[1,2] = M[1,2]
  Q[1,3] = M[1,2]
  Q[1,5] = q[1]
  
  Q[2,1] = M[2,1]
  Q[2,4] = M[1,2]
  
  Q[3,1] = M[2,1]
  Q[3,4] = M[1,2]
  
  Q[4,2] = M[2,1]
  Q[4,3] = M[2,1]
  Q[4,5] = q[2]
  diag(Q) = -rowSums(Q)
  return(Q)
}

temp <- function(u,t){
  return(exp(-2*u*t))
}

findtEnd <- function(u){
  thr = 1e-14
  t = seq(0, 1000, by=1)
  p = unlist(lapply(X=t, FUN=temp, u=u))
  return(max(which(p > thr)))
}

computeWeights55 <- function(d, L){
  
  x <- unlist(as.vector(read.table("xvals.txt")))
  w <- unlist(as.vector(read.table("weights.txt")))
  
  w = w/((d*L)/50)
  x = x/((d*L)/50)
  return(list(w=w, x=x))
  
}

computeWeightsMean <- function(u){
  x = c(0.07053988969199015, 0.3721268180016133, 0.9165821024832778, 1.707306531028346,
        2.749199255309431, 4.048925313850894, 5.615174970861625, 7.459017453671069,
        9.594392869581098, 12.03880254696431, 14.81429344263073, 17.94889552051938,
        21.47878824028502, 25.45170279318691, 29.9325546317006, 35.01343424047902,
        40.83305705672853, 47.61999404734653, 55.8107957500639, 66.52441652561578)
  w = c(0.1687468018511152, 0.2912543620060685, 0.266686102867001, 0.1660024532695055,
        0.07482606466879245, 0.02496441730928314, 0.006202550844572207, 0.001144962386476896,
        0.0001557417730278113, 1.54014408652249e-05, 1.086486366517979e-06, 5.3301209095567e-08,
        1.757981179050555e-09, 3.725502402512292e-11, 4.767529251578155e-13, 3.372844243362382e-15,
        1.155014339500436e-17, 1.539522140582338e-20, 5.28644272556909e-24, 1.656456612498994e-28)
  
  
  w = w*(1+x)/(2*u)
  x = x/(2*u)
  return(list(w=w, x=x))
}
  
CalculateMean <- function(u, Q, q){
  quad <- computeWeightsMean(u)
  # in morgans
  t = quad$x
  w = quad$w

  Q = Q[-5,]
  Q = Q[,-5]
  
  p = rep(0, length(t))
  
  for (i in 1:length(t)){
    # exact
    P = expm(Q*t[i])
    #p[i] = q[1]*P[1,1] + q[2]*P[1,4]
    p[i] = q[1]*P[2,1] + q[2]*P[2,4]
    
  }
  integral <- w%*%p
  return(integral)
}

calculateInnerIntegral <- function(d, L, Q, q){
  quad <- computeWeights55(d, L)
  # in morgans
  t = quad$x
  w = quad$w
  
  Q = Q[-5,]
  Q = Q[,-5]
  
  p = rep(0, length(t))

    for (i in 1:length(t)){
    # exact
    P = expm(Q*t[i])
    #p[i] = q[1]*P[1,1] + q[2]*P[1,4]
    p[i] = q[1]*P[2,1] + q[2]*P[2,4]
  }
  integral <- w%*%p
  return(integral)
}

CalculateDoubleIntegral <- function(m, L, Q, q){
  # in cM
  #m = 4
  # in cM
  #L = 278
  
  n = 1000
  dgrid <- seq(from = m/L, to =1, length.out = n)
  f <- rep(0, n)
  for (i in 1:n){
    f[i] <- (1-dgrid[i])*calculateInnerIntegral(dgrid[i], L, Q, q)
  }
  
  return(2*sintegral(dgrid, f, n =n)$value)
}

# diploid pop size
N = c(1000, 1000)
mrates = c(0.01,0.01)
M = matrix(nrow = 2, ncol =2, 0)
M[1,2] = mrates[1]
M[2,1] = mrates[2]
diag(M) = -rowSums(M)
q = 1/(2*N)
Q = makeQ(M,q)
L = 286
m = 4
#approx <- (log(L/m)-0.5) * (m/L) * CalculateMean(u = (1e-8*m*1e6), Q, q)
#int <- CalculateDoubleIntegral(m, L, Q, q)

mrates <- seq(0.01, 0.5, length.out = 10)
approxs <- rep(0, 10)
ints <- rep(0, 10)
for (i in 1:length(mrates)){
  M = matrix(nrow = 2, ncol =2, 0)
  M[1,2] = mrates[i]
  M[2,1] = mrates[i]
  Q = makeQ(M,q)
  approxs[i] <- (log(L/m)-1+(m/L)) * (m/L) * CalculateMean(u = (1e-8*m*1e6), Q, q)
  #approxs[i] <- log(L/m) * (m/L) * CalculateMean(u = (1e-8*m*1e6), Q, q)
  ints[i] <- CalculateDoubleIntegral(m, L, Q, q)
}

plot(ints, approxs, log = "xy", main = "N=1000", xlab = "exact", ylab = "approx")
abline(a = 0, b = 1, col = "red", lwd = 2)
