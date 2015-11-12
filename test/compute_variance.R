library(expm)
#source("myImagePlot.R")
rm(list = ls())

trapezoidal.integration = function(x, f)
{
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }
  
  # check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand "f" is not numeric.')
  }
  
  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(f))
  {
    stop('The lengths of the variable of integration and the integrand do not match.')
  }
  
  ### finish checks
  
  # obtain length of variable of integration and integrand
  n = length(x)
  
  # integrate using the trapezoidal rule
  integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))
  
  # print the definite integral
  return(integral)
}

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

computeWeights <- function(d, L){
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
  
  
  w = w*(d/x)
  x = x/((d*L)/50)
  return(list(w=w, x=x))
}



calculateInnerIntegral <- function(d, L, Q, M, q){
  quad <- computeWeights(d, L)
  # in morgans
  t = quad$x
  w = quad$w
  
  Q = Q[-5,]
  Q = Q[,-5]
  
  p = rep(0, length(t))

    for (i in 1:length(t)){
    # exact
    P = expm(Q*t[i])
    p[i] = q[1]*P[1,1] + q[2]*P[1,4]
 
  }
  integral <- w%*%p
  return(integral)
}

CalculateDoubleIntegral <- function(Q, M, q){
  # in cM
  m = 4
  # in cM
  L = 278
  
  n = 50000
  dgrid <- seq(from = m/L, to =1, length.out = n)
  f <- rep(0, n)
  for (i in 1:n){
    f[i] <- (1-dgrid[i])*calculateInnerIntegral(dgrid[i], L, Q, M, q)
  }
  
  integral <- trapezoidal.integration(dgrid, f)
  return(2*integral)
}

# diploid pop size
N = c(10000, 10000)
m = c(0.0,0.0)
M = matrix(nrow = 2, ncol =2, 0)
M[1,2] = m[1]
M[2,1] = m[2]
diag(M) = -rowSums(M)
q = 1/(2*N)
Q = makeQ(M,q)
#int <- calculateIntegral(r, L, Q, M, q)
int <- CalculateDoubleIntegral(Q, M, q)
