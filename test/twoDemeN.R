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

computeWeights <- function(u){
  x = c(0.118440697736960550688, 0.3973475034735802657556, 0.8365549141880933313119, 1.437175158191620443607,
        2.200789508440616292336, 3.129448303166859096349, 4.225699164493802071261, 5.492626704368934083587,
        6.933903364122364597039, 8.553853192793023779194, 10.35753137020864105106, 12.35082332811269876439,
        14.54056869943518703492, 16.93471724415800802837, 19.54252664684054185266, 22.37481610233449499411,
        25.44429563058376261798, 28.76600031447167014762, 32.35787326932856805551, 36.24156497875364752439,
        40.44355691460364227197, 44.99678841355200250088, 49.94309754094208987181, 55.33704611950810443499,
        61.25224904369593075136, 67.79260716731075303985, 75.11420274687672563149, 83.47405073153149030595,
        93.36359463048878316735, 106.0462505962874034422)
  w = c(0.02093564741472521761, 0.09585049298017654367, 0.18833296435057945936, 0.23281944819987904471,
        0.2060782293528492151, 0.138528960450616358, 0.07293919110208096649, 0.030605607903988887905,
        0.010333948458420042431, 0.002821608083735993584, 6.2402663742264620427E-4, 1.1168849922460852198E-4,
        1.6129719270580565631E-5, 1.87044426274856472768E-6, 1.72995513372709914535E-7, 1.26506996496773906645E-8,
        7.2352574135703022224E-10, 3.19320138447436406004E-11, 1.069761647687436460972E-12, 2.66597906070505518515E-14,
        4.82019019925788439097E-16, 6.12740480626441608041E-18, 5.26125812567892365789E-20, 2.89562589607893296815E-22,
        9.51695437836864011982E-25, 1.69046847745875738033E-27, 1.39738002075239812243E-30, 4.20697826929603166432E-34,
        2.89826026866498969507E-38, 1.411587124593531584E-43)
  
  w = w*(1/(L*2*r*L))
  x = x/(2*r*L)
  return(list(w=w, x=x))
}

#calculateIntegral2 <- function(L,r){
#  ret = computeWeights(u)
#  t = ret$x
#  w = ret$w
#  n = length(t)
#  P = rep(0, n)
#  print(t)
#  for (i in 1:n){
#    P[i] = expm(Q*t[i])[1,4]
#  }
#  p = c(0, diff(P)/(t[2:n]-t[1:(n-1)]))
#  int = w%*%p
#  return(int)
#}

calculateIntegral <- function(L, r, Q, M, q){
  n = 10000
  tEnd = findtEnd(r*L)
  t = seq(0, tEnd, length.out=n)
  # WITH THE CDF
  P = rep(0, n)
  for (i in 1:n){
    # (1,1) -> C
    P[i] = expm(Q*t[i])[1,5]
  }
  p = c(0, diff(P)/(t[2:n]-t[1:(n-1)]))

  p = p*exp(-2*r*t*L)*2*r*t
  w = w*(1/(L*2*r*L))
  integral = trapezoidal.integration(t, p)
  print(integral)
  
  # WITH THE DENSITY
  Q = Q[-5,]
  Q = Q[,-5]
  p = rep(0, n)
  papprox = rep(0, n)
  for (i in 1:n){
    # exact
    P = expm(Q*t[i])
    p[i] = q[1]*P[1,1] + q[2]*P[1,4]
    # approximation
    Papprox= expm(M*t[i])
    papprox[i] = q[1]*Papprox[1,1]*Papprox[1,1] + q[2]*Papprox[1,2]*Papprox[1,2]
  }
  p = p*exp(-2*r*t*L)*2*r*t
  papprox = papprox*exp(-2*r*t*L)*2*r*t
  
  integral = trapezoidal.integration(t, p)
  plot(t, p)
  print(integral)
  integral = trapezoidal.integration(t, papprox)
  print(integral)
  
  ret = computeWeights(r, L)
  x = ret$x
  w = ret$w
  papprox = rep(0, length(x))
  for (i in 1:length(x)){
    Papprox = expm(M*x[i])
    papprox[i] = q[1]*Papprox[1,1]*Papprox[1,1] + q[2]*Papprox[1,2]*Papprox[1,2] 
  }
  print(w%*%papprox)
}

N = c(10000, 10000)
m = c(0.0,0.0)
M = matrix(nrow = 2, ncol =2, 0)
M[1,2] = m[1]
M[2,1] = m[2]
diag(M) = -rowSums(M)
q = 1/(2*N)
r = 1e-8
Q = makeQ(M,q)
L = 3e6
number_constantsize_m = (2*2*N*r)/(1+2*2*N*L*r)^2
calculateIntegral(L, r, Q, M, q)
#int2 = calculateIntegral2(L, r)
