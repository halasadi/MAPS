library(expm)
source("myImagePlot.R")

makeQ <- function(m, q){
  Q = matrix(nrow=4, ncol=4, 0)
  Q[1,2] = 2*m
  Q[1,4] = q[1]
  Q[2,1] = m
  Q[2,3] = m
  Q[3,2] = 2*m
  Q[3,4] = q[2]
  diag(Q) = -rowSums(Q)
  return(Q)
}

computeWeights <- function(u){
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

# the probability (density) that an IBD segment greater than u coalesces at time t
erlang <- function(u,t){
  return(exp(-2*u*t)*(1+2*u))
}


findtEnd <- function(u){
  thr = 1e-8
  t = seq(0, 1000, by=1)
  p = unlist(lapply(X=t, FUN=erlang, u=u))
  return(max(which(p > thr)))
}


calculateIntegral <- function(u){
  #ret = computeWeights(u)
  #x = ret$x
  #w = ret$w
  #P = rep(0, length(x))
  #for (i in 1:length(x)){
  #  P[i] = expm(Q*x[i])[2,4]
  #}
  #p = c(0, diff(P)/(x[2:length(x)]-x[1:(length(x)-1)]))
  #int_g = w%*%p
  tEnd = findtEnd(u)
  
  n = 10000
  x = seq(0, tEnd, length.out=n)
  P = rep(0, n)
  for (i in 1:n){
    P[i] = expm(Q*x[i])[2,4]
  }
  p = c(0, diff(P)/(x[2:n]-x[1:(n-1)]))
  p = p*exp(-2*u*x)*(1+2*u*x) 
  int_t = 0.5*sum((x[2:n] - x[1:(n-1)]) * (p[2:n] + p[1:(n-1)]))
  
  return(int_t)
  
}

m = 0.01
N = c(1000, 1000)
q = 1/(2*N)
Q = makeQ(m,q)
u = 0.02
genomeSize = 2.99e9
int = calculateIntegral(u)


# between deme
b_F = ((m*(4*m+3*u))/(2*N[1]*u*(2*m+u)*(2*m+u)))/2
# within deme
w_F = (1/(2*N[1]*u) + (m+u)/(2*N[1]*((2*m+u)^2)))/2

file = "data/run_1/l_1.3e8_N_1000_m_0.01_mu_1.25e-8_r_1e-8_2D.ibd.2000000"
sim = read.csv(file, sep = " ", header=FALSE)
sim[,81] = NULL
sim = as.matrix(sim)
diag(sim) = NA
sim = sim/genomeSize
pop1 = sim[1:40, 1:40]
pop2 = sim[41:80, 41:80]
betw_pops = sim[1:40, 41:80]

m1 = mean(pop1[upper.tri(pop1)])
m2 = mean(pop2[upper.tri(pop2)])
m3 = mean(betw_pops)

myImagePlot(log(sim+1))



