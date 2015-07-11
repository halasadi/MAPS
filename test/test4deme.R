library(expm)

makeQ <- function(M, q){
  Q = matrix(nrow=11, ncol=11, 0)
  Q[1,2] = 2*M[1,2]
  Q[1,3] = 2*M[1,3]
  Q[1,11] = q[1]
  Q[2,1] = M[2,1]
  Q[2,4] = M[2,4]
  Q[2,5] = M[1,2]
  Q[2,6] = M[1,3]
  Q[3,1] = M[3,1]
  Q[3,4] = M[3,4]
  Q[3,6] = M[1,2]
  Q[3,8] = M[1,3]
  Q[4,2] = M[4,2]
  Q[4,3] = M[4,3]
  Q[4,7] = M[1,2]
  Q[4,9] = M[1,3]
  Q[5,2] = 2*M[2,1]
  Q[5,7] = 2*M[2,4]
  Q[5,11] = q[2]
  Q[6,2] = M[3,1]
  Q[6,3] = M[2,1]
  Q[6,7] = M[3,4]
  Q[6,9] = M[2,4]
  Q[7,2] = 0
  Q[7,4] = M[2,1]
  Q[7,5] = M[4,2]
  Q[7,6] = M[4,3]
  Q[7,10] = M[2,4]
  Q[8,3] = 2*M[3,1]
  Q[8,9] = 2*M[3,4]
  Q[8,11] = q[3]
  Q[9,4] = M[3,1]
  Q[9,6] = M[4,2]
  Q[9,8] = M[4,3]
  Q[9,10] = M[3,4]
  Q[10,7] = 2*M[4,2]
  Q[10,9] = 2*M[4,3]
  Q[10,11] = q[4]

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
norm_vec <- function(x) sqrt(sum(x^2))

krylov <- function(Qrate, n, mkrylov, t){
  Q = matrix(nrow=n, ncol = mkrylov, 0)
  Q[n,1] = 1
  H = matrix(nrow=mkrylov, ncol=mkrylov, 0)
  for (k in 2:mkrylov){
    q = Q[,k-1]
    z = Qrate%*%q
    for (i in 1:k){
      H[i, k-1] = Q[,i]%*%z
      z = z - H[i,k-1]*Q[,i]
    }
    H[k, k-1] = norm_vec(z)
    if (H[k,k-1] == 0){
      break
    }
    Q[,k] = z/H[k, k-1]
  }
  print(H)
  E = expm(H*t)
  return(Q%*%E%*%t(Q)%*%Q[,1])
}

calculateIntegral <- function(u){
  ret = computeWeights(u)
  x = ret$x
  w = ret$w
  P = matrix(0, nrow=length(x), ncol = 11)
  Pkrylov = P 
  v = rep(0, 11)
  v[11] = 1
  for (i in 1:length(x)){
    P[i,] = expm(Q*x[i])%*%v
    Pkrylov[i,] = krylov(Q, 11, 28, x[i])
  }
  
  int = rep(0, 11)
  for (i in 1:11){
    p = c(0, diff(P[,i])/(x[2:length(x)]-x[1:(length(x)-1)]))
    int[i] = w%*%p
  }
  
  intkrylov = rep(0, 11)
  for (i in 1:11){
    pkrylov = c(0, diff(Pkrylov[,i])/(x[2:length(x)]-x[1:(length(x)-1)]))
    intkrylov[i] = w%*%pkrylov
  }
  print(intkrylov-int)
  return(int)
}

#m = c(0.0001, 0.0001, 0.0001, 0.0001)
M = matrix(nrow=4, ncol=4, 0)
#for (i in 1:4){
#  for (j in 1:4){
#    M[i,j] = (m[1] + m[2])/2
#    M[j,i] = M[i,j]
#  }
#}
M[1,2] = M[2,1] = 0.143995
M[1,3] = M[3,1] = 0.213838
M[2,4] = M[4,2] = 0.213838
M[3,4] = M[4,3] = 0.28368
#q = c(0.0001, 0.0001, 0.0001, 0.0001)
q = c(0.00548318, 0.00945489, 0.00344319, 0.00353327)
Q = makeQ(M,q)
u = 0.02
int = calculateIntegral(u)
