library(expm)
#source("myImagePlot.R")

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


t = 1
m = c(0.00001, 0.005)
m = c(0,0)
N = c(10300, 20004)
q = 1/(2*N)
M = matrix(nrow = 2, ncol =2, 0)
M[1,2] = m[1]
M[2,1] = m[2]
M = M
Q = makeQ(M,q)
#Q = Q[-5,]
#Q = Q[,-5]
P = expm(Q*t)


diag(M) = -rowSums(M)
Papprox = expm(M*t)
print(Papprox[1,1]*Papprox[1,1])
print(P[1,1])
print("----------")

print(Papprox[1,1]*Papprox[2,2])
print(P[2,2])
print("----------")

print(Papprox[2,2]*Papprox[2,2])
print(P[4,4])
print("----------")

Qsym = makeQSym(M, q)
Psym = expm(Qsym*t)
