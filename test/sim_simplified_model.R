library(expm)

computeWeights <- function(r, L){
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


calculateIntegral <- function(M, q, r, L){
  ret = computeWeights(r, L)
  x = ret$x
  w = ret$w
  d = nrow(M)
  lambda = matrix(0, nrow = d, ncol = d, 0)
  for (t in 1:length(x)){
    Pm = expm(M*x[t])
    for (i in 1:d){
      for (j in i:d){
        lambda[i,j] = lambda[i,j] + w[t]*sum(Pm[i,]*Pm[j,]*q)
        lambda[j,i] = lambda[i,j]
      }
    }
  }

  return(3e9*lambda)
}

N = rep(10000, 12)
q = 1/(2*N)

M = matrix(nrow = 12, ncol = 12, 0)
const = 0.05;
M[1,2] = M[2,1] = const
M[2,3] = M[3,2] = const
M[1,5] = M[5,1] = const
M[2,5] = M[5,2] = const
M[2,6] = M[6,2] = const
M[3,6] = M[6,3] = const
M[3,4] = M[4,3] = const
M[3,7] = M[7,3] = const
M[7,4] = M[4,7] = const
M[7,8] = M[8,7] = const
M[4,8] = M[8,4] = const
M[5,6] = M[6,5] = const
M[6,7] = M[7,6] = const
M[5,9] = M[9,5] = const
M[10,5] = M[5,10] = const
M[9,10] = M[10,9] = const
M[10,6] = M[6,10] = const
M[10,11] = M[11,10] = const
M[11,6] = M[6,11] = const
M[11,7] = M[7,11] = const
M[11,12] = M[12,11] = const
M[7,12] = M[12,7] = const
M[12,8] = M[8,12] = const
 

# set these equations to zero for a barrier in the middle 

M[10, 11] = M[11,10] = 0
M[10, 6] = M[6,10] = 0
M[6,11] = M[11,6] = 0
M[11,7] = M[7,11] = 0
M[2,6] = M[6,2] = 0
M[6,7] = M[7,6] = 0
M[7,3] = M[3,7] = 0
M[6,3] = M[3,6] = 0
M[2,3] = M[3,2] = 0

diag(M) = -rowSums(M);
r = 1e-8;
L = 8e6;

ndemes = 20
int = calculateIntegral(M, q, r, L)
nsamp_per_deme = 50

fileName = paste0("data/sim_model/ndemes_20", L, "_M_", const, "_N_", N[1], "_nsamp_per_deme_", nsamp_per_deme)

lambda_row = matrix(rep(int, each = nsamp_per_deme), byrow = TRUE, ncol = ncol(int)*nsamp_per_deme)
lambda_m = matrix(rep(lambda_row, each = nsamp_per_deme), byrow = FALSE, nrow = nrow(int)*nsamp_per_deme)

# deme
nsamp = nsamp_per_deme*ndemes
#Sobs = matrix(0, nsamp, nsamp)
#for (i in 1:nsamp){
#  for (j in i:nsamp){
#    Sobs[i,j] = rpois(n = 1, lambda = lambda_m[i,j])
#    Sobs[j,i] = Sobs[i,j]
#  }
#}

#write.table(Sobs, file = paste0(fileName, ".sims"), quote = FALSE, row.names = FALSE, col.names = FALSE)
deme_coords = matrix(nrow = ndemes, ncol = 2)
#deme_coords[,1] = c(0, 1, 2, 3, 0.5, 1.5, 2.5, 3.50, 0, 1, 2, 3)
#deme_coords[,2] = c(0, 0, 0, 0, 0.86600, 0.86600, 0.86600, 0.86600, 1.73210, 1.73210, 1.73210, 1.73210)
deme_coords[,1] = c(0, 1, 2, 3, 4, 0.5, 1.5, 2.5, 3.50, 4.5, 0, 1, 2, 3, 4, 0.5, 1.5, 2.5, 3.50, 4.5)
deme_coords[,2] = c(0, 0, 0, 0, 0, 0.86600, 0.86600, 0.86600, 0.86600,  0.86600, 1.73210, 1.73210, 1.73210, 1.73210, 1.73210, 2.5981,2.5981,2.5981,2.5981,2.5981)
write.table(deme_coords, file = paste0(fileName, ".demes"), quote = FALSE, row.names = FALSE, col.names = FALSE)
coords = matrix(rep(deme_coords, each = nsamp_per_deme), byrow = FALSE, nrow = nsamp, ncol = 2)
write.table(coords, file = paste0(fileName, ".coord"), quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(matrix(sapply(1:ndemes, function(x) rep(x, nsamp_per_deme)), byrow=TRUE), file = paste0(fileName, ".ipmap"), quote=FALSE, row.names=FALSE, col.names=FALSE)

#outer_coords = matrix(nrow = 5, ncol = 2, c(0, 0, 0, 1.8, 3.5, 1.8, 3.5, 0, 0 ,0), byrow=TRUE)
outer_coords = matrix(nrow = 5, ncol = 2, c(0, 0, 0, 2.7, 5.5, 2.7, 5.5, 0, 0 ,0), byrow=TRUE)
write.table(outer_coords, file = paste0(fileName, ".outer"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#adjMatrix = matrix(0, nrow = 23, ncol = 2)
#adjMatrix[,1] = c(1, 2, 3, 5, 2, 5, 6, 6, 6, 7, 7, 8, 7, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12)
#adjMatrix[,2] = c(2, 3, 4, 1, 5, 6, 3, 7, 2, 3, 4, 4, 8, 5, 10, 5, 6, 11, 6, 7, 12, 7, 8)
#write.table(adjMatrix, file = paste0(fileName, ".edges"), quote = FALSE, row.names = FALSE, col.names = FALSE)

