estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
mean = 0.00001
var = 0.0000000001
params = estBetaParams(mean, var)
hist(rbeta(10000, shape1 = params$alpha, shape2 = params$beta), 1000, xlim = c(0, 10*mean), 200000)

