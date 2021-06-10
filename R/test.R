library(Rcpp)
Rcpp::sourceCpp('src/gmm.cpp')

d <- 10
pvec <- 1/1:d
pvec <- pvec/sum(pvec)
mus <- 1:d

n <- 10000
X <- rmultinom(1:d, n, pvec)
Y <- rep(0,n)
sigma = 0.3

tot = 0
for(i in 1:d){
  for(j in (tot+1):(tot+X[i])){
    Y[j] <- rnorm(1, mean = mus[i], sd = sigma)
  }
  tot = tot + X[i]
}
plot(density(Y, bw = 0.1), type = "l", col = "blue", lwd = 2)


fit = gmm(Y, d, max_iter = 10000, tol = 1e-5)
plot(density(Y, bw = 0.1), type = "l", col = "blue", lwd = 2)
for(i in 1:d){
  # overlay density of specific gaussian
  mu = fit$mu[i]
  sigma = fit$sd[i]
  p = fit$pi[i]
  x <- seq(-1,d+2, by = 0.01)
  points(x,p*dnorm(x, mean = mu, sd = sigma), type = "l", col = rainbow(d)[i])
}


