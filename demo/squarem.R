
require("SQUAREM")  # master funcion and the 4 different SQUAREM algorithms
require("RcppSqEM")
require("setRNG")   # only needed to reproduce the same results
require(profvis)
library(microbenchmark)

# data for Poisson mixture estimation  from Hasselblad (JASA 1969)
poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))

y <- poissmix.dat$freq
tol <- 1.e-08

# generate a random initial guess for 3 parameters
setRNG(list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=123))
p0 <- c(runif(1),runif(2,0,4))

# fixptfn = function that computes a single update of fixed-point iteration
# objfn = underlying objective function that needs to be maximized
#

# EM algorithm

poissmix.em <- function(p,y) {
  # The fixed point mapping giving a single E and M step of the EM algorithm
  #
  pnew <- rep(NA,3)
  i <- 0:(length(y)-1)
  zi <- p[1]*exp(-p[2])*p[2]^i / (p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i)
  pnew[1] <- sum(y*zi)/sum(y)
  pnew[2] <- sum(y*i*zi)/sum(y*zi)
  pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
  p <- pnew
  return(pnew)
}

poissmix.loglik <- function(p,y) {
  # Objective function whose local minimum is a fixed point \
  # negative log-likelihood of binary poisson mixture
  i <- 0:(length(y)-1)
  loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) +
                    (1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
  return ( -sum(loglik) )
}

# generate a random initial guess for 3 parameters
setRNG(list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=123))
p0 <- c(runif(1),runif(2,0,4))

# fixptfn = function that computes a single update of fixed-point iteration
# objfn = underlying objective function that needs to be minimized  (here it is the negative log-likelihood)
#

profvis({
for(i in 1:100){
  # EM algorithm
  em <- fpiter(p=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))

  # First-order SQUAREM algorithm with SqS3 method
  pureR_S3 <- SQUAREM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))

  # First-order SQUAREM algorithm with SqS3 method and Rcpp
  Rcpp_S3 <- RcppSqEM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))

  # First-order SQUAREM algorithm with SqS2 method
  pureR_S2 <- SQUAREM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(method=2, tol=tol))

  # First-order SQUAREM algorithm with SqS2 method and Rcpp
  Rcpp_S2 <- RcppSqEM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(method=2, tol=tol))

  # First-order SQUAREM algorithm with SqS3 method; non-monotone
  pureR_kr <- SQUAREM::squarem(par=p0,y=y, fixptfn=poissmix.em, control=list(tol=tol, kr=0.1))

  # First-order SQUAREM algorithm with SqS3 method and Rcpp; non-monotone
  Rcpp_kr <- RcppSqEM::squarem(par=p0,y=y, fixptfn=poissmix.em, control=list(tol=tol, kr=0.1))

  # Comparison of converged parameter estimates
  par.mat <- rbind(em$par, pureR_S3$par, pureR_S2$par, pureR_kr$par,
                   Rcpp_S3$par, Rcpp_S2$par, Rcpp_kr$par)
  par.mat
}
})

microbenchmark(
  # EM algorithm
  em = fpiter(p=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol)),

  # First-order SQUAREM algorithm with SqS3 method
  pureR_S3 = SQUAREM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol)),

  # First-order SQUAREM algorithm with SqS3 method and Rcpp
  Rcpp_S3 = RcppSqEM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol)),

  # First-order SQUAREM algorithm with SqS2 method
  pureR_S2 = SQUAREM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(method=2, tol=tol)),

  # First-order SQUAREM algorithm with SqS2 method and Rcpp
  Rcpp_S2 = RcppSqEM::squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(method=2, tol=tol)),

  # First-order SQUAREM algorithm with SqS3 method; non-monotone
  pureR_kr = SQUAREM::squarem(par=p0,y=y, fixptfn=poissmix.em, control=list(tol=tol, kr=0.1)),

  # First-order SQUAREM algorithm with SqS3 method and Rcpp; non-monotone
  Rcpp_kr = RcppSqEM::squarem(par=p0,y=y, fixptfn=poissmix.em, control=list(tol=tol, kr=0.1)),

  times = 50
)

