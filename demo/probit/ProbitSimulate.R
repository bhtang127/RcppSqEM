
ProbitSimulate <- function(beta.vec, X) {
    n <- nrow(X)
    lin.pred <- X%*%beta.vec
    pp <- pnorm(lin.pred)
    UU <- runif(n)
    yy <- rep(0, n)
    yy[UU < pp] <- 1
    return(yy)
}
