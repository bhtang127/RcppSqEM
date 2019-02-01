
ProbitUpdate <- function(beta.hat, X, y) {

  linear.pred <- as.vector(X%*%beta.hat)
  ZZ <- linear.pred + y*InvMillsRatio(linear.pred, lt=FALSE) - (1 - y)*InvMillsRatio(linear.pred, lt=TRUE)

  newbeta.hat <- solve(crossprod(X, X), crossprod(X, ZZ))
  return(newbeta.hat)
}
