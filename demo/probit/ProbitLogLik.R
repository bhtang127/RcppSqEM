
ProbitLogLik <- function(beta.hat, X, y) {

  fvals1 <- pnorm(X%*%beta.hat, log.p=TRUE)
  fvals2 <- pnorm(X%*%beta.hat, log.p=TRUE, lower.tail=FALSE)
  ans <- sum(y*fvals1) + sum((1 - y)*fvals2)
  return(ans)
}
