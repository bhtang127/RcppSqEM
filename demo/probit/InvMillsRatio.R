

InvMillsRatio <- function(x, lt) {

  ans <- exp(dnorm(x, log=TRUE) - pnorm(-x, lower.tail=lt, log.p=TRUE))
  return(ans)
}





