squarem <- function(par, fixptfn, objfn, ... , control=list()) {
  #
  # The Master function for SQUAREM acceleration of "ANY" fixed-point iteration
  # A wrapper that calls appropriate function below (squarem1, squarem2, cyclem1, cyclem2)
  #
  # Author:  Ravi Varadhan, Johns Hopkins University
  #
  # See Varadhan and Roland (Scandinavian J Statistics 2008)
  # See Roland, Varadhan and Frangakis (Applied Numerical Mathematics 2007)
  #
  # Last modified: June 25, 2010
  # Last modified: August 09, 2010
  # Last modified: December 31, 2010
  ###########################################################################
  #
  # par = starting value of parameter vector
  # fixptfn = fixed-point iteration F(x)
  # for which the solution: F(x*) = x* is sought
  # objfn = underlying objective function which is minimized at x*
  # method = 1, 2, or 3, indicating the type of steplength to be used
  #
  #
  control.default <- list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,
                          tol=1.e-07, maxiter=1500, trace=FALSE, intermed=FALSE)
  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])

  ctrl <- modifyList(control.default, control)

  if (ctrl$K > 1 & !(ctrl$method %in% c("rre", "mpe")) ) ctrl$method <- "rre"
  if (ctrl$K == 1 & !(ctrl$method %in% c(1,2,3)) ) ctrl$method <- 3

  if (!missing(objfn) ) {
    if (ctrl$K == 1) sqobj <- squarem1(par, fixptfn, objfn, ... , control=ctrl)
    else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe")) sqobj <- cyclem1(par, fixptfn, objfn, ... , control=ctrl)
  } else {
    if (ctrl$K == 1) sqobj <- squarem2(par, fixptfn, ... , control=ctrl)
    else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe"))  sqobj <- cyclem2(par, fixptfn, ... , control=ctrl)
  }
  return(sqobj)
}


###########################################################
# Wraps for cpp version of the SQUAREM acceleration
#
# Date: Jan 30, 2019
# Last modified: Jan 30 2019
###########################################################
squarem1 <- function(par, fixptfn, objfn, ... , control=list()) {
  # par = starting value of parameter vector
  # fixptfn = fixed-point iteration F(x)
  # for which the solution: F(x*) = x* is sought
  # objfn = underlying objective function which is minimized at x*
  #
  #
  control.default <- list(K=1, square=TRUE, method=3,
                          step.min0=1, step.max0=1, mstep=4,
                          kr=1, objfn.inc=1, tol=1.e-07,
                          maxiter=1500, trace=FALSE, intermed=FALSE)

  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)

  env <- list2env(list(...))

  outC <- squarem1_cpp(par, fixptfn, objfn, ctrl, env)
  if (!is.null(names(par)))
    names(outC$p) <- names(par)

  return(outC)
}

###################################################################
squarem2 <- function(par, fixptfn, ... , control=list() ) {

  # par = starting value of parameter vector
  # fixptfn = fixed-point iteration F(x)
  # for which the solution: F(x*) = x* is sought
  # method = 1, 2, or 3, indicating the type of steplength to be used
  #
  # See Varadhan and Roland (Scandinavian J Statistics 2008)
  #
  control.default <- list(K=1, square=TRUE, method=3,
                          step.min0=1, step.max0=1, mstep=4,
                          kr=1, objfn.inc=1, tol=1.e-07,
                          maxiter=1500, trace=FALSE, intermed=FALSE)

  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)

  env <- list2env(list(...))

  outC <- squarem2_cpp(par, fixptfn, ctrl, env)
  if (!is.null(names(par)))
    names(outC$p) <- names(par)

  return(outC)
}


#############################################
# Squared-Cycled wraps (still working)
#
# Latest modofications: Jan 30, 2019
#
#############################################
cyclem1 <- function(par, fixptfn, objfn, control=list(), ...) {
  "working"
}

##################################################################
cyclem2 <- function(par, fixptfn, control=list(), ...) {
  "working"
}
###################################################################



