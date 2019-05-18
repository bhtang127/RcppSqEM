\name{squarem}
\alias{squarem}
\title{Squared Extrapolation Methods for Accelerating Slowly-Convergent
Fixed-Point Iterations}
\description{Globally-convergent, partially monotone, acceleration schemes for
  accelerating the convergence of \emph{any} smooth, monotone, slowly-converging
  contraction mapping. It can be used to accelerate the convergence of a wide
  variety of iterations including the expectation-maximization (EM) algorithms
  and its variants, majorization-minimization (MM) algorithm, power method for
  dominant eigenvalue-eigenvector, Google's page-rank algorithm, and
  multi-dimensional scaling.
}
\usage{
  squarem(par, fixptfn, objfn, ... , control=list())
  }
\arguments{
    \item{par}{A vector of parameters denoting the initial guess for the
      fixed-point.}

    \item{fixptfn}{A vector function, $F$ that denotes the fixed-point
      mapping.  This function is the most essential input in the package.
      It should accept a parameter vector as input and should return a
      parameter vector of same length. This function defines the fixed-point
      iteration: \eqn{x_{k+1} = F(x_k)}{x[k+1] = F(x[k])}.
      In the case of EM algorithm, \eqn{F}{F} defines a single E and M step.}

    \item{objfn}{This is a scalar function, \eqn{L}{L}, that denotes
      a ''merit'' function which attains its local minimum at the fixed-point
      of \eqn{F}{F}.  This function should accept a parameter vector as
      input and should return a scalar value.  In the EM algorithm, the merit
      function \eqn{L}{L} is the negative of log-likelihood. In some problems, a natural
      merit function may not exist, in which case the algorithm works with
      only \code{fixptfn}. The merit function function \code{objfn} does not
      have to be specified, even when a natural merit function is available,
      especially when its computation is expensive.}

    \item{control}{A list of control parameters specifing any changes to
      default values of algorithm control parameters.  Full names of control list elements must be specified, otherwise, user-specifications are ignored.
	See *Details*.}
    \item{...}{Arguments passed to \code{fixptfn} and  \code{objfn}.}
}

\value{A list with the following components:
\item{par}{Parameter, \eqn{x*}{x*} that are the fixed-point of \eqn{F}{F}
  such that \eqn{x* = F(x*)}{x* = F(x*)}, if convergence is successful.}

\item{value.objfn}{The value of the objective function $L$ at termination.}

\item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}

\item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.  \code{0}
  indicates successful convergence, whereas \code{1} denotes failure to
  converge.}
}

\details{
 The function \code{squarem} is a general-purpose algorithm for accelerating
 the convergence of \emph{any} slowly-convergent (smooth) fixed-point iteration.  Full names of

Default values of \code{control} are:
    \code{K=1},
    \code{method=3},
    \code{square=TRUE},
    \code{step.min0=1},
    \code{step.max0=1},
    \code{mstep=4},
    \code{objfn.inc=1},
    \code{kr=1},
    \code{tol=1e-07},
    \code{maxiter=1500},
    \code{trace=FALSE}.

\describe{
\item{\code{K}}{An integer denoting the order of the SQUAREM scheme.
      Default is 1, which is a first-order scheme developed in Varadhan and
      Roland (2008). Our experience is that first-order schemes are adequate
      for most problems. \code{K=2,3} may provide greater speed in some
      problems, although they are less reliable than the first-order schemes.}

\item{\code{method}}{Either an integer or a character variable that denotes the
      particular SQUAREM scheme to be used.  When \code{K=1}, method should be
      an integer, either 1, 2, or 3.  These correspond to the 3 schemes
      discussed in Varadhan and Roland (2008).  Default is \code{method=3}.
      When K > 1, method should be a character string, either \code{''RRE''}
      or \code{''MPE''}.  These correspond to the reduced-rank extrapolation
      or squared minimal=polynomial extrapolation (See Roland, Varadhan, and
      Frangakis (2007)).  Default is ''RRE''.}
\describe{
	\item{\code{square}}{A logical variable indicating whether or not a squared
	  extrapolation scheme should be used.  Our experience is that the
	  squared extrapolation schemes are faster and more stable than the
	  unsquared schemes. Hence, we have set the default as TRUE.}
	\item{\code{step.min0}}{A scalar denoting the minimum steplength taken by a
	  SQUAREM algorithm.  Default is 1.  For contractive fixed-point
	  iterations (e.g. EM and MM), this defualt works well.  In problems
	  where an eigenvalue of the Jacobian of $F$ is outside of the
	  interval \eqn{(0,1)}{(0,1)}, \code{step.min0} should be less than 1
	  or even negative in some cases.}
	\item{\code{step.max0}}{A positive-valued scalar denoting the initial value
	  of the maximum steplength taken by a SQUAREM algorithm.
	  Default is 1.  When the steplength computed by SQUAREM exceeds
	  \code{step.max0}, the steplength is set equal to step.max0, but
	  then step.max0 is increased by a factor of mstep.}
	\item{\code{mstep}}{A scalar greater than 1.  When the steplength computed
	  by SQUAREM exceeds  \code{step.max0}, the steplength is set equal
	  to  \code{step.max0}, but  \code{step.max0} is increased by a factor
	  of  \code{mstep}. Default is 4.}
	\item{\code{objfn.inc}}{A non-negative scalar that dictates the degree of
	  non-montonicity.  Default is 1.  Set  \code{objfn.inc = 0} to
	  obtain monotone convergence. Setting  \code{objfn.inc = Inf} gives a
	  non-monotone scheme.  In-between values result in partially-monotone
	  convergence.}
	\item{\code{kr}}{A non-negative scalar that dictates the degree of
	  non-montonicity.  Default is 1.  Set  \code{kr = 0} to obtain
	  monotone convergence. Setting  \code{kr = Inf} gives a non-monotone
	  scheme.  In-between values result in partially-monotone convergence.  This parameter is only used when  \code{objfn} is not specified by user.}
	\item{\code{tol}}{A small, positive scalar that determines when iterations
	  should be terminated.  Iteration is terminated when
	  \eqn{||x_k - F(x_k)|| \leq tol}{abs(x[k] - F(x[k])) <= tol}.
	  Default is \code{1.e-07}.}
	\item{\code{maxiter}}{An integer denoting the maximum limit on the number
	  of evaluations of  \code{fixptfn}, \eqn{F}{F}.  Default is 1500.}
	}
\item{trace}{A logical variable denoting whether some of the intermediate
      results of iterations should be displayed to the user.
      Default is \code{FALSE}.}
}
}
\references{
R Varadhan and C Roland (2008), Simple and globally convergent numerical
  schemes for accelerating the convergence of any EM algorithm,
  \emph{Scandinavian Journal of Statistics}, 35:335-353.

C Roland, R Varadhan, and CE Frangakis (2007), Squared polynomial
  extrapolation methods with cycling: an application to the positron emission
  tomography problem, \emph{Numerical Algorithms}, 44:159-172.

}


\keyword{optimization, EM algorithm}
