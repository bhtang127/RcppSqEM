// Cpp Version of SQUAREM algorithm
// Author: Bohao Tang, Ravi Varadhan

#include <Rcpp.h>

// useful class for function evaluation
// motivated by RcppDE
class EvalBase {
public:
    EvalBase(const std::string fname_) : neval(0), fname(fname_) {};
    virtual double eval(SEXP par) = 0;
    unsigned long getNumEvals() { return neval; }
protected:
    unsigned long int neval;
    const std::string fname;
};

class EvalR : public EvalBase {
public:
    EvalR(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {} 
    Rcpp::NumericVector eval(SEXP par) {
        neval++;
        return defaultfn(par);
    }
private:
    SEXP fcall, env;
    Rcpp::NumericVector defaultfn(SEXP par) {               // essentialy same as the old evaluate
        try {
            SEXP fn = ::Rf_lang3(fcall, par, R_DotsSymbol); // this could be done with Rcpp
            SEXP res = ::Rf_eval(fn, env);                  // but is still a lot slower right now
            return(res); 
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...) { 
            const std::string err = fname + "Function Evaluation Error (Unknown Reason)";
            ::Rf_error(err.c_str()); 
        }
        return NA_REAL; 
    }
};

typedef double (*funcPtr)(SEXP, SEXP);
typedef double (*funcPtrTest)(SEXP);

// For compiled version of fixpt/obj fn, still working

class EvalCompiled : public EvalBase {
public:
    EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP env_) {
        funptr = *(xptr);
        env = env_;
    };
    EvalCompiled(SEXP xps, SEXP env_) {
        Rcpp::XPtr<funcPtr> xptr(xps);
        funptr = *(xptr);
        env = env_;
    };
    double eval(SEXP par) {
        neval++;
        return funptr(par, env);
    }
private:
    funcPtr funptr;
    SEXP env;
};

// [[Rcpp::export]]
Rcpp::List squarem1_cpp(Rcpp::NumericVector par, 
                        SEXP fixptfn, SEXP objfn,
                        const Rcpp::List & ctrl,
                        SEXP env) {

#if RCPP_DEV_VERSION >= RcppDevVersion(0,12,17,1)
    Rcpp::SuspendRNGSynchronizationScope rngScope;
#endif

    int method       = Rcpp::as<int>(ctrl["method"]);
    int maxiter      = Rcpp::as<int>(ctrl["maxiter"]);
    double tol       = Rcpp::as<double>(ctrl["tol"]);
    double step.min  = Rcpp::as<double>(ctrl["step.min0"]);
    double step.max0 = Rcpp::as<double>(ctrl["step.max0"]);
    double step.max  = Rcpp::as<double>(ctrl["step.max0"]);
    double mstep     = Rcpp::as<double>(ctrl["mstep"]);
    double objfn.inc = Rcpp::as<double>(ctrl["objfn.inc"]);
    int trace        = Rcpp::as<int>(ctrl["trace"]);
    int intermed     = Rcpp::as<int>(ctrl["intermed"]);
    int np           = par.size();

    EvalBase *fn = NULL;                   
    if (TYPEOF(fixptfn) == EXTPTRSXP) {                
        fn = new EvalCompiled(fixptfn, env); // working
    } else {                                         
        fn = new EvalR(fixptfn, env); 
    }

    EvalBase *obj = NULL;                   
    if (TYPEOF(objfn) == EXTPTRSXP) {                
        obj = new EvalCompiled(objfn, env); // working
    } else {                                         
        obj = new EvalR(objfn, env); 
    }

    int iter, converge = 1;
    Rcpp::NumericVector lold, p = Rcpp::clone(par);

    if (trace) Rcpp::Rcout << "Squarem-1" << endl;

    lold = obj->eval(p);
    if (trace) Rcpp::Rcout << "Objective fn:" << lold[0] << endl;

    


    return Rcpp::List::create(Rcpp::Named("bestmem")   = 123,   // and return a named list with results to R
                              ); 
}

// if (trace) cat("Squarem-1 \n")

// if (missing(objfn)) stop("\n squarem2 should be used if objective function is not available \n\n")

// iter <- 1
// objval <- rep(NA,1)
// p <- par

// lold <- objfn(p, ...)
// leval <- 1
// if (trace) cat("Objective fn: ", lold, "\n")
// feval <- 0
// conv <- TRUE
// p.inter <- c(p, lold)

// while (feval < maxiter) {
// 	extrap <- TRUE
// 	p1 <- try(fixptfn(p, ...),silent=TRUE)
// 	feval <- feval + 1
// 	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) stop("Error in function evaluation")
// 	q1 <- p1 - p
// 	sr2 <- crossprod(q1)
// 	if (sqrt(sr2) < tol) break

// 	p2 <- try(fixptfn(p1, ...),silent=TRUE)
// 	feval <- feval + 1
// 	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) stop("Error in function evaluation")

// 	q2 <- p2 - p1
// 	sq2 <- sqrt(crossprod(q2))
// 	if (sq2 < tol) break
// 	sv2 <- crossprod(q2-q1)
// 	srv <- crossprod(q1, q2-q1)

// 	alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2)) 

//  	alpha <- max(step.min, min(step.max, alpha))
// 	p.new <- p + 2*alpha*q1 + alpha^2*(q2-q1)
// 	if (abs(alpha - 1) > 0.01 ) {
// 		p.new <- try(fixptfn(p.new , ...),silent=TRUE)   # stabilization step
// 		feval <- feval + 1
// 	}

// 	if (class(p.new) == "try-error" | any(is.nan(p.new))) {
// 	p.new <- p2
// 	lnew <- try(objfn(p2, ...), silent=TRUE)
// 	leval <- leval + 1
// 	if (alpha == step.max) step.max <- max(step.max0, step.max/mstep)
// 	alpha <- 1
// 	extrap <- FALSE
// 	} else {
// 		if (is.finite(objfn.inc)) {
// 			lnew <- try(objfn(p.new, ...), silent=TRUE)
// 			leval <- leval + 1
// 		} else lnew <- lold
// 		if (class(lnew) == "try-error" | is.nan(lnew) | 
// 		(lnew > lold + objfn.inc)) {
// 			p.new <- p2
// 			lnew <- try(objfn(p2, ...), silent=TRUE)
// 			leval <- leval + 1
// 			if (alpha==step.max) step.max <- max(step.max0, step.max/mstep)
// 			alpha <- 1
// 			extrap <- FALSE
// 		}
// 	}	

// 	if (alpha == step.max) step.max <- mstep*step.max
// 	if (step.min < 0 & alpha == step.min) step.min <- mstep*step.min

// 	p <- p.new
// 	if (!is.nan(lnew)) lold <- lnew 
// 	if (trace) cat("Objective fn: ", lnew, "  Extrapolation: ", extrap, "  Steplength: ", alpha, "\n")
//   if(intermed) p.inter <- rbind(p.inter, c(p, lnew))

// 	iter <- iter+1	
// }
// 	if (feval >= maxiter) conv <- FALSE
// 	if (is.infinite(objfn.inc)) {
// 		lold <- objfn(p, ...)
// 		leval <- leval + 1
// 		}

// rownames(p.inter) <- NULL
// if (!intermed) return(list(par=p, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv))
//    else return(list(par=p, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv, p.inter=p.inter))