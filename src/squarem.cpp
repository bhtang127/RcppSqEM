// Cpp Version of SQUAREM algorithm
// Author: Bohao Tang, Ravi Varadhan

#include <Rcpp.h>
#include <cmath>

// useful class for function evaluation
// motivated by RcppDE
class EvalBase {
public:
    EvalBase() : neval(0) {};
    virtual Rcpp::NumericVector eval(SEXP par) = 0;
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
    Rcpp::NumericVector eval(SEXP par) {
        neval++;
        return funptr(par, env);
    }
private:
    funcPtr funptr;
    SEXP env;
};

// useful function
bool any_naC(Rcpp::NumericVector x) {
  return Rcpp::is_true(Rcpp::any(Rcpp::is_na(x)));
}

// [[Rcpp::export]]
Rcpp::List squarem1_cpp(Rcpp::NumericVector par,
                        SEXP fixptfn, SEXP objfn,
                        const Rcpp::List & ctrl,
                        SEXP env) {

  int method            = Rcpp::as<int>(ctrl["method"]);
  unsigned long maxiter = Rcpp::as<unsigned long>(ctrl["maxiter"]);
  double tol            = Rcpp::as<double>(ctrl["tol"]);
  double step_min       = Rcpp::as<double>(ctrl["step.min0"]);
  double step_max0      = Rcpp::as<double>(ctrl["step.max0"]);
  double step_max       = Rcpp::as<double>(ctrl["step.max0"]);
  double mstep          = Rcpp::as<double>(ctrl["mstep"]);
  double objfn_inc      = Rcpp::as<double>(ctrl["objfn.inc"]);
  int trace             = Rcpp::as<int>(ctrl["trace"]);
  int intermed          = Rcpp::as<int>(ctrl["intermed"]);
  int np                = par.size();

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

  int iter = 0, converge = 1, extrap, row_save=0;
  double sr2, sq2, sv2, srv, alpha;
  Rcpp::NumericVector lold, lnew(1), p = Rcpp::clone(par);
  Rcpp::NumericVector p1, q1, p2, q2, p_new, p_aux;
  Rcpp::NumericVector obj_save(maxiter);
  Rcpp::NumericMatrix par_save(maxiter, np);

  if (trace) Rcpp::Rcout << "Squarem-1 \n" ;

  lold = obj->eval(p);
  if (trace) Rcpp::Rcout << "Objective fn:" << lold[0] << "\n";

  obj_save[row_save] = lold[0];
  for (int i=0; i<np; i++) par_save(row_save, i) = p[i];
  row_save++;

  while (fn->getNumEvals() <= maxiter){
    extrap = 1;
    p1 = fn->eval(p);
    q1 = p1 - p;
    sr2 = Rcpp::sum(q1 * q1);
    if (std::sqrt(sr2) < tol) break;

    p2 = fn->eval(p1);
    q2 = p2 - p1;
    sq2 = Rcpp::sum(q2 * q2);
    if (std::sqrt(sq2) < tol) break;

    sv2 = Rcpp::sum((q2-q1) * (q2-q1));
    srv = Rcpp::sum(q1 * (q2-q1));

    if (method == 1) alpha = -srv/sv2;
    else if (method == 2) alpha = -sr2/srv;
    else alpha = std::sqrt(sr2/sv2);

    alpha = std::max(step_min, std::min(step_max, alpha));
    p_new = p + 2*alpha * q1 +  std::pow(alpha,2) * (q2-q1);
    if (std::abs(alpha - 1) > 0.01) {
      try{
        p_new = fn->eval(p_new);
        if (any_naC(p_new))
          throw "Value Error";
      } catch(...) {
        p_new = NA_REAL;
      }
    }

    if (any_naC(p_new)) {
      p_new = p2;
      try {lnew = obj->eval(p_new);}
      catch(...) {lnew(0) = NA_REAL;}
      if (alpha == step_max)
        step_max = std::max(step_max0, step_max/mstep);
      alpha = 1;
      extrap = 0;
    } else {
      if (std::isfinite(objfn_inc)) {
        try {lnew = obj->eval(p_new);}
        catch(...) {lnew(0) = NA_REAL;}
      } else lnew = lold;
      if (any_naC(lnew) || lnew[0] > lold[0] + objfn_inc){
        p_new = p2;
        try {lnew = obj->eval(p_new);}
        catch(...) {lnew(0) = NA_REAL;}
        if (alpha == step_max)
          step_max = std::max(step_max0, step_max/mstep);
        alpha = 1;
        extrap = 0;
      }
    }

    if (alpha == step_max) step_max = mstep * step_max;
    if (step_min < 0 && alpha == step_min) step_min = mstep*step_min;

    p = p_new;
    if (! any_naC(lnew)) lold = lnew;
    if (trace) {
      Rcpp::Rcout << "Objective fn:" << lnew[0];
      Rcpp::Rcout << " Extrapolation: " << extrap;
      Rcpp::Rcout << " Steplength: " << alpha << "\n";
    }
    if (intermed) {
      obj_save[row_save] = lnew[0];
      for (int i=0; i<np; i++) par_save(row_save, i) = p[i];
      row_save++;
    }

    iter++;
  }

  if (fn->getNumEvals() >= maxiter) converge = 0;
  if (std::isfinite(objfn_inc)) {
    lold = obj->eval(p);
  }

  if (!intermed)
    return Rcpp::List::create(Rcpp::Named("par")         = p,
                              Rcpp::Named("value.objfn") = lold,
                              Rcpp::Named("iter")        = iter,
                              Rcpp::Named("fpevals")     = fn->getNumEvals(),
                              Rcpp::Named("objfevals")   = obj->getNumEvals(),
                              Rcpp::Named("convergence") = converge
                              );
  else
    return Rcpp::List::create(Rcpp::Named("par")           = p,
                              Rcpp::Named("value.objfn")   = lold,
                              Rcpp::Named("iter")          = iter,
                              Rcpp::Named("fpevals")       = fn->getNumEvals(),
                              Rcpp::Named("objfevals")     = obj->getNumEvals(),
                              Rcpp::Named("convergence")   = converge,
                              Rcpp::Named("stored_objval") = obj_save[Rcpp::Range(0,row_save)],
                              Rcpp::Named("stored_par")    = par_save(Rcpp::Range(0,row_save),
                                                                      Rcpp::Range(0,np))
                              );
}

// [[Rcpp::export]]
Rcpp::List squarem2_cpp(Rcpp::NumericVector par,
                        SEXP fixptfn,
                        const Rcpp::List & ctrl,
                        SEXP env) {

  int method            = Rcpp::as<int>(ctrl["method"]);
  unsigned long maxiter = Rcpp::as<unsigned long>(ctrl["maxiter"]);
  double tol            = Rcpp::as<double>(ctrl["tol"]);
  double step_min       = Rcpp::as<double>(ctrl["step.min0"]);
  double step_max0      = Rcpp::as<double>(ctrl["step.max0"]);
  double step_max       = Rcpp::as<double>(ctrl["step.max0"]);
  double mstep          = Rcpp::as<double>(ctrl["mstep"]);
  double kr             = Rcpp::as<double>(ctrl["kr"]);
  int trace             = Rcpp::as<int>(ctrl["trace"]);

  EvalBase *fn = NULL;
  if (TYPEOF(fixptfn) == EXTPTRSXP) {
    fn = new EvalCompiled(fixptfn, env); // working
  } else {
    fn = new EvalR(fixptfn, env);
  }

  int iter = 0, converge = 1, extrap;
  double sr2, sq2, sv2, srv, alpha, res, parnorm, kres;
  Rcpp::NumericVector p = Rcpp::clone(par);
  Rcpp::NumericVector p1, q1, p2, q2, p_new, p_aux;

  if (trace) Rcpp::Rcout << "Squarem-2 \n" ;

  while (fn->getNumEvals() <= maxiter){
    extrap = 1;
    p1 = fn->eval(p);
    if (any_naC(p1)) break;
    q1 = p1 - p;
    sr2 = Rcpp::sum(q1 * q1);
    if (std::sqrt(sr2) < tol) break;

    p2 = fn->eval(p1);
    if (any_naC(p2)) break;
    q2 = p2 - p1;
    sq2 = Rcpp::sum(q2 * q2);
    res = sq2;
    if (std::sqrt(sq2) < tol) break;

    sv2 = Rcpp::sum((q2-q1) * (q2-q1));
    srv = Rcpp::sum(q1 * (q2-q1));

    if (method == 1) alpha = -srv/sv2;
    else if (method == 2) alpha = -sr2/srv;
    else alpha = std::sqrt(sr2/sv2);

    alpha = std::max(step_min, std::min(step_max, alpha));
    p_new = p + 2*alpha * q1 +  std::pow(alpha,2) * (q2-q1);

    if (std::abs(alpha - 1) > 0.01) {
      try{
        p_aux = fn->eval(p_new);
        if (any_naC(p_aux))
          throw "Value Error";
      } catch(...) {
        p_aux = NA_REAL;
      }

      if (any_naC(p_aux)) {
        p_new = p2;
        if (alpha == step_max)
          step_max = std::max(step_max0, step_max/mstep);
        alpha = 1;
        extrap = 0;
      } else {
        res = std::sqrt(Rcpp::sum((p_aux - p_new) * (p_aux - p_new)));
        parnorm = std::sqrt(Rcpp::sum(p2 * p2) / p2.size());
        kres = kr * (1 + parnorm) + sq2;
        if (res <= kres) p_new = p_aux;
        else p_new = p2;
        if (res > kres) {
          if (alpha == step_max)
            step_max = std::max(step_max0, step_max/mstep);
          alpha = 1;
          extrap = 0;
        }
      }
    }

    if (alpha == step_max) step_max = mstep * step_max;
    if (step_min < 0 && alpha == step_min) step_min = mstep*step_min;

    p = p_new;
    if (trace) {
      Rcpp::Rcout << "Residual:" << res;
      Rcpp::Rcout << " Extrapolation: " << extrap;
      Rcpp::Rcout << " Steplength: " << alpha << "\n";
    }

    iter++;
  }

  return Rcpp::List::create(Rcpp::Named("par")           = p,
                            Rcpp::Named("value.objfn")   = NA_REAL,
                            Rcpp::Named("iter")          = iter,
                            Rcpp::Named("fpevals")       = fn->getNumEvals(),
                            Rcpp::Named("objfevals")     = 0,
                            Rcpp::Named("convergence")   = converge
                            );
}
