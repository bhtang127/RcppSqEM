sourceCpp(code='
  #include <Rcpp.h>

  // [[Rcpp::export]]
  SEXP test_C(SEXP fn, Rcpp::NumericVector par, SEXP rho) {
    SEXP f = ::Rf_lang3(fn, par, R_DotsSymbol); // this could be done with Rcpp
    SEXP sf = ::Rf_eval(f, rho);
    double f_result = Rcpp::NumericVector(sf)[0];
    return Rcpp::wrap(f_result + par.size()*10000);
  }'
)

sourceCpp(file = "src/testC.cpp")

f = function(p, X=0, y=1){
  p*y + X
}

test = function(p, f, ...){
  env = list2env(list(...))
  test_C(f, p, env)
}
