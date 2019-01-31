#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

void testerror_C(const std::string name) {
  const std::string err = name + " haha\n";
  Rcpp::Rcout << "haha" <<std::endl;
  ::Rf_error(err.c_str());
}
// [[Rcpp::export]]
void open(){
  testerror_C("haha");
}
