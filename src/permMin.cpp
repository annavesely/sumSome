#include <Rcpp.h>
using namespace Rcpp;

//' @title Comparing Vector Elements to a Value
//' @description Internal function. It determines whether all the elements of a vector, eventually except the first one, are equal to a given value.
//' @usage permMin (X, B, truncTo)
//' @param X numeric vector
//' @param B length of \code{X}
//' @param truncTo numeric value
//' @author Anna Vesely
//' @return \code{permMin} returns \code{TRUE} if all elements of \code{X}, excluding the first one, are equal to \code{truncTo}, and returns \code{FALSE} otherwise.
//' @export
// [[Rcpp::export]]



bool permMin (const NumericVector &X, const int &B, const double &truncTo){
  int i = 0;
  bool cond = TRUE;
  while(i < B-1 && cond){
    ++i;
    cond = (X[i] == truncTo);
  }
  return cond;
}