#include <Rcpp.h>
using namespace Rcpp;


//' @title Compare Vector Elements to Value
//' @description Internal function, called in \code{sum.internal}.
//' It determines whether all the elements of a vector, eventually except the first one, are equal to a given value.
//' @usage permMin(X, B, truncTo)
//' @param X numeric vector.
//' @param B length of \code{X}.
//' @param truncTo numeric value.
//' @author Anna Vesely.
//' @return It returns \code{TRUE} if all elements of \code{X}, excluding the first one, are equal to \code{truncTo},
//' and returns \code{FALSE} otherwise.
//' @noRd
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





//' @title First Selected Elements
//' @description Internal function, called in \code{sortSum}.
//' It selects the first elements of a vector that do not exceed a given threshold.
//' @usage firstInS(z, s, X, f)
//' @param z number of selected elements.
//' @param s integer threshold.
//' @param X integer vector.
//' @param f length of \code{X}.
//' @author Anna Vesely.
//' @return It returns a ligical vector of length \code{f}, where \code{TRUE} values correspond
//' to the first \code{z} elements of \code{X} that are not greater than \code{s}.
//' @noRd

LogicalVector firstInS (const int &z, const int &s, const IntegerVector &X, const int &f){
  LogicalVector sel (f, FALSE);
  int found = 0;
  int i = 0;
  while(found < z){
    if(X[i] <= s){
      sel[i] = TRUE;
      ++found;
    }
    ++i;
  }
  return sel;
}





//' @title Element Selection
//' @description Internal function, called in \code{buildMatrices}.
//' It selects the elements of a vector that appear in a second vector.
//' @usage selIndices (indices, X, H, f)
//' @param indices integer vector.
//' @param X integer vector.
//' @param H length of \code{indices}.
//' @param f length of \code{X}.
//' @author Anna Vesely.
//' @return It returns a logical vector of length \code{f}, where \code{TRUE} values correspond
//' to elements of \code{X} that appear in \code{indices}.
//' @noRd

LogicalVector selIndices (const IntegerVector &indices, const IntegerVector &X,
                          const int &H, const int &f){
  LogicalVector out (f, FALSE); // TRUE if X[i] in indices
  int found = 0; // number of indices found in X
  int i = 0;
  while(i < f && found < H){
    for(int h = 0; h <= H-1; ++h){
      if(X[i] == indices[h]){
        out[i] = TRUE;
        ++ found;
        break;
      }
    }
    ++i;
  }
  return out;
}




//' @title Sign of Vector Elements
//' @description Internal function, called in \code{findCol}.
//' It determines whether all the elements of a vector are non-negative.
//' @usage allPos(X, B)
//' @param X numeric vector.
//' @param B length of \code{X}.
//' @author Anna Vesely.
//' @return It returns \code{TRUE} if all elements of \code{X} are non-negative, and \code{FALSE} otherwise.
//' @noRd

bool allPos(const NumericVector &X, const int &B){
  bool out = TRUE;
  for (int i = 0; i <= B-1; ++i) {
    if (X[i] < 0){
      out = FALSE;
      break;
    }
  }
  return out;
}





//' @title Upper Bound Behavior
//' @description Internal function, called in \code{cumulativeMatrices}.
//' It characterizes the behavior of the upper bound in a subspace.
//' @usage findCol(j1, j2, R, c, B)
//' @param j1 integer value, not smaller than the index of the last column of \code{R} having no negative elements.
//' @param j2 integer value, not smaller than the index of the last column of \code{R} having at least one positive element.
//' @param R matrix for the upper bound.
//' @param c number of columns in \code{R}.
//' @param B number of rows in \code{R}.
//' @author Anna Vesely.
//' @return It updates the values of \code{j1} and \code{j2}, as the indices of the last columns of \code{R}
//' having no negative elements (\code{j1}) and containing at least one positive element (\code{j2}).
//' If such a column does not exist, the corresponding index is set to 0.
//' As a result, the upper bound is increasing up to size \code{j1}, and decreasing after size \code{j2}.
//' @examples
//' R <- matrix(c(rep(0,5), rep(1,4), rep(-1,1), rep(1,3), rep(-1,2), rep(1,2), rep(-1,3)), ncol=5, byrow=TRUE)
//' findCol(5, 5, R, 5, 4) # j1=2, j2=4
//' @noRd

void findCol (int &j1, int &j2, const NumericMatrix &R, const int &c, const int &B){
  int j = std::min(c, j1);
  
  // we go back to the left until we find a column which is entirely non-negative
  // if there is no such column, j1=0
  bool cond = FALSE;
  while (!cond && j > 0){
    cond = allPos(R(_,j-1), B);
    if(!cond){--j;}
  }
  j1 = j;
  
  // we go back to the left until we find a column which is not entirely non-positive
  // if there is no such column, j1=0
  j = std::min(c, j2);
  cond = TRUE;
  while (cond && j > 0){
    cond = allPos(-R(_,j-1), B);
    if(cond){--j;}
  }
  j2 = j;
}





//' @title Matrices for a Bound
//' @description Internal function, called in \code{buildMatrices}.
//' It updates the matrices for the lower/upper bound in a subspace.
//' @usage cumulativeMatrices(Asum, A, fixed, j1, j2, z, f, B, getJ)
//' @param Asum matrix of cumulative sums for the lower/upper bound.
//' @param A matrix for the lower/upper bound.
//' @param fixed numeric vector, sum by row of fixed columns.
//' @param j1 integer value.
//' @param j2 integer value.
//' @param z number of selected columns, equal to \code{s-TD+1},
//' where \code{TD} is a candidate value for the number of true discoveries.
//' @param f total number of non-fixed variables.
//' @param B number of transformations.
//' @param getJ logical, \code{TRUE} to update \code{j1} and \code{j2}.
//' @author Anna Vesely.
//' @return It updates the matrices \code{Asum} and \code{A}
//' in the original space or a subspace generated by the branch and bound procedure.
//' If \code{getJ} is \code{TRUE}, it updates the indices \code{j1} and \code{j2}
//' by applying \code{findCol} to \code{A}.
//' @noRd

void cumulativeMatrices(NumericMatrix &Asum, NumericMatrix &A, const NumericVector &fixed,
                        int &j1, int &j2, const int &z, const int &f, const int &B, const bool getJ){
  
  NumericVector firstCol = clone(fixed);
  if(z == 1){firstCol += A(_,0);}
  if(z > 1){firstCol += rowSums(A(_,Range(0, z-1)));}
  
  NumericMatrix temp (B, f-z+1);
  Asum = temp;
  
  if(z == f){
    Asum(_,0) = firstCol;
    if(getJ){j1 = 0; j2 = 0;}
  }
  else{
    NumericMatrix rest (B, f-z);
    if(z == f-1){rest(_,0) = A(_,z);}else{rest = A(_,Range(z, f-1));}
    Asum = cbind(firstCol, rest);
    if(getJ){findCol(j1, j2, rest, f-z, B);}
    for(int j = 1; j <= f-z; ++j){Asum(_,j) = Asum(_,j) + Asum(_,j-1);} // Asum = cumulative sums of its own columns
  }
  
}





//' @title Sorting and Cumulative Sums
//' @description Internal function, called in \code{buildMatrices} and \code{checkTD}.
//' It sorts the matrices that will be used to apply the shortcut in a subspace.
//' @usage sortSum(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, B, getDsum)
//' @param Dsum matrix of cumulative sums for the lower bound.
//' @param Rsum matrix of cumulative sums for the upper bound.
//' @param D matrix for the lower bound.
//' @param I matrix of indices corresponding to elements in \code{R}.
//' @param R matrix for the upper bound.
//' @param fixed numeric vector, sum by row of fixed columns.
//' @param j1 integer value.
//' @param j2 integer value.
//' @param z number of selected columns, equal to \code{s-TD+1},
//' where \code{TD} is a candidate value for the number of true discoveries.
//' @param s size of the subset of interest.
//' @param f total number of non-fixed variables.
//' @param B number of transformations.
//' @param getDsum logical, \code{TRUE} to update \code{Dsum}.
//' @author Anna Vesely.
//' @return It returns a list with the matrices \code{Dtemp}, \code{Itemp} and \code{Rtemp},
//' obtained from \code{D}, \code{I} and \code{R} by moving at the beginning
//' the first \code{z} columns from the subset of interest. It updates \code{Rsum} and,
//' \code{getDsum} is \code{TRUE}, \code{Dsum}.
//' It updates the indices \code{j1} and \code{j2} by applying \code{findCol} to \code{R}.
//' @noRd

List sortSum(NumericMatrix &Dsum, NumericMatrix &Rsum, const NumericMatrix &D, const IntegerMatrix &I,
                   const NumericMatrix &R, const NumericVector &fixed, int &j1, int &j2, const int &z,
                   const int &s, const int &f, const int &B,
                   const bool getDsum){
  
  LogicalVector sel (f);
  IntegerVector ind = Range(0, f-1); // indices
  IntegerVector indS;
  IntegerVector indC;
  NumericMatrix Dtemp (B, f);
  
  // Dtemp = D where the first z elements from S are at the beginning
  if(getDsum){
    sel = firstInS (z, s, I(0,_), f);
    indS = ind[sel];
    indC = ind[!sel];
    for(int i = 0; i <= z-1; ++i){Dtemp(_,i) = D(_,indS[i]);}
    for(int i = 0; i <= f-z-1; ++i){Dtemp(_,z+i) = D(_,indC[i]);}
    cumulativeMatrices(Dsum, Dtemp, fixed, j1, j2, z, f, B, FALSE);
  }
  
  // update R, Rsum and I
  IntegerMatrix Itemp = clone(I);
  NumericMatrix Rtemp = clone(R);
  
  // move the first z indices in S at the beginning
  for(int b = 0; b <= B-1; ++b){
    sel = firstInS (z, s, I(b,_), f);
    indS = ind[sel];
    indC = ind[!sel];
    for(int i = 0; i <= z-1; ++i){
      Itemp(b,i) = I(b,indS[i]);
      Rtemp(b,i) = R(b,indS[i]);
    }
    for(int i = 0; i <= f-z-1; ++i){
      Itemp(b,z+i) = I(b,indC[i]);
      Rtemp(b,z+i) = R(b,indC[i]);
    }
  }
  
  cumulativeMatrices(Rsum, Rtemp, fixed, j1, j2, z, f, B, TRUE);
  
  List out = List::create(Named("Dtemp") = Dtemp, Named("Itemp") = Itemp, Named("Rtemp") = Rtemp);
  return out;
}




//' @title Building Matrices
//' @description Internal function, called in \code{goLeft} and \code{checkTD}.
//' It updates the matrices that will be used to apply the shortcut in a subspace.
//' @usage buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f0, B, getDsum)
//' @param Dsum matrix of cumulative sums for the lower bound.
//' @param Rsum matrix of cumulative sums for the upper bound.
//' @param D matrix for the lower bound.
//' @param I matrix of indices corresponding to elements in \code{R}.
//' @param R matrix for the upper bound.
//' @param fixed numeric vector, sum by row of fixed columns.
//' @param j1 integer value.
//' @param j2 integer value.
//' @param z number of selected columns, equal to \code{s-TD+1},
//' where \code{TD} is a candidate value for the number of true discoveries.
//' @param s size of the subset of interest.
//' @param f total number of non-fixed variables in the new subspace.
//' @param f0 total number of non-fixed variables in the current space.
//' @param B number of transformations.
//' @param getDsum logical, \code{TRUE} to update \code{Dsum} and \code{D}.
//' @author Anna Vesely.
//' @return It updates the matrices \code{Rsum}, \code{I} and \code{R}
//' in the original space or a subspace generated by the branch and bound procedure.
//' It updates the indices \code{j1} and \code{j2} by applying \code{findCol} to \code{R}.
//' If \code{getDsum} is \code{TRUE}, it updates the matrices \code{Dsum} and \code{D}.
//' @noRd


void buildMatrices(NumericMatrix &Dsum, NumericMatrix &Rsum, NumericMatrix &D, IntegerMatrix &I,
                   NumericMatrix &R, const NumericVector &fixed, int &j1, int &j2, const int &z,
                   const int &s, const int &f, const int &f0, const int &B,
                   const bool getDsum){
  
  NumericMatrix temp (B, 1);
  List res;
  
  if(f == 0){
    temp(_,0) = clone(fixed);
    if(getDsum){Dsum = clone(temp);}
    Rsum = clone(temp);
    j1 = 0;
    j2 = 0;
  }
  
  // update D and Dsum, if required
  if(getDsum){
    if(f == 1){
      temp(_,0) = D(_,0);
      D = as<NumericMatrix>(clone(temp));
    }
    else{D = D(_,Range(0, f-1));}
  }
  
  // update R, Rsum and I
  if(f < f0){
    IntegerMatrix Itemp (B, f);
    NumericMatrix Rtemp(B, f);
    
    IntegerVector i = I(0,_); // row in I0
    IntegerVector itemp (f); // row in Itemp
    NumericVector r (f0); // row in R
    NumericVector rtemp (f); // row in Rtemp
    
    IntegerVector indices (f0 - f);
    if(f0-f == 1){indices = i[f];}else{indices = i[Range(f, f0-1)];}// last f0-f indices in I
    LogicalVector sel (f0);
    
    // remove last f0-f indices in I
    for(int b = 0; b <= B-1; ++b){
      i = I(b,_);
      r = R(b,_);
      sel = selIndices(indices, i, f0-f, f0);
      itemp = i[!sel];
      rtemp = r[!sel];
      Itemp(b,_) = itemp;
      Rtemp(b,_) = rtemp;
    }
    I = as<IntegerMatrix>(clone(Itemp));
    R = as<NumericMatrix>(clone(Rtemp));
  }
  
  res = sortSum(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, B, getDsum);
}





//' @title Left Node by Removing Variable
//' @description Internal function, called in \code{checkTD}.
//' It generates a left node, i.e. a subspace obtained by removing a variable in the branch and bound procedure.
//' Moreover, it saves a list of elements characterizing the corresponding right node, obtained by fixing the variable.
//' @usage goLeft(vMin, vMax, Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B)
//' @param vMin minimum size to be checked.
//' @param vMax maximum size to be checked.
//' @param Dsum matrix of cumulative sums for the lower bound.
//' @param Rsum matrix of cumulative sums for the upper bound.
//' @param D matrix for the lower bound.
//' @param I matrix of indices corresponding to elements in \code{R}.
//' @param R matrix for the upper bound.
//' @param fixed numeric vector, sum by row of fixed columns.
//' @param fixedLast numeric vector, column fixed in the last step.
//' @param j1 integer value.
//' @param j2 integer value.
//' @param z number of selected columns, equal to \code{s-TD+1}, where \code{TD} is a candidate value for the number of true discoveries.
//' @param s size of the subset of interest.
//' @param f total number of non-fixed variables in the current space.
//' @param B number of transformations.
//' @author Anna Vesely.
//' @return It updates the matrices \code{Rsum}, \code{I} and \code{R},
//' and the indices \code{j1} and \code{j2} in the new left subspace.
//' It returns a list containing the values of \code{z}, \code{vMin}, \code{vMax}, \code{fixed} and \code{f}
//' in the corresponding right subspace. 
//' @noRd

List goLeft (const int &vMin, const int &vMax, NumericMatrix &Dsum, NumericMatrix &Rsum,
             NumericMatrix &D, IntegerMatrix &I, NumericMatrix &R, const NumericVector &fixed,
             NumericVector &fixedLast, int &j1, int &j2, const int &z, const int &s,
             int &f, const int &B){
  
  bool lastFromS = (I(0,f-1) <= s); // true if the index removed is in S
  --f;
  buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f+1, B, FALSE); // getDsum = FALSE
  
  fixedLast = D(_,f);
  int zNew = z;
  int vMinNew = vMin;
  int vMaxNew = vMax;
  if(lastFromS){zNew = std::max(zNew - 1, 0);}
  else{
    vMinNew = std::max(vMin - 1, 0);
    vMaxNew = std::max(vMax - 1, 0);
  }
  
  List out = List::create(Named("z") = zNew, Named("vMin") = vMinNew,
                          Named("vMax") = vMaxNew, Named("fixed") = fixed + fixedLast,
                          Named("f") = f);
  return out;
}





//' @title Sign of Quantile
//' @description Internal function, called in \code{computeBounds}.
//' It determines whether the quantile of a vector is negative.
//' @usage Q(X, k, B)
//' @param X numeric vector.
//' @param k integer value.
//' @param B length of \code{X}.
//' @author Anna Vesely.
//' @return It returns \code{TRUE} if the \code{k}-th smallest element of \code{X} is negative, and \code{FALSE} otherwise.
//' @noRd

bool Q(const NumericVector &X, const int &k, const int &B){
  bool out = TRUE; // TRUE if q < 0
  int t = B-k+1; // threshold, minimum number of non-negative elements needed to have q >= 0
  int n = 0; // number of non-negative elements
  
  for (int i = 0; i <= B-1; ++i) {
    if (X[i] >= 0){
      ++n;
      // if n reaches the threshold, then q >= 0
      if(n >= t){out = FALSE; break;}
      // if the negative elements found are already k or more, then q < 0
      else if(i + 1 - n >= k) {break;}
    }
  }
  return out;
}





//' @title Shortcut
//' @description Internal function, called in \code{checkTD}.
//' It applies the shortcut, checking the bounds for a given collection of sizes.
//' @usage computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, both)
//' @param vMin minimum size to be checked.
//' @param vMax maximum size to be checked.
//' @param rej logical, \code{TRUE} if no non-rejection has been found.
//' @param indSizes logical, \code{TRUE} if the shortcut is unsure for some sizes.
//' @param Dsum matrix of cumulative sums for the lower bound.
//' @param Rsum matrix of cumulative sums for the upper bound.
//' @param k integer value that determines the quantiles.
//' @param B number of transformations.
//' @param j1 integer value.
//' @param j2 integer value.
//' @param both logical, \code{TRUE} to check both bounds, \code{FALSE} to check only the upper bound.
//' @author Anna Vesely.
//' @return It updates \code{rej} and \code{indSizes}, and eventually \code{vMin} and \code{vMax}.
//' A non-rejection corresponds to \code{rej=FALSE} and \code{indSizes=FALSE},
//' a rejection to \code{rej=TRUE} and \code{indSizes=FALSE}, and
//' an unsure outcome to \code{rej=TRUE} and \code{indSizes=TRUE}.
//' @noRd

void computeBounds(int &vMin, int &vMax, bool &rej, bool &indSizes,
                   NumericMatrix &Dsum, NumericMatrix &Rsum,
                   const int &k, const int &B, const int &j1, const int &j2, const bool both){
  
  int H = vMax - vMin + 1;
  IntegerVector ind (H);
  if(H == 1){ind(0) = vMin;}
  else{ind = Range(vMin, vMax);}
  bool low = TRUE;
  LogicalVector up (H, FALSE); // keeps track of unsure sizes
  
  // bounds for v = ind[h0], where h0 is the first index h such that ind[h] >= j1 (h0 = H-1 if it does not exist)
  int h0 = 0;
  if(j1 >= vMin){h0 = std::min(j1 - vMin, H-1);}
  int h = h0;
  if(both){low = Q(Dsum(_,ind[h]), k, B);}
  up[h] = !Q(Rsum(_,ind[h]), k, B);
  
  // bounds for v in [j1, j2] (the behavior of the upper bound is not known)
  while (low && ind[h] < j2 && h < H-1){
    ++h;
    if(both){low = Q(Dsum(_,ind[h]), k, B);}
    up[h] = !Q(Rsum(_,ind[h]), k, B);
  }
  
  // bounds for v > j2 (the upper bound is decreasing)
  while(low && up[h] && h < H-1){
    ++h;
    if(both){low = Q(Dsum(_,ind[h]), k, B);}
    up[h] = !Q(Rsum(_,ind[h]), k, B);
  }
  
  // bounds for v < j1 (the upper bound is increasing)
  h = h0;
  while(low && up[h] && h > 0){
    --h;
    if(both){low = Q(Dsum(_,ind[h]), k, B);}
    up[h] = !Q(Rsum(_,ind[h]), k, B);
  }
  
  rej = low;
  ind = ind[up];
  H = ind.length();
  if(!rej || H == 0){indSizes = FALSE;}
  else{
    vMin = ind[0];
    vMax = ind[H-1];
  }
}








//' @title Check Candidate TD Value
//' @description Internal function, called in \code{bisectionTD}.
//' It checks whether a candidate value is a valid lower confidence bound for the number of true discoveries
//' within a subset of interest.
//' @usage checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, nMax)
//' @param TD candidate value for the number of true discoveries.
//' @param D0 matrix for the lower bound.
//' @param I0 matrix of indices corresponding to elements in \code{R0}.
//' @param R0 matrix for the upper bound.
//' @param s size of the subset of interest.
//' @param f0 total number of variables.
//' @param k integer value that determines the quantiles.
//' @param B number of transformations.
//' @param BAB number of iterations.
//' @param nMax maximum number of iterations.
//' @author Anna Vesely.
//' @return It returns a list containing \code{rej} (\code{TRUE} if no non-rejection was found),
//' \code{indSizes} (\code{TRUE} if the outcome is unsure) and \code{BAB} (number of iterations).
//' The candidate value \code{TD} is a valid confidence lower bound when \code{rej=TRUE} and \code{indSizes=FALSE},
//' it is not when \code{rej=FALSE} and \code{indSizes=FALSE}, and
//' the outcome is unsure when \code{rej=TRUE} and \code{indSizes=TRUE}.
//' The number of iterations \code{BAB} is updated.
//' @noRd

List checkTD(const int &TD, const NumericMatrix &D0, const IntegerMatrix &I0, const NumericMatrix &R0,
             const int &s, const int &f0, const int &k, const int &B,
             int &BAB, const int &nMax){
  
  ++BAB;
  int z = s - TD + 1;
  NumericMatrix Dsum (B, f0 - z + 1);
  NumericMatrix Rsum (B, f0 - z + 1);
  int j1 = f0 - z;
  int j2 = f0 - z;
  NumericVector fixed (B);
  
  List mat = sortSum(Dsum, Rsum, D0, I0, R0, fixed, j1, j2, z, s, f0, B, TRUE); // getDsum = TRUE
  NumericMatrix D = as<NumericMatrix>(mat["Dtemp"]);
  IntegerMatrix I = as<IntegerMatrix>(mat["Itemp"]);
  NumericMatrix R = as<NumericMatrix>(mat["Rtemp"]);
  
  NumericMatrix D0sort = clone(D);
  IntegerMatrix I0sort = clone(I0);
  I0sort(0,_) = I(0,_);

  int vMin = 0;
  int vMax = f0 - z + 1;
  bool rej = TRUE;
  bool indSizes = TRUE;
  computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, TRUE);
  List out = List::create(Named("rej") = rej, Named("indSizes") = indSizes, Named("BAB") = BAB);
  if(!indSizes || BAB >= nMax){return out;} // if the outcome is sure or iterations reach the maximum
  
  // BRANCH AND BOUND  ------------------------------------------------------------------------------------
  
  NumericVector fixedLast (B);
  int f = f0;
  
  List rightNodes (nMax - BAB); // list of right nodes (z, vMin, vMax, fixed, f)
  int n = 0; // number of right nodes that need to be checked
  List node (5); // right node to be added
  
  while(rej && indSizes && BAB < nMax){
    
    // while the outcome is unsure, remove variables and explore the upper bound
    while(indSizes && BAB < nMax){
      ++BAB;
      node = goLeft(vMin, vMax, Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B);
      computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, FALSE);
      rightNodes[n] = clone(node); // add node to the list
      ++n;
    }
    
    // if the first right node is closed, explore right nodes from the list
    // until either an unsure outcome is found or the list becomes empty
    while(rej && !indSizes && n > 0 && BAB < nMax){
      ++BAB;
      node = as<List>(rightNodes[n-1]);
      z = node["z"];
      vMin = node["vMin"];
      vMax = node["vMax"];
      fixed = node["fixed"];
      f = node["f"];
      D = clone(D0sort);
      I = clone(I0sort);
      R = clone(R0);
      j1 = f0 - z;
      j2 = f0 - z;
      buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f0, B, TRUE); // getDsum = TRUE
      computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, TRUE);
      --n;
      // note: not necessary to remove rightNodes[n], they can we overwritten if necessary
    }
  }
  
  out["rej"] = rej; out["indSizes"] = (rej && (indSizes || n>0)); out["BAB"] = BAB;
  return out;
}





//' @title Binary Search for the Number of True Discoveries
//' @description Internal function, called in \code{sum.internal}.
//' It employs a binary search to determine a lower confidence bound for the number of true discoveries
//' within a subset of interest.
//' @usage bisectionTD(D0, I0, R0, s, f0, k, B, nMax)
//' @param D0 matrix for the lower bound.
//' @param I0 matrix of indices corresponding to elements in \code{R0}.
//' @param R0 matrix for the upper bound.
//' @param s size of the subset of interest.
//' @param f0 total number of variables.
//' @param k integer value that determines the quantiles.
//' @param B number of transformations.
//' @param nMax maximum number of iterations.
//' @author Anna Vesely.
//' @return It returns a list containing \code{TDmin} (a valid lower confidence bound for the number of true discoveries),
//' \code{TDmax} (the maximum bound that could be found under convergence of the algorithm),
//' and \code{BAB} (number of iterations).
//' @noRd
// [[Rcpp::export]]

List bisectionTD(const NumericMatrix &D0, const IntegerMatrix &I0, const NumericMatrix &R0,
                 const int &s, const int &f0, const int &k, const int &B, const int &nMax){
  
  // NOTE: in the code, TDmax = smallest size with sure non-reject
  // but in the outcome, TDmax = smallest size with sure non-reject - 1 = maximum possible TD value
  List out = List::create(Named("TDmin") = 0, Named("TDmax") = s, Named("BAB") = 0); // [0,s]
  int TDmin = 0;
  int TDmax = s + 1; 
  int BAB = 0;
  
  // check TD = 1
  List res = checkTD(1, D0, I0, R0, s, f0, k, B, BAB, BAB+1);
  if(!res["rej"]){out["TDmax"] = 0; out["BAB"] = BAB; return out;} // sure non-reject: [0,0]
  else if(!res["indSizes"]){TDmin = 1;} // sure reject: [1,s]
  
  // check TD = s
  res = checkTD(s, D0, I0, R0, s, f0, k, B, BAB, BAB+1);
  if(!res["rej"]){TDmax = s;} // sure non-reject: [0,s-1]
  else if(!res["indSizes"]){out["TDmin"] = s; out["BAB"] = BAB; return out;} // sure reject: [s,s]
  
  int TD;
  int iMin = 0; // minimum size that leads to unsure
  int iMax = 0; // maximum size that leads to unsure
  bool cond = TRUE; // FALSE if we need to explore outside [iMin,iMax]
  
  while(TDmax - TDmin > 1 && BAB < nMax){
    TD = (TDmax + TDmin)/2;
    res = checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, BAB+1);
    if(!res["rej"]){TDmax = TD;} // sure non-reject: [TDmin,TD-1]
    else if(!res["indSizes"]){TDmin = TD;} // sure reject: [TD,TDmax-1]
    else{ // unsure outcome
      iMin = TD;
      iMax = TD;
      cond = TRUE;
      
      // TDmin = maximum value that is surely rejected without BAB
      while(iMin - TDmin > 1 && cond && BAB < nMax){
        TD = (TDmin + iMin)/2;
        res = checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, BAB+1);
        if(!res["rej"]){TDmax = TD; cond = FALSE;} // will explore [TDmin,TD-1]
        else if(!res["indSizes"]){TDmin = TD;}
        else{iMin = TD;}
      }
      
      // TDmax = minimum value that is surely not rejected without BAB
      while(TDmax - iMax > 1 && cond && BAB < nMax){
        TD = (iMax + TDmax)/2;
        res = checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, BAB+1);
        if(!res["rej"]){TDmax = TD;}
        else if(!res["indSizes"]){TDmin = TD; cond = FALSE;} // will explore [TD,TDmax-1]
        else{iMax = TD;}
      }
      
      if(cond){
        while(TDmax - TDmin > 1 && BAB < nMax){
          TD = (TDmax + TDmin)/2;
          res = checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, nMax);
          if(!res["rej"]){TDmax = TD;}
          else if(!res["indSizes"]){TDmin = TD;}
          else{break;}
        }
      }
    }
  }
  out["TDmin"] = TDmin; out["TDmax"] = TDmax - 1; out["BAB"] = BAB;
  return out;
}












