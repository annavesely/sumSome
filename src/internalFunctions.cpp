#include <Rcpp.h>
using namespace Rcpp;


//' @title Compare Vector Elements to Value
//' @description Internal function, called in \code{sumCt.internal}.
//' It determines whether all the elements of a vector, eventually except the first one, are equal to a given value.
//' @usage permMin(X, B, truncTo)
//' @param X numeric vector
//' @param B length of \code{X}
//' @param truncTo numeric value
//' @author Anna Vesely
//' @return It returns \code{TRUE} if all elements of \code{X}, excluding the first one, are equal to \code{truncTo}, and returns \code{FALSE} otherwise.
//' @keywords internal
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




//' @title Sign of Vector Elements
//' @description Internal function, called in \code{findCol}.
//' It determines whether all the elements of a vector are non-negative.
//' @usage allPos(X, B)
//' @param X numeric vector
//' @param B length of \code{X}
//' @author Anna Vesely
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





//' @title Sign of Matrix Columns
//' @description Internal function, called in \code{cumulativeMatrices}.
//' Given a matrix where columns are in decreasing order, it finds the last column having no negative elements,
//' and the last column having at least one positive element.
//' @usage findCol(j1, j2, R, c, B)
//' @param j1 integer value, not smaller than the index of the last column of \code{R} having no negative elements
//' @param j2 integer value, not smaller than the index of the last column of \code{R} having at least one positive element
//' @param R numeric matrix
//' @param c number of columns in \code{R}
//' @param B number of rows in \code{R}
//' @author Anna Vesely
//' @return It updates the values of \code{j1} and \code{j2}, as the indices of the last columns of \code{R}
//' having no negative elements (\code{j1}) and containing at least one positive element (\code{j2}).
//' If such a column does not exist, the corresponding index is set to 0.
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





//' @title Matrix of Cumulative Sums
//' @description Internal function, called in \code{buildMatrices}.
//' It updates a matrix of comulative sums.
//' @usage cumulativeMatrices(Asum, A, fixed, j1, j2, z, f, B, getJ)
//' @param Asum numeric matrix
//' @param A numeric matrix
//' @param fixed numeric vector
//' @param j1 integer value
//' @param j2 integer value
//' @param z number of selected columns
//' @param f number of columns in \code{A}
//' @param B number of rows in \code{A}
//' @param getJ boolean, \code{TRUE} to update indices
//' @author Anna Vesely
//' @return It sums by row \code{fixed} and the first \code{z} columns of \code{A}, then combines by row
//' the result and the remaining \code{f-z} columns of \code{A}.
//' \code{Asum} is updated as the corresponding matrix of cumulative sums.
//' If \code{getJ} is \code{TRUE}, then \code{j1} and \code{j2} are updated by applying \code{findCol}
//' to the remaining \code{f-z} columns of \code{A}.
//' @noRd

void cumulativeMatrices(NumericMatrix &Asum, NumericMatrix &A, const NumericVector &fixed,
                        int &j1, int &j2, const int &z, const int &f, const int &B, const bool getJ){
  
  NumericVector firstCol = clone(fixed);
  if(z > 0){firstCol += rowSums(A(_,Range(0, z-1)));}
  NumericMatrix rest = A(_,Range(z, f-1));
  Asum = cbind(firstCol, rest);
  if(getJ){findCol(j1, j2, rest, f-z, B);}
  
  // Asum = cumulative sums of its own columns
  for(int j = 1; j <= f-z; ++j){
    Asum(_,j) = Asum(_,j) + Asum(_,j-1);
  }
}





//' @title First Selected Elements
//' @description Internal function, called in \code{buildMatrices}.
//' It selects the first elements of a vector that do not exceed a given threshold.
//' @usage firstInS(z, s, X, f)
//' @param z number of selected elements
//' @param s integer threshold
//' @param X integer vector
//' @param f length of \code{X}
//' @author Anna Vesely
//' @return It returns a boolean vector of length \code{f}, where \code{TRUE} values correspond
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
//' @param indices integer vector
//' @param X integer vector
//' @param H length of \code{indices}
//' @param f length of \code{X}
//' @author Anna Vesely
//' @return It returns a boolean vector of length \code{f}, where \code{TRUE} values correspond
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





//' @title Building Matrices
//' @description Internal function, called in \code{goLeft}, \code{goRight} and \code{checkTD}.
//' It updates the matrices that will be used in \code{checkTD} to apply the shortcut in a subspace.
//' @usage buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f0, B, getDsum)
//' @param Dsum matrix of cumulative sums for the lower bound
//' @param Rsum matrix of cumulative sums for the upper bound
//' @param D matrix for the lower bound
//' @param I matrix of indices corresponding to elements in \code{R}
//' @param R matrix for the upper bound
//' @param fixed numeric vector, sum by row of fixed columns
//' @param j1 integer value
//' @param j2 integer value
//' @param z number of selected columns, equal to \code{s-TD+1}, where \code{TD} is a candidate value for the number of true discoveries
//' @param s size of the subset under closed testing
//' @param f number of non-fixed variables in the new subspace
//' @param f0 number of non-fixed variables in the current space
//' @param B number of transformations
//' @param getDsum boolean, \code{TRUE} to update \code{Dsum} and \code{D} 
//' @author Anna Vesely
//' @return It updates the matrices \code{Rsum}, \code{I} and \code{R} used to compute the upper bound
//' in the original space or a subspace generated by the branch and bound procedure.
//' It updates the indices \code{j1} and \code{j2} by applying \code{findCol} to \code{R}.
//' If \code{getDsum} is \code{TRUE}, it updates the matrices \code{Dsum} and \code{D} used to compute the lower bound.
//' @noRd

void buildMatrices(NumericMatrix &Dsum, NumericMatrix &Rsum, NumericMatrix &D, IntegerMatrix &I,
                   NumericMatrix &R, const NumericVector &fixed, int &j1, int &j2, const int &z,
                   const int &s, const int &f, const int &f0, const int &B,
                   const bool getDsum){
  
  LogicalVector sel (f0);
  IntegerVector ind = Range(0, f-1); // indices
  IntegerVector indS;
  IntegerVector indC;
  
  // update D and Dsum, if required
  if(getDsum){
    NumericMatrix Dtemp = D(_,Range(0, f-1));
    if (z > 0){
      sel = firstInS (z, s, I(0,_), f);
      indS = ind[sel];
      indC = ind[!sel];
      for(int i = 0; i <= z-1; ++i){Dtemp(_,i) = D(_,indS[i]);}
      for(int i = 0; i <= f-z-1; ++i){Dtemp(_,z+i) = D(_,indC[i]);}
    }
    D = clone(Dtemp);
    cumulativeMatrices(Dsum, D, fixed, j1, j2, z, f, B, FALSE);
  }
  
  // update R, Rsum and I
  IntegerMatrix Itemp (B, f);
  NumericMatrix Rtemp(B, f);
  
  // if necessary, reduce dimension
  if(f == f0){
    Itemp = clone(I);
    Rtemp = clone(R);
  }else{
    IntegerVector i = I(0,_); // row in I
    IntegerVector itemp (f); // row in Itemp
    NumericVector r (f0); // row in R
    NumericVector rtemp (f); // row in Rtemp
    IntegerVector indices = i[Range(f, f0-1)]; // last f0 - f indices in I
    
    for(int b = 0; b <= B-1; ++b){
      i = I(b,_);
      r = R(b,_);
      sel = selIndices(indices, i, f0-f, f0);
      itemp = i[!sel];
      rtemp = r[!sel];
      Itemp(b,_) = itemp;
      Rtemp(b,_) = rtemp;
    }
    I = I(_,Range(0,f-1));
    R = R(_,Range(0,f-1));
  }
  
  // move the first z indices in S at the beginning
  if(z > 0){
    for(int b = 0; b <= B-1; ++b){
      sel = firstInS (z, s, Itemp(b,_), f);
      indS = ind[sel];
      indC = ind[!sel];
      for(int i = 0; i <= z-1; ++i){
        I(b,i) = Itemp(b,indS[i]);
        R(b,i) = Rtemp(b,indS[i]);
      }
      for(int i = 0; i <= f-z-1; ++i){
        I(b,z+i) = Itemp(b,indC[i]);
        R(b,z+i) = Rtemp(b,indC[i]);
      }
    }
  }else{
    I <- clone(Itemp);
    R <- clone(Rtemp);
  }
  
  // compute cumulative sums and indices j1, j2
  cumulativeMatrices(Rsum, R, fixed, j1, j2, z, f, B, TRUE);
}





//' @title Creating Left Node by Removing Variable
//' @description Internal function, called in \code{checkTD}.
//' It generates a left node, i.e. a subspace obtained by removing a variable in the branch and bound procedure.
//' Moreover, it saves a list of elements characterizing the corresponding right node, obtained by fixing the variable.
//' @usage goLeft(vMin, vMax, Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B, lastFromS)
//' @param vMin minimum size to be checked
//' @param vMax maximum size to be checked
//' @param Dsum matrix of cumulative sums for the lower bound
//' @param Rsum matrix of cumulative sums for the upper bound
//' @param D matrix for the lower bound
//' @param I matrix of indices corresponding to elements in \code{R}
//' @param R matrix for the upper bound
//' @param fixed numeric vector, sum by row of fixed columns
//' @param fixedLast numeric vector, column fixed in the last step
//' @param j1 integer value
//' @param j2 integer value
//' @param z number of selected columns, equal to \code{s-TD+1}, where \code{TD} is a candidate value for the number of true discoveries
//' @param s size of the subset under closed testing
//' @param f number of non-fixed variables in the current space
//' @param B number of transformations
//' @param lastFromS boolean, \code{TRUE} if the removed variable is in the subset under closed testing
//' @author Anna Vesely
//' @return It updates \code{lastFromS}, the matrices \code{Rsum}, \code{I} and \code{R},
//' and the indices \code{j1} and \code{j2} in the new left subspace.
//' It returns a list containing the values of \code{z}, \code{vMin}, \code{vMax}, \code{fixed} and \code{f}
//' in the corresponding right subspace. 
//' @noRd

List goLeft (const int &vMin, const int &vMax, NumericMatrix &Dsum, NumericMatrix &Rsum,
             NumericMatrix &D, IntegerMatrix &I, NumericMatrix &R, const NumericVector &fixed,
             NumericVector &fixedLast, int &j1, int &j2, const int &z, const int &s,
             int &f, const int &B, bool &lastFromS){
  
  lastFromS = (I(0,f-1) <= s); // true if the index removed is in S
  --f;
  buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f+1, B, FALSE); // getDsum = FALSE
  
  fixedLast = D(_,f);
  int zNew = z;
  int vMinNew = vMin;
  int vMaxNew = vMax;
  if(lastFromS){--zNew;}
  else{
    vMinNew = std::max(vMin - 1, 0);
    vMaxNew = std::max(vMax - 1, 0);
  }
  
  List out = List::create(Named("z") = zNew, Named("vMin") = vMinNew,
                          Named("vMax") = vMaxNew, Named("fixed") = fixed + fixedLast,
                          Named("f") = f);
  return out;
}





//' @title Creating Right Node by Fixing Variable
//' @description Internal function, called in \code{checkTD}.
//' It generates a right node, i.e. a subspace obtained by fixing a variable in the branch and bound procedure,
//' after closing the corresponding left node.
//' @usage goRight(Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B, lastFromS)
//' @param Dsum matrix of cumulative sums for the lower bound
//' @param Rsum matrix of cumulative sums for the upper bound
//' @param D matrix for the lower bound
//' @param I matrix of indices corresponding to elements in \code{R}
//' @param R matrix for the upper bound
//' @param fixed numeric vector, sum by row of fixed columns
//' @param fixedLast numeric vector, column fixed in the last step
//' @param j1 integer value
//' @param j2 integer value
//' @param z number of selected columns, equal to \code{s-TD+1}, where \code{TD} is a candidate value for the number of true discoveries
//' @param s size of the subset under closed testing
//' @param f number of non-fixed variables in the current space
//' @param B number of transformations
//' @param lastFromS boolean, \code{TRUE} if the fixed variable is in the subset under closed testing
//' @author Anna Vesely
//' @return It updates the matrices \code{Dsum}, \code{Rsum}, \code{D}, \code{I} and \code{R},
//' and the indices \code{j1} and \code{j2} in the new right subspace.
//' @noRd

void goRight (NumericMatrix &Dsum, NumericMatrix &Rsum,
              NumericMatrix &D, IntegerMatrix &I, NumericMatrix &R, const NumericVector &fixed,
              const NumericVector &fixedLast, int &j1, int &j2, const int &z, const int &s,
              const int &f, const int &B, const bool &lastFromS){
  
  if(lastFromS){buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f, B, TRUE);}
  else{
    D = D(_,Range(0,f-1));
    Dsum = Dsum(_,Range(0,f-z));
    for (int b = 0; b <= B-1; ++b){ // add fixedLast to Dsum and Rsum
      for (int i = 0; i <= f; ++i){
        Dsum(b,i) += fixedLast[b];
        Rsum(b,i) += fixedLast[b];
      }
    }
  }
}





//' @title Sign of Quantile
//' @description Internal function, called in \code{computeBounds}.
//' It determines whether the quantile of a vector is negative.
//' @usage Q(X, k, B)
//' @param X numeric vector
//' @param k integer value
//' @param B length of \code{X}
//' @author Anna Vesely
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
//' @param vMin minimum size to be checked
//' @param vMax maximum size to be checked
//' @param rej boolean, \code{TRUE} if no non-rejection has been found
//' @param indSizes boolean, \code{TRUE} if the shortcut is unsure for some sizes
//' @param Dsum matrix of cumulative sums for the lower bound
//' @param Rsum matrix of cumulative sums for the upper bound
//' @param k integer value that determines the quantiles
//' @param B number of transformations
//' @param j1 integer value
//' @param j2 integer value
//' @param both boolean, \code{TRUE} to check both bounds, \code{FALSE} to check only the upper bound
//' @author Anna Vesely
//' @return It updates \code{rej} and \code{indSizes}, and eventually \code{vMin} and \code{vMax}.
//' A non-rejection corresponds to \code{rej=FALSE} and \code{indSizes=FALSE},
//' a rejection to \code{rej=TRUE} and \code{indSizes=FALSE}, and
//' an unsure outcome to \code{rej=TRUE} and \code{indSizes=TRUE}.
//' @noRd

void computeBounds(int &vMin, int &vMax, bool &rej, bool &indSizes,
                   NumericMatrix &Dsum, NumericMatrix &Rsum,
                   const int &k, const int &B, const int &j1, const int &j2, const bool both){
  
  IntegerVector ind = Range(vMin, vMax);
  int H = ind.length();
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
//' It checks whether a candidate value is a valid lower confidence bound for the number of true discoveries within a subset.
//' @usage checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, nMax)
//' @param TD candidate value for the number of true discoveries
//' @param D0 matrix for the lower bound
//' @param I0 matrix of indices corresponding to elements in \code{R0}
//' @param R0 matrix for the upper bound
//' @param s size of the subset under closed testing
//' @param f0 number of variables
//' @param k integer value that determines the quantiles
//' @param B number of transformations
//' @param BAB number of iterations
//' @param nMax maximum number of iterations
//' @author Anna Vesely
//' @return It returns a list containing \code{rej} (\code{TRUE} if no non-rejection was found),
//' \code{indSizes} (\code{TRUE} if the outcome is unsure) and \code{BAB} (number of iterations).
//' \code{TD} is a valid confidence lower bound when \code{rej=TRUE} and \code{indSizes=FALSE},
//' it is not when \code{rej=FALSE} and \code{indSizes=FALSE}, and
//' the outcome is unsure when \code{rej=TRUE} and \code{indSizes=TRUE}.
//' @keywords internal
// [[Rcpp::export]]

List checkTD(const int &TD, const NumericMatrix &D0, const IntegerMatrix &I0, const NumericMatrix &R0,
             const int &s, const int &f0, const int &k, const int &B,
             int &BAB, const int &nMax){
  
  ++BAB;
  int z = s - TD + 1;
  NumericMatrix D = clone(D0);
  IntegerMatrix I = clone(I0);
  NumericMatrix R = clone(R0);
  NumericMatrix Dsum (B, f0 - z + 1);
  NumericMatrix Rsum (B, f0 - z + 1);
  int j1 = f0 - z;
  int j2 = f0 - z;
  NumericVector fixed (B);
  
  buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f0, f0, B, TRUE); // getDsum = TRUE
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
  bool lastFromS = FALSE;
  
  List rightNodes (nMax - BAB); // list of right nodes (z, vMin, vMax, fixed, f)
  int n = 0; // number of right nodes that need to be checked
  List node (5); // right node to be added
  
  while(rej && indSizes && BAB < nMax){
    
    // while the outcome is unsure, remove variables and explore the upper bound
    while(indSizes && BAB < nMax){
      ++BAB;
      node = goLeft(vMin, vMax, Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B, lastFromS);
      computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, FALSE);
      ++n;
      if(indSizes){rightNodes[n-1] = clone(node);} // add node to the list
    }
    
    // proceed with the node on the right of the last closed node
    // everything is the same, eventually except z, vMin, vMax, fixed
    if(BAB < nMax){
      ++BAB;
      z = node["z"];
      vMin = node["vMin"];
      vMax = node["vMax"];
      fixed = node["fixed"];
      goRight(Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B, lastFromS);
      computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, TRUE);
      --n;
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
      D = clone(D0);
      I = clone(I0);
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





