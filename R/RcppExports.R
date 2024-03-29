# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

findDiscSum <- function(s, m, u, v, cs) {
    .Call(`_sumSome_findDiscSum`, s, m, u, v, cs)
}

#' @title First Selected Elements
#' @description Internal function, called in \code{sortSum}.
#' It selects the first elements of a vector that do not exceed a given threshold.
#' @usage firstInS(z, s, X, f)
#' @param z number of selected elements.
#' @param s integer threshold.
#' @param X integer vector.
#' @param f length of \code{X}.
#' @author Anna Vesely.
#' @return It returns a ligical vector of length \code{f}, where \code{TRUE} values correspond
#' to the first \code{z} elements of \code{X} that are not greater than \code{s}.
#' @noRd
NULL

#' @title Element Selection
#' @description Internal function, called in \code{buildMatrices}.
#' It selects the elements of a vector that appear in a second vector.
#' @usage selIndices (indices, X, H, f)
#' @param indices integer vector.
#' @param X integer vector.
#' @param H length of \code{indices}.
#' @param f length of \code{X}.
#' @author Anna Vesely.
#' @return It returns a logical vector of length \code{f}, where \code{TRUE} values correspond
#' to elements of \code{X} that appear in \code{indices}.
#' @noRd
NULL

#' @title Sign of Vector Elements
#' @description Internal function, called in \code{findCol}.
#' It determines whether all the elements of a vector are non-negative.
#' @usage allPos(X, B)
#' @param X numeric vector.
#' @param B length of \code{X}.
#' @author Anna Vesely.
#' @return It returns \code{TRUE} if all elements of \code{X} are non-negative, and \code{FALSE} otherwise.
#' @noRd
NULL

#' @title Upper Bound Behavior
#' @description Internal function, called in \code{cumulativeMatrices}.
#' It characterizes the behavior of the upper bound in a subspace.
#' @usage findCol(j1, j2, R, c, B)
#' @param j1 integer value, not smaller than the index of the last column of \code{R} having no negative elements.
#' @param j2 integer value, not smaller than the index of the last column of \code{R} having at least one positive element.
#' @param R matrix for the upper bound.
#' @param c number of columns in \code{R}.
#' @param B number of rows in \code{R}.
#' @author Anna Vesely.
#' @return It updates the values of \code{j1} and \code{j2}, as the indices of the last columns of \code{R}
#' having no negative elements (\code{j1}) and containing at least one positive element (\code{j2}).
#' If such a column does not exist, the corresponding index is set to 0.
#' As a result, the upper bound is increasing up to size \code{j1}, and decreasing after size \code{j2}.
#' @examples
#' R <- matrix(c(rep(0,5), rep(1,4), rep(-1,1), rep(1,3), rep(-1,2), rep(1,2), rep(-1,3)), ncol=5, byrow=TRUE)
#' findCol(5, 5, R, 5, 4) # j1=2, j2=4
#' @noRd
NULL

#' @title Matrices for a Bound
#' @description Internal function, called in \code{buildMatrices}.
#' It updates the matrices for the lower/upper bound in a subspace.
#' @usage cumulativeMatrices(Asum, A, fixed, j1, j2, z, f, B, getJ)
#' @param Asum matrix of cumulative sums for the lower/upper bound.
#' @param A matrix for the lower/upper bound.
#' @param fixed numeric vector, sum by row of fixed columns.
#' @param j1 integer value.
#' @param j2 integer value.
#' @param z number of selected columns, equal to \code{s-TD+1},
#' where \code{TD} is a candidate value for the number of true discoveries.
#' @param f total number of non-fixed variables.
#' @param B number of transformations.
#' @param getJ logical, \code{TRUE} to update \code{j1} and \code{j2}.
#' @author Anna Vesely.
#' @return It updates the matrices \code{Asum} and \code{A}
#' in the original space or a subspace generated by the branch and bound procedure.
#' If \code{getJ} is \code{TRUE}, it updates the indices \code{j1} and \code{j2}
#' by applying \code{findCol} to \code{A}.
#' @noRd
NULL

#' @title Sorting and Cumulative Sums
#' @description Internal function, called in \code{buildMatrices} and \code{checkTD}.
#' It sorts the matrices that will be used to apply the shortcut in a subspace.
#' @usage sortSum(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, B, getDsum)
#' @param Dsum matrix of cumulative sums for the lower bound.
#' @param Rsum matrix of cumulative sums for the upper bound.
#' @param D matrix for the lower bound.
#' @param I matrix of indices corresponding to elements in \code{R}.
#' @param R matrix for the upper bound.
#' @param fixed numeric vector, sum by row of fixed columns.
#' @param j1 integer value.
#' @param j2 integer value.
#' @param z number of selected columns, equal to \code{s-TD+1},
#' where \code{TD} is a candidate value for the number of true discoveries.
#' @param s size of the subset of interest.
#' @param f total number of non-fixed variables.
#' @param B number of transformations.
#' @param getDsum logical, \code{TRUE} to update \code{Dsum}.
#' @author Anna Vesely.
#' @return It returns a list with the matrices \code{Dtemp}, \code{Itemp} and \code{Rtemp},
#' obtained from \code{D}, \code{I} and \code{R} by moving at the beginning
#' the first \code{z} columns from the subset of interest. It updates \code{Rsum} and,
#' \code{getDsum} is \code{TRUE}, \code{Dsum}.
#' It updates the indices \code{j1} and \code{j2} by applying \code{findCol} to \code{R}.
#' @noRd
NULL

#' @title Building Matrices
#' @description Internal function, called in \code{goLeft} and \code{checkTD}.
#' It updates the matrices that will be used to apply the shortcut in a subspace.
#' @usage buildMatrices(Dsum, Rsum, D, I, R, fixed, j1, j2, z, s, f, f0, B, getDsum)
#' @param Dsum matrix of cumulative sums for the lower bound.
#' @param Rsum matrix of cumulative sums for the upper bound.
#' @param D matrix for the lower bound.
#' @param I matrix of indices corresponding to elements in \code{R}.
#' @param R matrix for the upper bound.
#' @param fixed numeric vector, sum by row of fixed columns.
#' @param j1 integer value.
#' @param j2 integer value.
#' @param z number of selected columns, equal to \code{s-TD+1},
#' where \code{TD} is a candidate value for the number of true discoveries.
#' @param s size of the subset of interest.
#' @param f total number of non-fixed variables in the new subspace.
#' @param f0 total number of non-fixed variables in the current space.
#' @param B number of transformations.
#' @param getDsum logical, \code{TRUE} to update \code{Dsum} and \code{D}.
#' @author Anna Vesely.
#' @return It updates the matrices \code{Rsum}, \code{I} and \code{R}
#' in the original space or a subspace generated by the branch and bound procedure.
#' It updates the indices \code{j1} and \code{j2} by applying \code{findCol} to \code{R}.
#' If \code{getDsum} is \code{TRUE}, it updates the matrices \code{Dsum} and \code{D}.
#' @noRd
NULL

#' @title Left Node by Removing Variable
#' @description Internal function, called in \code{checkTD}.
#' It generates a left node, i.e. a subspace obtained by removing a variable in the branch and bound procedure.
#' Moreover, it saves a list of elements characterizing the corresponding right node, obtained by fixing the variable.
#' @usage goLeft(vMin, vMax, Dsum, Rsum, D, I, R, fixed, fixedLast, j1, j2, z, s, f, B)
#' @param vMin minimum size to be checked.
#' @param vMax maximum size to be checked.
#' @param Dsum matrix of cumulative sums for the lower bound.
#' @param Rsum matrix of cumulative sums for the upper bound.
#' @param D matrix for the lower bound.
#' @param I matrix of indices corresponding to elements in \code{R}.
#' @param R matrix for the upper bound.
#' @param fixed numeric vector, sum by row of fixed columns.
#' @param fixedLast numeric vector, column fixed in the last step.
#' @param j1 integer value.
#' @param j2 integer value.
#' @param z number of selected columns, equal to \code{s-TD+1}, where \code{TD} is a candidate value for the number of true discoveries.
#' @param s size of the subset of interest.
#' @param f total number of non-fixed variables in the current space.
#' @param B number of transformations.
#' @author Anna Vesely.
#' @return It updates the matrices \code{Rsum}, \code{I} and \code{R},
#' and the indices \code{j1} and \code{j2} in the new left subspace.
#' It returns a list containing the values of \code{z}, \code{vMin}, \code{vMax}, \code{fixed} and \code{f}
#' in the corresponding right subspace. 
#' @noRd
NULL

#' @title Sign of Quantile
#' @description Internal function, called in \code{computeBounds}.
#' It determines whether the quantile of a vector is negative.
#' @usage Q(X, k, B)
#' @param X numeric vector.
#' @param k integer value.
#' @param B length of \code{X}.
#' @author Anna Vesely.
#' @return It returns \code{TRUE} if the \code{k}-th smallest element of \code{X} is negative, and \code{FALSE} otherwise.
#' @noRd
NULL

#' @title Shortcut
#' @description Internal function, called in \code{checkTD}.
#' It applies the shortcut, checking the bounds for a given collection of sizes.
#' @usage computeBounds(vMin, vMax, rej, indSizes, Dsum, Rsum, k, B, j1, j2, both)
#' @param vMin minimum size to be checked.
#' @param vMax maximum size to be checked.
#' @param rej logical, \code{TRUE} if no non-rejection has been found.
#' @param indSizes logical, \code{TRUE} if the shortcut is unsure for some sizes.
#' @param Dsum matrix of cumulative sums for the lower bound.
#' @param Rsum matrix of cumulative sums for the upper bound.
#' @param k integer value that determines the quantiles.
#' @param B number of transformations.
#' @param j1 integer value.
#' @param j2 integer value.
#' @param both logical, \code{TRUE} to check both bounds, \code{FALSE} to check only the upper bound.
#' @author Anna Vesely.
#' @return It updates \code{rej} and \code{indSizes}, and eventually \code{vMin} and \code{vMax}.
#' A non-rejection corresponds to \code{rej=FALSE} and \code{indSizes=FALSE},
#' a rejection to \code{rej=TRUE} and \code{indSizes=FALSE}, and
#' an unsure outcome to \code{rej=TRUE} and \code{indSizes=TRUE}.
#' @noRd
NULL

#' @title Check Candidate TD Value
#' @description Internal function, called in \code{bisectionTD}.
#' It checks whether a candidate value is a valid lower confidence bound for the number of true discoveries
#' within a subset of interest.
#' @usage checkTD(TD, D0, I0, R0, s, f0, k, B, BAB, nMax)
#' @param TD candidate value for the number of true discoveries.
#' @param D0 matrix for the lower bound.
#' @param I0 matrix of indices corresponding to elements in \code{R0}.
#' @param R0 matrix for the upper bound.
#' @param s size of the subset of interest.
#' @param f0 total number of variables.
#' @param k integer value that determines the quantiles.
#' @param B number of transformations.
#' @param BAB number of iterations.
#' @param nMax maximum number of iterations.
#' @author Anna Vesely.
#' @return It returns a list containing \code{rej} (\code{TRUE} if no non-rejection was found),
#' \code{indSizes} (\code{TRUE} if the outcome is unsure) and \code{BAB} (number of iterations).
#' The candidate value \code{TD} is a valid confidence lower bound when \code{rej=TRUE} and \code{indSizes=FALSE},
#' it is not when \code{rej=FALSE} and \code{indSizes=FALSE}, and
#' the outcome is unsure when \code{rej=TRUE} and \code{indSizes=TRUE}.
#' The number of iterations \code{BAB} is updated.
#' @noRd
NULL

#' @title Compare Vector Elements to Value
#' @description Internal function, called in \code{sum.internal}.
#' It determines whether all the elements of a vector, eventually except the first one, are equal to a given value.
#' @usage permMin(X, B, truncTo)
#' @param X numeric vector.
#' @param B length of \code{X}.
#' @param truncTo numeric value.
#' @author Anna Vesely.
#' @return It returns \code{TRUE} if all elements of \code{X}, excluding the first one, are equal to \code{truncTo},
#' and returns \code{FALSE} otherwise.
#' @noRd
permMin <- function(X, B, truncTo) {
    .Call(`_sumSome_permMin`, X, B, truncTo)
}

#' @title Binary Search for the Number of True Discoveries
#' @description Internal function, called in \code{sum.internal}.
#' It employs a binary search to determine a lower confidence bound for the number of true discoveries
#' within a subset of interest.
#' @usage bisectionTD(D0, I0, R0, s, f0, k, B, nMax)
#' @param D0 matrix for the lower bound.
#' @param I0 matrix of indices corresponding to elements in \code{R0}.
#' @param R0 matrix for the upper bound.
#' @param s size of the subset of interest.
#' @param f0 total number of variables.
#' @param k integer value that determines the quantiles.
#' @param B number of transformations.
#' @param nMax maximum number of iterations.
#' @author Anna Vesely.
#' @return It returns a list containing \code{TDmin} (a valid lower confidence bound for the number of true discoveries),
#' \code{TDmax} (the maximum bound that could be found under convergence of the algorithm),
#' and \code{BAB} (number of iterations).
#' @noRd
bisectionTD <- function(D0, I0, R0, s, f0, k, B, nMax) {
    .Call(`_sumSome_bisectionTD`, D0, I0, R0, s, f0, k, B, nMax)
}

