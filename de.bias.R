#' This function calculates the bias term given the following inputs:
#'
#' @param X0 This is a (k x n0) matrix of feature values for the control group.
#' @param X1 This is a (k x n1) matrix of feature values for the treatment group.
#' @param Ypre This is a (T0 x n1+n0) matrix of pre-treatment outcome values.
#' @param D This is a binary vector of length (n1+n0) indicating membership in one of the two groups.
#' @param W This is a (n1 x n0) matrix of weights that will be used to adjust the predicted values.
#'
#' @return It returns the bias of the aggregated synthetic control (see Abadie and L'Hour (2021) for more details)
#'
#' @author Onur Düzgün

de.bias <- function(X0, X1, Ypre, D, W){
  
  # Adding some checks
  if (!is.matrix(X0) || !is.matrix(X1) || !is.matrix(Ypre) || !is.vector(D) || !is.matrix(W)) {
    stop("X0, X1, Ypre, and W must be matrices and D must be a vector")
  }
  if (!is.numeric(X0) || !is.numeric(X1) || !is.numeric(Ypre) || !is.numeric(W)) {
    stop("X0, X1, Ypre, and W must contain only numeric values")
  }
  if (any(!is.finite(X0)) || any(!is.finite(X1)) ||
      any(!is.finite(Ypre)) || any(!is.finite(W))) {
    stop("Inputs must contain only finite values")
  }
  if (nrow(X0) != nrow(X1)) {
    stop("X0 and X1 must have the same number of rows")
  }
  if (length(D) != ncol(X0) + ncol(X1) || ncol(Ypre) != length(D)) {
    stop("D and Ypre must be conformable with X0 and X1")
  }
  if (any(is.na(D)) || !all(D %in% c(0, 1))) {
    stop("D must be a binary vector containing only 0 and 1")
  }
  if (sum(D == 0) != ncol(X0) || sum(D == 1) != ncol(X1)) {
    stop("The group counts in D must match the columns of X0 and X1")
  }
  if (!identical(dim(W), c(ncol(X1), ncol(X0)))) {
    stop("W must be an n1 x n0 matrix")
  }
  
  # Creating augmented feature matrices by setting the last element of each column to 1
  augmented_X0 <- rbind(X0, rep(1, ncol(X0)))
  augmented_X1 <- rbind(X1, rep(1, ncol(X1)))
  
  # Calculating the regression matrix for the untreated units
  weight_matrix <- tryCatch(
    solve(augmented_X0 %*% t(augmented_X0)) %*% augmented_X0,
    error = function(e) {
      stop("Error calculating weight matrix: augmented_X0 must have full row rank")
    }
  )
  
  # Pre-allocating dbias
  dbias <- matrix(nrow = nrow(Ypre), ncol = ncol(X1))
  
  # Calculating dbias given the weights
  for (t in 1:nrow(Ypre)) {
    mean0 <- weight_matrix %*% Ypre[t, D == 0]
    dbias[t, ] <- t(augmented_X1) %*% mean0 -
      W %*% (t(augmented_X0) %*% mean0)
  }
  
  return(dbias)
}
