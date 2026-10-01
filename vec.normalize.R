#' This function takes X0 and X1 and divides the elements by the square root of their row variances.
#' The row variances are calculated jointly across the control and treatment units.
#' The function also includes checks and error handling.
#'
#' @param X0 A (k x n0) matrix of observables of the control group.
#' @param X1 A (k x n1) matrix of observables of the treatment group.
#'
#' @return It returns a list containing the two normalized matrices.
#'
#' @author Onur Düzgün


vec.normalize <- function(X0, X1){
  # Adding some checks
  if (!is.matrix(X0) && !is.data.frame(X0)) {
    stop("X0 must be a matrix or data frame")
  }
  if (!is.matrix(X1) && !is.data.frame(X1)) {
    stop("X1 must be a matrix or data frame")
  }
  
  X0 <- as.matrix(X0)
  X1 <- as.matrix(X1)
  
  if (!is.numeric(X0) || !is.numeric(X1)) {
    stop("X0 and X1 must contain only numeric values")
  }
  if (nrow(X0) != nrow(X1)) {
    stop("X0 and X1 must have the same number of rows")
  }
  if (ncol(X0) + ncol(X1) < 2) {
    stop("At least two units are required to calculate row variances")
  }
  if (any(!is.finite(X0)) || any(!is.finite(X1))) {
    stop("X0 and X1 must contain only finite values")
  }

  # Normalizing X
  combined_X <- cbind(X0, X1)
  row_sd <- sqrt(diag(cov(t(combined_X))))
  
  # Handling cases where the variance is zero
  if (any(!is.finite(row_sd)) || any(row_sd == 0)) {
    stop("Cannot normalize data. The variance of one or more rows is zero or undefined. Note: all variables must vary across units.")
  }

  X0_norm <- X0 / row_sd
  X1_norm <- X1 / row_sd
  
  return(list(X0 = X0_norm, X1 = X1_norm))
}
