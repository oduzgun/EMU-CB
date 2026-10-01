#' The function allows the user to choose from three different options (crossproduct, square, or combined)
#' to specify the type of covariance matrix to be returned. Default is "crossproduct".
#'
#' The function now includes input validation to ensure that the covariate matrix is valid and suitable for analysis.
#' It is also able to remove unnecessary columns with only zero values that may be generated through two-way interactions.
#'
#' @author Onur Düzgün

cov.specifier <- function(cov_matrix, method = "crossproduct") {
    
  # Input validation
  if (!is.matrix(cov_matrix) || nrow(cov_matrix) == 0 || ncol(cov_matrix) == 0) {
    stop("Input must be a non-empty matrix")
  }
  if (!is.numeric(cov_matrix)) {
    stop("Covariate matrix must contain only numeric values")
  }
  if (is.null(colnames(cov_matrix)) || any(colnames(cov_matrix) == "") ||
      length(unique(colnames(cov_matrix))) != ncol(cov_matrix)) {
    stop("Covariate matrix must have unique, non-empty column names")
  }
  if (nrow(cov_matrix) < 2) {
    stop("Covariate matrix must have at least two rows")
  }
  if (any(is.na(cov_matrix))) {
    stop("Covariate matrix must not contain missing values")
  }
  if (any(!is.finite(cov_matrix))) {
    stop("Covariate matrix must contain only finite values")
  }
  if (any(apply(cov_matrix, 2, var) == 0)) {
    stop("Covariate matrix must not have zero variance columns")
  }
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("Invalid option")
  }
    
  # Option 1: Returns original matrix and cross product matrix
  if (method == "crossproduct") {
    # Calculating the product matrix of each combination of columns
    k <- ncol(cov_matrix)
    n <- nrow(cov_matrix)
    n_products <- k * (k - 1) / 2
    if (as.double(n) * n_products > .Machine$integer.max) {
      stop("Output matrix is too large to be computed")
    }
    product_matrix <- matrix(NA_real_, n, n_products)
    
    count <- 1
    if (k > 1) {
      for (i in 1:(k - 1)) {
        index <- count:(count + k - i - 1)
        product_matrix[, index] <- cov_matrix[, i] *
          cov_matrix[, (i + 1):k, drop = FALSE]
        colnames(product_matrix)[index] <- paste0(
          colnames(cov_matrix)[i], ".", colnames(cov_matrix)[(i + 1):k]
        )
        count <- count + k - i
      }
    }
     
    # Combining the original matrix with the new matrix
    out <- cbind(cov_matrix, product_matrix)

    # Removing columns with all zero values
    zero_cols <- apply(out, 2, function(x) all(x == 0))
    out <- out[, !zero_cols, drop = FALSE]
                       
    return(out)
    
  # Option 2: Returns original matrix and squares of each column
  } else if (method == "square") {
    # Adding squared terms to output matrix
    cov_matrix.sq <- cov_matrix^2
    cov_matrix.sq.names <- paste0(colnames(cov_matrix), "^2")
    colnames(cov_matrix.sq) <- cov_matrix.sq.names
   
    # Combining the original matrix with the new matrix
    out <- cbind(cov_matrix, cov_matrix.sq)
                           
    return(out)

  # Option 3: Calculate combined matrix
  } else if (method == "combined") {
    # Calculating the product matrix of each combination of columns
    k <- ncol(cov_matrix)
    n <- nrow(cov_matrix)
    n_products <- k * (k - 1) / 2
    if (as.double(n) * n_products > .Machine$integer.max) {
      stop("Output matrix is too large to be computed")
    }
    product_matrix <- matrix(NA_real_, n, n_products)
    
    count <- 1
    if (k > 1) {
      for (i in 1:(k - 1)) {
        index <- count:(count + k - i - 1)
        product_matrix[, index] <- cov_matrix[, i] *
          cov_matrix[, (i + 1):k, drop = FALSE]
        colnames(product_matrix)[index] <- paste0(
          colnames(cov_matrix)[i], ".", colnames(cov_matrix)[(i + 1):k]
        )
        count <- count + k - i
      }
    }
      
    # Adding squared terms to output matrix
    cov_matrix.sq <- cov_matrix^2
    cov_matrix.sq.names <- paste0(colnames(cov_matrix), "^2")
    colnames(cov_matrix.sq) <- cov_matrix.sq.names
    
    # Combining the original matrix with the new matrices
    out <- cbind(cov_matrix, product_matrix, cov_matrix.sq)

    # Removing columns with all zero values
    zero_cols <- apply(out, 2, function(x) all(x == 0))
    out <- out[, !zero_cols, drop = FALSE]
    
    return(out)
    
  # Lastly, handling invalid options
  } else {
    stop("Invalid option")
  }
}
