#' Display the summary statistics of covariates before and after matching.
#'
#' This function takes three arguments:
#' @param MatchBalance_output A list containing data about the results of the MatchBalance() function before and after matching.
#' @param variable_names A character vector of variable names being analyzed.
#' @param after_matching A logical value indicating whether to use data from before or after matching.
#'
#' The function returns a matrix with balance statistics for the specified variables.
#' @return A matrix containing the mean treatment, mean control, standard deviation difference, standard deviation difference pooled,
#' variance ratio, T p-value, KS p-value, mean difference, median difference, and maximum difference for each variable.
#'
#' @author Onur Düzgün


cov.balcheck <- function(MatchBalance_output, variable_names, after_matching = TRUE) {
  
  # Adding some checks
  if (!is.list(MatchBalance_output) ||
      !all(c("BeforeMatching", "AfterMatching") %in% names(MatchBalance_output))) {
    stop("Error: MatchBalance_output must contain BeforeMatching and AfterMatching.")
  }
  if (!is.character(variable_names)) {
    stop("Error: variable_names must be a character vector.")
  }
  if (!is.logical(after_matching) || length(after_matching) != 1 || is.na(after_matching)) {
    stop("Error: after_matching must be a single logical value.")
  }
  
  # Selecting the balance results needed in the output matrix
  balance_output <- if (after_matching) {
    MatchBalance_output$AfterMatching
  } else {
    MatchBalance_output$BeforeMatching
  }
  
  if (is.null(balance_output) || !is.list(balance_output)) {
    stop("Error: The requested MatchBalance() output is not available.")
  }
  if (length(variable_names) != length(balance_output)) {
    stop("Error: variable_names must have the same number of elements as the selected MatchBalance() output.")
  }
  
  # Checking that the necessary elements exist in each balanceUV object
  required_elements <- c("mean.Tr", "mean.Co", "sdiff", "sdiff.pooled",
                         "var.ratio", "p.value", "qqsummary")
  
  if (!all(vapply(balance_output, function(x) {
    is.list(x) &&
      all(required_elements %in% names(x)) &&
      is.list(x$qqsummary) &&
      all(c("meandiff", "mediandiff", "maxdiff") %in% names(x$qqsummary))
  }, logical(1)))) {
    stop("Error: The necessary elements do not exist in the selected MatchBalance() output.")
  }
  
  # Creating the output matrix
  output_matrix <- t(vapply(balance_output, function(x) {
    ks_pvalue <- if (!is.null(x$ks) && !is.null(x$ks$ks.boot.pvalue)) {
      x$ks$ks.boot.pvalue
    } else {
      NA_real_
    }
    
    c(mean_treatment = x$mean.Tr,
      mean_control = x$mean.Co,
      sdiff = x$sdiff,
      sdiff_pooled = x$sdiff.pooled,
      var_ratio = x$var.ratio,
      T_pval = x$p.value,
      KS_pval = ks_pvalue,
      qqmeandiff = x$qqsummary$meandiff,
      qqmediandiff = x$qqsummary$mediandiff,
      qqmaxdiff = x$qqsummary$maxdiff)
  }, numeric(10)))
  
  # Setting the row names
  rownames(output_matrix) <- variable_names
  
  return(output_matrix)
}
