#' Calculate the balance of covariates before and after matching
#'
#' This function takes three arguments:
#' @param matchbal_output: A list containing data about the results of the entropy balancing method before and after matching.
#' @param variable_names: A character vector of the names of the variables being analyzed.
#' @param after_matching: A logical value indicating whether to use the data from before or after matching. 
#'
#' The function returns a matrix with the balance statistics for the specified variables.
#' @return A matrix containing the mean treatment, mean control, standard deviation difference, standard deviation difference pooled, 
#' variance ratio, T p-value, KS p-value, mean difference, median difference, and maximum difference for each variable.


cov.balcheck <- function(matchbal_output, variable_names, after_matching = TRUE) {
  # Use the select function to select only the columns needed in the output matrix
  output_matrix <- ifelse(after_matching == FALSE, 
                         select(matchbal_output$BeforeMatching, mean.Tr, mean.Co, sdiff, sdiff.pooled, 
                         var.ratio, p.value, ks.boot.pvalue = ks$ks.boot.pvalue, meandiff = qqsummary$meandiff, 
                         mediandiff = qqsummary$mediandiff, maxdiff = qqsummary$maxdiff),
                         select(matchbal_output$AfterMatching, mean.Tr, mean.Co, sdiff, sdiff.pooled, 
                         var.ratio, p.value, ks.boot.pvalue = ks$ks.boot.pvalue, meandiff = qqsummary$meandiff, 
                         mediandiff = qqsummary$mediandiff, maxdiff = qqsummary$maxdiff))
  # Set the column names and row names
  colnames(output_matrix) <- c("mean_treatment", "mean_control", "sdiff", "sdiff_pooled", 
                               "var_ratio", "T_pval", "KS_pval", "qqmeandiff", "qqmediandiff", "qqmaxdiff")
  rownames(output_matrix) <- variable_names
  
  return(output_matrix)
}