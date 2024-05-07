#' Translate P-values into significance stars
#'
#' @param p_val A numeric value
#'
#' @return A string of asterisks representing the significance level of the input values
#' @export
#'
#' @examples sig_stars <- add_sig_stars(c(0.01, 0.0001, 5)) # returns c("*", "***", "")
add_sig_stars <- function(p_val) {
  sig_stars <- case_when(p_val < 0.001 ~ "***",
                         p_val < 0.01 ~ "**",
                         p_val <0.05 ~ "*",
                         TRUE ~ "")
  return(sig_stars)
}


#' Calculate correlation statistics between sets of features
#'
#' @param data a data frame with columns as features and rows as observations
#' @param target_features a vector of column names (to be correlated with the comparison features)
#' @param comparison_features a vector of column names (to be correlated with the target features)
#'
#' @return A data frame of target-comparison feature names,
#'  with the correlation coefficient, p-value and FDR adjusted p-value
#' @export
#'
#' @examples
run_serial_correlation_tests <- function(data, target_features, comparison_features){
  cor_df <- data.frame()
  for (target_feature in target_features){
    # Use lapply and cor.test to calculate correlations and p-values
    results <- lapply(comparison_features, function(comparison_feature)
      cor.test(data[[target_feature]], data[[comparison_feature]]))
    
    # Extract correlation coefficients and p-values from the results
    cor_coef <- sapply(results, function(x) x$estimate)
    p_vals <- sapply(results, function(x) x$p.value)
    cor_df <- data.frame(target_feature = target_feature,
                         comparison_feature = comparison_features,
                         correlation = cor_coef,
                         p_value = p_vals)
    
  }
  # Adjust the significace to multiple testing
  cor_df <- cor_df %>% 
    mutate(p_value_adj = p.adjust(p_value, method = "fdr"),
           p_value_adj_stars = add_sig_stars(p_value_adj))
  return(cor_df)
}