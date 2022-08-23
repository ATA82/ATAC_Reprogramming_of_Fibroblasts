format_fe_results <- function(fe_results){
  fe_results <- fe_results[, c("source", "direction", "term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size", "intersection", "precision", "recall")]
  colnames(fe_results) <-    c("source", "direction", "term_id", "term_name", "p_value", "term_size", "query_size", "number_of_shared_genes", "shared_genes", "precision", "recall")
  return(fe_results)
}