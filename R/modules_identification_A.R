

#' @title modules_identification
#' @description  function to get modules identification
#' @param data refer to time series of committed data
#' @return   gene modules
#' @export





modules_identification_A <- function(data){
  
  data_a <- data
  expr_matrix_a <- t(data_a[2:12,])
  
  expr_matrix <- expr_matrix_a
  
  cor_matrix <- cor(expr_matrix,method = 'pearson')
  
  
  
  library(pheatmap)
  
  pheatmap(cor_matrix, 
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean", 
           clustering_method = "complete", 
           breaks = seq(-0.8, 1, length.out = 51),
           color = colorRampPalette(c('#0000B5',"white", "#FF0000"))(50)
  )
  
  
}


