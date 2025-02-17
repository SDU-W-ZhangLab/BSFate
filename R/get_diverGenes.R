
#' @title get_diverGenes
#' @description  function to get diver genes
#' @param ODE_pair_results refer to bistable circuitGenes
#' @return  diver genes
#' @export






get_diverGenes <- function(ODE_pair_results){
  
  results <- ODE_pair_results
  results <- as.matrix(results)
  sorted_scores <- results
  filtered_score <- sorted_scores 
  colnames(filtered_score) <- c("RSS")
  filtered_score <- as.data.frame(filtered_score)
  
  
  
  
  expanded_data <- filtered_score %>%
    rownames_to_column("Gene") %>%  
    separate(Gene, into = c("Gene1", "Gene2"), sep = "_") %>%  
    pivot_longer(cols = c("Gene1", "Gene2"), names_to = "GeneType", values_to = "Gene") %>%
    dplyr::select(Gene, RSS)
  
  
  scores <- expanded_data %>%
    group_by(Gene) %>%
    summarise(min_RSS = min(RSS)) %>%
    arrange(min_RSS)  # 按 t-value 降序排序
  
  
  return(scores)
}
