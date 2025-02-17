

#' @title modules_identification
#' @description  function to get modules identification
#' @param data refer to time series of committed data
#' @param ODE_pair_results refer to bistable circuitGenes
#' @return   gene modules
#' @export



modules_identification <- function(data,ODE_single_results){
  
  
  string_db <- STRINGdb$new( version="11.0b", species=10090,                  
                             score_threshold=400, input_directory="")
  
  gene=ODE_single_results[,1]
  
  
  
  gene <- gene %>% bitr(fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db", 
                        drop = T)
  
  
  data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                        removeUnmappedRows = TRUE)
  
  
  data_links <- data_mapped$STRING_id[1:100] %>% string_db$get_interactions()
  
  links <- data_links %>%
    mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
    mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
    dplyr::select(from, to , last_col()) %>% 
    dplyr::rename(weight = combined_score)
  
  nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  
  net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
  
  igraph::V(net)$deg <- igraph::degree(net) 
  igraph::V(net)$size <- igraph::degree(net)/5 #
  igraph::E(net)$width <- igraph::E(net)$weight/10
  
  
  from_counts <- table(links$from)
  to_counts <- table(links$to)
  
  
  from_counts_df <- as.data.frame(from_counts)
  to_counts_df <- as.data.frame(to_counts)
  
  
  colnames(from_counts_df) <- c("from", "from_c")
  colnames(to_counts_df) <- c("to", "to_c")
  
  links <- as_tibble(links)
  
  
  links_2 <- links %>%
    left_join(from_counts_df, by = "from") %>%
    left_join(to_counts_df, by = "to") %>%
    filter(!(from_c == 1 & to_c == 1)) %>%
    dplyr::select(1, 2, 3)  # 使用 dplyr::select() 显式调用
  
  
  
  nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  
  net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
  
  igraph::V(net_2)$deg <- igraph::degree(net_2)
  igraph::V(net_2)$size <- igraph::degree(net_2)/5
  igraph::E(net_2)$width <- igraph::E(net_2)$weight/10
  
  
  nodes_2 <- nodes_2 %>%
    mutate(module = case_when(
      gene %in% module_A ~ "Module_A",
      gene %in% module_B ~ "Module_B",
      TRUE ~ "Other"
    ))
  
  
  net_2 <- igraph::graph_from_data_frame(d = links_2, vertices = nodes_2, directed = F)
  
  
  igraph::V(net_2)$deg <- igraph::degree(net_2)
  igraph::V(net_2)$size <- igraph::degree(net_2) / 5
  igraph::E(net_2)$width <- igraph::E(net_2)$weight / 10
  
  
  
  
  g <- net_2  
  
  ######################################################
  
  ########################### PPI _plot
  
  ################################################## 
  
  ggraph(g, layout = layout_components(net_2)) + 
    geom_edge_fan(aes(edge_width = width), color = "#C0C0C0", show.legend = F) +
    geom_node_point(aes(size = size, color = module), alpha = 0.7) +
    geom_node_text(aes(filter = deg > 5, label = name), size = 3.5, repel = TRUE) +
    scale_edge_width(range = c(0.2, 1)) +
    scale_size_continuous(range = c(1, 10)) +
    scale_color_manual(values = c("Module_A" = "#E41A1C", "Module_B" = "#1E90FF", "Other" = "gray")) + 
    guides(size = FALSE, color = guide_legend(title = "Module")) +
    theme_graph() +
    theme(legend.position = "right")
  
  
  our_net <- net_2
  
  BSfate_all <- names(degree(our_net)[order(degree(our_net),decreasing=T)])
  
  module_A_PPI <-   intersect(BSfate_all,A)
  module_B_PPI <-   intersect(BSfate_all,B)
  
  PPI_gene <- cbind(module_A_PPI,module_B_PPI)
  
  
  
  ######################################################
  
  ########################### heatmap plot
  
  ################################################## 
  
  
  data_n <- data
  data_n <- t(apply(data_n,1,function(x) {
    (x - min(x))/(max(x) - min (x))
  }))
  cell_names <- colnames(data_n)
  time_point <- sub(".*_([^_]+)_[^_]+$", "\\1", cell_names)
  
  A <- unique(PPI_plot[,1])
  B <- unique(PPI_plot[,2])
  
  genes <- c(A,B)
  
  expr_matrix <- t(data_n[genes,])
  
  cor_matrix <- cor(expr_matrix,method = 'spearman')
  
  pheatmap(cor_matrix, 
           # cluster_rows = F,
           # cluster_cols = F,
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean", 
           clustering_method = "complete", 
           breaks = seq(-1, 1, length.out = 51),
           color = colorRampPalette(c('#000080',"white", "red"))(50)
  )
  
  
  return(PPI_gene)
}



