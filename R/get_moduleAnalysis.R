
#' @title get_moduleAnalysis
#' @description  function to get module Analysis
#' @param data refer to time series of committed data
#' @param PPI_gene refer to bistable circuitGenes
#' @return   gene modules Analysis
#' @export



get_moduleAnalysis <- function(data,PPI_gene){
  
  A=intersect(PPI_gene[,1],Switch_TF)
  B=intersect(PPI_gene[,2],Transient_TF)
  expr_df <- data
  
  cell_names <- colnames(expr_df)
  time_point <- sub(".*_([^_]+)_[^_]+$", "\\1", cell_names)
  
  
  library(circlize)
  library(ComplexHeatmap)
  
  
  module1 <- B
  module2 <- A
  exp_matrix <- expr_df
  
  table(time_point)
  
  exp_matrix00 <- exp_matrix[c(Switch_TF, Transient_TF),c(1:70)]
  exp_matrix12 <- exp_matrix[c(Switch_TF, Transient_TF),c(71:140)]
  exp_matrix24 <- exp_matrix[c(Switch_TF, Transient_TF),c(141:210)]
  exp_matrix48 <- exp_matrix[c(Switch_TF, Transient_TF),c(211:280)]
  exp_matrix72 <- exp_matrix[c(Switch_TF, Transient_TF),c(281:350)]
  exp_matrix96 <- exp_matrix[c(Switch_TF, Transient_TF),c(351:421)]
  
  dim(exp_matrix72)
  
  expression_matrix <- list(
    "Time 1" = exp_matrix00,
    "Time 2" = exp_matrix12,
    "Time 3" = exp_matrix24,
    "Time 4" = exp_matrix48,
    "Time 5" = exp_matrix72,
    "Time 6" = exp_matrix96
  )
  
  
  threshold <- 2
  co_expression_matrices <- list()
  
  for (time in names(expression_matrix)) {
    expr <- expression_matrix[[time]]
    
    co_expression_matrix <- matrix(0, nrow = length(module1), ncol = length(module2),
                                   dimnames = list(module1, module2))
    for (gene1 in module1) {
      for (gene2 in module2) {
        co_expression_matrix[gene1, gene2] <- sum(expr[gene1, ] > threshold & expr[gene2, ] > threshold)
      }
    }
    co_expression_matrices[[time]] <- co_expression_matrix
  }
  
  ######################################################
  
  ###########################  plot
  
  ################################################## 
  par(mfrow = c(1, 6))  
  
  for (time in names(co_expression_matrices)[1:6]) {
    circos.clear()  
    
    filtered_connections <- as.data.frame(as.table(co_expression_matrices[[time]]))
    filtered_connections[which(filtered_connections[,3]==0),3]=0.0001
    filtered_connections <- filtered_connections[filtered_connections$Freq > 0, ]
    
    max_value <- max(unlist(co_expression_matrices))
    min_value <- min(unlist(co_expression_matrices))
    
    connection_colors <- colorRamp2(c( 0,max_value/2,max_value), c('white', '#ffa17a','#fb575f'))
    
    grid.col <- c(rep('#74aed4',length(module1)),rep('#F3b3c6',length(module2)))
    
    chordDiagram(
      filtered_connections, grid.col = grid.col,
      col = connection_colors(filtered_connections$Freq),
      transparency = 0.3
    )
    
    title(main = paste("Time:", time))
  }
  
  gene_expr <- data
  module1 <- intersect(PPI_gene[,1],Switch_TF)
  module2 <- intersect(PPI_gene[,2],Transient_TF)
  
  
  
  expr1 <- gene_expr[module1, ]
  expr2 <- gene_expr[module2, ]
  
  
  E_M1 <- colMeans(expr1)
  E_M2 <- colMeans(expr2)
  
  
  
  
  alpha <- 1.0     
  beta <- 0.5     
  
  
  max_sum <- max(E_M1 + E_M2)
  max_diff <- max(abs(E_M1 - E_M2))
  
  
  scores <- sapply(seq_along(E_M1), function(i) {
    a <- E_M1[i]
    b <- E_M2[i]
    total_sum <- a + b
    diff <- abs(a - b)
    alpha * (total_sum / max_sum) - beta * (diff / max_diff)
  })
  
  
  max_index <- which.max(scores)
  max_score <- scores[max_index]
  
  
  cat(sprintf("The highest score for module collaborative expression: %.3f\n", max_score))
  cat(sprintf("Time point: %d\n", max_index))
  
  ######################################################
  
  ########################### Module Synergy Score.
  
  ################################################## 
  plot(scores,ylim=c(0,1), type = "l", col = "#ff948a", lwd = 1,
       xlab = "cells", ylab = "Synergy Score", frame.plot = FALSE,las=1)
  points(max_index, max_score, col = "#c82423", pch = 19, cex = 1.2)  # 标记最大点
  
  abline(v = max_index, col = "#708090", lty = 2, lwd = 2)
  
  text(200, 0.4, labels=c("x=",max_index), cex=1, col="#708090")
  
  
  
  antagonistic_scores <- sapply(seq_along(E_M1), function(i) {
    a <- E_M1[i]
    b <- E_M2[i]
    total_sum <- a + b
    diff <- abs(a - b)
    antagonistic_score <- alpha * ((max_diff - diff) / max_diff) - beta * (total_sum / max_sum)
    return(antagonistic_score)
  })
  
  
  max_antagonistic_index <- which.max(antagonistic_scores)
  max_antagonistic_score <- antagonistic_scores[max_antagonistic_index]
  
  
  cat(sprintf("The highest score for module antagonism: %.3f\n", max_antagonistic_score))
  cat(sprintf("Time point:%d\n", max_antagonistic_index))
  
  ######################################################
  
  ########################### Module Antagonism Score.
  
  ################################################## 
  plot(antagonistic_scores, type = "l", col = "#ff948a", lwd = 1,
       xlab = "cells", ylab = "Antagonistic Score",frame.plot = FALSE,las=1)
  
  text(100, 0, labels=c("x=",antagonistic_scores), cex=1, col="#708090")
}



