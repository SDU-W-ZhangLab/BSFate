#' @title modules_identification
#' @description  function to get modules identification
#' @param data refer to time series of committed data
#' @param Switch_TF Pro-screening candidate genes of Switch_TF
#' @param Transient_TF Pro-screening candidate genes of Transient_TF
#' @return   gene modules
#' @export


get_moduleAnalysis_A <- function(data,Switch_TF,Transient_TF){
  
  
  set.seed(123)
  module1 <- Transient_TF 
  module2 <- Switch_TF
  exp_matrix  <- data
  
  
  exp_matrix00 <- exp_matrix[c(Switch_TF,Transient_TF),c(1:100)]
  exp_matrix06 <- exp_matrix[c(Switch_TF,Transient_TF),c(101:200)]
  exp_matrix12 <- exp_matrix[c(Switch_TF,Transient_TF),c(201:500)]
  exp_matrix24 <- exp_matrix[c(Switch_TF,Transient_TF),c(501:700)]
  exp_matrix36 <- exp_matrix[c(Switch_TF,Transient_TF),c(701:1000)]

  expression_matrix <- list(
    "Time 1" = exp_matrix00,
    "Time 2" = exp_matrix06,
    "Time 3" = exp_matrix12,
    "Time 4" = exp_matrix24,
    "Time 5" = exp_matrix36
  )
  

  threshold <- 1
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
  
  # 验证共表达矩阵
  for (time in names(co_expression_matrices)) {
    print(paste("共表达矩阵 - 时间点:", time))
    print(co_expression_matrices[[time]])
  }
  
  ######################################################
  
  ########################### Module plot
  
  ################################################## 
  
  par(mfrow = c(1, 5)) 
  for (time in names(co_expression_matrices)[1:5]) {
    circos.clear()  # 清除绘图状态
    # time <- c("Time 5")
    filtered_connections <- as.data.frame(as.table(co_expression_matrices[[time]]))
    filtered_connections[which(filtered_connections[,3]==0),3]=0.0001
    filtered_connections <- filtered_connections[filtered_connections$Freq > 0, ]
    
    max_value <- max(unlist(co_expression_matrices))
    min_value <- min(unlist(co_expression_matrices))
    # connection_colors <- colorRamp2(c(min_value, 0, max_value), c("blue", "white", "red"))
    # connection_colors <- colorRamp2(c(min_value, -5, 0, 5, max_value), 
    #                                 c("darkblue", "blue", "white", "pink", "red"))
    connection_colors <- colorRamp2(c( 0,max_value/2,max_value), c('white', '#ffa17a','#fb575f'))
    
    grid.col <- c(rep('#74aed4',length(module1)),rep('#F3b3c6',length(module2)))
    
    chordDiagram(
      filtered_connections, grid.col = grid.col,
      col = connection_colors(filtered_connections$Freq),
      transparency = 0.3
    )
    
    title(main = paste("Time:", time))
  }
  
  
  
  ######################################################
  
  ########################### Module Synergy Score.
  
  ################################################## 
  module1 <- Transient_TF  
  module2 <- Switch_TF
  exp_matrix  <- data
  
  
  exp_matrix00 <- exp_matrix[c(Switch_TF,Transient_TF),c(1:100)]
  exp_matrix06 <- exp_matrix[c(Switch_TF,Transient_TF),c(101:200)]
  exp_matrix12 <- exp_matrix[c(Switch_TF,Transient_TF),c(201:500)]
  exp_matrix24 <- exp_matrix[c(Switch_TF,Transient_TF),c(501:700)]
  exp_matrix36 <- exp_matrix[c(Switch_TF,Transient_TF),c(701:1000)]
  
  exp_matrix_all <- list(exp_matrix00,exp_matrix06,exp_matrix12,exp_matrix24,exp_matrix36)
  scores_all <- list()
  for(i in 1:5){
    gene_expr <- exp_matrix_all[[i]]
    module1 <- Switch_TF
    module2 <- Transient_TF

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
    scores <- list(scores)
    scores_all <-c(scores_all,scores)
  }
  
  
  
  
  par(mfrow = c(1,2), mar = c(2, 2, 1, 1), oma = c(1, 1, 1, 1))
  
  boxplot(scores_all,
          col = c("#8ccdbf", "#cde0a5", "#f9db95", "#ef8476", "#c5a8ce"),
          names = c("1-200", "201-400", "401-600", "601-800", "801-1000"),ylab="Synergy Score")
  
  
  
  ######################################################
  
  ########################### Module Antagonism Score.
  
  ################################################## 
  
  antagonistic_scores_all <- c()
  for( i in 1:5){
    
    gene_expr <- exp_matrix_all[[i]]
    module1 <- Switch_TF
    module2 <- Transient_TF
    

    expr1 <- gene_expr[module1, ]
    expr2 <- gene_expr[module2, ]
    

    E_M1 <- colMeans(expr1)
    E_M2 <- colMeans(expr2)
    

    antagonistic_scores <- sapply(seq_along(E_M1), function(i) {
      a <- E_M1[i]
      b <- E_M2[i]
      total_sum <- a + b
      diff <- abs(a - b)
      antagonistic_score <- alpha * ((max_diff - diff) / max_diff) - beta * (total_sum / max_sum)
      return(antagonistic_score)
    })
    antagonistic_scores <- list(antagonistic_scores)
    antagonistic_scores_all <- c(antagonistic_scores_all,antagonistic_scores)
  }
  
  
  

  boxplot(antagonistic_scores_all, col = c("#8ccdbf", "#cde0a5", "#f9db95", "#ef8476", "#c5a8ce"),
          names = c("1-200", "201-400", "401-600", "601-800", "801-1000"),ylab = "antagonistic_scores")
  
  
  
}

