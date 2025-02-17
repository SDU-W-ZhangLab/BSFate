
#' @title get_bistable_circuitGenes
#' @description  function to get bistable circuit genes
#' @param data refer to time series of committed data
#' @param Switch_TF Pro-screening candidate genes of Switch_TF
#' @param Transient_TF Pro-screening candidate genes of Transient_TF
#' @return  bistable circuitGenes
#' @export






get_bistable_circuitGenes <- function(data,Switch_TF,Transient_TF){
  
  
  gene_pairs <-data.frame(expand.grid(Switch_TF, Transient_TF), stringsAsFactors = FALSE)
  colnames(gene_pairs) <- c("GeneSet_1", "GeneSet_2")
  
  
  normalize_matrix <- function(matrix) {
    min_val <- min(matrix)
    max_val <- max(matrix)
    (matrix - min_val) / (max_val - min_val)
  }
  
  
  
  gene_model_full <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- a * x^n / (S^n + x^n) + b * S^n / (S^n + y^n) - c * x
      dy <- d * y^m / (S^m + y^m) + e * S^m / (S^m + x^m) - f * y
      list(c(dx, dy))
    })
  }
  initial_params_full <- c(a = 0.5, b = 0.5, c = 1, d = 0.5, e = 0.5, f = 1, S = 0.5, n = 4, m = 4)
  
  
  resid_func <- function(params, t, x, y, model_func) {
    state0 <- c(x = 0.05, y = 0.01)
    out <- ode(y = state0, times = t, func = model_func, parms = params)
    out <- as.data.frame(out)
    residuals <- c(x - out$x, y - out$y)
    return(residuals)
  }
  
  Rss=c()
  
  for ( i in 1:nrow(gene_pairs)){
    
    exp_matrix_genepair=t(data[c(as.character(gene_pairs[i,1]),as.character(gene_pairs[i,2])),])
    exp_matrix <- normalize_matrix(exp_matrix_genepair)
    initial_params_full <- c(a = 0.5, b = 0.5, c = 1, d = 0.5, e = 0.5, f = 1, S = 0.5, n = 4, m = 4)
    
    df <- data.frame(
      time = seq(0, 1, length.out = nrow(exp_matrix)),
      x = as.numeric(exp_matrix[, 1]),
      y = as.numeric(exp_matrix[, 2])
    )
    fit_full <- nls.lm(
      par = initial_params_full, 
      fn = function(params) resid_func(params, df$time, df$x, df$y, gene_model_full), 
      control = nls.lm.control(maxiter = 50)
    )
    RSS_full <- sum(fit_full$fvec^2)
    names(RSS_full)=paste(c(as.character(gene_pairs[i,1]),as.character(gene_pairs[i,2])),collapse  = "_")
    Rss=c(Rss,RSS_full)
  }
  
  results <- Rss[order(Rss,decreasing=F)]
  
  return(results)
  
}
