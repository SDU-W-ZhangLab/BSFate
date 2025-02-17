
#' @title get_bistable_circuitTest
#' @description  function to get bistable circuit genes
#' @param data refer to time series of committed data
#' @param Switch_TF Pro-screening candidate genes of Switch_TF
#' @param Transient_TF Pro-screening candidate genes of Transient_TF
#' @return  bistable circuitGenes
#' @export
#' @name get_bistable_circuitTest


library(ggplot2) #library for plotting
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
require(FME)
library(patchwork)



get_bistable_circuitTest <- function(data,Switch_TF,Transient_TF){
  
  
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
  
  gene_model_simple <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- a * x^n / (S^n + x^n) - c * x
      dy <- d * y^m / (S^m + y^m) - f * y
      list(c(dx, dy))
    })
  }
  
 
  resid_func <- function(params, t, x, y, model_func) {
    state0 <- c(x = 0.05, y = 0.01)
    out <- ode(y = state0, times = t, func = model_func, parms = params)
    out <- as.data.frame(out)
    residuals <- c(x - out$x, y - out$y)
    return(residuals)
  }
  
  
  get_LPT=function(df){
    
    initial_params_full <- c(a = 0.5, b = 0.5, c = 1, d = 0.5, e = 0.5, f = 1, S = 0.5, n = 4, m = 4)
    initial_params_simple <- c(a = 0.5, c = 1, d = 0.5, f = 0.5, S = 0.5, n = 4, m = 4)
    
    
 
    fit_full <- nls.lm(
      par = initial_params_full, 
      fn = function(params) resid_func(params, df$time, df$x, df$y, gene_model_full), 
      control = nls.lm.control(maxiter = 50)
    )
    
    fit_simple <- nls.lm(
      par = initial_params_simple, 
      fn = function(params) resid_func(params, df$time, df$x, df$y, gene_model_simple), 
      control = nls.lm.control(maxiter = 100)
    )
    
  
    params_full <- fit_full$par
    params_simple <- fit_simple$par
    RSS_full <- sum(fit_full$fvec^2)
    RSS_simple <- sum(fit_simple$fvec^2)
    
    
   
    N <- nrow(df) 
    k_full <- length(coef(fit_full))  
    k_simple <- length(coef(fit_simple))  
    delta_k <- k_full - k_simple  
    

    chi2 <- N * log(RSS_simple / RSS_full)
    print(chi2)

    p_value <- pchisq(chi2, df = delta_k, lower.tail = FALSE)
    
    
    predictions_full <- ode(y = c(x = 0.05, y = 0.01), 
                            times = df$time, 
                            func = gene_model_full, 
                            parms = params_full)
    
  
    LPT_test=list(RSS_full,RSS_simple,chi2,p_value, df,predictions_full)
    names(LPT_test)=c("RSS_full","RSS_simple","chi2","p_value","true_value","predict_value")
    
    return(LPT_test)
  }
  

  
  Fit_result=list()
  
  for ( i in 1:nrow(gene_pairs)){
    
    exp_matrix_genepair=t(data[c(as.character(gene_pairs[i,1]),as.character(gene_pairs[i,2])),])
    exp_matrix <- normalize_matrix(exp_matrix_genepair)

    df <- data.frame(
      time = seq(0, 1, length.out = nrow(exp_matrix)),
      x = as.numeric(exp_matrix[, 1]),
      y = as.numeric(exp_matrix[, 2])
    )

    Fit_result[[i]]=get_LPT(df)

  }
  
  names(Fit_result)=apply(gene_pairs,1,function (x) paste(x,collapse  = "_"))
  
  Fit_full_RSS=unlist(sapply(Fit_result,function(x) x[1]))
  names(Fit_full_RSS)=names(Fit_result)
  ranked_bistable=Fit_full_RSS[order(Fit_full_RSS,decreasing = F)]
  
  return(ranked_bistable)
  
}


