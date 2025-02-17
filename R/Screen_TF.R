
#' @title screen switch-like and transient TFs
#' @description  function to Pro-screening candidate genes
#' @param scExp refer to time series of committed data
#' @return  Pro-screening candidate genes
#' @export







Switch_nls<-function(scExp){
  
  nls_fit<-function(Data){
    # Data <- as.numeric(neuron_exprs_mat[c("Tuj1"),])
    data <- list(
      y = as.vector(Data),
      x = (1:length(Data)/length(Data))
    )
    
    a <- round(length(which(Data<(max(Data)/2)))/length(Data),1)
    
    nls_model <- try(nls(y ~ k / (1 + exp(-10*σ * (x - μ))),data = data, start = list( k = max(Data), σ = 1, μ = a)),silent = T)
    
    if(!'try-error' %in% class(nls_model))
    {
      
      params <- coef(nls_model)
      fitted_values=predict(nls_model)
      RSS <- sum((data[["y"]] - fitted_values)^2)
      TSS <- sum((data[["y"]] - mean(data[["y"]]))^2)
      rsquared <- 1 - RSS / TSS
      
    }else{
      params=rep(NA,3)
      rsquared=NA
      fitted_values=rep(NA,length(Data))
    }
    model_summary=c(rsquared,params,fitted_values)
    names(model_summary)=c("R-squared","k", "σ","μ")
    return(model_summary)
  }
  
  switch_test=t(apply(scExp,1,nls_fit))
  switch_test=as.data.frame(na.omit(switch_test))
  switch_test <- na.omit(switch_test)
  
  return(switch_test)}



Transient_nls<-function(scExp){
  
  nls_fit<-function(Data){
    
    # Data <- neuron_exprs_mat[c("Scl"),]
    # plot(neuron_exprs_mat[c("Scl"),])
    data <- list(
      y = as.vector(Data),
      x = (1:length(Data)/length(Data))
    )
    
    
    a <- round(which(Data==(max(Data)))/length(Data),1)
    
    
    nls_model <- try(nls(y ~ k *exp(-10*σ * (x - μ)^2),data = data, start = list(k = max(Data), σ = 1, μ = a)),silent = T)
    
    
    
    if(!'try-error' %in% class(nls_model))           
    {
      
      params <- coef(nls_model)
      fitted_values=predict(nls_model)
      
      RSS <- sum((data[["y"]] - fitted_values)^2)
      TSS <- sum((data[["y"]] - mean(data[["y"]]))^2)
      rsquared <- 1 - RSS / TSS
      
    }else{
      params=rep(NA,3)
      rsquared=NA
      fitted_values=rep(NA,length(Data))
    }
    model_summary=c(rsquared,params,fitted_values)
    names(model_summary)=c("R-squared","k", "σ","μ")  ###a b c
    return(model_summary)
  }
  
  
  transient_test=t(apply(scExp,1,nls_fit))
  transient_test=as.data.frame(na.omit(transient_test))
  transient_test <- na.omit(transient_test)
  return(transient_test)}






Screen_TF<-function(Switch_test, Tansient_test,top_n){
  
  Switch_up_test=Switch_test[Switch_test[,"σ"]>0.1&Switch_test[,"σ"]<15&Switch_test[,"μ"]<0.8&Switch_test[,"μ"]>0.2&Switch_test[,"k"]>0.2,]
  Switch_up_test=Switch_up_test[order(Switch_up_test[,"R-squared"],decreasing = T),]
  Switch_up_test_top_n=Switch_up_test[1:min(top_n,nrow(Switch_up_test)),]
  
  
  Tansient_test=Tansient_test[which(Tansient_test[,"μ"]<0.7&Tansient_test[,"μ"]>0.3&Tansient_test[,"σ"]>0.1&Tansient_test[,"σ"]<15&Tansient_test[,"k"]>0.2),]
  Tansient_test=Tansient_test[order(Tansient_test[,"R-squared"],decreasing = T),]
  Tansient_test_top_n=Tansient_test[1:min(top_n,nrow(Tansient_test)),]
  
  TF_candidate=cbind(rownames(Switch_up_test_top_n),rownames(Tansient_test_top_n))
  
  return(TF_candidate)}
