

library(minpack.lm)  # 加载 minpack.lm 包
library(dplyr)
library(tidyr)
library(tidyverse) 
library(dyno)
library(R.matlab)
library(deSolve)
library(sindyr)
library(nleqslv)##求解非线性微分方程组
library(numDeriv)  ###计算雅可比
library(magicfor)
library(car) 
library(mgcv)  # GAM函数

library(ggraph)
library(deSolve)

library(circlize)
library(ComplexHeatmap)

###################
########数据读取
###################

setwd("")

astrocyte_exprs_mat <- t(read.table("astrocyte_exprs_mat.txt"))

###################
########基因预筛选
###################
Switch_TF=c("Scl","Aldh1L","Hes5","Stat3")
Transient_TF=c("Olig2","Sox8","Myt1L","Mash1","Brn2","Zic1","Tuj1")



###################################################
########  cell fate driver gene identification
###########################################



ODE_pair_results=get_bistable_circuitTest(astrocyte_exprs_mat,Switch_TF,Transient_TF)

driver_genes <- get_diverGenes(ODE_pair_results)


###################
########Gene modules identification
###################
Gene_modules <- modules_identification_A(astrocyte_exprs_mat)

###################
########Gene modules Analysis
###################
module_analysis <- get_moduleAnalysis_A(astrocyte_exprs_mat,Switch_TF,Transient_TF)








