
library(minpack.lm)  
library(dplyr)
library(tidyr)
library(tidyverse) 
library(dyno)
library(R.matlab)
library(deSolve)
library(sindyr)
library(nleqslv)
library(numDeriv)  
library(magicfor)
library(car) 
library(mgcv)  
library(ggraph)
library(deSolve)

library(circlize)
library(ComplexHeatmap)




astrocyte_exprs_mat <- t(read.table("data/astrocyte_exprs_mat.txt"))

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








