
library(minpack.lm)  # 

library(dplyr)
library(tidyr)
library(tidyverse) 
library(dyno)
library(tidyverse)
library(R.matlab)
library(deSolve)
library(sindyr)
library(nleqslv)##
library(numDeriv)  ###
library(magicfor)
library(car) 
library(mgcv)  # 
library(tidyverse)
library(clusterProfiler) #

library(STRINGdb)
library(igraph)
library(ggraph)


###################
########read data 
###################

load("scExp_mESC.RData")
load("pesudo_mESC.RData")

pesudo_mESC=as.matrix(pesudo_mESC)
scExp_mESC=as.matrix(scExp_mESC)

scExp_mESC <- rbind(pesudo_mESC[,1],scExp_mESC)
scExp_mESC=scExp_mESC[,order(scExp_mESC[1,])]


TF_mouse= read.table("TF_mouse.txt")
TF_mouse=(TF_mouse)[,1]

scExp_mESC = scExp_mESC[intersect(TF_mouse,rownames(scExp_mESC)),]

###################
######## data preprocessing
###################

row_variances <- apply(scExp_mESC, 1, var)
mean(row_variances)

scExp_mESC <- scExp_mESC[row_variances >= 0.03, ]
dim(scExp_mESC)  ###row gene 


###################
########Pro-screening candidate genes
###################
Switch_test_mESC=Switch_nls(scExp_mESC)
Transient_test_mESC=Transient_nls(scExp_mESC)

top_n=30
Pro_screen_TFs=Screen_TF(Switch_test_mESC, Transient_test_mESC,top_n)

Switch_TF=unique(Pro_screen_TFs[,1])
Transient_TF=unique(Pro_screen_TFs[,2])


####################################################
########  cell fate driver gene identification
####################################################

bistable_circuitGenes <- get_bistable_circuitGenes(scExp_mESC,Switch_TF,Transient_TF)

driver_genes <- get_diverGenes(bistable_circuitGenes)

###################
########Gene modules identification
###################
Gene_modules  <- modulePlot(scExp_mESC,driver_genes)


###################
########Gene modules Analysis
###################
module_analysis <- get_moduleAnalysis(scExp_mESC,Gene_modules)


