# BSFate: A Dynamical ODE Framework for Uncovering Regulatory Modules of Cell Fate Determination from Single-Cell RNA-Seq Data

BSFate is a computational framework designed to uncover the regulatory mechanisms governing cell fate decisions by integrating **Ordinary Differential Equation (ODE) modeling** with **single-cell RNA sequencing (scRNA-seq) data**. It identifies **bistable switch circuits** that regulate binary cell fate transitions and systematically determines **gene modules** associated with both target and alternative fates. By leveraging directed differentiation data, BSFate **pinpoints key driver regulators**, providing a quantitative assessment of **synergistic and antagonistic interactions** between gene modules. This ena
bles a deeper understanding of the dynamic regulatory landscape that drives cell fate commitment. Through applications to both simulated and experimental datasets, BSFate demonstrates high accuracy in identifying crucial regulators and offers new insights into the mechanistic basis of differentiation processes.

## Modeling a two-gene bistable switch governing cell fate transitions.
![Figure1-01](https://github.com/user-attachments/assets/cdf11b33-e1d1-40bd-a038-01430f410388)

**A. Waddington’s epigenetic landscape**
Illustrating the differentiation process. A **multipotent progenitor cell** (represented by a ball) traverses a dynamic landscape with bifurcating valleys, representing alternative cell fate trajectories.  
**B. A two-gene bistable switch circuit**
Regulating cell fate decisions. The circuit consists of two **mutually inhibitory genes**, \( g_1 \) and \( g_2 \), where \( x_1 \) and \( x_2 \) represent their expression levels.  
**C. Stable point analysis of the bistable system**
The parameters are set as follows:  
\[
a_1 = a_2 = b_1 = b_2 = k_1 = k_2 = 1
\]
Moderate sensitivity thresholds:  
\[
S_{a_1} = S_{b_1} = S_{a_2} = S_{b_2} = 0.5
\]
Sharp **bistable switching behavior** (Hill coefficients \( m = n = 4 \)). **Nullclines** for \( dx_1/dt = 0 \) and \( dx_2/dt = 0 \) are plotted, with their intersections indicating equilibrium points.  
**D. Basins of attraction**
This plot shows how different initial conditions converge toward distinct stable states. The **separatrix** defines the **critical threshold**, beyond which cells commit to one of the two alternative fates.




## **Installation**
We recommend installing the CytoTRACE 2 package using the devtools package from the R console. If you do not have devtools installed, you can install it by running install.packages("devtools") in the R console.
```
devtools::install_github("SDU-W-Zhanglab/BSFate")

library(BSFate)
```

## **Running** **BSfate**
Running BSfate is very simple, after loading the library, it only takes four steps to complete the prediction of cell fate determinants:

1.Obtain TF gene expression matrix and pseudo-time

2.Screen TFs that display switch-like or transient activation patterns

3.Calculate the significance scores of the candidate cell fate determinants

4.Order cell fate determinants by significance scores

The following shows specific applications on simulated data and two sets of real data.

- Simulated DATASET

The switch gene and transient gene of simulated data are used as known data here because of the simple branching of simulated data. Below, we will take Astrocyte lineage as an example.

![image](https://github.com/user-attachments/assets/07227edc-ca28-406e-a11e-c8e8a60189b7)
![image](https://github.com/user-attachments/assets/22846047-128b-45bc-a402-8862fcbd4654)




### Step 1 and Step 2: Known Data Acquisition

```r
# Load required libraries
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
library(BSFate)
library(circlize)
library(ComplexHeatmap)

# Load the astrocyte expression data
astrocyte_exprs_mat <- t(read.table("data/astrocyte_exprs_mat.txt"))

# Define transcription factors
Switch_TF = c("Scl", "Aldh1L", "Hes5", "Stat3")
Transient_TF = c("Olig2", "Sox8", "Myt1L", "Mash1", "Brn2", "Zic1", "Tuj1")
```

### Step 3: Cell Fate Driver Gene Identification

```r
# Perform bistable circuit test to identify potential driver genes
ODE_pair_results = get_bistable_circuitTest(astrocyte_exprs_mat, Switch_TF, Transient_TF)

# Identify driver genes from the results
driver_genes <- get_diverGenes(ODE_pair_results)
```

### Step 4: Gene Modules Identification

```r
# Identify gene modules associated with cell fate
Gene_modules <- modules_identification_A(astrocyte_exprs_mat)
```

### Step 5: Gene Modules Analysis

```r
# Perform analysis on gene modules with respect to transcription factors
module_analysis <- get_moduleAnalysis_A(astrocyte_exprs_mat, Switch_TF, Transient_TF)
```


- hESC DATASET
```
###step 1 Known data acquisition
library(stats)
library(entropy)
library("RColorBrewer")
data(scExp_hESC)
data(pesudo_hESC)
data(TF_human)
TF_Human=(TF_human)[,1]
scExp_hESC_TF=scExp_hESC[intersect(rownames(scExp_hESC),TF_Human),order(pesudo_hESC[,1],decreasing = F)]

###step 2 Screem candidate cell fate determinants
Switch_test_hESC=Switch_nls(scExp_hESC_TF)
Tansient_test_hESC=Tansient_nls(scExp_hESC_TF)
Pro_screen_TFs=Screen_TF(Switch_test_hESC,Tansient_test_hESC,top_n=20)
```
![tupian5_1](./image/hotplot_hESC.jpg)
![tupian5_2](./image/switch_transient_TF_hESC.jpg)

```
####step 3 Calculate the significance scores
Hs_switch=Pro_screen_TFs[,1]
Hs_transient=Pro_screen_TFs[,2]
SignificanceScore=get_SignificanceScore(scExp_hESC_TF,Hs_switch,Hs_transient,0.1)

####step 4 Order cell fate determinants
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)
```
![tupian6_1](./image/score_hESC.jpg)
![tupian6_2](./image/results_hESC.jpg)
```
Following are comparison of BSFate with Monocle3-DE, tradeSeq, ImpulseDE2, and scFates, on astrocyte branch. The precision is plotted as a function of the number of top-ranked genes involved in the four ground truth gene sets. And the prioritization of driver regulators predicted by BSFate. TFs are sorted in descending order according to their significance sores. When TFs are annotated within the ground truth sets, they are connected.
```
![tupian7_1](./image/other_DE_compare_hESC.jpg)
![tupian7_2](./image/GO_results_hESC.jpg)

## **Input** **files**

BSfate requires a single-cell RNA-sequencing gene expression object as input. This can include raw counts, as well as normalized counts, as long as normalization preserves ranking of input gene values within a cell. The gene expression matrix takes the following form.

![tupian8](./image/Input.jpg)

## **Output**

The output of BSfate is TFs and significance score ranking.

![tupian9](./image/Output.jpg)

## **Authors**

BSFate was developed in the Zhang lab by Naiqian Zhang,Yumiao Hou, Wenchao Xiu, Junchao Wang, Ling Sun, Yisheng Huang, Yong Zhang, Naiqian Zhang
## **Contact**

If you have any questions, please contact the BSfate team at nqzhang@email.sdu.edu.cn

## **Funding**

This work has been supported by the National Natural Science Foundation of China [62072277].
