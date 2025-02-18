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

### Simulated DATASET

The switch gene and transient gene of simulated data are used as known data here because of the simple branching of simulated data. Below, we will take Astrocyte lineage as an example.

![S_1_画板 1](https://github.com/user-attachments/assets/1e1c4c9c-ec1c-4cdd-99fe-8965ffd00013)




#### Step 1 and Step 2: Known Data Acquisition

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

#### Step 3: Cell Fate Driver Gene Identification

```r
# Perform bistable circuit test to identify potential driver genes
ODE_pair_results = get_bistable_circuitTest(astrocyte_exprs_mat, Switch_TF, Transient_TF)

# Identify driver genes from the results
driver_genes <- get_diverGenes(ODE_pair_results)
```

#### Step 4: Gene Modules Identification

```r
# Identify gene modules associated with cell fate
Gene_modules <- modules_identification_A(astrocyte_exprs_mat)
```
![image](https://github.com/user-attachments/assets/fb3903b5-ee8b-4ee8-aba6-d76742584415)

#### Step 5: Gene Modules Analysis

```r
# Perform analysis on gene modules with respect to transcription factors
module_analysis <- get_moduleAnalysis_A(astrocyte_exprs_mat, Switch_TF, Transient_TF)
```
![image](https://github.com/user-attachments/assets/974eb456-c1a4-450c-b528-99aebab426b5)



### hESC DATASET
  
#### Step 1: Load Required Libraries

```r
# Load the required libraries
library(BSFate)
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
library(clusterProfiler)
library(STRINGdb)
library(igraph)
library(ggraph)
```

#### Step 2: Read Data

```r
# Load data
load("scExp_mESC.RData")
load("pesudo_mESC.RData")

# Convert to matrix format
pesudo_mESC = as.matrix(pesudo_mESC)
scExp_mESC = as.matrix(scExp_mESC)

# Combine pseudo-time data with expression data
scExp_mESC <- rbind(pesudo_mESC[, 1], scExp_mESC)

# Order the columns by the pseudo-time
scExp_mESC = scExp_mESC[, order(scExp_mESC[1, ])]

# Load transcription factor data
TF_mouse = read.table("TF_mouse.txt")
TF_mouse = TF_mouse[, 1]

# Subset the expression data based on common transcription factors
scExp_mESC = scExp_mESC[intersect(TF_mouse, rownames(scExp_mESC)), ]
```

#### Step 3: Data Preprocessing

```r
# Calculate row variances
row_variances <- apply(scExp_mESC, 1, var)

# Filter out rows with low variance
scExp_mESC <- scExp_mESC[row_variances >= 0.03, ]

# Display the dimensions of the filtered data
dim(scExp_mESC)
```

#### Step 4: Pre-screening Candidate Genes

```r
# Perform switch and transient tests
Switch_test_mESC = Switch_nls(scExp_mESC)
Transient_test_mESC = Transient_nls(scExp_mESC)

# Select top N candidate genes
top_n = 30
Pro_screen_TFs = Screen_TF(Switch_test_mESC, Transient_test_mESC, top_n)

# Extract the unique switch and transient transcription factors
Switch_TF = unique(Pro_screen_TFs[, 1])
Transient_TF = unique(Pro_screen_TFs[, 2])
```

#### Step 5: Cell Fate Driver Gene Identification

```r
# Identify bistable circuit genes
bistable_circuitGenes = get_bistable_circuitGenes(scExp_mESC, Switch_TF, Transient_TF)

# Identify driver genes from the bistable circuit genes
driver_genes = get_diverGenes(bistable_circuitGenes)
```

#### Step 6: Gene Modules Identification

```r
# Identify gene modules associated with cell fate
Gene_modules = modulePlot(scExp_mESC, driver_genes)
```

#### Step 7: Gene Modules Analysis

```r
# Analyze the gene modules
module_analysis = get_moduleAnalysis(scExp_mESC, Gene_modules)
```


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
