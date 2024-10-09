---
title: "16wk Tgfbr1^M318R/+ vs WT murine aorta scRNAseq CoGAPs Analysis"
#adapted from R vignette for CoGAPS https://www.biorxiv.org/content/10.1101/2022.07.09.499398v1 
#https://satijalab.org/signac/articles/monocle.html 
author: "Emily Bramel, Wendy A. Espinoza Camejo, Jacob T. Mitchell"
date: "last updated 9/11/23"
---
##install necessary packages (only once)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.14")

#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils',), force=TRUE)

#install.packages("devtools")
#devtools::install_github("FertigLab/CoGAPS")
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("MarioniLab/DropletUtils")
#devtools::install_github("satijalab/seurat-wrappers")
#if (!require("BiocManager", quietly = TRUE))

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Nebulosa")

#Load Libraries
library(projectR)
library(CoGAPS)
library(ggpubr)
library(monocle3)
library(CoGAPS)
library(Matrix)
library(dplyr)
library(DropletUtils)
library(SeuratWrappers)

#visualization packages 
library("fgsea", quietly = TRUE)
library("msigdbr", quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(forcats)
library(ggplot2)
library("Nebulosa")

# Load data
#first run through line 202 of Bramel_16wkControlandLDS_MouseAorta_scRNAseqScript.R 
aorta.combined.sct.VSMCsonly <- subset(aorta.combined.sct_named, idents = c("VSMCs"))
aorta.combined.sct.VSMCsforCoGAPs<- subset(aorta.combined.sct.VSMCsonly, idents = c("VSMCs"))
aorta.combined.sct.VSMCsforCoGAPs <- RunPCA(aorta.combined.sct.VSMCsforCoGAPs, verbose = FALSE)
aorta.combined.sct.VSMCsforCoGAPs <- RunUMAP(aorta.combined.sct.VSMCsforCoGAPs, reduction = "pca", dims = 1:30)
aorta.combined.sct.VSMCsforCoGAPs <- FindNeighbors(aorta.combined.sct.VSMCsforCoGAPs, dims = 1:30, verbose = FALSE)
aorta.combined.sct.VSMCsforCoGAPs <- FindClusters(aorta.combined.sct.VSMCsforCoGAPs, verbose = FALSE, resolution = 0.15)

#after clustering set the default assay to RNA and find the most varibale genes 
DefaultAssay(aorta.combined.sct.VSMCsforCoGAPs)<- "RNA"
dat <- GetAssayData(object = aorta.combined.sct.VSMCsforCoGAPs, assay="RNA")
dat <- as.data.frame(dat)
write.csv(dat, "./aortamatrix.csv")


#dat <- t(dat) #to transform data if necessary this is not needed

#Set parameters before running CoGAPS. Parameters are managed with a CogapsParams object. This object will store all parameters needed to run CoGAPS and provides a simple interface for viewing and setting the parameter values.
params <- CogapsParams(nPatterns=8, nIterations=30000, seed=42,
                       distributed="Single-Cell", # or distributed="Single-Cell","Genome-wide"
                       sparseOptimization=TRUE)

#To run distributed CoGAPS, which is strongly recommended for most large datasets, you must also call the setDistributedParams function.
params <- setDistributedParams(params, nSets=20, minNS=8, maxNS=23, cut=8)

#run coGAPs
result<- CoGAPS(
  data = "./aortamatrix.csv",
  params = params)

saveRDS(result, "CoGAPS_8patterns_16wkLDSmice_result.rds")
