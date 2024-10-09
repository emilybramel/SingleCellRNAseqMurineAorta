---
title: "Re-analysis of LDS Human Aortic Samples from Pedroza et al. 2023 scRNAseq Analysis includes ProjectR and CoGAPS"
#adapted from https://satijalab.org/seurat/ 
#reference for Pedroza et al. 2023 https://doi.org/10.1016/j.jtcvs.2023.07.023 
authors: "Emily Bramel, Wendy A. Espinoza Camejo, Jacob T. Mitchell"
date: "last updated 5/7/24" 
---
## Load required packages
library(Seurat)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 8000 * 1024^2)
library(SeuratData)
library(SeuratWrappers)
library(plyr)
library(dplyr)
library(ggplot2)
library(sctransform)
library(ggrepel)
library(umap)
library(knitr)
library(cowplot)
library(glmGamPoi)
library(viridis)

## Lists of Interesting Genes
Kalluri.markers <- c('MYH11','TPM2','MYL9','ACTA2','TAGLN',
                     'VWA1','PRNP','CNP','GPM6B','MBP',
                     'PF4','C1QB','RETNLA','C1QA','LYZ2',
                     'BTG1','LTB','CORO1A','RAC2','CD52',
                     'SERPINF1','LUM','CLEC3B','GSN','DCN',
                     'GPIHBP1','PECAM1','CCL21A','CYTL1','FABP4') #cell type markers from published scRNAseq dataset of mouse aorta: https://pubmed.ncbi.nlm.nih.gov/31146585/ 

canonical.markers <- c('MYH11','TPM2','MYL9','ACTA2','TAGLN','CNN1','SMTN','VWA1','PRNP','CNP','GPM6B','MBP',
                       'C1QB','C1QA','PTPRC','CD52',
                       'SERPINF1','LUM','CLEC3B',
                       'PECAM1','GPIHBP1','FABP4') #aortic cell type markers, adapted from Kalluri markers 

## Load single cell data set from appropriate file path using Read10x function (this should be a non-normalized matrix)
AllSamples_data <-Read10X(data.dir="~/Desktop/HumanLDS/filtered_feature_bc_matrix/")

## Initialize the Seurat object with the raw (non-normalized data)
AllSamples<- CreateSeuratObject(counts = AllSamples_data, project = "humanLDS", min.cells = 3, min.features = 200)

## Create percent mitochondirial derived reads as metadata, visualize QC metrics as Violin plot  
AllSamples[["percent.MT"]] <- PercentageFeatureSet(AllSamples, pattern = "^MT-")
# Save plots as objects
AllSamples_RNAQC_plot1 <- VlnPlot(AllSamples, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
AllSamples_plot1 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "percent.MT")
AllSamples_plot2 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
AllSamples_RNAQC_plot2 <- AllSamples_plot1 + AllSamples_plot2
# Visualize plots to analyze
AllSamples_RNAQC_plot1
AllSamples_RNAQC_plot2

## Filtering cells
length(colnames(AllSamples))
length(colnames(subset(AllSamples, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 30000 & percent.MT < 20)))
# Filter All cells 
AllSamples <- subset(AllSamples, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 30000  & percent.MT < 20)
# Create RNAQCplots, post filtering
AllSamples_filteredRNAQC_plot1 <- VlnPlot(AllSamples, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
AllSamples_filteredplot1 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "percent.MT")
AllSamples_filteredplot2 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
AllSamples_filteredRNAQC_plot2 <- AllSamples_filteredplot1+ AllSamples_filteredplot2
# Visualize plots
AllSamples_filteredRNAQC_plot1
AllSamples_filteredRNAQC_plot2

##Add metadata for each mouse ID, CRE, genotype, and batch within the dataset 
length(grep("*-1",colnames(AllSamples)))
length(grep("*-2",colnames(AllSamples)))
length(grep("*-3",colnames(AllSamples)))
length(grep("*-4",colnames(AllSamples)))
length(grep("*-5",colnames(AllSamples)))
length(grep("*-6",colnames(AllSamples)))
length(grep("*-7",colnames(AllSamples)))
length(grep("*-8",colnames(AllSamples)))
length(grep("*-9",colnames(AllSamples)))
length(grep("*-10",colnames(AllSamples)))
length(grep("*-11",colnames(AllSamples)))
length(grep("*-12",colnames(AllSamples)))
length(grep("*-13",colnames(AllSamples)))
length(grep("*-14",colnames(AllSamples)))
length(grep("*-15",colnames(AllSamples)))
length(grep("*-16",colnames(AllSamples)))

#sampleID
AllSamples <- AddMetaData(object = AllSamples, col.name = "sampleID", metadata = 
                            c(rep(x = "g140root",length(grep("*-1$",colnames(AllSamples)))),
                              rep("g140asc",length(grep("*-2",colnames(AllSamples)))),
                              rep("g129root",length(grep("*-3",colnames(AllSamples)))),
                              rep("g129asc",length(grep("*-4",colnames(AllSamples)))),
                              rep("donor1root",length(grep("*-5",colnames(AllSamples)))), 
                              rep("donor1asc",length(grep("*-6",colnames(AllSamples)))),
                              rep("g150root",length(grep("*-7",colnames(AllSamples)))),
                              rep("g150asc",length(grep("*-8",colnames(AllSamples)))),
                              rep("g163root",length(grep("*-9",colnames(AllSamples)))),
                              rep("g163asc",length(grep("*-10",colnames(AllSamples)))),
                              rep("donor2root",length(grep("*-11",colnames(AllSamples)))),
                              rep("donor2asc",length(grep("*-12",colnames(AllSamples)))),
                              rep("g165root",length(grep("*-13",colnames(AllSamples)))),
                              rep("g165asc",length(grep("*-14",colnames(AllSamples)))),
                              rep("m11root",length(grep("*-15",colnames(AllSamples)))),
                              rep("m11asc",length(grep("*-16",colnames(AllSamples))))))
head(AllSamples$sampleID)
tail(AllSamples$sampleID)

#aorticregion
AllSamples <- AddMetaData(object = AllSamples, col.name = "aorticregion", metadata = 
                            c(rep(x = "root",length(grep("*-1$",colnames(AllSamples)))),
                              rep("asc",length(grep("*-2",colnames(AllSamples)))),
                              rep("root",length(grep("*-3",colnames(AllSamples)))),
                              rep("asc",length(grep("*-4",colnames(AllSamples)))),
                              rep("root",length(grep("*-5",colnames(AllSamples)))), 
                              rep("asc",length(grep("*-6",colnames(AllSamples)))),
                              rep("root",length(grep("*-7",colnames(AllSamples)))),
                              rep("asc",length(grep("*-8",colnames(AllSamples)))),
                              rep("root",length(grep("*-9",colnames(AllSamples)))),
                              rep("asc",length(grep("*-10",colnames(AllSamples)))),
                              rep("root",length(grep("*-11",colnames(AllSamples)))),
                              rep("asc",length(grep("*-12",colnames(AllSamples)))),
                              rep("root",length(grep("*-13",colnames(AllSamples)))),
                              rep("asc",length(grep("*-14",colnames(AllSamples)))),
                              rep("root",length(grep("*-15",colnames(AllSamples)))),
                              rep("asc",length(grep("*-16",colnames(AllSamples))))))
head(AllSamples$aorticregion)
tail(AllSamples$aorticregion)

#geno 
AllSamples <- AddMetaData(object = AllSamples, col.name = "geno", metadata = 
                            c(rep(x = "LDSroot",length(grep("*-1$",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-2",colnames(AllSamples)))),
                              rep("LDSroot",length(grep("*-3",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-4",colnames(AllSamples)))),
                              rep("CTLroot",length(grep("*-5",colnames(AllSamples)))), 
                              rep("CTLasc",length(grep("*-6",colnames(AllSamples)))),
                              rep("LDSroot",length(grep("*-7",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-8",colnames(AllSamples)))),
                              rep("LDSroot",length(grep("*-9",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-10",colnames(AllSamples)))),
                              rep("CTLroot",length(grep("*-11",colnames(AllSamples)))),
                              rep("CTLasc",length(grep("*-12",colnames(AllSamples)))),
                              rep("LDSroot",length(grep("*-13",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-14",colnames(AllSamples)))),
                              rep("LDSroot",length(grep("*-15",colnames(AllSamples)))),
                              rep("LDSasc",length(grep("*-16",colnames(AllSamples))))))
head(AllSamples$geno)
tail(AllSamples$geno)

#genoregioncombined
AllSamples <- AddMetaData(object = AllSamples, col.name = "genoregioncombined", metadata = 
                            c(rep(x = "LDS",length(grep("*-1$",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-2",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-3",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-4",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-5",colnames(AllSamples)))), 
                              rep("CTL",length(grep("*-6",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-7",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-8",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-9",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-10",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-11",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-12",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-13",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-14",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-15",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-16",colnames(AllSamples))))))


#Normalize data using SCTransform
AllSamples.sct<- SCTransform(AllSamples, verbose = FALSE)

#Run non-linear dimensional reduction 
aorta.combined.sct<- AllSamples.sct
aorta.combined.sct <- RunPCA(aorta.combined.sct, verbose = FALSE)
aorta.combined.sct <- RunUMAP(aorta.combined.sct, reduction = "pca", dims = 1:30)
aorta.combined.sct <- FindNeighbors(aorta.combined.sct, dims = 1:30, verbose = FALSE)
aorta.combined.sct <- FindClusters(aorta.combined.sct, verbose = FALSE, resolution = 0.10) #resolution used was 0.35 
DimPlot(aorta.combined.sct, label = TRUE) + NoLegend()

#Also normalize and scale just the RNA assay (to use for differential expression analysis)
aorta.combined.sct <-NormalizeData(aorta.combined.sct, assay='RNA') 
aorta.combined.sct <-ScaleData(aorta.combined.sct, assay='RNA') 

#Visualize UMAPs of data set split by metadata 
DimPlot(aorta.combined.sct, reduction = "umap") 
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "aorticregion")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "aorticregion")
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "geno")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "geno")
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "sampleID")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "sampleID")

# Identify major cell types using canonical markers/those validated by earlier aortic scRNAseq expereiments 
DotPlot(aorta.combined.sct, features = canonical.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(7.5,12.5,16.5,19.5))
DotPlot(aorta.combined.sct, features = Kalluri.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))

## Name clusters ***this needs to be fixed cluster number is slightly differnt due to updates*** 8/4/22
aorta.combined.sct_named <- RenameIdents(aorta.combined.sct, `0` = "VSMC1", `1` = "Fibroblasts", `2` = "Endothelial", `3` = "Immune", `4` = "VSMC2", `5` = "Immune", `6` = "Endothelial", `7` = "Neuron")
aorta.combined.sct_named2 <- RenameIdents(aorta.combined.sct, `0` = "VSMCs", `1` = "Fibroblasts", `2` = "Endothelial", `3` = "Immune", `4` = "VSMCs", `5` = "Immune", `6` = "Endothelial", `7` = "Neuron")

DimPlot(aorta.combined.sct_named, reduction = "umap", split.by = "sampleID")
DimPlot(aorta.combined.sct_named, reduction = "umap", split.by = "geno")
DimPlot(aorta.combined.sct_named, reduction = "umap", split.by = "genoregioncombined")
genoregioncombined

#Subset samples from each aortic region 
#aorta.combined.sct_named_Rootonly <- subset(aorta.combined.sct_named, aorticregion=="root")
#aorta.combined.sct_named_Asconly <- subset(aorta.combined.sct_named, aorticregion=="asc")

##### Analysis of data collapsed into major clusters

#create a version of the file where idents are split by genotype for downstream analysis 
aorta.combined.sct.named.splitbygeno <- aorta.combined.sct_named
aorta.combined.sct.named.splitbygeno$geno <- paste(Idents(aorta.combined.sct.named.splitbygeno), aorta.combined.sct.named.splitbygeno$geno, sep = "_")
aorta.combined.sct.named.splitbygeno$celltype <- Idents(aorta.combined.sct.named.splitbygeno)
Idents(aorta.combined.sct.named.splitbygeno) <- "geno"

#UMAP split by relevant metadata saved as objects for PDF export 
p1<- DimPlot(aorta.combined.sct_named, label = TRUE) 
p2<- DimPlot(aorta.combined.sct_named, label = TRUE, split.by = "geno") 
p3<- DimPlot(aorta.combined.sct_named, reduction = "umap", split.by = "sampleID")

FeaturePlot(aorta.combined.sct_named, features = canonical.markers)

##### Analysis of VSMCs major subclusters (not re-clustered after subsetting)
#subset 
aorta.combined.sct.VSMCsonly <- subset(aorta.combined.sct_named , idents = c("VSMC1","VSMC2"))
#aorta.combined.sct.VSMCsonly <- subset(aorta.combined.sct_named2 , idents = c("VSMCs")) #use if you want all VSMCs combined 
#create a version of the file where idents are split by genotype for downstream analysis 

aorta.combined.sct.VSMCsonly.splitbygeno <- aorta.combined.sct.VSMCsonly
aorta.combined.sct.VSMCsonly.splitbygeno$geno <- paste(Idents(aorta.combined.sct.VSMCsonly.splitbygeno), aorta.combined.sct.VSMCsonly.splitbygeno$geno, sep = "_")
aorta.combined.sct.VSMCsonly.splitbygeno$celltype <- Idents(aorta.combined.sct.VSMCsonly.splitbygeno)
Idents(aorta.combined.sct.VSMCsonly.splitbygeno) <- "geno"


aorta.combined.sct.VSMCsonly.splitbygeno <- aorta.combined.sct.VSMCsonly
aorta.combined.sct.VSMCsonly.splitbygeno$genoregioncombined <- paste(Idents(aorta.combined.sct.VSMCsonly.splitbygeno), aorta.combined.sct.VSMCsonly.splitbygeno$genoregioncombined, sep = "_")
aorta.combined.sct.VSMCsonly.splitbygeno$celltype <- Idents(aorta.combined.sct.VSMCsonly.splitbygeno)
Idents(aorta.combined.sct.VSMCsonly.splitbygeno) <- "genoregioncombined"

aorta.combined.sct.VSMCsonly.splitbyaorticregion <- aorta.combined.sct.VSMCsonly
aorta.combined.sct.VSMCsonly.splitbyaorticregion$aorticregion <- paste(Idents(aorta.combined.sct.VSMCsonly.splitbyaorticregion), aorta.combined.sct.VSMCsonly.splitbyaorticregion$aorticregion, sep = "_")
aorta.combined.sct.VSMCsonly.splitbyaorticregion$celltype <- Idents(aorta.combined.sct.VSMCsonly.splitbyaorticregion)
Idents(aorta.combined.sct.VSMCsonly.splitbyaorticregion) <- "aorticregion"

#Visualize cell type defining markers in the collapsed clusters 
r5<-DotPlot(aorta.combined.sct.VSMCsonly, assay="RNA", features = canonical.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))
r5.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = canonical.markers)
r6<-DotPlot(aorta.combined.sct.VSMCsonly, assay="RNA", features = RegionSpecific, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))
r6.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = RegionSpecific)
r6.2<- DotPlot(aorta.combined.sct.VSMCsonly,  assay="RNA", features = c('PTPRZ1','POSTN','TES','GATA4','ENPEP','NOTCH3','RGS5'), dot.scale = 6,scale.min=0) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(3.5,6.5))
r9<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Gata4", "Ptprz1"), min.cutoff=0,max.cutoff=10)
r10<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("GATA4", "PTPRZ1"),min.cutoff=0,max.cutoff=5, blend = TRUE, cols = c("lightgrey", "red", "cyan3")) & xlim(-5,10) &ylim(0,10)
r10.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("ENPEP", "TES"),min.cutoff=0,max.cutoff=5, blend = TRUE, cols = c("lightgrey", "red", "cyan3")) & xlim(-5,10) &ylim(0,10)
r10.2<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Ptprz1","Postn"),min.cutoff=0,max.cutoff=5, blend = TRUE,cols = c("lightgrey", "red", "cyan3")) & xlim(-5,10) &ylim(0,10)
r11<- DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = c('TGFB1','TGFB2','TGFB3','CTGF','SERPINE1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 100) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS","VSMC3_CTL","VSMC3_LDS" )) + geom_hline(yintercept = c(seq(2.5,20.5,2)))
r12<- DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = c('Col3a1','Col1a1','Col1a2','Mfap4','Mfap5','Eln','Rock1','Fbln1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 100) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS","VSMC3_CTL","VSMC3_LDS" )) + geom_hline(yintercept = c(seq(2.5,20.5,2)))
leda <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('Meg3','Acta2','Bax','F2r','Fn1','Gata3','Serpine1', 'Tgfbr1', 'Id2', 'Id4', 'Egr1'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis()
leda <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('Meg3','Tgfbr1','Acta2','Bax','F2r','Fn1','Gata3','Serpine1','Id2','Id4','Egr1','Ppp1r15a','Cebpb','Btg2','Zfp36','Sgk1','Ier5','Atf3','Ier3'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis()
meg3<- DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('MEG3'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 50) + RotatedAxis()
DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = c('GATA4',), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC2_CTLasc","VSMC2_LDSasc","VSMC2_CTLroot","VSMC2_LDSroot" ))
DotPlot(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay="RNA", features = c('GATA4'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_asc","VSMC1_root","VSMC2_asc","VSMC2_root"))
Gata4regulated <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('GATA4','JUN','PPP1R15A','CEBPB','BTG2','PIM1','HES1','FOS','IER5','KLF2','IER2'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_LDS","VSMC2_LDS"))
Autophagy <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('CDKN1A', 'CXCL1', 'GATA4', 'THBS1', 'BMP2', 'IGFBP5','IRF1'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_LDS","VSMC2_LDS"))


#plots for paper
#vsmc1 blue, vsmc2red 

q1<- DimPlot(aorta.combined.sct_named, label = TRUE)
q2<- DimPlot(aorta.combined.sct_named, label = TRUE, split.by = "geno") 
r1<- DimPlot(aorta.combined.sct.VSMCsonly, label = FALSE, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))
p1<-DotPlot(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay="RNA", features = c('GATA4','ENPEP','NOTCH3','TES','PTPRZ1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_asc","VSMC1_root","VSMC2_asc","VSMC2_root"))
p2<-DotPlot(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay="RNA", features = c('GATA4','ENPEP','NOTCH3','TES','PTPRZ1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMCs_asc","VSMCs_root"))
r2<- DimPlot(aorta.combined.sct.VSMCsonly, label = TRUE, split.by = "genoregioncombined") 
r6.2<- DotPlot(aorta.combined.sct.VSMCsonly,  assay="RNA", features = c('GATA4','ENPEP','NOTCH3','TES','PTPRZ1'), dot.scale = 6,scale.min=0,scale.max = 50) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(3.5,6.5))
DotPlot(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay="RNA", features = c('GATA4','ENPEP','NOTCH3','TES','PTPRZ1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_asc","VSMC1_root","VSMC2_asc","VSMC2_root"))


ENRICHR<-c('BTG2', 'CSF1', 'CEBPD','CEBPB', 'TNC', 'SLC2A3', 'SOCS3', 'ZFP36', 'NFIL3', 'PLAU', 'MYC', 'CCL2', 'PHLDA2', 'JUNB', 'PHLDA1', 'EGR1', 'JUN', 'GADD45B', 'FOS', 'NR4A2', 'NR4A1', 'IL6', 'BMP2', 'NR4A3', 'FOSB', 'ATF3','RGS4', 'PRRX1', 'LUM', 'IGFBP4', 'CXCL1', 'PCOLCE', 'MSX1', 'THBS1', 'DCN', 'IL32','LAMA3', 'THY1', 'CDH6', 'CXCL12', 'COL5A3', 'SNAI2', 'TIMP3', 'SLIT3', 'CDKN1A', 'GSN', 'IRF1', 'HMGB2', 'KRT18', 'HGF', 'PLPPR4', 'AVPR1A', 'PLK2', 'HES1', 'SRPX', 'MT2A', 'ANGPTL4', 'PGF', 'MT1E', 'LXN', 'TGFB3', 'CITED2', 'PGAM2', 'NDRG1','CP' )
ENRICHRsubset<-c('CEBPD','CEBPB', 'MYC', 'EGR1','CDKN1A', 'IRF1', 'JUN','FOS','IL6', 'BMP2', 'MSX1', 'THBS1', 'DCN' , 'TGFB3', 'CITED2')


DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = ENRICHRsubset, dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))

	
Subset<-c('MYH11','CNN1','TET2','KLF4','OLFM2','SOX9','TCF21','MALAT1','TWIST1','CEBPB','DCN')
Fig2H <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = Subset, dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0, scale.max=20) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))

#Generate tables of cell number per cluster, and how these numbers are split by genotype (are there changes in relative proportions of cells in disease?)
r7<- table(Idents(aorta.combined.sct.VSMCsonly))
r8<- table(Idents(aorta.combined.sct.VSMCsonly), aorta.combined.sct.VSMCsonly$genoregioncombined)


VSMC1.LDSRootvsDonorRoot.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDSroot", ident.2 = "VSMC1_CTLroot", verbose = TRUE, logfc.threshold = 0.00)
VSMC2.LDSRootvsDonorRoot.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC2_LDSroot", ident.2 = "VSMC2_CTLroot", verbose = TRUE,  logfc.threshold = 0.00)
VSMC1.LDSAscvsDonorAsc.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDSasc", ident.2 = "VSMC1_CTLasc", verbose = TRUE, logfc.threshold = 0.00)
VSMC2.LDSAscvsDonorAsc.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC2_LDSasc", ident.2 = "VSMC2_CTLasc", verbose = TRUE,  logfc.threshold = 0.00)
VSMC1.LDSvsVSMC2.LDS.ROOT.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDSroot", ident.2 = "VSMC2_LDSroot", verbose = FALSE, logfc.threshold = 0.00)
VSMC1.DonorvsVSMC2.Donor.ROOT.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_CTLroot", ident.2 = "VSMC2_CTLroot", verbose = FALSE, logfc.threshold = 0.00)

VSMC1.LDSvsVSMC2.LDS.ASC.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDSasc", ident.2 = "VSMC2_LDSasc", verbose = FALSE, logfc.threshold = 0.00)
VSMC1.DonorvsVSMC2.Donor.ASC.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_CTLasc", ident.2 = "VSMC2_CTLasc", verbose = FALSE, logfc.threshold = 0.00)

VSMC1.LDSvsVSMC2.LDS.human.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDS", ident.2 = "VSMC2_LDS", verbose = FALSE, logfc.threshold = 0.00)


write.csv(VSMC1.LDSRootvsDonorRoot.diffexp, file = '~/Desktop/VSMC1.LDSRootvsDonorRoot.diffexp_RNAassay.csv') 
write.csv(VSMC2.LDSRootvsDonorRoot.diffexp, file = '~/Desktop/VSMC2.LDSRootvsDonorRoot.diffexp_RNAassay.csv') 
write.csv(VSMC1.LDSAscvsDonorAsc.diffexp, file = '~/Desktop/VSMC1.LDSAscvsDonorAsc.diffexp_RNAassay.csv') 
write.csv(VSMC2.LDSAscvsDonorAsc.diffexp, file = '~/Desktop/VSMC2.LDSAscvsDonorAsc.diffexp_RNAassay.csv') 
write.csv(VSMC1.LDSvsVSMC2.LDS.ROOT.diffexp , file = '~/Desktop/VSMC1.LDSvsVSMC2.LDS.ROOT.diffexp.csv') 
write.csv(VSMC1.DonorvsVSMC2.Donor.ROOT.diffexp , file = '~/Desktop/VSMC1.DonorvsVSMC2.Donor.ROOT.diffexp.csv') 
write.csv(VSMC1.LDSvsVSMC2.LDS.ASC.diffexp , file = '~/Desktop/VSMC1.LDSvsVSMC2.LDS.ASC.diffexp.csv') 
write.csv(VSMC1.DonorvsVSMC2.Donor.ASC.diffexp , file = '~/Desktop/VSMC1.DonorvsVSMC2.Donor.ASC.diffexp.csv') 
write.csv(VSMC1.LDSvsVSMC2.LDS.human.diffexp , file = '~/Desktop/VSMC1.LDSvsVSMC2.LDS.human.diffexp.csv') 

aorta.combined.sct.VSMCsonly.splitbyaorticregion
r7<- table(Idents(aorta.combined.sct.VSMCsonly))
r8<- table(Idents(aorta.combined.sct.VSMCsonly), aorta.combined.sct.VSMCsonly$aorticregion)
VSMC1.RootvsAsc.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay='RNA', ident.1 = "VSMC1_root", ident.2 = "VSMC1_asc", verbose = TRUE, logfc.threshold = 0.00)
VSMC2.RootvsAsc.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay='RNA', ident.1 = "VSMC2_root", ident.2 = "VSMC2_asc", verbose = TRUE,  logfc.threshold = 0.00)

VSMC1.RootvsVSMC2.Root.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay='RNA', ident.1 = "VSMC1_root", ident.2 = "VSMC2_root", verbose = FALSE, logfc.threshold = 0.00)
VSMC1.AscvsVSMC2.Asc.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbyaorticregion,assay='RNA', ident.1 = "VSMC1_asc", ident.2 = "VSMC2_asc", verbose = FALSE, logfc.threshold = 0.00)


write.csv(VSMC1.RootvsAsc.diffexp, file = '~/Desktop/VSMC1.RootvsAsc.diffexp.csv') 
write.csv(VSMC2.RootvsAsc.diffexp, file = '~/Desktop/VSMC2.RootvsAsc.diffexp.csv') 
write.csv(VSMC1.RootvsVSMC2.Root.diffexp, file = '~/Desktop/VSMC1.RootvsVSMC2.Root.diffexp.csv') 
write.csv(VSMC1.AscvsVSMC2.Asc.diffexp, file = '~/Desktop/VSMC1.AscvsVSMC2.Asc.diffexp.csv') 


pdf(file = '~/Desktop/HumanscRNAseq_forFig2.pdf')  
r1
p1 
p2
dev.off()



#projectR 

#replace human genes in the Seurat object's RNA assay to mouse gene names 
test<- aorta.combined.sct.VSMCsonly

test@assays$SCT@counts@Dimnames[[1]] <-genes$MusGenes
test@assays$SCT@data@Dimnames[[1]] <- genes$MusGenes
test@assays$SCT@meta.features <- genes

scaled_mat <- test@assays$SCT@scale.data
head(rownames(scaled_mat))
dim(genes)


#updated with Jacob T. Mitchell, converting human gene names to mouse so they match the patterns 

test<- aorta.combined.sct.VSMCsonly
library(biomaRt)
huGenes <- rownames(aorta.combined.sct.VSMCsonly)
useEnsembl(biomart = "ensembl")
# Basic function to convert mouse to human gene names
convertHumantoMouseGeneList <- function(x){
  require("biomaRt")
  m = useMart(host="https://dec2021.archive.ensembl.org/", 
              biomart="ENSEMBL_MART_ENSEMBL", 
              dataset="mmusculus_gene_ensembl")
  h = useMart(host="https://dec2021.archive.ensembl.org/", 
              biomart="ENSEMBL_MART_ENSEMBL",  dataset = "hsapiens_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                   values = x , mart = h, attributesL = c("mgi_symbol"),
                   martL = m, uniqueRows=T)
  
  mouse_ortho <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                        values = x,
                        mart = h,
                        attributesL = c("mgi_symbol"), martL = m,
                        uniqueRows=T)
  return(mouse_ortho)
}
mouse_ortho <- convertHumantoMouseGeneList(huGenes)

# pull out human gene names from scale data
scaled_mat <- aorta.combined.sct.VSMCsonly@assays$SCT@scale.data
mat_hu_genes <- rownames(scaled_mat)

mat_mus_genes <- plyr::mapvalues(x = mat_hu_genes,
                                 from = mouse_ortho$HGNC.symbol,
                                 to = mouse_ortho$MGI.symbol)
head(mat_hu_genes)
head(mat_mus_genes)

# convert row names to mouse orthologs
rownames(scaled_mat) <- mat_mus_genes
library(CoGAPS)
library(projectR)

result <- readRDS("~/Downloads/CoGAPS_8patterns_16wkLDSmice_result.rds")
coGAPsMatrix<-result@featureLoadings #to turn coGAPs result into a usable matrix 

projectR_result_t <-projectR(data= scaled_mat, loadings= coGAPsMatrix)


pattern_weights <- as.data.frame(t(projectR_result_t))

for(pattern in colnames(pattern_weights)) {
  test <- AddMetaData(object = test,metadata = pattern_weights[[pattern]],
                      col.name = pattern)
}


DimPlot(samples.combined.sct.VSMCsonly, reduction = "umap",label = "True") & xlim(-6.5,5) &ylim(-10,0)
DotPlot(samples.combined.sct.VSMCsonly, features= c('Gata4',"Notch3",'Enpep','Ptprz1','Tes'), dot.scale = 6, scale.min = 0) + RotatedAxis()

VlnPlot(test , features = c("Pattern_4"),pt.size = 0.05, alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))
VlnPlot(test , features = c("Pattern_5"),pt.size = 0.05,alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))
FeaturePlot(test, c("Pattern_1","Pattern_2",'Pattern_3',"Pattern_4","Pattern_5","Pattern_6", "Pattern_7","Pattern_8"), pt.size = 0.1, min.cutoff=0, max.cutoff=2, label=TRUE) & scale_color_viridis_c(option="plasma") & xlim(-10,15) &ylim(-8,15)
FeaturePlot(test, c("Pattern_4","Pattern_5"), pt.size = 0.1, min.cutoff=0, max.cutoff=2) & scale_color_viridis_c(option="plasma") & xlim(-10,15) &ylim(-8,15)
FeaturePlot(test, c("Pattern_5"), pt.size = 0.1, label=TRUE) & scale_color_viridis_c(option="magma") & xlim(-10,15) &ylim(-8,15)


test$cluster <- Idents(aorta.combined.sct.VSMCsonly)
vsmc1_pat4 <- test@meta.data[test@meta.data$cluster == "VSMC1", "Pattern_4"]
vsmc2_pat4 <- test@meta.data[test@meta.data$cluster == "VSMC2", "Pattern_4"]
wilcox.test(x = vsmc1_pat4, y = vsmc2_pat4)

vsmc1_pat5 <- test@meta.data[test@meta.data$cluster == "VSMC1", "Pattern_5"]
vsmc2_pat5 <- test@meta.data[test@meta.data$cluster == "VSMC2", "Pattern_5"]
wilcox.test(x = vsmc1_pat5, y = vsmc2_pat5)



#####testing pattern list 
pm<- patternMarkers(result)
write.csv(pm, file = '../macfarlanelab/Desktop/For GATA4 paper/PatternMarkers.csv')

pm4<- (pm$PatternMarkers[4])
pm5<- (pm$PatternMarkers[5])
write.csv(pm4, file = '../macfarlanelab/Desktop/For GATA4 paper/PatternMarkers4.csv')
write.csv(pm5, file = '../macfarlanelab/Desktop/For GATA4 paper/PatternMarkers5.csv')

