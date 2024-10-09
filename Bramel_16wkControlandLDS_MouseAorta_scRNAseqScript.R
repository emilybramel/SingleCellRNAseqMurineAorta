---
title: "16wk Tgfbr1^M318R/+ vs WT Seurat scRNAseq Analysis (includes RPCA batch integration)"
#adapted from https://satijalab.org/seurat/ 
author: "Emily Bramel"
date last updated: "7/24/24" 
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
Kalluri.markers <- c('Myh11','Tpm2','Myl9','Acta2','Tagln',
                     'Vwa1','Prnp','Cnp','Gpm6b','Mbp',
                     'Pf4','C1qb','Retnla','C1qa','Lyz2',
                     'Btg1','Ltb','Coro1a','Rac2','Cd52',
                     'Serpinf1','Lum','Clec3b','Gsn','Dcn',
                     'Gpihbp1','Pecam1','Ccl21a','Cytl1','Fabp4') #cell type markers from published scRNAseq dataset of mouse aorta: https://pubmed.ncbi.nlm.nih.gov/31146585/ 

canonical.markers <- c('Myh11','Tpm2','Myl9','Acta2','Tagln','Cnn1','Smtn',
                       'C1qb','C1qa','Ptprc','Cd52',
                       'Serpinf1','Lum','Clec3b',
                       'Pecam1','Gpihbp1','Fabp4') #aortic cell type markers, adapted from Kalluri markers 

FApanel<-c('Fblim1', 'Svil','Tln1','Macf1', 'Trip6','Zrp1','Tgfb1i1')
NotchPathway<-c('Jag1','Heyl','Hes1','Nov')
WntSignaling<-c('Wif1','Dact3','Sfrp2','Nxn','Lgr','Sost','Wisp2','Wnt16')
FAandCytoskeleton<-c('Syne1','Fblim1','Filip1l','Svil','Synpo2','Zyx','Flna','Tln1','Macf1','Pdlim1','Pdlim2','Pdlim5','Pfn1','Palld','Tgfb1i1','Msn')
AngII_TGFBsignaling<-c('Ctgf','Serpine1','Serpinh1','Dcn','Tgfb2','Inhba','Pmepa1','Skil','Nr4a1','Ace','Thbs1','Ppp1r14a','Ppp1r15a','Ppp1r12b','Rgs5','Nox4','Ltbp1','Ltbp2','Ltbp3','Ltbp4','Postn','Map3k4')
MatrixProteins<- c('Mfap4','Mfap5','Fn1','Eln','Bgn','Aebp1','Fbln1','Fbln5','Col8a1','Col1a1','Col5a2','Col3a1','Col6a2','Col15a1','Col4a6','Loxl1','Lox','Col1a2')

SHF<- c('Nkx2','Gata4','Isl1','Tbx1','Foxc1','Foxc2','Fgf8','Fgf10','Foxh1')
CNC<- c('Tfap2','Msx1','Msx2','Dlx3','Dlx5','Pax3','Pax7','Zic1','Snai2','Foxd3','Sox9','Sox10')
RegionSpecific <- c('Gata4', 'Enpep','Notch3','Tes','Ptprz1', 'Rgs5')

Focaladhesioncomplexcomponents <- c('Fblim1','Svil','Tes','Tln1','Vcl','Zyx','Filip1l','Flna','Lims2','Nedd9','Pdlim1','Pdlim3','Pdlim4','Pdlim5','Pdlim7','Synpo2','Actn1','Actn2','Actn4','Abra','Rnd1','Rnd3','Rock1','Rasgrf2','Pkp4','Itsn1','Itga1','Itga4','Itga5','Itga10','Itga6','Itga7','Itga8','Itga9','Itgb1','Itgb5','Itgb7','Itgb8','Itgbl1','Enah','Actc1','Acta2','Actb','Actg1','Actg2','Actr3','Gsn','Dstn','Pfn1','Palld','Vasp')
Focaladhesioncomplexcomponents1 <- c('Fblim1','Svil','Tes','Tln1','Vcl','Zyx','Filip1l','Flna','Pdlim3','Rnd1','Rnd3','Rock1','Pkp4','Itga1','Itga4','Itga5','Itga8','Itga9','Itgb1','Itgb5','Enah','Acta2','Actg1','Actg2')
Matrixcomponents <- c('Lox','Loxl1','Loxl3','Postn','Col1a1','Col1a2','Col3a1','Col5a2','Adam15','Adamts1','Adamts15','Adamts9','Htra1','Mmp14','Mmp2','Acan','Bgn','Vcan','Cspg4','Sdc4','Sdc1','Sdc3','Fbn1','Eln','Mfap4','Fbln5','Npnt','Fn1','Thbs1','Aebp1','Ecm1','Comp','Igfbp2','Spp1')
Matrixcomponents1 <- c('Lox','Loxl1','Loxl3','Postn','Col1a1','Col1a2','Col3a1','Adam15','Adamts1','Adamts9','Htra1','Mmp2','Acan','Bgn','Vcan','Cspg4','Sdc4','Sdc1','Sdc3','Fbn1','Eln','Mfap4','Fbln5','Npnt','Fn1','Thbs1','Aebp1','Ecm1','Igfbp2')
CellSignaling<- c('Acvrl1','Bmp2','Ctgf','Inhba','Ltbp1','Ltbp2','Ltbp3','Ltbp4','Serpine1','Pmepa1','Nr4a1','Scube3','Tgfb1','Tgfb1i1','Tgfb2','Tgfb3','Tgfbr3','Smad6','Smad7','Vasn','Jag1','Jag2','Notch3','Notch4','Hes1','Heyl','Nov','Fgf1','Fgf2','Fos','Hgf','Id3','Pdgfa','Pdgfrb','Nrp1','Ldb2','Cxcl2','Edn1','Ednra','Ednrb','Enpp2','Yes1','Jund','Jun','Junb','Myc','Glul','Gria1','Lmcd1','Rcan1')
CellSignaling1<- c('Ctgf','Inhba','Ltbp1','Ltbp2','Ltbp3','Ltbp4','Serpine1','Pmepa1','Nr4a1','Tgfb1','Tgfb2','Tgfb3','Tgfbr3','Smad7','Tgfb1i1','Vasn','Jag1','Jag2','Notch3','Hes1','Heyl','Nov','Fgf1','Id3','Pdgfa','Pdgfrb','Fos','Jund','Jun','Junb','Myc','Rcan1')
Epigeneticregulation <- c('Hdac7','Cited4','Jarid2','Prdm6','Syne1','Sap30','Fhl2','Sirt1','Csrp2','Smtn','Tagln2','Myh11','Mylk','Myocd','Mustn1','Cnn1','Vim','Tcf21','Fbxo30','Hlf','Fmo3','Maob','Sod3','Nox4','Lrp1')
Epigeneticregulation1 <-c('Prdm6','Syne1','Sap30','Fhl2','Csrp2','Smtn','Myh11','Mylk','Myocd','Cnn1','Sod3','Nox4','Lrp1')

Anchoring<-c('Fblim1','Filip1l','Homer1','Lims2','Nedd9','Pdlim1','Pdlim5','Pdlim7','Sh3pxd2a','Sntb1','Svil','Tln1','Vcl','Sgcg','Tgfb1i1','Zyx')
Actin<- c('Actc1','Actg2','Actn1','Actn4','Arc','Flna','Wdr1','Actg1','Actr3','Cfl2','Enah','Pdlim3','Pfn1','Synpo2','Fbxl22','Gsn')
Rho<- c('Arhgap22','Arhgap4','Dock4','Pak1','Pkd1','Rock1','Net1','Abra','Elmo1','Eps8','Pak7','Pkp4','Prex1','Rnd2','Rras','Tec','Rnd1','Rnd3','Nuak1')
Integrins<- c('Itga1','Itga7','Itga8','Itga9','Itga5','Itgb1','Hmmr','Itga10','Itga6','Itgb5','Sema5a','Itgbl1')
Collagen<- c('Col14a1','Col1a2','Col27a1','Col4a1','Col4a2','Col4a5','Col4a6','Col5a1','Col6a1','Col6a2','Col15a1','Col18a1','Col23a1','Col6a3','Loxl1','Loxl3','Lox')
Proteases<- c('Adamts10','Adamtsl1','Htra1','Mmp15','Mmp17','Mmp23','Timp3','Wfikkn2','Adamts8','Adamts9','Adamts2','Fgl2')
Proteoglycans<- c('Acan','Bgn','Vcan','Cspg4','Sdc1','Hspg2','Gpc3','Sdc4')
ElasticFiber<- c('Eln','Fbln5','Mfap4','Fbn1','Efemp2')
IntegrinLigands<- c('Thbs1','Lgals3','Npnt','Fn1','Tnc','Lgals3bp')
OtherMatrix<- c('Aebp1','Postn','Sparc','Sulf1','Ecm1','Igfbp2')
TGFB<- c('Ltbp1','Ltbp2','Ltbp3','Ltbp4','Tgfb1','Tgfb2','Tgfb3','Pmepa1','Skil','Nr4a1','Nrep','Wwp2','Smad7','Tgfbr3','Acvrl1','Ctgf','Fbxo2','Scube3','Aspn','Bmp3','Gdf6','Aurka','Inhba','Bmp1','Vasn','Bambi','Cav1','Cav2','Ptrf')
Wnt <- c('Amotl2','Ccdc88c','Frat2','Frzb','Fzd8','Jade1','Wisp1','Wnt9a','Sost','Dkk3','Fzd1','Caprin2','Lgr5','Lgr6','Mdfic','Rspo1','Sfrp2','Tle4','Wnt16')
Notch<- c('Notch3','Hey1','Hey2','Heyl','Nov','Rbp1','Jag1','Rbpj','Notch2')
NFAT<- c('Lmcd1','Rcan1','Rcan2','Myoz2')
Growth <- c('Bdnf','Fgf16','Flt1','Npr1','Npr3','Ntrk2','Ntrk3','Pdgfb','Pdgfrb','Spry2','Vegfc','Pdgfa','Nrp2','Fgf1','Fgf2')
Gprotein<-c('Edn1','Gng2','Paqr6','Rgs2','Rgs5','Sgpp1','Gnb4','Gng8','Grb14','Ace','Arrdc3','Cav3','Ednrb','Gng11','Gpr45','Rgs3','Gnai1')
Epigenetic1<- c('Auts2','Baz1a','Brd2','Cited4','Fhl2','H2afx','Jade1','Jarid2','Sap30','Arid5a','Arid5b','Asf1b','H2afz','Prdm6','Sirt1','Syne1')
Epigenetic2<- c('Bmyc','Egr1','Fosb','Id4','Jun','Osr1','Sox4','Id3','Gata4','Relb','Scx','Yes1','Id2','Fos','Junb')
Smoothmuscle<-c('Cnn1','Csrp2','Mrvi1','Myh11','Mylk','Mylk3','Mylk4','Myocd','Nos3','Tnnc1','Cnn2','Mustn1','Smtn','Tagln','Vim','Tagln2','Acta2')
Miscellaneous<- c('Sod3','Nox4','Lrp1','Ldb2','Cxcl12','Gria3','Wsb1','Fbxl22','Fbxo2','Fbxo30','Fbxo32','Fbxo33','Myd88','Neat1','Malat1','Meg3')

## Load single cell data set from appropriate file path using Read10x function (this should be a non-normalized matrix)
AllSamples_data <-Read10X(data.dir="~/Downloads/16wk-MR_AHF_CNC_nonNormalized_11.19.20/AllSamples/")

## Initialize the Seurat object with the raw (non-normalized data)
AllSamples<- CreateSeuratObject(counts = AllSamples_data, project = "16wkMRandWT", min.cells = 3, min.features = 200)

## Create percent mitochondirial derived reads as metadata, visualize QC metrics as Violin plot  
AllSamples[["percent.mt"]] <- PercentageFeatureSet(AllSamples, pattern = "^mt-")
# Save plots as objects
AllSamples_RNAQC_plot1 <- VlnPlot(AllSamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
AllSamples_plot1 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "percent.mt")
AllSamples_plot2 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
AllSamples_RNAQC_plot2 <- AllSamples_plot1 + AllSamples_plot2
# Visualize plots to analyze
AllSamples_RNAQC_plot1
AllSamples_RNAQC_plot2

## Filtering cells
# Reduces AllSamples data set from 30,704 cells to 24,971 {~19% reduction}
length(colnames(AllSamples))
length(colnames(subset(AllSamples, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 1500 & nCount_RNA < 25000 & percent.mt < 20)))
# Filter All cells 
AllSamples <- subset(AllSamples, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 1500 & nCount_RNA < 25000 & percent.mt < 20)
# Create RNAQCplots, post filtering
AllSamples_filteredRNAQC_plot1 <- VlnPlot(AllSamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
AllSamples_filteredplot1 <- FeatureScatter(AllSamples, feature1 = "nCount_RNA", feature2 = "percent.mt")
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

#mouseID
AllSamples <- AddMetaData(object = AllSamples, col.name = "mouseID", metadata = 
                            c(rep(x = "B530",length(grep("*-1",colnames(AllSamples)))),
                              rep("B388",length(grep("*-2",colnames(AllSamples)))),
                              rep("B492",length(grep("*-3",colnames(AllSamples)))),
                              rep("B493",length(grep("*-4",colnames(AllSamples)))),
                              rep("B531",length(grep("*-5",colnames(AllSamples)))), 
                              rep("B389",length(grep("*-6",colnames(AllSamples))))))
head(AllSamples$mouseID)
tail(AllSamples$mouseID)

#CRE
AllSamples <- AddMetaData(object = AllSamples, col.name = "CRE", metadata = 
                            c(rep(x = "CNC",length(grep("*-1",colnames(AllSamples)))),
                              rep("CNC",length(grep("*-2",colnames(AllSamples)))),
                              rep("AHF",length(grep("*-3",colnames(AllSamples)))),
                              rep("AHF",length(grep("*-4",colnames(AllSamples)))),
                              rep("None",length(grep("*-5",colnames(AllSamples)))), 
                              rep("None",length(grep("*-6",colnames(AllSamples))))))
head(AllSamples$CRE)
tail(AllSamples$CRE)

#geno
AllSamples <- AddMetaData(object = AllSamples, col.name = "geno", metadata = 
                            c(rep(x = "CTL",length(grep("*-1",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-2",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-3",colnames(AllSamples)))),
                              rep("CTL",length(grep("*-4",colnames(AllSamples)))),
                              rep("LDS",length(grep("*-5",colnames(AllSamples)))), 
                              rep("LDS",length(grep("*-6",colnames(AllSamples))))))
head(AllSamples$geno)
tail(AllSamples$geno)

#batch
AllSamples <- AddMetaData(object = AllSamples, col.name = "batch", metadata = 
                            c(rep(x = "2",length(grep("*-1",colnames(AllSamples)))),
                              rep("1",length(grep("*-2",colnames(AllSamples)))),
                              rep("2",length(grep("*-3",colnames(AllSamples)))),
                              rep("2",length(grep("*-4",colnames(AllSamples)))),
                              rep("2",length(grep("*-5",colnames(AllSamples)))), 
                              rep("1",length(grep("*-6",colnames(AllSamples))))))
head(AllSamples$batch)
tail(AllSamples$batch)

#Normalize data using SCTransform
AllSamples.sct<- SCTransform(AllSamples, verbose = FALSE)

#Integrate batches using RPCA 
AllSamples.list <- SplitObject(AllSamples.sct, split.by = "batch")
AllSamples.list <- lapply(X = AllSamples.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = AllSamples.list, nfeatures = 3000)
AllSamples.list <- PrepSCTIntegration(object.list = AllSamples.list, anchor.features = features)
AllSamples.list <- lapply(X = AllSamples.list, FUN = RunPCA, features = features)
aorta.anchors <- FindIntegrationAnchors(object.list = AllSamples.list, normalization.method = "SCT", 
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
aorta.combined.sct <- IntegrateData(anchorset = aorta.anchors, normalization.method = "SCT", dims = 1:30)

#Run non-linear dimensional reduction 
aorta.combined.sct <- RunPCA(aorta.combined.sct, verbose = FALSE)
aorta.combined.sct <- RunUMAP(aorta.combined.sct, reduction = "pca", dims = 1:30)
aorta.combined.sct <- FindNeighbors(aorta.combined.sct, dims = 1:30, verbose = FALSE)
aorta.combined.sct <- FindClusters(aorta.combined.sct, verbose = FALSE, resolution = 0.25) #resolution used was 0.25 
DimPlot(aorta.combined.sct, label = TRUE) + NoLegend()

#Also normalize and scale just the RNA assay (to use for differential expression analysis)
aorta.combined.sct <-NormalizeData(aorta.combined.sct, assay='RNA') 
aorta.combined.sct <-ScaleData(aorta.combined.sct, assay='RNA') 

#Visualize UMAPs of data set split by metadata 
DimPlot(aorta.combined.sct, reduction = "umap") 
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "CRE")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "CRE")
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "geno")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "geno")
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "mouseID")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "mouseID")
DimPlot(aorta.combined.sct, reduction = "umap", group.by = "batch")
DimPlot(aorta.combined.sct, reduction = "umap", split.by = "batch")

# Identify major cell types using canonical markers/those validated by earlier aortic scRNAseq expereiments 
DotPlot(aorta.combined.sct, features = canonical.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))

## Name clusters 
aorta.combined.sct_named <- RenameIdents(aorta.combined.sct, `0` = "VSMCs", `1` = "Fibroblasts", `2` = "Fibroblasts", `3` = "Fibroblasts", `4` = "VSMCs", `5` = "Fibroblasts", `6` = "Immune", `7` = "Endothelial", `8` = "Endothelial",`9`="Immune")
aorta.combined.sct_named2 <- RenameIdents(aorta.combined.sct, `0` = "VSMC1", `1` = "Fibro1", `2` = "Fibro2", `3` = "Fibro3", `4` = "VSMC2", `5` = "Fibro4", `6` = "Immune1", `7` = "Endo1", `8` = "Endo2",`9`="Immune2")

#Subset just samples with CNC (Wnt1v2) CRE 
aorta.combined.sct_named2_CNConly <- subset(aorta.combined.sct_named2, CRE=="CNC")

c1<-DotPlot(aorta.combined.sct_named2_CNConly, assay="RNA", features = c('GFP'), dot.scale = 6,scale.min=0,scale.max = 10, col.min = -1,col.max=1) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1","VSMC2" )) +geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))


##### Analysis of data collapsed into four major clusters: VSMCs, Fibroblasts, Endothelial Cells, and Leukocytes 

#create a version of the file where idents are split by genotype for downstream analysis 
aorta.combined.sct.named.splitbygeno <- aorta.combined.sct_named
aorta.combined.sct.named.splitbygeno$geno <- paste(Idents(aorta.combined.sct.named.splitbygeno), aorta.combined.sct.named.splitbygeno$geno, sep = "_")
aorta.combined.sct.named.splitbygeno$celltype <- Idents(aorta.combined.sct.named.splitbygeno)
Idents(aorta.combined.sct.named.splitbygeno) <- "geno"

#UMAP split by relevant metadata saved as objects for PDF export 
p1<- DimPlot(aorta.combined.sct_named, label = FALSE) 
p2<- DimPlot(aorta.combined.sct_named, label = TRUE, split.by = "geno") 
p3<- DimPlot(aorta.combined.sct_named2, reduction = "umap", split.by = "mouseID")
p4<- DimPlot(aorta.combined.sct_named, reduction = "umap", split.by = "batch")

#Visualize cell type defining markers in the collapsed clusters 
p5<-DotPlot(aorta.combined.sct_named , assay="RNA", features = canonical.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))
p6<- FeaturePlot(aorta.combined.sct_named, features = canonical.markers)
p6.1 <- FeaturePlot(aorta.combined.sct_named, features = c("Myh11","Lum","Pecam1","Ptprc"))

#Generate tables of cell number per cluster, and how these numbers are split by genotype 
p7<- table(Idents(aorta.combined.sct_named))
p8<- table(Idents(aorta.combined.sct_named), aorta.combined.sct_named$geno)

#Pseudo bulk excel sheet of average expression of all transcripts per cluster split by genotype 
Collapsed4MajorClusters_CTLvsMR_AverageExpression<- AverageExpression(aorta.combined.sct.named.splitbygeno, assay='RNA', verbose = TRUE)
write.csv(Collapsed4MajorClusters_CTLvsMR_AverageExpression, file = '../macfarlanelab/Desktop/Collapsed4MajorClusters_CTLvsMR_AverageExpression.csv')


##### Analysis of VSMCs major subclusters (not re-clustered after subsetting)
#subset 
aorta.combined.sct.VSMCsonly <- subset(aorta.combined.sct_named2 , idents = c("VSMC1","VSMC2"))


#create a version of the file where idents are split by genotype for downstream analysis 
aorta.combined.sct.VSMCsonly.splitbygeno <- aorta.combined.sct.VSMCsonly
aorta.combined.sct.VSMCsonly.splitbygeno$geno <- paste(Idents(aorta.combined.sct.VSMCsonly.splitbygeno), aorta.combined.sct.VSMCsonly.splitbygeno$geno, sep = "_")
aorta.combined.sct.VSMCsonly.splitbygeno$celltype <- Idents(aorta.combined.sct.VSMCsonly.splitbygeno)
Idents(aorta.combined.sct.VSMCsonly.splitbygeno) <- "geno"

#UMAP split by relevant metadata saved as objects for PDF export 
r1<- DimPlot(aorta.combined.sct.VSMCsonly, label = TRUE) +xlim(-5,7)+ylim(-12,0)
r2<- DimPlot(aorta.combined.sct.VSMCsonly, label = FALSE , split.by = "geno") +xlim(-5,7)+ylim(-12,-1)
r3<- DimPlot(aorta.combined.sct.VSMCsonly, reduction = "umap", split.by = "mouseID") +xlim(-5,7)+ylim(-12,0)
r4<- DimPlot(aorta.combined.sct.VSMCsonly, reduction = "umap", split.by = "batch") +xlim(-5,7)+ylim(-12,0)

#Visualize cell type defining markers in the collapsed clusters 
r5<-DotPlot(aorta.combined.sct.VSMCsonly, assay="RNA", features = canonical.markers, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))
r5.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = canonical.markers)
r6<-DotPlot(aorta.combined.sct.VSMCsonly, assay="RNA", features = RegionSpecific, dot.scale = 6) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,24.5,1))) + geom_vline(xintercept = c(7.5,11.5,14.5))
r6.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = RegionSpecific)
r6.2<- DotPlot(aorta.combined.sct.VSMCsonly,  assay="RNA", features = c('Ptprz1','Postn','Tes','Gata4','Enpep','Notch3'), dot.scale = 6,scale.min=0) + RotatedAxis() + geom_hline(yintercept = c(seq(1.5,20.5,1))) + geom_vline(xintercept = c(3.5,6.5))
r9<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Gata4", "Ptprz1"), min.cutoff=0,max.cutoff=10) +xlim(-5,7)+ylim(-12,0)
r10<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Gata4", "Ptprz1"),min.cutoff=0,max.cutoff=5, blend = TRUE, cols = c("lightgrey", "red", "cyan3")) & xlim(-5,7) &ylim(-12,0)
r10.1<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Enpep", "Tes"),min.cutoff=0,max.cutoff=5, blend = TRUE, cols = c("lightgrey", "red", "cyan3")) & xlim(-5,7) &ylim(-12,0)
r10.2<- FeaturePlot(aorta.combined.sct.VSMCsonly, features = c("Ptprz1","Postn"),min.cutoff=0,max.cutoff=5, blend = TRUE,cols = c("lightgrey", "red", "cyan3")) & xlim(-5,7) &ylim(-12,0)
r11<- DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = c('Tgfb1','Tgfb2','Tgfb3','Ctgf','Serpine1','Pmepa1','Nr4a1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 100) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS","VSMC3_CTL","VSMC3_LDS" )) + geom_hline(yintercept = c(seq(2.5,20.5,2)))
r12<- DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = c('Col3a1','Col1a1','Col1a2','Mfap4','Mfap5','Eln','Rock1','Fbln1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 100) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS","VSMC3_CTL","VSMC3_LDS" )) + geom_hline(yintercept = c(seq(2.5,20.5,2)))
Gata4regulated <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('Gata4','Jun','Ppp1r15a','Cebpb','Btg2','Pim1','Hes1','Fos','Ier5','Klf2','Ier2'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))
Autophagy <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('Cdkn1a', 'Cxcl1', 'Gata4', 'Thbs1', 'Bmp2', 'Igfbp5', 'Irf1'), dot.scale = 6,col.min=-1.5, col.max=1.5) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))
r13 <-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = c('Efemp2','Adam10','Scarb2','Socs2','Rgs4','Gja1','Fn1','Gas1','Cebpb','Cebpd','Mef2c','Klf2','Il6','Gadd45g','Nfkbia','Malat1'), dot.scale = 6, scale.min=0, scale.max=50) + RotatedAxis() + scale_y_discrete(limits = c("VSMCs_CTL","VSMCs_LDS"))


#VSMC1 vs VSMC2 in ctl and LDS 
Subset<-c('Myh11','Cnn1','Tet2','Klf4','Olfm2','Sox9','Tcf21','Malat1','Twist1','Cebpb','Dcn')
r14<-DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno, assay="RNA", features = Subset, dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0, scale.max=20) + RotatedAxis() + scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))


#Generate tables of cell number per cluster, and how these numbers are split by genotype 
r7<- table(Idents(aorta.combined.sct.VSMCsonly))
r8<- table(Idents(aorta.combined.sct.VSMCsonly), aorta.combined.sct.VSMCsonly$geno)

DefaultAssay( aorta.combined.sct.VSMCsonly.splitbygeno) <- "RNA"
aorta.combined.sct.VSMCsonly.splitbygeno<-JoinLayers(aorta.combined.sct.VSMCsonly.splitbygeno)

VSMC1.MRvsWT.VSMCsOnly.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDS", ident.2 = "VSMC1_CTL", verbose = FALSE, logfc.threshold = 0.00)
VSMC2.MRvsWT.VSMCsOnly.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC2_LDS", ident.2 = "VSMC2_CTL", verbose = FALSE,logfc.threshold = 0.00)
VSMC1.MRvsVSMC2.MR.VSMCsOnly.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_LDS", ident.2 = "VSMC2_LDS", verbose = FALSE, logfc.threshold = 0.00)
VSMC1.CTLvsVSMC2.CTL.VSMCsOnly.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMC1_CTL", ident.2 = "VSMC2_CTL", verbose = FALSE, logfc.threshold = 0.00)

AllVSMCs.MRvsWT.VSMCsOnly.diffexp <- FindMarkers(aorta.combined.sct.VSMCsonly.splitbygeno,assay='RNA', ident.1 = "VSMCs_LDS", ident.2 = "VSMCs_CTL", verbose = FALSE,   logfc.threshold = 0.00)

write.csv(VSMC1.MRvsWT.VSMCsOnly.diffexp, file = '~/Desktop/VSMC1.MRvsWT.VSMCsOnly.diffexp.csv') 
write.csv(VSMC2.MRvsWT.VSMCsOnly.diffexp, file = '~/Desktop/VSMC2.MRvsWT.VSMCsOnly.diffexp.csv') 
write.csv(VSMC1.MRvsVSMC2.MR.VSMCsOnly.diffexp, file = '~/Desktop/VSMC1.MRvsVSMC2.MR.VSMCsOnly.diffexp.csv') 
write.csv(VSMC1.CTLvsVSMC2.CTL.VSMCsOnly.diffexp, file = '~/Desktop/VSMC1.CTLvsVSMC2.CTL.VSMCsOnly.diffexp.csv') 

write.csv(AllVSMCs.MRvsWT.VSMCsOnly.diffexp, file = '~/Desktop/AllVSMCs.MRvsWT.diffexp.csv') 

DefaultAssay(aorta.combined.sct.VSMCsonly) <- "RNA"
aorta.combined.sct.VSMCsonly<-JoinLayers(aorta.combined.sct.VSMCsonly)

VSMC1.markers <- FindMarkers(aorta.combined.sct.VSMCsonly,assay='RNA', ident.1 = "VSMC1", verbose = FALSE)
VSMC2.markers <- FindMarkers(aorta.combined.sct.VSMCsonly,assay='RNA', ident.1 = "VSMC2", verbose = FALSE)

write.csv(VSMC1.markers, file = '~/Desktop/VSMC1.markers.csv') 
write.csv(VSMC2.markers, file = '~/Desktop/VSMC2.markers.csv') 


#Gata4 expressing cells 
table(Idents(aorta.combined.sct.VSMCsonly), aorta.combined.sct.VSMCsonly$geno)
sum(GetAssayData(object = aorta.combined.sct.VSMCsonly.splitbygeno, layer ="data")["Gata4",]<0.1) #1165 VSMCS have non-zero Gata4 expression 
Gata4.VSMCsonly<-subset(x = aorta.combined.sct.VSMCsonly.splitbygeno, subset = Gata4 > 0)
Gata4neg.VSMCsonly<-subset(x = aorta.combined.sct.VSMCsonly.splitbygeno, subset = Gata4 < 0.1)
table(Idents(Gata4neg.VSMCsonly)) # ~9% of VSMCs in cluster 1 are GATA4 positive while 22% are GATA4 positive in cluster 2 
# create a data table for stacked bar plot 
# create a dataset
condition <- c(rep("Gata4+" , 2) , rep("Gata4-" , 2))
celltype <- rep(c("VSMC1","VSMC2"),2)
data <- data.frame(condition,celltype))
data$percentage<- c(0.0894,0.9106, 0.2240, 0.7760)

ggplot(data, aes(fill=condition, y=percentage, x=celltype)) + geom_bar(position="stack", stat="identity")
#scale_fill_manual(values=c("#8B93FF","#D575FE","#E18A00","#BE9C00","#24B700","#00C1AB", "#00BBDA" ,"#00ACFC" ,"#FF65AC", "#F8766D" ,"#8CAB00","#00BE70","#F962DD"

#must run all VSMCs together (not split to generate this file)
Gata4pos.MRvsWT.VSMCsOnly.diffexp <- FindMarkers(Gata4.VSMCsonly,assay='RNA', ident.1 = "VSMCs_LDS", ident.2 = "VSMCs_CTL", verbose = FALSE,logfc.threshold = 0.00)
write.csv(Gata4pos.MRvsWT.VSMCsOnly.diffexp, file = '~/Desktop/Gata4pos.MRvsWT.VSMCsOnly.diffexp.csv') 

Gata4neg.MRvsWT.VSMCsOnly.diffexp <- FindMarkers(Gata4neg.VSMCsonly,assay='RNA', ident.1 = "VSMCs_LDS", ident.2 = "VSMCs_CTL", verbose = FALSE,logfc.threshold = 0.00)
write.csv(Gata4neg.MRvsWT.VSMCsOnly.diffexp, file = '~/Desktop/Gata4neg.MRvsWT.VSMCsOnly.diffexp.csv') 

#Pseudo bulk excel sheet of average expression of all transcripts per cluster split by genotype 
VSMCSubClusters_CTLvsMR_AverageExpression<- AverageExpression(aorta.combined.sct.VSMCsonly.splitbygeno, assay='RNA', verbose = TRUE)
write.csv(VSMCSubClusters_CTLvsMR_AverageExpression, file = '../macfarlanelab/Desktop/VSMC1and2_CTLvsMR_pseudobulk.csv')

#https://satijalab.org/seurat/reference/pseudobulkexpression

VSMCsCTLvsMR_pseudobulk<-PseudobulkExpression(aorta.combined.sct.VSMCsonly.splitbygeno, assay='RNA', verbose = TRUE)
write.csv(VSMCs_CTLvsMR_AverageExpression, file = '../macfarlanelab/Desktop/VSMCs_CTLvsMR_pseudobulk.csv')


p1<-DimPlot(aorta.combined.sct.VSMCsonly, label = FALSE, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))
p2<-DotPlot(aorta.combined.sct.VSMCsonly,  assay="RNA", features = c('Gata4','Enpep','Notch3','Tes','Ptprz1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis() 

#save plots 
pdf(file = '~/Documents/Emily Bramel/16wkCTLandMR_2VSMCclusters.pdf')  
p1
p5
q1
r1
r2
r6.2
r10 
a1
c1
dev.off()


p1<-DimPlot(aorta.combined.sct.VSMCsonly, label = FALSE, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D')) & xlim(-5,7) &ylim(-12,-2)
p1.2<-DimPlot(aorta.combined.sct.VSMCsonly, label = FALSE, split.by="geno", cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D')) & xlim(-5,7) &ylim(-12,-2)

SHFvsCNC<- c('Tnnt2','Dcn','Des','Lum','Jund','Junb','Cebpd','Meg3','Klf4','Ptprz1','Tes','Sost','Ncam1','Mylk4','Mylk','Mt1','Gm26771','Filip1l','Efhd1','Cped1','Abcc9', 'Rgs4')
DotPlot(aorta.combined.sct.VSMCsonly,assay="RNA", features = SHFvsCNC, dot.scale = 6, scale.min=0) + RotatedAxis()


pp2<-DotPlot(aorta.combined.sct.VSMCsonly,  assay="RNA", features = c('Gata4','Enpep','Notch3','Tes','Ptprz1'), dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis() 
ENRICHR<-c('Btg2', 'Csf1', 'Cebpd','Cebpb', 'Tnc', 'Slc2a3', 'Socs3', 'Zfp36', 'Nfil3', 'Plau', 'Myc', 'Ccl2', 'Phlda2', 'Junb', 'Phlda1', 'Egr1', 'Jun', 'Gadd45b', 'Fos', 'Nr4a2', 'Nr4a1', 'Il6', 'Bmp2', 'Nr4a3', 'Fosb', 'Atf3','Rgs4', 'Prrx1', 'Lum', 'Igfbp4', 'Cxcl1', 'Pcolce', 'Msx1', 'Thbs1', 'Dcn', 'Il32','Lama3', 'Thy1', 'Cdh6', 'Cxcl12', 'Col5a3', 'Snai2', 'Timp3', 'Slit3', 'Cdkn1a', 'Gsn', 'Irf1', 'Hmgb2', 'Krt18', 'Hgf', 'Plppr4', 'Avpr1a', 'Plk2', 'Hes1', 'Srpx', 'Mt2a', 'Angptl4', 'Pgf', 'Mt1e', 'Lxn', 'Tgfb3', 'Cited2', 'Pgam2', 'Ndrg1','Cp' )
ENRICHRsubset<-c('Cebpd','Cebpb', 'Myc', 'Egr1','Cdkn1a', 'Irf1', 'Jun','Fos','Il6', 'Bmp2', 'Msx1', 'Thbs1', 'Dcn' , 'Tgfb3', 'Cited2')

DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = ENRICHRsubset, dot.scale = 6,col.min=-1.5, col.max=1.5, scale.min=0,scale.max = 20) + RotatedAxis()+ scale_y_discrete(limits = c("VSMC1_CTL","VSMC1_LDS","VSMC2_CTL","VSMC2_LDS"))

AllVSMCsUpandDown<-c('Atf3', 'Btg1', 'Btg2', 'Cdkn1a', 'Ddit3', 'Fos', 'Hbegf', 'Ier3', 'Ier5', 'Jun', 'Klf4', 'Ndrg1', 'Ppp1r15a', 'Tob1', 'Txnip', 'Zbtb16', 'Bcl2l1', 'Cebpb', 'Cebpd', 'Irf1', 'Junb', 'Mcl1', 'Prkcd', 'Socs3', 'Efemp2', 'Eln', 'Fbln1', 'Fbln5', 'Fbn1', 'Fn1', 'Itga5', 'Itgb3', 'Loxl1', 'Loxl3', 'Ltbp1', 'Ltbp3', 'Mfap4', 'Mfap5', 'Cav3', 'Col1a1', 'Col1a2', 'Col4a6', 'Col5a2', 'Col6a2', 'Ctnnb1', 'Flnc', 'Itga9', 'Lamc1', 'Mylk4', 'Pdgfrb', 'Tnxb', 'Ap3m1', 'Cd164', 'Ctsb', 'Ctsc', 'Ctsh', 'Gaa', 'Gns', 'Gusb', 'Hexa', 'Hexb', 'Lamp1', 'Laptm4a', 'Laptm4b', 'Lgmn', 'Psap', 'Scarb2', 'Adamts2', 'Bmp1', 'Col14a1', 'Col15a1', 'Col3a1', 'Col4a5', 'Col5a1', 'Col6a1', 'Col6a3', 'Colgalt1', 'Crtap', 'Mmp17', 'P3h3', 'Pcolce2', 'Serpinh1', 'Timp2', 'Ankh', 'Dcn', 'Diablo', 'F2r', 'Gadd45B', 'Gsn', 'H10', 'Hmgb2', 'Hspb1', 'Lum', 'Rock1', 'Sc5d', 'Smad7', 'Timp3', 'Ccn1', 'Dusp1', 'Egr1', 'Fosb', 'Jund', 'Mt2a')


AllVSMCsUpandDownSubset<-c('Atf3', 'Btg1', 'Btg2', 'Cdkn1a', 'Ddit3', 'Fos', 'Hbegf', 'Ier3', 'Ier5', 'Klf4', 'Ndrg1', 'Ppp1r15a', 'Tob1', 'Txnip', 'Zbtb16', 'Cebpb', 'Cebpd', 'Irf1', 'Junb', 'Mcl1', 'Prkcd', 'Socs3', 'Efemp2', 'Eln', 'Fbln1', 'Fbln5', 'Fbn1', 'Fn1', 'Itga5', 'Itgb3', 'Loxl1', 'Loxl3', 'Ltbp1', 'Ltbp3', 'Mfap4', 'Mfap5', 'Cav3', 'Col1a1', 'Col1a2', 'Col4a6', 'Col5a2', 'Col6a2', 'Ctnnb1', 'Flnc', 'Itga9', 'Lamc1', 'Mylk4', 'Pdgfrb', 'Tnxb', 'Ap3m1', 'Cd164', 'Ctsb', 'Ctsc', 'Ctsh', 'Gaa', 'Gns', 'Gusb', 'Hexa', 'Hexb', 'Lamp1', 'Laptm4a', 'Laptm4b', 'Lgmn', 'Psap', 'Scarb2', 'Adamts2', 'Bmp1', 'Col14a1', 'Col15a1', 'Col3a1', 'Col4a5', 'Col5a1', 'Col6a1', 'Col6a3', 'Colgalt1', 'Crtap', 'Mmp17', 'P3h3', 'Pcolce2', 'Serpinh1', 'Timp2', 'Ankh', 'Dcn', 'Diablo', 'F2r', 'Gadd45B', 'Gsn', 'H10', 'Hmgb2', 'Hspb1', 'Lum', 'Rock1', 'Sc5d', 'Smad7', 'Timp3', 'Ccn1', 'Dusp1', 'Egr1', 'Fosb', 'Jund', 'Mt2a')

AllVSMCsUpandDownSubset<-c('Eln','Carmn','Smtn', 'Efemp2', 'Fbn1', 'Fn1', 'Itga5', 'Itga9', 'Loxl1', 'Loxl3', 'Ltbp1', 'Ltbp3', 'Mfap4', 'Mfap5', 'Col1a1', 'Col1a2', 'Col4a6', 'Col5a2', 'Col6a2', 'Flnc' , 'Lamc1', 'Mylk4','Lamp1', 'Laptm4a', 'Laptm4b', 'Lgmn','Cnn1', 'Rock1', 'Egr1', 'Smad7', 'Klf4', 'Ndrg1', 'Ppp1r15a', 'Zbtb16', 'Cebpb', 'Cebpd', 'Irf1', 'Junb', 'Meg3','Ier2','Ier3','Klf2','Klf7','Cyr61','Ctgf','Gata4','Thbs1')
DotPlot(aorta.combined.sct.VSMCsonly.splitbygeno,assay="RNA", features = AllVSMCsUpandDownSubset, dot.scale = 6, scale.min=0, col.min=-0, col.max=2.5) + RotatedAxis()+ scale_y_discrete(limits = c("VSMCs_CTL","VSMCs_LDS"))

pdf(file = '~/Desktop/Mouse_plotsforFig2.pdf')  
p1
p1.2
p2
dev.off()



#CoGAPS 
library(projectR)
library(CoGAPS)
library(ggpubr)
test<- aorta.combined.sct.VSMCsonly
scaled_mat <- aorta.combined.sct.VSMCsonly@assays$SCT@scale.data

result <- readRDS("~/Downloads/CoGAPS_8patterns_16wkLDSmice_result.rds")
coGAPsMatrix<-result@featureLoadings #to turn coGAPs result into a usable matrix 

projectR_result_t <-projectR(data= scaled_mat, loadings= coGAPsMatrix)

pattern_weights <- as.data.frame(t(projectR_result_t))

for(pattern in colnames(pattern_weights)) {
  test <- AddMetaData(object = test,metadata = pattern_weights[[pattern]],
                      col.name = pattern)
}

#Visualization of CoGAPS patterns 
test$cluster <- Idents(aorta.combined.sct.VSMCsonly)

DimPlot(samples.combined.sct.VSMCsonly, reduction = "umap",label = "True") & xlim(-6.5,5) &ylim(-10,0)
DotPlot(samples.combined.sct.VSMCsonly, features= c('Gata4',"Notch3",'Enpep','Ptprz1','Tes'), dot.scale = 6, scale.min = 0) + RotatedAxis()

VlnPlot(test, features = c("Pattern_4"),pt.size = 0.05, alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))

p+ stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2 )
VlnPlot(test , features = c("Pattern_5"),pt.size = 0.05,alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))

FeaturePlot(test, c("Pattern_1","Pattern_2",'Pattern_3',"Pattern_4","Pattern_5","Pattern_6", "Pattern_7","Pattern_8"), pt.size = 0.1, min.cutoff=0, max.cutoff=2, label=TRUE) & scale_color_viridis_c(option="plasma") & xlim(-5,7) &ylim(-12,-1)
FeaturePlot(test, c("Pattern_4","Pattern_5"), pt.size = 0.1, min.cutoff=0, max.cutoff=2) & scale_color_viridis_c(option="plasma") &  xlim(-5,7) &ylim(-12,-1)
FeaturePlot(test, c("Pattern_5"), pt.size = 0.1, label=TRUE) & scale_color_viridis_c(option="magma") & xlim(-10,15) &ylim(-8,15)


test$cluster <- Idents(aorta.combined.sct.VSMCsonly)
vsmc1_pat4 <- test@meta.data[test@meta.data$cluster == "VSMC1", "Pattern_4"]
vsmc2_pat4 <- test@meta.data[test@meta.data$cluster == "VSMC2", "Pattern_4"]
wilcox.test(x = vsmc1_pat4, y = vsmc2_pat4)

vsmc1_pat5 <- test@meta.data[test@meta.data$cluster == "VSMC1", "Pattern_5"]
vsmc2_pat5 <- test@meta.data[test@meta.data$cluster == "VSMC2", "Pattern_5"]
wilcox.test(x = vsmc1_pat5, y = vsmc2_pat5)


VlnPlot(test , features = c("Pattern_4"),pt.size = 0.05, alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))+stat_pvalue_manual(
  , label = "wilcox.test, p = {p}")


VlnPlot(test , features = c("Pattern_4"),pt.size = 0.05, alpha=0.6, cols=c('VSMC1'='#00BFC4','VSMC2'='#F8766D'))+ stat_compare_means(comparisons =test$Pattern_4, label = "p.signif")


