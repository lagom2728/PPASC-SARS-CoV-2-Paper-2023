##### Code for scRNA merge & Integration
##### Author: Hyundong Yoon (lagom2728@snu.ac.kr) / Logan S. Dean (Ls.Dean@colostate.edu) 

## STEP 1 : Calling packages 

library(Seurat)
library(harmony)
library(scuttle)
library(scater)
library(scran)
library(pheatmap)
library(ggplot2)
library(dplyr)

## STEP 2 : Read integrated data - Integ_PPASC_Normal / Integ_SEVERE_Normal

Integ_object1 = readRDS("~/003_Single_cell_RNA_Integration/003_Integ_PPASC_Normal")
Integ_object2 = readRDS("~/003_Single_cell_RNA_Integration/003_Integ_SEVERE_Normal")

## STEP 3 : UMAP visualizationi for cell cluster check

DimPlot(Integ_object1, label = T, group.by = "seurat_clusters", raster = F)
DimPlot(Integ_object2, label = T, group.by = "seurat_clusters", raster = F)

## STEP 4 : Cannonical marker gene define for cell type profile 

Erythrocyte = c("HBB","HBA2","HBA1")
Bcell <- c("MS4A1","CD22","CD19","CD79A")
NKcell = c("NCAM1","KLRD1","GNLY","NKG7","GZMH","CD247","CD160","GZMB")
NKTcell = c("ZNF683","TYROBP")
Platelet = c("PPBP","CD9","PF4")
CD8Tcell = c("CD8A","CD8B","CD3D","CD3E","CD3G")
CD4Tcell = c("CD40LG","CTLA4","IL7R","LTB")
CD16_Monocyte = c("FCGR3A","MS4A7","HES4","TCF7L2","MTSS1")
CD14_Monocyte = c("CD14","CD68","S100A12")
DC = c("CD1C","CD83","CD86","CLEC10A","CST3","FCER1A")
HSPC <- c("CD34","MYB","FLT3")
Neutrophil = c("FUT4","MPO","ANPEP","TNFSF13B")

markers = c(Erythrocyte,Bcell,NKcell,NKTcell,Platelet,
            CD8Tcell,CD4Tcell,CD16_Monocyte,CD14_Monocyte,
            DC,HSPC,Neutrophil)

## STEP 5 : Cell type define - Integ_PPASC_Normal 

Integ_object1 <- RenameIdents(object = Integ_object1,
                              "0" = "NKcell",
                              "1" = "CD14+Monocyte",
                              "2" = "CD4+Tcell",
                              "3" = "CD4+Tcell",
                              "4" = "CD4+Tcell",
                              "5" = "CD8+Tcell",
                              "6" = "CD4+Tcell",
                              "7" = "NKTcell",
                              "8" = "NKcell",
                              "9" = "Bcell",
                              "10" = "CD4+Tcell",
                              "11" = "NKcell",
                              "12" = "NKcell",
                              "13" = "CD8+Tcell",
                              "14" = "CD16+Monocyte",
                              "15" = "CD16+Monocyte",
                              "16" = "Erythrocyte",
                              "17" = "Dendritic_cell",
                              "18" = "NKcell",
                              "19" = "Bcell",
                              "20" = "NKcell",
                              "21" = "NKcell",
                              "22" = "CD14+Monocyte",
                              "23" = "CD14+Monocyte",
                              "24" = "HSC",
                              "25" = "Platelet_cell")

Integ_object$ident <- Integ_object@active.ident
Integ_object$celltypes <- Integ_object@active.ident
Integ_object = SetIdent(Integ_object, value = Integ_object@active.ident)

DimPlot(Integ_object, label = T, raster = F)

## STEP 6 : Cell type define - Integ_SEVERE_Normal 

Integ_object2 <- RenameIdents(object = Integ_object2,
                             "0" = "NKTcell",
                             "1" = "CD8Tcell",
                             "2" = "Platelet",
                             "3" = "Dendritic_cell",
                             "4" = "CD14_Monocyte",
                             "5" = "CD8Tcell",
                             "6" = "Neutrophil",
                             "7" = "Bcell",
                             "8" = "CD16_Monocyte",
                             "9" = "CD8Tcell",
                             "10" = "Platelet",
                             "11" = "Erythrocyte",
                             "12" = "NKTcell")

Integ_object2$ident <- Integ_object2@active.ident
Integ_object2$celltypes <- Integ_object2@active.ident

DimPlot(Integ_object2, label = T, raster = F)


## STEP 7 : Save cell type profiled data - Integ_PPASC_Normal / Integ_SEVERE_Normal 

saveRDS(Integ_object1,"~/003_Single_cell_RNA_Integration/003_Integ_PPASC_Normal")
saveRDS(Integ_object2,"~/003_Single_cell_RNA_Integration/003_Integ_SEVERE_Normal")



















