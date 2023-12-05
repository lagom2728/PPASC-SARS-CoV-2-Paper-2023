##### Code for scRNA merge & Integration
##### Author: Hyundong Yoon (lagom2728@snu.ac.kr) / Logan S. Dean (Ls.Dean@colostate.edu) 

## STEP 1 : Calling packages 

library(Seurat)
library(harmony)
library(scuttle)
library(scater)
library(scran)

## STEP 2 : Read filtered data - removed Low QC cells

sce.1 = readRDS("~/002_Single_cell_RNA_QC/002_PPASC1_QC")
sce.2 = readRDS("~/002_Single_cell_RNA_QC/002_PPASC2_QC")
sce.3 = readRDS("~/002_Single_cell_RNA_QC/002_Hawaii_Normal_QC")
sce.4 = readRDS("~/002_Single_cell_RNA_QC/002_KAIST_Normal_QC")
sce.5 = readRDS("~/002_Single_cell_RNA_QC/002_KAIST_S215_QC")
sce.6 = readRDS("~/002_Single_cell_RNA_QC/002_KAIST_S216_QC")
sce.7 = readRDS("~/002_Single_cell_RNA_QC/002_KAIST_S221_QC")
sce.8 = readRDS("~/002_Single_cell_RNA_QC/002_KAIST_S222_QC")

## STEP 3 : Assuming the list 

sce.list <- list(sce.1, sce.2, sce.3, sce.4, sce.5, sce.6, sce.7, sce.8)
projects <- c("case1", "case2", "ctrl2", "Normal220", "KAIST_S215", "KAIST_S216", "KAIST_S221", "KAIST_S222")
conditions <- c("PASC", "PASC", "Control", "Control", "Severe", "Severe", "Severe", "Severe")

seurat.list <- list()

## STEP 4 : Create New Seurat object

for (i in seq_along(sce.list)) {
  seurat.list[[i]] <- CreateSeuratObject(
    counts = counts(sce.list[[i]]),
    min.cells = 0, min.features = 0,
    project = projects[i],
    assay = "RNA"
  )
  seurat.list[[i]]$condition <- conditions[i]
}

## STEP 5 : Merge each data 

merged.list <- list()
merged.list[["Integ_PPASC_Normal"]] <- merge(x = seurat.list[[1]], y = seurat.list[2:4], project = "Hawaii", add.cell.ids = projects[1:4])
merged.list[["Integ_SEVERE_Normal"]] <- merge(x = seurat.list[[3]], y = seurat.list[4:8], project = "Hawaii", add.cell.ids = projects[3:8])

## STEP 6 : Merge data normalization

for (i in seq_along(merged.list)) {
  merged.list[[i]] <- as.SingleCellExperiment(merged.list[[i]])
  clusters <- quickCluster(merged.list[[i]])
  merged.list[[i]] <- computeSumFactors(merged.list[[i]], clusters = clusters)
  merged.list[[i]] <- logNormCounts(merged.list[[i]], pseudo_count = 1)
  gene_variance <- modelGeneVar(merged.list[[i]])
  df <- data.frame(mean = gene_variance$mean, total = gene_variance$total)
  ggplot(df) +
    geom_point(aes(x = mean, y = total), shape = 1, size = 2) +
    stat_function(fun = metadata(gene_variance)$trend, color = "blue") +
    labs(x = "mean_log-expression", y = "variance") +
    scale_x_continuous(breaks = seq(0, 8, 2)) +
    theme_classic()
}

## STEP 7 : Merge data HVG selection & Scaling

for (i in seq_along(merged.list)) {
  hvg.norm <- getTopHVGs(modelGeneVar(merged.list[[i]]), fdr.threshold = 0.05)
  seurat <- as.Seurat(logNormCounts(merged.list[[i]], pseudo_count = 1), counts = "counts", data = "logcounts")
  VariableFeatures(seurat) <- hvg.norm
  all.genes <- rownames(seurat)
  seurat <- ScaleData(seurat, features = all.genes)
  seurat.list[[i]] <- seurat
}

## STEP 8 : Run PCA

for (i in seq_along(seurat.list)) {
  seurat.list[[i]] <- RunPCA(seurat.list[[i]], features = VariableFeatures(seurat.list[[i]]))
}

## STEP 9 : Dimension reduction

for (i in seq_along(seurat.list)) {
  seurat.list[[i]] <- seurat.list[[i]] %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.7) %>%
    RunUMAP(dims = 1:20) %>%
    identity()
}

## STEP 10 : UMAP visualization for check Batch effects

for (i in seq_along(seurat.list)) {
  DimPlot(seurat.list[[i]], group.by = "orig.ident", label = T, raster = F)
  if (i %% 2 == 1) {  # Only plot condition for odd-indexed datasets
    DimPlot(seurat.list[[i]], group.by = "condition", label = T, raster = F)
  }
}

## STEP 11 : Harmony Integration for remove Batch effects

for (i in seq_along(seurat.list)) {
  seurat.list[[i]] <- seurat.list[[i]] %>% RunHarmony("orig.ident", plot_convergence = F, max.iter.harmony = 100)
  seurat.list[[i]] <- seurat.list[[i]] %>%
    RunUMAP(reduction = "harmony", dims = 1:50) %>%
    FindNeighbors(reduction = "harmony", dims = 1:50) %>%
    FindClusters(resolution = 0.7) %>%
    identity()
}

## STEP 12 : UMAP visualization for check Batch effects

for (i in seq_along(seurat.list)) {
  DimPlot(seurat.list[[i]], group.by = "orig.ident", label = T, raster = F)
  if (i %% 2 == 1) {  # Only plot condition for odd-indexed datasets
    DimPlot(seurat.list[[i]], group.by = "condition", label = T, raster = F)
  }
}

## STEP 13 : Save Integrated data - Integ_PPASC_Normal / Integ_SEVERE_Normal 

saveRDS(seurat.list[[1]],"~/003_Single_cell_RNA_Integration/003_Integ_PPASC_Normal")
saveRDS(seurat.list[[2]],"~/003_Single_cell_RNA_Integration/003_Integ_SEVERE_Normal")
