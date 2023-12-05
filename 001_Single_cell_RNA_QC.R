##### Code for scRNA QC 
##### Author: Hyundong Yoon (lagom2728@snu.ac.kr) / Logan S. Dean (Ls.Dean@colostate.edu) 

## STEP 1. Calling packages 

library(Seurat)
library(DropletUtils)
library(scuttle)
library(scater)

## STEP 2 : Set directory 

# PPASC data  

Dir1 = "~/PPASC1"
Dir2 = "~/PPASC2"

# Healthy Normal data 

Dir3 = "~/Hawaii_Normal"
Dir4 = "~/KAIST_Normal_220"

# Severe data 

Dir5 = "~/KAIST_S215"
Dir6 = "~/KAIST_S216"
Dir7 = "~/KAIST_S221"
Dir8 = "~/KAIST_S222"

## STEP 3 : Remove empty droplets 

empty_droplet_list <- list()  

empty_droplet_list[[1]] <- Dir1
empty_droplet_list[[2]] <- Dir2
empty_droplet_list[[3]] <- Dir3
empty_droplet_list[[4]] <- Dir4
empty_droplet_list[[5]] <- Dir5
empty_droplet_list[[6]] <- Dir6
empty_droplet_list[[7]] <- Dir7
empty_droplet_list[[8]] <- Dir8

for (i in 1:8) {
  sce <- empty_droplet_list[[i]]
  
  sce.1 <- read10xCounts(Dir1)
  rownames(sce.1) = uniquifyFeatureNames(rowData(sce.1)$ID, rowData(sce.1)$Symbol)
  colnames(sce.1) = sce.1$Barcode
  my.counts.1=counts(sce.1)
  e.out.1 <- emptyDrops(my.counts.1)
  is.cell.1 <- e.out.1$FDR <= 0.01
  sum(is.cell.1, na.rm=TRUE)
  is.cell.1[is.na(is.cell.1)] <- FALSE
  names(is.cell.1) <- colnames(sce.1)
  sce.1$cells_kept <- is.cell.1
  sce.1 <- sce.1[, sce.1$cells_kept == T]
  empty_droplet_list[[i]] <- sce.1
}

for (i in 1:8) {
  sce <- empty_droplet_list[[i]]
  base::print(dim(sce))}

## STEP 4 : Save filtered data - removed empty droplet

saveRDS(empty_droplet_list[[1]],"~/001_Removed_Empty/001_PPASC1_Empty")
saveRDS(empty_droplet_list[[2]],"~/001_Removed_Empty/001_PPASC2_Empty")
saveRDS(empty_droplet_list[[3]],"~/001_Removed_Empty/001_Hawaii_Normal_Empty")
saveRDS(empty_droplet_list[[4]],"~/001_Removed_Empty/001_KAIST_Normal_Empty")
saveRDS(empty_droplet_list[[5]],"~/001_Removed_Empty/001_KAIST_S215_Empty")
saveRDS(empty_droplet_list[[6]],"~/001_Removed_Empty/001_KAIST_S216_Empty")
saveRDS(empty_droplet_list[[7]],"~/001_Removed_Empty/001_KAIST_S221_Empty")
saveRDS(empty_droplet_list[[8]],"~/001_Removed_Empty/001_KAIST_S222_Empty")

## STEP 5 : Read filtered data - removed empty droplet

sce.1 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_PPASC1_Empty")
sce.2 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_PPASC2_Empty")
sce.3 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_Hawaii_Normal_Empty")
sce.4 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_KAIST_Normal_Empty")
sce.5 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_KAIST_S215_Empty")
sce.6 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_KAIST_S216_Empty")
sce.7 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_KAIST_S221_Empty")
sce.8 = readRDS("E:/Hawaii_Git_Hub/001_Single_cell_RNA_QC/001_KAIST_S222_Empty")

## STEP 6 : Cell QC calculation 

sce_list <- list()  

sce_list[[1]] <- sce.1
sce_list[[2]] <- sce.2
sce_list[[3]] <- sce.3
sce_list[[4]] <- sce.4
sce_list[[5]] <- sce.5
sce_list[[6]] <- sce.6
sce_list[[7]] <- sce.7
sce_list[[8]] <- sce.8

for (i in 1:8) {
  sce <- sce_list[[i]]
  
  mtgenes <- rowData(sce)[grep("MT-", rowData(sce)$Symbol), ]$Symbol
  is.mito <- rownames(sce) %in% mtgenes
  table(is.mito)
  
  sce <- addPerCellQC(sce, subsets = list(MT = is.mito), percent_top = c(50, 100, 200, 500), use_altexps = TRUE)
  sce$log10_sum <- log10(sce$sum + 1)
  sce$log10_detected <- log10(sce$detected + 1)
  variable_list <- c("sum", "detected", "log10_sum", "log10_detected",
                     "subsets_MT_percent", "percent.top_50", "percent.top_100",
                     "percent.top_200", "percent.top_500", "total",
                     "subsets_MT_percent", "subsets_MT_detected", "subsets_MT_sum")
  sce <- runColDataPCA(sce, variables = variable_list)
  
  sce_list[[i]] <- sce
}

for (i in 1:8) {
  sce <- sce_list[[i]]
  base::print(dim(sce))}

## STEP 7 : Low QC cell filtration - based on UMI count 

UMI_thresholds <- c(2.3, 2.5, 2.5, 3.4, 2.6, 2.4, 2.2, 2.3)
Doublet_thresholds <- c(4.6, 4.3, 4.3, 4.3, 4.5, 4.5, 4.5, 4.6)

for (num in 1:8) {
  h <- hist(sce_list[[num]]$log10_sum, xlab="log10(Library sizes)", main="", 
            breaks=125, col="grey80", ylab="Number of cells")
  rug(jitter(sce_list[[num]]$log10_sum))
  
  cuts <- cut(h$breaks, c(0, UMI_thresholds[num], Doublet_thresholds[num], Inf))
  
  plot(h, col=c("red","blue","orange")[cuts], xlab = "Total_UMI", main = paste("Total Cell distribution KAIST-", num))
  rug(jitter(sce_list[[num]]$log10_sum))
  abline(v=UMI_thresholds[num], col="red")
  abline(v=Doublet_thresholds[num], col="orange")
  
  sce_discard_counts <- sce_list[[num]]$log10_sum < UMI_thresholds[num] 
  sce_discard_Doublet = sce_list[[num]]$log10_sum > Doublet_thresholds[num] 
  sce_list[[num]]$discard_total <- (sce_discard_counts | sce_discard_Doublet)
  
  plotReducedDim(sce_list[[num]], dimred="PCA_coldata", colour_by="discard_total")
  
  sce_list[[num]] <- sce_list[[num]][, !sce_list[[num]]$discard_total]
  cat("Sample", num, "dimensions after cell QC:", toString(dim(sce_list[[num]])), "\n")
}

## STEP 8 : Low QC cell filtration - based on Feature count 

Feature_thresholds <- c(2.4, 2.4, 2.4, 3.0, 2.3, 2.3, 2.1, 2.15)
Double_Feature_thresholds <- c(4.4, 4.4, 4.4, 4.4, 4.5, 4.5, 4.5, 4.5)

for (num in 1:8) {
  h <- hist(sce_list[[num]]$log10_detected, xlab="log10(Library sizes)", main="", 
            breaks=125, col="grey80", ylab="Number of cells")
  rug(jitter(sce_list[[num]]$log10_detected))
  
  cuts <- cut(h$breaks, c(0, Feature_thresholds[num], Double_Feature_thresholds[num], Inf))
  
  plot(h, col=c("red","blue")[cuts], xlab = "Total_Feature", main = paste("Total Feature distribution KAIST-", num))
  abline(v=Feature_thresholds[num], col="red")
  abline(v=Double_Feature_thresholds[num], col="orange")
  
  sce_discard_counts <- sce_list[[num]]$log10_sum < UMI_thresholds[num]
  sce_discard_Doublet = sce_list[[num]]$log10_sum > Doublet_thresholds[num]
  sce_discard_Feature = sce_list[[num]]$log10_detected < Feature_thresholds[num]
  sce_list[[num]]$discard_total <- (sce_discard_counts | sce_discard_Doublet | sce_discard_Feature)
  
  plotReducedDim(sce_list[[num]], dimred="PCA_coldata", colour_by="discard_total")
  
  sce_list[[num]] <- sce_list[[num]][, !sce_list[[num]]$discard_total]
  cat("Sample", num, "dimensions after feature QC:", toString(dim(sce_list[[num]])), "\n")
}

## STEP 9 : Low QC cell filtration - based on mitochondria % 

Mito_thresholds <- c(30, 30, 30, 30, 20, 20, 20, 20)

for (num in 1:8) {
  h <- hist(sce_list[[num]]$subsets_MT_percent, xlab="log10(Library sizes)", main="", 
            breaks=125, col="grey80", ylab="Number of Genes")
  rug(jitter(sce_list[[num]]$subsets_MT_percent))
  
  cuts <- cut(h$breaks, c(0, Mito_thresholds[num], Inf))
  
  plot(h, col=c("blue","red")[cuts], xlab = "Total_Mitochondria", main = paste("Total Mitochondria distribution KAIST-", num))
  abline(v=Mito_thresholds[num], col="red")
  
  sce_discard_counts <- sce_list[[num]]$log10_sum < UMI_thresholds[num]
  sce_discard_Doublet = sce_list[[num]]$log10_sum > Doublet_thresholds[num]
  sce_discard_Feature = sce_list[[num]]$log10_detected < Feature_thresholds[num]
  sce_discard_mt = sce_list[[num]]$subsets_MT_percent > Mito_thresholds[num]
  sce_list[[num]]$discard_total <- (sce_discard_counts | sce_discard_Doublet | sce_discard_Feature | sce_discard_mt)
  
  plotReducedDim(sce_list[[num]], dimred="PCA_coldata", colour_by="discard_total")
  
  sce_list[[num]] <- sce_list[[num]][, !sce_list[[num]]$discard_total]
  cat("Sample", num, "dimensions after mitochondrial QC:", toString(dim(sce_list[[num]])), "\n")
}

## STEP 10 : Save filtered data - removed Low QC cells

saveRDS(sce_list[1],"~/002_Single_cell_RNA_QC/002_PPASC1_QC")
saveRDS(sce_list[2],"~/002_Single_cell_RNA_QC/002_PPASC2_QC")
saveRDS(sce_list[3],"~/002_Single_cell_RNA_QC/002_Hawaii_Normal_QC")
saveRDS(sce_list[4],"~/002_Single_cell_RNA_QC/002_KAIST_Normal_QC")
saveRDS(sce_list[5],"~/002_Single_cell_RNA_QC/002_KAIST_S215_QC")
saveRDS(sce_list[6],"~/002_Single_cell_RNA_QC/002_KAIST_S216_QC")
saveRDS(sce_list[7],"~/002_Single_cell_RNA_QC/002_KAIST_S221_QC")
saveRDS(sce_list[8],"~/002_Single_cell_RNA_QC/002_KAIST_S222_QC")
