######################################################################################
## We process the single cell data from Puram et al in accordance with Arora et al. ##
######################################################################################

library(Seurat)
library(reshape2)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)

# Load in and process single-cell data 
scData <- read.delim2("Data/GSE103322_HNSCC_all_data.txt")
scData$X  <- as.list(sapply(scData$X , function(x) gsub("\'", "", x)))
rownames(scData) <- scData$X
scData$X <- NULL

metaData <- scData[0:5,]
metaData <- as.data.frame(t(metaData))
expressionData <- scData[6:23691,]
sData <- CreateSeuratObject(expressionData, project = "puram_data")
sData <- AddMetaData(sData,metadata = metaData)

sData <- FindVariableFeatures(sData)
sData <- ScaleData(sData)
sData <- Seurat::RunPCA(sData, verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(., dims = 1:30)%>%
  FindClusters(., resolution = 0.5)

sData@meta.data$non.cancer.cell.type <- gsub("-","",as.character(sData@meta.data$non.cancer.cell.type))

# Rename cancer cells
sData@meta.data$non.cancer.cell.type <- gsub("0","Cancer cell",as.character(sData@meta.data$non.cancer.cell.type))
Seurat::DimPlot(sData,
                group.by = "RNA_snn_res.0.5",
                label = TRUE) + Seurat::NoLegend()

# Identify celltypes based on published markers
cluster_ids <- read.csv(file = "Data/cluster_identification.csv")

Idents(sData) = "RNA_snn_res.0.5"
new.cluster.ids <- cluster_ids$cell_type
names(new.cluster.ids) <- levels(sData)
sData <- RenameIdents(sData, new.cluster.ids)
sData@meta.data$cellype_fine <- Idents(sData)

# Make factor with cell names
sData@meta.data$cellype_fine <- factor(sData@meta.data$cellype_fine, levels = 
                                               c("myofibroblast","cancer cell","B cell","ecm-myCAF",
                                                 "Intermediate fibroblast","detox-iCAF","macrophage",
                                                 "endothelial","dendritic ","mast","conventional CD4+ T-helper cells",
                                                 "cytotoxic CD8+ T ","Tregs","cytotoxic CD8+ T exhausted"))

Idents(sData) = "cellype_fine"
puram_data <- sData
save(puram_data, file = "./Data/puram_data.Robj")                    
