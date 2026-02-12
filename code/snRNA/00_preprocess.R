#----------------------------------------------------------------------------
# 00_preprocess.R
#----------------------------------------------------------------------------
library(Seurat)
library(dplyr)

#---- make SeuratObject -----#
sampleName <- list.files("data")
d <- lapply(sampleName, function(x) Read10X(paste0("data/", x)))
names(d) <- sampleName

pre_seu <- lapply(names(d), function(x) CreateSeuratObject(d[[x]], project = x))
pre_seu <- merge(pre_seu[[1]], pre_seu[-1])
pre_seu <- JoinLayers(pre_seu)

pre_seu[["percent.mt"]] <- PercentageFeatureSet(pre_seu, pattern = "^MT-")
p1 <- VlnPlot(pre_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
p2 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("output/Plots/00_QCplots.pdf")
p1
p2
p3
dev.off()

saveRDS(pre_seu, "rds/00_pre_seu.rds")

#---- Quality filter -----#
seu <- subset(pre_seu, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 10) #21523 cells

saveRDS(seu, "rds/00_seu.rds")
