#----------------------------------------------------------------------------
# 04_clustering_sample.R
#----------------------------------------------------------------------------
#Seurat5
library(Seurat)
library(parallel)
library(ggplot2)
library(ggrastr)
library(gghalves)
library(scales)
library(dplyr)
library(alphahull)
library(viridisLite)
library(stringr)
library(GeneOverlap)
library(ComplexHeatmap)
library(circlize)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

seu <- readRDS("rds/01_seu.rds")

#--------- clustering per sample ----------#
seu$Sample <- str_split(seu$orig.ident, "_", simplify = T)[,1]
sample_list <- unique(seu$Sample)

seu_ls <- lapply(sample_list, function(x){
  out <- seu[,which(seu$Sample == x)]
  out <- NormalizeData(out, normalization.method = "LogNormalize", scale.factor = 10000)
  out <- FindVariableFeatures(out, selection.method = "vst", nfeatures = 2000)
  out <- ScaleData(out)
  out <- RunPCA(out)
  out <- FindNeighbors(out, dims = 1:30)
  out <- RunUMAP(out, dims = 1:30)
  out <- FindClusters(out, resolution = 0.6)
  return(out)
})
names(seu_ls) <- sample_list


cluster_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
                    "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D",
                    "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A")

p1 <- lapply(sample_list, function(x) DimPlot(seu_ls[[x]], cols = cluster_colors) + theme_cb() + ggtitle(x))
pdf("output/Plots/04_UMAP_sample.pdf", width = 4, height = 4)
p1
dev.off()

RP_genes <- readRDS("rds/02_RP_genes.rds")

seu_ls <- lapply(seu_ls, function(x) AddModuleScore(x, features = RP_genes, name = "RP"))

p2 <-lapply(sample_list, function(x){
  out <- lapply(paste0("RP", c(1:8)), function(y){
    FeaturePlot(seu_ls[[x]], features = y) + 
      scale_color_gradientn(colours = c(rev(viridis(6, option = "G")))) + 
      ggtitle(paste0(y, ": ", x)) + theme_cb()
  })
  return(out)
})
pdf("output/Plots/04_UMAPOL_RP_sample.pdf", width = 4, height = 4)
p2
dev.off()

#--------- clustering per sample ----------#
Hwang_list <- read.csv("ref/Hwang.csv")
tmp1 <- colnames(Hwang_list)
Hwang_list <- lapply(tmp1, function(x) Hwang_list[,x])
names(Hwang_list) <- tmp1

JI <- lapply(RP_genes, function(x){
  out <- lapply(Hwang_list, function(y){
    gov_tmp <- newGeneOverlap(x, y, genome.size=36601)
    gov_tmp <- testGeneOverlap(gov_tmp)
    return(gov_tmp@Jaccard)
  })
  out <- unlist(out)
  return(out)
})
JI <- do.call(rbind, JI)

FD <- lapply(RP_genes, function(x){
  out <- lapply(Hwang_list, function(y){
    gov_tmp <- newGeneOverlap(x, y, genome.size=36601)
    gov_tmp <- testGeneOverlap(gov_tmp)
    return(gov_tmp@pval)
  })
  out <- unlist(out)
  out <- p.adjust(out, method = "fdr")
  return(out)
})
FD <- do.call(rbind, FD)
FD <- -log10(FD)
FD[is.infinite(FD)] <- max(FD[which(is.finite(FD)==T)])

## visualization ##
col_fun1 <- colorRamp2(seq(0,0.5,0.1), c("white", rev(viridis(5, option = "A"))))
col_fun2 <- colorRamp2(c(0,1.9999,2,50,100,150,200), c("gray", "gray", rev(viridis(5, option = "E"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")

ht1 <- Heatmap(JI, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh, 
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
p3 <- draw(ht1)
ht2 <- Heatmap(FD, name = "-log10(FDR)", col = col_fun2, 
               cluster_rows = fh, cluster_columns = fh, 
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
               heatmap_legend_param = list(at = c(0,2,50,100,150,200)))
p4 <- draw(ht2)

ht3 <- Heatmap(FD, name = "-log10(FDR)", col = col_fun2, 
               cluster_rows = F, cluster_columns = F, 
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
               heatmap_legend_param = list(at = c(0,2,50,100,150,200)))
p5 <- draw(ht3)

pdf("output/Plots/04_Heatmap_RPvsHwang.pdf",width = 5, height = 4)
p3
p4
p5
dev.off()
