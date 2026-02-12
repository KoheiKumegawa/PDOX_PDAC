#-------------------------------------------------------------------------------
# 02_RNAanalysis.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggrastr)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(stringr)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}
"%ni%" <- Negate("%in%")

se <- readRDS("rds/01_se.rds")

#----- PCA -----#
se_pc <- se[which(rowData(se)$gene_type == "protein_coding")]

row_vars1 <- rowVars(assays(se_pc)$log2tpm)
idx1 <- order(row_vars1, decreasing = T)[c(1:4000)]
pca1 <- prcomp(t(assays(se_pc)$log2tpm[idx1,]))

df <- pca1$x[, c(1:2)] %>% data.frame
df$Sample <- colData(se)$SampleID
df$SampleType <- se$SampleType

summary(pca1)
### PCA plot ###
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = SampleType, label = Sample)) + geom_point(alpha = 0.5, size = 4) + 
  theme_cb() + scale_color_manual(values = c("Primary"="black","Metastasis"="red")) + 
  geom_vline(xintercept = 0, lty = "dashed", col = "gray") + geom_hline(yintercept = 0, lty = "dashed", col = "gray") +
  labs(x = "PC1 (27.4% variance explained)", y = "PC2 (16.9% variance explained)") +
  theme(legend.key.size = unit(0.4, 'cm')) + geom_text_repel()
pdf("output/Plots/02_PCA.pdf", width = 4, height = 3)
p1
dev.off()

### PC loading ###
pca1_loading <- t(t(pca1$rotation)*pca1$sdev)[,c(1:2)] %>% as.data.frame()
pca1_loading$gene <- rowData(se)[rownames(pca1_loading),]$gene_name
write.csv(pca1_loading, "output/Tables/01_PCA_loading.csv")

#----- k-means heatmap -----#
mtx1 <- assays(se_pc)$log2tpm[idx1,]
mtx1_z <- t(scale(t(mtx1)))
colnames(mtx1_z) <- se$SampleID

col_fun1 <- colorRamp2(c(-2,-1,0,1,2), viridis(5, option = "A"))
km1 <- kmeans(mtx1_z, 5)
ha1 <- HeatmapAnnotation(SampleType = colData(se)$SampleType, col = list(SampleType = c("Primary"="black","Metastasis"="red")))
ht1 <- Heatmap(mtx1_z, col = col_fun1, name = "z-score", row_split = km1$cluster, bottom_annotation = ha1, 
               show_row_names = F, cluster_columns = F, cluster_rows = F, use_raster = T)

p2 <- draw(ht1)
df1 <- data.frame(gene_id = names(km1$cluster), kmeans_clst = as.numeric(km1$cluster))
df1$gene_name <- rowData(se)[df1$gene_id,]$gene_name

pdf("output/Plots/02_Heatmap_zscore_RNA.pdf", width = 4, height = 5)
p2
dev.off()
write.csv(df1, "output/Tables/02_Heatmap_zscore_RNA.csv")

#----- marker expression -----#
classical_markers <- c(
  "GATA6", "KRT19", "TFF1", "TFF2", "TFF3", 
  "MUC1", "MUC5AC", "AGR2", "FOXA2", "SPDEF",
  "HNF1A", "HNF4A", "EPCAM"
)

basal_like_markers <- c(
  "KRT5", "KRT6A", "KRT14", "TP63", "S100A2", 
  "VIM", "ZEB1", "JUN", "FOSL1", "MMP7", 
  "CD44", "ITGA6", "CXCL1"
)
idx2 <- which(rowData(se_pc)$gene_name %in% c(classical_markers, basal_like_markers))

mtx2 <- assays(se_pc)$log2tpm[idx2,]
rownames(mtx2) <- rowData(se_pc)$gene_name[idx2]
mtx2_z <- t(scale(t(mtx2)))
colnames(mtx2_z) <- se$SampleID

fh <- function(x) hclust(dist(x), method = "ward.D2")
ht2 <- Heatmap(mtx2_z, col = col_fun1, name = "z-score", bottom_annotation = ha1, 
               show_row_names = T, cluster_columns = fh, cluster_rows = fh)

p3 <- draw(ht2)
pdf("output/Plots/02_Heatmap_zscore_markers.pdf", width = 6, height = 6)
p3
dev.off()

colnames(mtx2) <- se$SampleID
mtx2_melt <- reshape2::melt(mtx2)
mtx2_melt$SampleType <- "Primary"
mtx2_melt$SampleType[mtx2_melt$Var2 %in% c("OP025A","OP026A","OP033A","OP038P","OP039A")] <- "Metastasis"
mtx2_melt$SampleType <- factor(mtx2_melt$SampleType, levels = c("Primary","Metastasis"))
mtx2_melt$Var1 <- factor(mtx2_melt$Var1, levels = c(classical_markers, basal_like_markers))

p4 <- ggplot(mtx2_melt, aes(x=SampleType,y=value)) + geom_jitter(height = 0, width = 0.1) + 
  theme(axis.text.x = element_text(angle=90,hjust = 1)) + 
  geom_signif(comparisons = list(c("Primary","Metastasis")), test = "wilcox.test") +
  facet_grid(~Var1) + labs(x="", y="log2(TPM+1)")
pdf("output/Plots/02_Expression_Log2TPM_markers.pdf", width = 14, height = 4)
p4
dev.off()

oxiphos <- c("NDUFA4", "NDUFB9", "NDUFA8", "NDUFC2", "NDUFA6", "NDUFB4", "NDUFB3", "NDUFA2", "NDUFB11",
  "UQCR10",
  "COX4I1", "COX6C", "COX7C", "COX6A1", "COX8A", "COX7B",
  "ATP5F1A", "ATP5F1E", "ATP5PB", "ATP5PD", "ATP5MD", "ATP5MPL",
  "VDAC1",
  "TIMM8B", "TOMM22", "TMEM70", "GHITM", "MPV17L", "MPV17L2",
  "CYCS")
idx3 <- which(rowData(se_pc)$gene_name %in% oxiphos)

mtx3 <- assays(se_pc)$log2tpm[idx3,]
rownames(mtx3) <- rowData(se_pc)$gene_name[idx3]
colnames(mtx3) <- se$SampleID
mtx3_melt <- reshape2::melt(mtx3)
mtx3_melt$SampleType <- "Primary"
mtx3_melt$SampleType[mtx3_melt$Var2 %in% c("OP025A","OP026A","OP033A","OP038P","OP039A")] <- "Metastasis"
mtx3_melt$SampleType <- factor(mtx3_melt$SampleType, levels = c("Primary","Metastasis"))
mtx3_melt$Var1 <- factor(mtx3_melt$Var1, levels = oxiphos)

p5 <- ggplot(mtx3_melt, aes(x=SampleType,y=value)) + geom_jitter(height = 0, width = 0.1) + 
  theme(axis.text.x = element_text(angle=90,hjust = 1)) + 
  geom_signif(comparisons = list(c("Primary","Metastasis")), test = "wilcox.test") +
  facet_grid(~Var1) + labs(x="", y="log2(TPM+1)")
pdf("output/Plots/02_Expression_Log2TPM_oxiphos.pdf", width = 20, height = 4)
p5
dev.off()

#----- mofitt genes -----#
genes <- c(
  "VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1", 
  "KRT15", "CTSV", "DHRS9", "AREG", "CST6", "SERPINB3", "KRT6C", 
  "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7", "KRT17", 
  "GPR87", "TNS4", "SLC2A1", "ANXA8L1", "BTNL8", "FAM3D", "PRR15L", 
  "AGR3", "CTSE", "TMEM238L", "LYZ", "TFF2", "TFF1", "ANXA10", 
  "LGALS4", "PLA2G10", "CEACAM6", "VSIG2", "TSPAN8", "ST6GALNAC1", 
  "AGR2", "TFF3", "CYP3A7", "MYO1A", "CLRN3", "KRT20", "CDH17", 
  "SPINK4", "REG4"
)

idx4 <- which(rowData(se)$gene_name %in% genes)
mtx4 <- assays(se)$log2tpm[idx4,]
rownames(mtx4) <- rowData(se)$gene_name[idx4]
colnames(mtx4) <- se$SampleID
mtx4 <- mtx4[genes,]
mtx4_z <- t(scale(t(mtx4)))

col_fun2 <- colorRamp2(c(0,2,4,6,8,10,12,14), viridis(8, option = "A"))

ht6 <- Heatmap(mtx4, name = "log2(TPM+1)", bottom_annotation = ha1, col = col_fun2,
               row_split = c(rep("basal-like",25),rep("classical",25)),
               show_row_names = T, cluster_columns = fh, cluster_rows = fh)
p6 <- draw(ht6)
pdf("output/Plots/02_Heatmap_Log2TPM_MoffittSignature.pdf", width = 6, height = 10)
p6
dev.off()

#--------- Difference between primary vs meta ---------#
source("code/edgeR_PairwiseFunction.R")
DiffTest <- edgeR_pairwise(se, compareCol = "SampleType", topGroup = "Metastasis", bottomGroup = "Primary")

df7 <- assay(DiffTest) %>% as.data.frame
df7$logFDR <- -log10(df7$FDR)
df7$SYMBOL <- rowData(DiffTest)$gene_sym

df7$Signif <- "NS"
df7$Signif[which(df7$log2FC > 0 & df7$FDR < 0.05)] <- "Meta_UP"
df7$Signif[which(df7$log2FC < 0 & df7$FDR < 0.05)] <- "Meta_DN"
df7$Signif <- factor(df7$Signif, levels = c("NS","Meta_UP","Meta_DN"))
df7 <- df7[order(df7$Signif),]

p7 <- ggplot(df7, aes(x=log2FC,y=logFDR,color=Signif)) + geom_point_rast() + theme_cb() +
  scale_color_manual(values = c("NS"="gray","Meta_UP"="red","Meta_DN"="blue")) +
  geom_vline(xintercept = 0) + geom_vline(xintercept = c(-1,1), lty="dashed")
pdf("output/Plots/02_Volcano_Meta_vs_Ctrl.pdf", width = 6, height = 4)
p7
dev.off()

write.csv(df7, "output/Tables/02_Volcano_Meta_vs_Ctrl.csv")
