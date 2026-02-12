#----------------------------------------------------------------------------
# 01_clustering.R
#----------------------------------------------------------------------------
#Seurat5
library(Seurat)
library(parallel)
library(ggplot2)
library(ggrastr)
library(scales)
library(dplyr)
library(alphahull)
library(viridisLite)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

seu <- readRDS("rds/02_seu.rds")
seu$Sample <- factor(seu$orig.ident, levels=c("xp026A-69P","xp026A-83P","xp026A-75P","xp026A-75A"))

#--------- UMAP plot ----------#
cluster_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
                    "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D",
                    "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A")
sample_colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% `names<-`(.,c("xp026A-69P","xp026A-83P","xp026A-75P","xp026A-75A"))

umap_data <- as.data.frame(seu@reductions$umap@cell.embeddings)
umap_data$cluster <- seu$seurat_clusters
umap_data$sample <- seu$Sample

#shaping
ashape_obj <- ashape(umap_data$umap_1, umap_data$umap_2, alpha = 0.16)
edges <- as.data.frame(ashape_obj$edges)
edges_df <- data.frame(
  x = umap_data$umap_1[edges$ind1],
  y = umap_data$umap_2[edges$ind1],
  xend = umap_data$umap_1[edges$ind2],
  yend = umap_data$umap_2[edges$ind2]
)

#centroids
centroids <- umap_data %>%
  group_by(cluster) %>%
  summarise(
    cx = mean(umap_1),
    cy = mean(umap_2)
  )

#plotting
p1 <- ggplot() +
  geom_point_rast(data = umap_data, aes(x = umap_1, y = umap_2, color = cluster), alpha = 0.5, size = 2) +
  scale_color_manual(values = cluster_colors) +
  geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") +
  geom_text(data = centroids, aes(x = cx, y = cy, label = cluster), 
            color = "black", size = 4, fontface = "bold") + theme_cb()
p2 <- ggplot() +
  geom_point_rast(data = umap_data, aes(x = umap_1, y = umap_2, color = sample), alpha = 0.3, size = 2) +
  scale_color_manual(values = sample_colors) +
  geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") + theme_cb()

pdf("output/Plots/01b_UMAP.pdf", width = 6, height = 5)
p1
p2
dev.off()
# 
# saveRDS(seu, "rds/01_seu.rds")
saveRDS(ashape_obj, "rds/01b_ashape_obj.rds")

#--------- UMAP Overlay ----------#
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

p3.1 <- mclapply(classical_markers, function(x){
  df_expr <- umap_data
  df_expr$expr <- seu[["RNA"]]$data[x,]
  df_expr <- df_expr[order(df_expr$expr),]
  out <- ggplot() +
    geom_point_rast(data = df_expr, aes(x = umap_1, y = umap_2, color = expr), alpha = 1, size = 0.1) +
    scale_color_gradientn(colours = c("white", rev(viridis(6, option = "A")))) +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") + theme_cb() + ggtitle(x)
  return(out)
}, mc.cores = 8)
p3.2 <- mclapply(basal_like_markers, function(x){
  df_expr <- umap_data
  df_expr$expr <- seu[["RNA"]]$data[x,]
  df_expr <- df_expr[order(df_expr$expr),]
  out <- ggplot() +
    geom_point_rast(data = df_expr, aes(x = umap_1, y = umap_2, color = expr), alpha = 1, size = 0.1) +
    scale_color_gradientn(colours = c("white", rev(viridis(6, option = "A")))) +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") + theme_cb() + ggtitle(x)
  return(out)
}, mc.cores = 8)


pdf("output/Plots/01b_UMAPOL_classicalMarkers.pdf", width = 4, height = 4)
p3.1
dev.off()
pdf("output/Plots/01b_UMAPOL_basalMarkers.pdf", width = 4, height = 4)
p3.2
dev.off()

#--------- Sample Mixture ----------#
df4 <- table(seu$seurat_clusters, seu$orig.ident)
df4 <- reshape2::melt(df4)
df4$Var1 <- factor(df4$Var1)

p4.1 <- ggplot(df4, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Cluster", y = "%Sample") + theme_cb() +
  theme(legend.key.size = unit(0.4, 'cm'))
p4.2 <- ggplot(df4, aes(x = Var2, y = value, fill = Var1)) + 
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = cluster_colors) + labs(x = "Sample", y = "%Cluster") + theme_cb() +
  theme(legend.key.size = unit(0.4, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/01b_Barplot_SampleCluster.pdf", width = 8, height = 4)
p4.1
p4.2
dev.off()

write.csv(table(seu$seurat_clusters, seu$orig.ident), "output/Tables/01b_Barplot_SampleCluster.csv")

#------------ DEG ------------#
DEG <- FindAllMarkers(seu, logfc.threshold = 0.25, only.pos = T)
write.csv(DEG, "output/Tables/01b_DEG.csv")

DEG_ls <- lapply(c(0:10), function(x) DEG$gene[DEG$avg_log2FC > 1 & DEG$p_val_adj < 0.01 & DEG$cluster == x])
names(DEG_ls) <- paste0("cluster", c(0:10))

source("code/enrichment_func.R")
gmt_files <- paste0("../scRNA_v1/ref/AnnoDB/", list.files("../scRNA_v1/ref/AnnoDB/", pattern = ".gmt"))
DIR <- "output/clusterProfiler/DEG"

for (i in 1:length(gmt_files)){
  skip_to_next <- FALSE
  gmt_file <- gmt_files[i]
  tryCatch(EnrichmentAnalysis(gmt_file, DEG_ls, DIR),
           error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}


oxphos <- c("NDUFA4", "NDUFB9", "NDUFA8", "NDUFC2", "NDUFA6", "NDUFB4", "NDUFB3", "NDUFA2", "NDUFB11",
  "UQCR10",
  "COX4I1", "COX6C", "COX7C", "COX6A1", "COX8A", "COX7B",
  "ATP5F1A", "ATP5F1E", "ATP5PB", "ATP5PD", "ATP5MD", "ATP5MPL",
  "VDAC1",
  "TIMM8B", "TOMM22", "TMEM70", "GHITM", "MPV17L", "MPV17L2",
  "CYCS")
oxphos <- intersect(oxphos, DEG_ls$cluster1)

p5 <- mclapply(oxphos, function(x){
  df_expr <- umap_data
  df_expr$expr <- seu[["RNA"]]$data[x,]
  df_expr <- df_expr[order(df_expr$expr),]
  out <- ggplot() +
    geom_point_rast(data = df_expr, aes(x = umap_1, y = umap_2, color = expr), alpha = 1, size = 0.1) +
    scale_color_gradientn(colours = c("white", rev(viridis(6, option = "A")))) +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") + theme_cb() + ggtitle(x)
  return(out)
}, mc.cores = 8)

pdf("output/Plots/01b_UMAPOL_OXPHOS.pdf", width = 4, height = 4)
p5
dev.off()

# DEG <- read.csv("output/Tables/01b_DEG.csv")
# DEG_ls <- lapply(c(0:10), function(x) DEG$gene[DEG$avg_log2FC > 1 & DEG$p_val_adj < 0.01 & DEG$cluster == x])
# names(DEG_ls) <- paste0("cluster", c(0:10))

df6 <- data.frame(cluster = factor(names(DEG_ls),levels = paste0("cluster",c(0:10))), DEGNO = sapply(DEG_ls, length))
p6 <- ggplot(df6, aes(x=cluster,y=DEGNO)) + geom_bar(stat="identity") + theme_cb() + 
  theme(axis.text.x = element_text(angle=90,hjust = 1)) + labs(x="",y="# DEG")
pdf("output/Plots/01b_Barplot_DEG_Numbers.pdf", height = 3, width = 3.5)
p6
dev.off()


p7.1 <- mclapply(DEG_ls$cluster0, function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
         scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)
p7.2 <- mclapply(DEG_ls$cluster1, function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
                 scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)
p7.3 <- mclapply(DEG_ls$cluster2, function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
                 scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)
p7.4 <- mclapply(DEG_ls$cluster3, function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
                 scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)
p7.5 <- mclapply(DEG_ls$cluster4, function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
                  scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)

pdf("output/Plots/01b_UMAPOL_Expr_DEG_Cluster0.pdf", height = 3, width = 3.5)
p7.1
dev.off()
pdf("output/Plots/01b_UMAPOL_Expr_DEG_Cluster1.pdf", height = 3, width = 3.5)
p7.2
dev.off()
pdf("output/Plots/01b_UMAPOL_Expr_DEG_Cluster2.pdf", height = 3, width = 3.5)
p7.3
dev.off()
pdf("output/Plots/01b_UMAPOL_Expr_DEG_Cluster3.pdf", height = 3, width = 3.5)
p7.4
dev.off()
pdf("output/Plots/01b_UMAPOL_Expr_DEG_Cluster4.pdf", height = 3, width = 3.5)
p7.5
dev.off()

#------------ EMT markers ------------#
p8 <- mclapply(c("SOX2","ZEB1","ZEB2","SNAI1","SNAI2","FOXC1","VIM"), function(x) FeaturePlot(seu, features = x, pt.size = 1, order = T) + 
                   scale_colour_gradientn(colours = rev(viridis(256, option = "A"))) + theme_cb(), mc.cores = 8)
pdf("output/Plots/01b_UMAPOL_Expr_EMT.pdf", height = 3, width = 3.5)
p8
dev.off()

p9 <- VlnPlot(seu, features = c("SOX2","ZEB1","ZEB2","SNAI1","SNAI2","FOXC1","VIM"), group.by = "Sample", raster = F)
pdf("output/Plots/01b_VlnPlot_Expr_EMT.pdf", height = 8, width = 6)
p9
dev.off()

avg <- AverageExpression(seu, return.seurat = T, group.by = "Sample")
mtx <- avg@assays$RNA$scale.data[c("ZEB1","ZEB2","SNAI1","SNAI2","FOXC1","VIM"),]

library(ComplexHeatmap)
library(circlize)
col_fun1 <- colorRamp2(c(-2,-1,0,1,2), rev(viridis(5, option = "A")))
fh = function(x) hclust(dist(x), method="ward.D2")
ht1 <- Heatmap(mtx, name = "Relative expression", cluster_columns = F, cluster_rows = fh, 
               show_row_dend = T, show_row_names = T, 
               col = col_fun1)
p10 <- draw(ht1)
pdf("output/Plots/01b_Heatmap_Expr_EMT.pdf", height = 4, width = 5)
p10
dev.off()

#------------ QC metrics ------------#
p11.1 <- VlnPlot(seu, group.by = "Sample", cols = rep("gray", 4), pt.size = 0, features = "nFeature_RNA") + theme(axis.text.x = element_text(angle=90, hjust = 1)) + NoLegend()
p11.2 <- VlnPlot(seu, group.by = "Sample", cols = rep("gray", 4), pt.size = 0, features = "nCount_RNA") + theme(axis.text.x = element_text(angle=90, hjust = 1)) + NoLegend()
p11.3 <- VlnPlot(seu, group.by = "Sample", cols = rep("gray", 4), pt.size = 0, features = "percent.mt") + theme(axis.text.x = element_text(angle=90, hjust = 1)) + NoLegend()
pdf("output/Plots/01b_VlnPlot_QC.pdf", width = 3, height = 3)
p11.1
p11.2
p11.3
dev.off()


