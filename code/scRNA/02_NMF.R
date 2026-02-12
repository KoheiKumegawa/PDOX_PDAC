#----------------------------------------------------------------------------
# 02_NMF.R
#----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(stringr)
library(ggplot2)
library(ggsignif)
library(GeneOverlap)
library(Seurat)
library(parallel)
library(ggrastr)
"%ni%" <- Negate("%in%")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#------------ Summarizing NMF results ------------#
OptimalK <- c(
  "OP001"=6,
  "OP012"=6,
  "OP014"=8,
  "OP016"=6,
  "OP018"=6,
  "OP019"=8,
  "OP021"=5,
  "OP022"=7,
  "OP025A"=8,
  "OP026A"=5,
  "OP030"=5,
  "OP031A"=6,
  "OP033A"=7,
  "OP037"=5,
  "OP037H"=6,
  "OP038P"=7,
  "OP039A"=7
)

spectra <- lapply(names(OptimalK), function(x){
  i <- as.numeric(OptimalK[x])
  tmp <- fread(paste0("output/cNMF/result1/",x,"/", x, ".gene_spectra_score.k_",i,".dt_0_2.txt")) %>% as.matrix
  rownames(tmp) <- tmp[,1]
  tmp <- tmp[,-1]
  tmp <- t(tmp)
  out <- apply(tmp, 2, function(n) rownames(tmp)[order(n, decreasing = T)[c(1:200)]])
  colnames(out) <- paste0(x, "_state", c(1:ncol(out)))
  return(out)
})
spectra_mtx <- do.call(cbind,spectra)

jaccard_mtx <- apply(spectra_mtx, 2, function(x) apply(spectra_mtx,2,function(y) jaccard(x,y)))
col_fun1 <- colorRamp2(seq(0,0.4,0.05), c("white", rev(viridis(8, option = "A"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht1 <- Heatmap(jaccard_mtx, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh, 
               #bottom_annotation = ha1, right_annotation = ra1,
               row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
               use_raster = T)
p1 <- draw(ht1)
write.csv(spectra_mtx, "output/Tables/02_spectra_all.csv")
write.csv(spectra_mtx[,row_order(p1)], "output/Tables/02_spectra_all_ordered.csv")

pdf("output/Plots/02_Heatmap_NMF_Jaccard.pdf",width = 9, height = 8)
p1
dev.off()

#compare Hwang et al
Hwang_list <- read.csv("ref/Hwang.csv")
jaccard_mtx2 <- apply(spectra_mtx, 2, function(x) apply(Hwang_list,2,function(y) jaccard(x,y)))

ht2 <- Heatmap(jaccard_mtx2, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh, 
               #bottom_annotation = ha1, right_annotation = ra1,
               row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
               use_raster = T)
p2 <- draw(ht2)
pdf("output/Plots/02_Heatmap_NMF_Jaccard_vsHwang.pdf",width = 9, height = 8)
p2
dev.off()

#------------ Recurrent programs ------------#
RP <- list(c(1:19), c(20:36), c(37:53), c(54:58), c(73:79), c(80:85), c(86:93), c(94:108))
names(RP) <- paste0("RP", c(1:8))

spectra_mtx_ordered <- spectra_mtx[,row_order(p1)]

RP_genes <- lapply(RP, function(x) names(which(table(spectra_mtx_ordered[,x]) >= length(x)*0.33)))

source("code/enrichment_func.R")
gmt_files <- paste0("ref/AnnoDB/", list.files("ref/AnnoDB/", pattern = ".gmt"))
DIR <- "output/clusterProfiler/RP"

for (i in 1:length(gmt_files)){
  skip_to_next <- FALSE
  gmt_file <- gmt_files[i]
  tryCatch(EnrichmentAnalysis(gmt_file, RP_genes, DIR),
           error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}

saveRDS(RP_genes, "rds/02_RP_genes.rds")

jaccard_mtx3 <- lapply(RP_genes, function(x) apply(Hwang_list,2,function(y) jaccard(x,y)))
jaccard_mtx3 <- do.call(rbind, jaccard_mtx3)

ht3 <- Heatmap(jaccard_mtx3, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh)
p3 <- draw(ht3)
pdf("output/Plots/02_Heatmap_NMF_Jaccard_RPvsHwang.pdf",width = 5, height = 4)
p3
dev.off()

#Visualize
seu <- readRDS("rds/01_seu.rds")
RP_genes <- readRDS("rds/02_RP_genes.rds")

seu <- AddModuleScore(seu, features = RP_genes, name = "RP")

umap_data <- as.data.frame(seu@reductions$umap@cell.embeddings)
umap_data$cluster <- seu$seurat_clusters
umap_data$sample <- seu$orig.ident
ashape_obj <- readRDS("rds/01_ashape_obj.rds")
edges <- as.data.frame(ashape_obj$edges)
edges_df <- data.frame(
  x = umap_data$umap_1[edges$ind1],
  y = umap_data$umap_2[edges$ind1],
  xend = umap_data$umap_1[edges$ind2],
  yend = umap_data$umap_2[edges$ind2]
)

p4 <- mclapply(c(1:8), function(i){
  x <- paste0("RP",i)
  df_expr <- umap_data
  df_expr$expr <- as.numeric(seu[[]][,x])
  df_expr <- df_expr[order(df_expr$expr),]
  out <- ggplot() +
    geom_point_rast(data = df_expr, aes(x = umap_1, y = umap_2, color = expr), alpha = 1, size = 0.1) +
    scale_color_gradientn(colours = c("white", rev(viridis(6, option = "G")))) +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black") + theme_cb() + 
    ggtitle(names(RP_genes)[i])
  return(out)
}, mc.cores = 8)
pdf("output/Plots/02_UMAPOL_ModuleScore_RP.pdf", width = 6, height = 6)
p4
dev.off()

cluster_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
                    "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D",
                    "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A")

seu$Sample <- stringr::str_split(seu$orig.ident, "_", simplify = T)[,1]
seu$Sample <- factor(seu$Sample, levels　=　rev(c("OP001","OP012","OP014","OP016","OP018","OP019","OP021","OP022","OP030","OP037",
                                                    "OP025A","OP026A","OP031A","OP033A","OP037H","OP038P","OP039A")))

p5 <- lapply(c(1:8), function(i)
  RidgePlot(seu, features = paste0("RP",i), group.by = "Sample", cols = c(rep("red",7), rep("gray", 10))) +
  labs(x="Module Score", y="Sample"))
pdf("output/Plots/02_Ridge_ModuleScore_RP_v2.pdf", width = 6, height = 6)
p5
dev.off()
lapply(names(RP_genes), function(x){
  write.table(data.frame(gene=RP_genes[[x]]), paste0("output/Tables/02_geneList_",x), col.names = F, row.names = F, quote = F)
})

#------------ Recurrent programs Violin ------------#
df6 <- seu[[]][, c("Sample", paste0("RP",c(1:8)))]
df6_melt <- reshape2::melt(df6)

df6_melt$Sample <- factor(df6_melt$Sample,
                          levels　=　c("OP001","OP012","OP014","OP016","OP018","OP019","OP021","OP022","OP030","OP037",
                                     "OP025A","OP026A","OP031A","OP033A","OP037H","OP038P","OP039A"))
ggplot(df6_melt, aes(x=Sample,y=value,fill=variable)) + 
  geom_violin() + geom_boxplot(outlier.shape = NA, width = 0.25) +
  facet_grid(~variable)
  

p7 <-lapply(c(1:8), function(i){
  tmp <- sapply(c("OP001","OP012","OP014","OP016","OP018","OP019","OP021","OP022","OP030","OP037",
                  "OP025A","OP026A","OP031A","OP033A","OP037H","OP038P","OP039A"), function(x) median(df6[,paste0("RP",i)][which(df6$Sample==x)]))
  df7 <- data.frame(score=tmp, SampleType=factor(c(rep("Primary",10),rep("Metastasis",7)), levels = c("Primary", "Metastasis")))
  out <- ggplot(df7, aes(x=SampleType,y=score)) + geom_jitter(height = 0, width = 0.2, size = 4, alpha=0.5) +
    geom_signif(comparisons = list(c("Primary","Metastasis")), test = "wilcox.test") +
    stat_summary(fun = median, geom = "crossbar", width = 1, color = "black") + theme_cb() +
    theme(axis.text.x = element_text(angle=90,hjust = 1)) + ggtitle(paste0("RP",i))
  return(out)
})
pdf("output/Plots/02_Jitter_ModuleScore_RP.pdf", width = 2, height = 4)
p7
dev.off()

ggplot(df7, aes(x=SampleType,y=score)) + geom_jitter(height = 0, width = 0.2, size = 4, alpha=0.5) +
  geom_signif(comparisons = list(c("Primary","Metastasis")), test = "wilcox.test") +
  stat_summary(fun = median, geom = "crossbar", width = 1, color = "black") + theme_cb() +
  theme(axis.text.x = element_text(angle=90,hjust = 1))

