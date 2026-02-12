#----------------------------------------------------------------------------
# 01_clustering.R
#----------------------------------------------------------------------------
library(ArchR)
library(ggrastr)
addArchRThreads(threads = 16)
addArchRGenome("hg19")

#---------- clustering ----------#
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI") 
arc <- addClusters(arc, reducedDims = "IterativeLSI", force = T)
arc <- addUMAP(arc, reducedDims = "IterativeLSI", force = T)

sample_colors <- RColorBrewer::brewer.pal(4, "Dark2")[c(3:4)] %>% `names<-`(.,c("P","A"))
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 3) 
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", pal = sample_colors, size = 3)
plotPDF(p1, p2, name = "01_UMAP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#Sample x Cluster matrix
cM <- table(arc$Sample, arc$Clusters)[, paste0("C", seq_along(unique(arc$Clusters)))]
write.csv(cM, "output/Tables/01_SampleClusterMatrix.csv")

#---------- gene activity ----------#
markers1 <- c(
  "VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1", 
  "KRT15", "CTSV", "DHRS9", "AREG", "CST6", "SERPINB3", "KRT6C", 
  "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7", "KRT17", 
  "GPR87", "TNS4", "SLC2A1", "ANXA8L1", "BTNL8", "FAM3D", "PRR15L", 
  "AGR3", "CTSE", "TMEM238L", "LYZ", "TFF2", "TFF1", "ANXA10", 
  "LGALS4", "PLA2G10", "CEACAM6", "VSIG2", "TSPAN8", "ST6GALNAC1", 
  "AGR2", "TFF3", "CYP3A7", "MYO1A", "CLRN3", "KRT20", "CDH17", 
  "SPINK4", "REG4"
) %>% intersect(., getFeatures(arc))

arc <- addImputeWeights(arc)
p3 <- plotEmbedding(arc,
                    colorBy = "GeneScoreMatrix", 
                    name = markers1, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
plotPDF(p3, name = "01_UMAPOL_GS_marker1.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

p4 <- plotEmbedding(arc,
                    colorBy = "GeneScoreMatrix", 
                    name = c("SOX2","ZEB1","ZEB2","SNAI1","SNAI2","VIM"), 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
plotPDF(p4, name = "01_UMAPOL_GS_marker2.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#GEP module score
GEPs <- read.csv("../mouse_scRNA_v2/output/Tables/02_spectra_T_k6.csv", header = T, row.names = 1)
GEPs <- lapply(c(1:6), function(i) intersect(GEPs[,i], getFeatures(arc))[c(1:100)])

arc <- addModuleScore(arc, useMatrix = "GeneScoreMatrix", features = GEPs)
p5 <- plotEmbedding(arc,
                    colorBy = "cellColData", 
                    name = paste0("Module", c(1:6)), 
                    embedding = "UMAP", pal = ArchRPalettes$comet,
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
plotPDF(p5, name = "01_UMAPOL_ModuleScore_GEP_v2.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

p6 <- plotEmbedding(arc,
                    colorBy = "cellColData", 
                    name = paste0("Module", c(1:6)), 
                    embedding = "UMAP", pal = ArchRPalettes$comet,
                    imputeWeights = NULL,
                    plotAs = "points", size = 3)
plotPDF(p6, name = "01_UMAPOL_ModuleScore_GEP_v2_NoImpute.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#---------- peak call ----------#
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters") 
arc <- addReproduciblePeakSet(arc, groupBy = "Clusters", pathToMacs2 = pathToMacs2, method = "q", cutOff = 0.05)
arc <- addPeakMatrix(arc)

#calculate motif deviation
arc <- addMotifAnnotations(arc, motifSet = "homer", annoName = "homer")
arc <- addMotifAnnotations(arc, motifSet = "cisbp", annoName = "cisbp")
arc <- addBgdPeaks(arc) %>% addDeviationsMatrix(., peakAnnotation = "homer") %>% addDeviationsMatrix(., peakAnnotation = "cisbp")

saveRDS(arc, "rds/01_arc.rds")

#---------- Motif score ----------#
sort(names(arc@peakAnnotation$cisbp$motifs))
sort(names(arc@peakAnnotation$homer$motifs))

p7 <- plotEmbedding(arc, 
                    colorBy = "cisbpMatrix", 
                    name = paste0("z:", c("ZEB1_157","SOX2_769","SNAI1_199","SNAI2_161","CTCF_177","MYC_55")), 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 4)
p8 <- plotEmbedding(arc,
                    colorBy = "cisbpMatrix",
                    name = paste0("z:", c("SOX2_769","STAT3_777","SMAD3_743")),
                    embedding = "UMAP",
                    imputeWeights = getImputeWeights(arc),
                    plotAs = "points", size = 4)
plotPDF(p7, name = "01_UMAPOL_MotifScore_v1.pdf", ArchRProj = arc, addDOC = FALSE, width = 4, height = 4)
plotPDF(p8, name = "01_UMAPOL_MotifScore_v2.pdf", ArchRProj = arc, addDOC = FALSE, width = 4, height = 4)

#---------- Correlation ----------#
ModuleScoreMat <- t(as.matrix(arc@cellColData[,paste0("Module", c(1:6))]))
ModuleScoreMat <- imputeMatrix(mat = ModuleScoreMat, imputeWeights = getImputeWeights(arc))
ModuleScoreMat <- t(ModuleScoreMat)

cisbp_se <- getMatrixFromProject(arc, useMatrix = "cisbpMatrix")
homer_se <- getMatrixFromProject(arc, useMatrix = "homerMatrix")

cisbp_mat <- imputeMatrix(assays(cisbp_se)$z, getImputeWeights(arc))
cisbp_mat <- t(cisbp_mat)
homer_mat <- imputeMatrix(assays(homer_se)$z, getImputeWeights(arc))
homer_mat <- t(homer_mat)

correlation_result1 <- lapply(c(1:6), function(i){
  value1 <- ModuleScoreMat[,i]
  tmp1 <- sapply(colnames(homer_mat), function(x) cor(value1, homer_mat[,x]))
  tmp2 <- sapply(colnames(cisbp_mat), function(x) cor(value1, cisbp_mat[,x]))
  
  names(tmp1) <- paste0("homer:", names(tmp1))
  names(tmp2) <- paste0("cisbp:", names(tmp2))
  
  out <- c(tmp1, tmp2)
  out <- out[order(names(out))]
  return(out)
})
correlation_result1 <- do.call(cbind, correlation_result1)
colnames(correlation_result1) <- paste0("GEP", c(1:6))

write.csv(correlation_result1, "output/Tables/01_correlation_ModuleScore_vs_MotifScore.csv")

#---------- DEG Correlation ----------#
# arc <- readRDS("rds/01_arc.rds")
DEG <- read.csv("../mouse_scRNA_v1/output/Tables/01b_DEG.csv")
DEG_ls <- lapply(c(0:10), function(x) DEG$gene[DEG$avg_log2FC > 1 & DEG$p_val_adj < 0.01 & DEG$cluster == x])
names(DEG_ls) <- paste0("cluster", c(0:10))
saveRDS(DEG_ls, "rds/01_scRNA_DEG_ls.rds")

DEG_ls <- lapply(DEG_ls, function(x) intersect(x, getFeatures(arc)))

arc <- addModuleScore(arc, useMatrix = "GeneScoreMatrix", features = DEG_ls)
p9 <- plotEmbedding(arc,
                    colorBy = "cellColData", 
                    name = paste0("Module.cluster", c(0:10)), 
                    embedding = "UMAP", pal = ArchRPalettes$comet,
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
plotPDF(p9, name = "01_UMAPOL_ModuleScore_DEG.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)


ModuleScoreMat2 <- t(as.matrix(arc@cellColData[,paste0("Module.cluster", c(0:10))]))
ModuleScoreMat2 <- imputeMatrix(mat = ModuleScoreMat2, imputeWeights = getImputeWeights(arc))
ModuleScoreMat2 <- t(ModuleScoreMat2)

correlation_result2 <- lapply(c(1:11), function(i){
  value1 <- ModuleScoreMat2[,i]
  tmp1 <- sapply(colnames(homer_mat), function(x) cor(value1, homer_mat[,x]))
  tmp2 <- sapply(colnames(cisbp_mat), function(x) cor(value1, cisbp_mat[,x]))
  
  names(tmp1) <- paste0("homer:", names(tmp1))
  names(tmp2) <- paste0("cisbp:", names(tmp2))
  
  out <- c(tmp1, tmp2)
  out <- out[order(names(out))]
  return(out)
})
correlation_result2 <- do.call(cbind, correlation_result2)
colnames(correlation_result2) <- paste0("cluster", c(0:10))

write.csv(correlation_result2, "output/Tables/01_correlation_ModuleScoreDEG_vs_MotifScore.csv")


library(ComplexHeatmap)
correlation_result2_cisbp <- correlation_result2[c(1:870),]

fh <- function(x) hclust(dist(x), method = "ward.D2")
row_labels <- ifelse(rownames(correlation_result2_cisbp) == "cisbp:SOX2_769", rownames(correlation_result2_cisbp), "")

ht1 <- Heatmap(correlation_result2_cisbp, cluster_rows = fh, cluster_columns = fh, name="Correlation",
               row_names_gp = gpar(fontsize = 8), row_labels = row_labels)
p10 <- draw(ht1)
pdf("output/Plots/01_Heatmap_Correlation_MSDEG_CISBP.pdf", width = 5, height = 6)
p10
dev.off()

row_labels[row_order(p10)]
write.csv(correlation_result2_cisbp[row_order(p10),column_order(p10)], "output/Tables/01_correlation_ModuleScoreDEG_vs_MotifCisbp_ordered.csv")

df11 <- data.frame(module4=ModuleScoreMat2[,"Module.cluster4"], sox2motif=cisbp_mat[,"SOX2_769"])
plot(df11)
p11 <- ggplot(df11, aes(x=module4, y=sox2motif)) + geom_point_rast(alpha=0.5, size=1.5) + theme_ArchR() +
  geom_smooth(se=F, method = "lm", color="red") + labs(x="Cluster4 DEG ModuleScore", y="SOX2 MotifScore")
pdf("output/Plots/01_Scatterplot_Correlation_Module4_SOX2.pdf", width = 4, height = 4)
p11
dev.off()

#all motifs vs cluster4 DEG
p12 <- lapply(colnames(cisbp_mat), function(x){
  df_tmp <- data.frame(modulescore=ModuleScoreMat2[,"Module.cluster4"], motifscore=cisbp_mat[,x])
  out <- ggplot(df_tmp, aes(x=modulescore, y=motifscore)) + geom_point_rast(alpha=0.5, size=1.5) + theme_ArchR() +
    geom_smooth(se=F, method = "lm", color="red") + labs(x="Cluster4 DEG ModuleScore", y="MotifScore") +
    ggtitle(x)
  return(out)
})
pdf("output/Plots/01_Scatterplot_Correlation_Module4_Allmotifs.pdf", width = 4, height = 4)
p12
dev.off()


#### gene score
arc <- addImputeWeights(arc)
p13 <- plotEmbedding(arc,
                    colorBy = "GeneScoreMatrix", 
                    name = DEG_ls$cluster1, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
p14 <- plotEmbedding(arc,
                    colorBy = "GeneScoreMatrix", 
                    name = DEG_ls$cluster4, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 3)
plotPDF(p14, name = "01_UMAPOL_GS_DEG_cluster4.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
plotPDF(p13, name = "01_UMAPOL_GS_DEG_cluster1.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#---------- QC metrics ----------#
p15.1 <- plotGroups(arc, groupBy = "Sample", name = "TSSEnrichment", plotAs = "violin", pal = rep("gray",2)) 
p15.2 <- plotGroups(arc, groupBy = "Sample", name = "nFrags", plotAs = "violin", pal = rep("gray",2))
plotPDF(p15.1,p15.2, name = "01_VlnPlot_QC.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
