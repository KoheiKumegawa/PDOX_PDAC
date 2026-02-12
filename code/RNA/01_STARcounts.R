#-------------------------------------------------------------------------------
# 01_STARcounts.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(stringr)
countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#----- generate count matrix (STAR counts) -----#
sampleName <- list.files("data/")
tmp <- str_split(sampleName, "-", simplify = T)[,2]
tmp <- str_split(tmp, "_S", simplify = T)[,1]
names(sampleName) <- tmp

counts_ls <- mclapply(sampleName, function(i){
  out <- fread(paste0("data/", i))
  out <- data.frame(gene = out$V1, count = out$V2)[-c(1:4),]
  colnames(out) <- c("gene", gsub("_ReadsPerGene.out.tab","",i))
  return(out)
}, mc.cores = 12)
counts <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), counts_ls)

#checking data
counts[c(1:5), c(1:5)]
which(is.na(counts))

rownames(counts) <- counts$gene
counts <- counts[,-1]

#----- Sample Info (colData) curation -----#
SampleData <- data.frame(row.names=colnames(counts),Sample=colnames(counts),SampleID=names(counts_ls))
SampleData$SampleType  <- c(rep("Primary",7), rep("Metastasis",5))

#change order
counts <- counts[,rownames(SampleData)]

#----- Feature Info (rowData) curation -----#
gff_anno <- rtracklayer::import.gff("../../ref/gencode.v36.annotation.gtf")

geneData <- DataFrame(row.names = gff_anno$gene_id,
                      gene_type = gff_anno$gene_type,
                      gene_name = gff_anno$gene_name)
geneData <- geneData[!duplicated(rownames(geneData)),]
geneData <- geneData[rownames(counts),]

#----- calculate TPM -----#
#gene length: https://www.biostars.org/p/83901/
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("../../ref/gencode.v36.annotation.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene") # then collect the exons per gene id
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene))) # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then

intersect(names(exonic.gene.sizes), rownames(counts)) %>% length #60660, all genes
len <- exonic.gene.sizes[rownames(counts)] %>% as.numeric()
tpms <- countToTpm(counts, len)

#----- finishing RNA-seq SummarizedExperiment -----#
se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           colData = DataFrame(SampleData),
                           rowData = geneData)
saveRDS(se, "rds/01_se.rds")
