#----------------------------------------------------------------------------
# 00_preprocess.R
#----------------------------------------------------------------------------
library(Seurat)
library(parallel)

#function
makeSeuratObj <- function(
    Dir_10x = NULL,
    SampleName = NULL
){
  data <- Read10X(Dir_10x)
  out <- CreateSeuratObject(counts = data$`Gene Expression`)
  
  #add %mitochondrial RNA 
  out[["percent.mt"]] <- PercentageFeatureSet(out, pattern = "^MT-")
  #rename samples
  out$orig.ident <- SampleName
  out <- RenameCells(out, new.names = paste0(SampleName, "#", colnames(out)))
  
  return(out)
}

#sample assignment
sampleName <- list.files("data/")
names(sampleName) <- sampleName

#make Seurat object
pre_seu_ls <- mclapply(seq_along(sampleName), function(x){
  a <- as.character(sampleName[x])
  b <- names(sampleName)[x]
  out <- makeSeuratObj(Dir_10x = paste0("data/", a), SampleName = b)
  return(out)
}, mc.cores = 4)

#quality filter
pre_seu_ls <- mclapply(pre_seu_ls, function(x) x[,which(x$nFeature_RNA > 500 & x$nFeature_RNA < 6000 & x$percent.mt < 10)])
seu <- merge(pre_seu_ls[[1]], pre_seu_ls[-1])
seu <- JoinLayers(seu)


saveRDS(seu, "rds/00_seu.rds")
