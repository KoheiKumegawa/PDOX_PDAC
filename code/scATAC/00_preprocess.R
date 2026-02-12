#----------------------------------------------------------------------------
# 00_preprocess.R
#----------------------------------------------------------------------------
library(ArchR)
addArchRThreads(threads = 16)
addArchRGenome("hg19")

#assign sample tsv files
sampleName <- list.files(path = "data/", pattern = "tsv.gz")
names(sampleName) <- c("A","P")
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("data//", outFile), 
                               sampleNames = names(sampleName),
                               minTSS = 4, 
                               minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               force = TRUE)

#infer doublets
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#make pre-filtered ArchRProject
pre_arc <- ArchRProject(ArrowFiles, outputDirectory = "output", copyArrows = F)
pre_arc <- filterDoublets(pre_arc)

#quality metrics
arc <- pre_arc[which(pre_arc$TSSEnrichment > 4 & pre_arc$nFrags > 1500)]

table(pre_arc$Sample)
# A    P 
# 1262  749 
table(arc$Sample)
# A   P 
# 789 508

#save rds
saveRDS(pre_arc, "rds/00_pre_arc.rds")
saveRDS(arc, "rds/00_arc.rds")
