library(broman)

setwd("/groups/umcg-fg/tmp04/projects/genenetwork/recount3")

load("/groups/umcg-fg/tmp04/projects/genenetwork/recount3/DataForPredictions.RData")

str(pcsAndMeta$excludeBasedOnPredictionCancer)


table(pcsAndMeta$Tissue == "Breast", pcsAndMeta$excludeBasedOnPredictionCancer)
table(pcsAndMeta$Tissue == "Brain", pcsAndMeta$excludeBasedOnPredictionCancer)

table(pcsAndMeta$Tissue == "Blood", pcsAndMeta$excludeBasedOnPredictionCancer)
table(pcsAndMeta$Tissue == "Bone Marrow", pcsAndMeta$excludeBasedOnPredictionCancer)
table(pcsAndMeta$Tissue == "Blood Vessel", pcsAndMeta$excludeBasedOnPredictionCancer)
      
brainCancer <- pcsAndMeta$Row.names[pcsAndMeta$Tissue == "Brain" & pcsAndMeta$excludeBasedOnPredictionCancer]

str(brainCancer)





sraFiles <- list.files(path="rse-sra", pattern="*.rda", full.names=TRUE, recursive=FALSE)
gtexFiles <- list.files(path="rse-gtex", pattern="*.rda", full.names=TRUE, recursive=FALSE)
allFiles <- c(sraFiles, gtexFiles, "rse-tcga/rseTCGA.rda", "rse-tcga/rse_ESCA_TCGA.rda")
str(allFiles)

perChunkExp <- sapply(allFiles, function(file){
  
  loadedObject <- load(file)
  
  sreObjects <- get(loadedObject[1])
  
  #sometimes single RSE is not in list. Put in list of one to make code uniform
  if(!is.list(sreObjects)){
    sreObjects <- list(sreObjects)
  }
  
  #sreObject <- sreObjects[[1]]
  
  perStudyExp <- lapply(sreObjects, function(sreObject){
    studyExp <- sreObject@assays@data@listData$raw_counts
    return(studyExp[,colnames(studyExp) %in% brainCancer, drop = F])
  })
  
  return(do.call(cbind, perStudyExp))
  
})

str(perChunkExp)

selectedSamplesExp <- do.call(cbind, perChunkExp)
str(selectedSamplesExp)

uniqueSamplesIndex <- match(brainCancer, colnames(selectedSamplesExp))
selectedSamplesExp <- selectedSamplesExp[,uniqueSamplesIndex]

str(selectedSamplesExp)

#save.image("brainCancer.RData")



load("brainCancer.RData")

totalSamples <- ncol(selectedSamplesExp)


hist(apply(selectedSamplesExp,1,function(a){sum(a==0)}))
hist(apply(selectedSamplesExp,2,function(a){sum(a==0)}))

selectedSamplesExp2 <- selectedSamplesExp[apply(selectedSamplesExp, 1, function(a){
    return((sum(a==0) / totalSamples) < 0.5)
  }),]
str(selectedSamplesExp2)

selectedSamplesExpNorm <- log2(selectedSamplesExp2 + 1)

selectedSamplesExpNorm <- normalize(selectedSamplesExpNorm)
dimnames(selectedSamplesExpNorm) <- dimnames(selectedSamplesExp2)

geneMean <- rowMeans(selectedSamplesExpNorm)


load(file = "Metadata/combinedMeta_2022_09_15.RData", verbose = T)

#TCGA samples don't have sra.sample_spots instead use recount_qc.bc_frag.count
missingSampleSpots <- is.na(combinedMeta[,"sra.sample_spots"])
combinedMeta[missingSampleSpots,"sra.sample_spots"] <- combinedMeta[missingSampleSpots,"recount_qc.bc_frag.count"]

#TCGA and Gtex don't report layout but is all paired
missingLayout <- is.na(combinedMeta[,"sra.library_layout"])
combinedMeta[missingLayout,"sra.library_layout"] <- "paired"

covariatesToCorrectFor <- read.delim("CovariateNames.txt", header = F)$V1
combinedMetaSelection <- combinedMeta[brainCancer,covariatesToCorrectFor]
combinedMetaSelection$sra.library_layout <- as.factor(combinedMetaSelection$sra.library_layout)

selectedSamplesExpNormCovCor <- apply(selectedSamplesExpNorm, 1 ,function(geneExp, combinedMetaSelection){
  return(residuals(lm(geneExp ~ . ,data = combinedMetaSelection)))
}, combinedMetaSelection = combinedMetaSelection)

selectedSamplesExpNormCovCor <- t(selectedSamplesExpNormCovCor)

selectedSamplesExpNormCovCor <- selectedSamplesExpNormCovCor + geneMean


str(selectedSamplesExpNormCovCor)

boxplot(selectedSamplesExpNormCovCor[,1:10])
boxplot(t(selectedSamplesExpNormCovCor)[,1:10])

write.table(selectedSamplesExpNormCovCor, file = "brainCancerQcNormCovCor.txt", sep = "\t", quote = F, col.names = NA)
write.table(selectedSamplesExp, file = "brainCancerTpm.txt", sep = "\t", quote = F, col.names = NA)
