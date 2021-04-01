#!/bin/env/Rscript

# Lab Research - HLA-DRB1 ASE
# 20210328 Sehwan Chun at Corestem, Inc.
# 3.1. HTSeq reads to log2-scale TMM normalized HLA-DRB1 table

#### 1. Library Loading ####
library(edgeR)

#### 2. Files Loading ####
sampleFiles = grep("HTseq_output.txt", list.files(getwd()), value = TRUE) #HTSeq reads
TotalCount = read.table("~~HLA-DRB1 counts file~~", header = T) #pre-defined HLA-DRB1 counts

#### 3.setup HTSeq reads as table and replace HLA-DRB1 reads ####
y = read.table(sampleFiles[1], sep = '\t', header = FALSE, skip = 1)
for (i in (2:length(sampleFiles))) {
  tmp = read.table(sampleFiles[i], sep = '\t' , header = FALSE, skip = 1)
  y = cbind(y, tmp[,2])
}

colnames(y) = c("genes",substr(sampleFiles,1,7))
rownames(y) = y$genes
y = subset(y, select = - 1)

sizeGeneCounts = dim(y)
y = y[1:(sizeGeneCounts[1]-5),]
HtseqDRB = y["HLA-DRB1",]

y["HLA-DRB1",] = TotalCount$Pair
y["HLA-DRB1",] = TotalCount$All

#### 4.TMM normalization with log2-scale for HLA-DRB1 counts ####
yEdged = DGEList(counts = y)
yEdged = calcNormFactors(yEdged, method = "TMM")
tmmEXP = cpm(yEdged, log = TRUE)
NormDRBCounts = as.numeric(tmmEXP["HLA-DRB1",])
TotalCount[,4] = NormDRBCounts

yEdged = DGEList(counts = y)
yEdged = calcNormFactors(yEdged, method = "TMM")
tmmEXP = cpm(yEdged)
NormDRBCounts = as.numeric(tmmEXP["HLA-DRB1",])
TotalCount[,5] = NormDRBCounts

colnames(TotalCount)[4] = "TMMPair"
colnames(TotalCount)[5] = "TMMAll"

TotalCount[,6] = TotalCount[,4]/median(TotalCount$TMMPair)
TotalCount[,7] = TotalCount[,5]/median(TotalCount$TMMAll)

colnames(TotalCount)[6] = "TMMMedianPair"
colnames(TotalCount)[7] = "TMMMedianAll"

#### 5. Save results files --> to use numeric traits of cis-eQTL ####
write.table(TotalCount,"TMMCount_LogCPM_KGP",sep = "\t", quote = F, row.names = F)
