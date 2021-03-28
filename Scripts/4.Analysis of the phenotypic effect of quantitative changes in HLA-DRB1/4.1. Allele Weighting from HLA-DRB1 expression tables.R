# Lab Research - HLA-DRB1 ASE
# 20210328 Sehwan Chun at Corestem, Inc.
# 4.1. Allele Weighting from HLA-DRB1 expression tables 

#### 0.Func setup ####
#First Inference Missing Values
InferFunc = function(completeTable,threshold){
  for (i in 1:nrow(completeTable)){
    for (j in 1:nrow(completeTable)){
      if(i < j){
        if(is.na(completeTable[i,j])){
          inferValues = c()
          for (k in 1:nrow(completeTable)){
            inferValuesCandidate = (completeTable[k,j]*completeTable[i,k])
            if (is.na(inferValuesCandidate) == F){
              inferValues = c(inferValues,inferValuesCandidate)
            }
            if(length(inferValues) >= threshold){
              completeTable[i,j] = mean(inferValues)
            }
          }
        } 
      } 
    }
  }
  for (i in 1:nrow(completeTable)){
    for (j in 1:nrow(completeTable)){
      completeTable[j,i] = 1/completeTable[i,j]
    }
  }
  
  return(completeTable)
}

#### 1. Library Loading ####
library(pmr)

#### 2. Files Loading ####
files = grep("table",list.files(), value = TRUE)

#### 3.Making Ratio Table ####
tmpTable = data.frame()
for (i in 1:length(files)){
  fileRead = read.table(files[i], header = T)
  allelePair = sort(c(as.character(fileRead[1,1]),as.character(fileRead[2,1])))
  if(allelePair[1] == as.character(fileRead[1,1])){
    allelePairRatio = fileRead[1,15]/fileRead[2,15]
    tmpTable[i,1] = allelePair[1]
    tmpTable[i,2] = allelePair[2]
    tmpTable[i,3] = paste0(allelePair[1],allelePair[2])
    tmpTable[i,4] = allelePairRatio
    tmpTable[i,5] = files[i]
    tmpTable[i,6] = fileRead[1,15]
    tmpTable[i,7] = fileRead[2,15]
  }else{
    allelePairRatio = fileRead[2,15]/fileRead[1,15]
    tmpTable[i,1] = allelePair[1]
    tmpTable[i,2] = allelePair[2]
    tmpTable[i,3] = paste0(allelePair[1],allelePair[2])
    tmpTable[i,4] = allelePairRatio
    tmpTable[i,5] = files[i]
    tmpTable[i,6] = fileRead[2,15]
    tmpTable[i,7] = fileRead[1,15]
    print(files[i])
  }
}

tmpTable = tmpTable[complete.cases(tmpTable),]
row.names(tmpTable) = 1:nrow(tmpTable) 
colnames(tmpTable) = c("allele1","allele2","pair","allele1/allele2","Sample","allele1Count","allele2Count")
write.table(tmpTable,"KORratio", sep = "\t", row.names = F, quote = F)

#### 4.Create Pairwise table ####
alleleLevels = sort(unique(as.factor(c(tmpTable$allele1,tmpTable$allele2))))
SafeAlleleLevels = as.character(alleleLevels[table(rbind(tmpTable$allele1,tmpTable$allele2)) >= 3])
completeTable = matrix(data = NA, nrow = length(SafeAlleleLevels), ncol = length(SafeAlleleLevels), dimnames = list(SafeAlleleLevels))
colnames(completeTable) = SafeAlleleLevels

for (i in 1:nrow(completeTable)){
  completeTable[i,i] = 1
}
safePair = c()
for (i in 1:nrow(tmpTable)){
  if(length(which(SafeAlleleLevels == tmpTable[i,1])) != 0 && length(which(SafeAlleleLevels == tmpTable[i,2])) != 0){
    if(length(which(tmpTable$pair == tmpTable[i,3])) >= 1){
      safePair = c(safePair, tmpTable[i,3])
    }
  }
}
safePair = sort(unique(safePair))

for (i in 1:length(safePair)){
  pairSamples = which(tmpTable$pair == safePair[i])
  
  pairMean = mean(tmpTable[pairSamples,4])
  
  pairRowLocation = tmpTable[pairSamples,1]
  pairColLocation = tmpTable[pairSamples,2]
    
  completeTable[pairRowLocation,pairColLocation] = pairMean
  completeTable[pairColLocation,pairRowLocation] = 1/pairMean
}


#### 5.Create completed Pairwise table ####
for (m in nrow(completeTable):0){
  completeTable = InferFunc(completeTable,m)
}

#### 6.Create AHP results from completed pairwise matrix ####
Results = ahp(completeTable)
resultTable = matrix(ncol = 2, nrow = nrow(completeTable))
resultTable[,1] = rownames(completeTable)
resultTable[,2] = Results$weighting
colnames(resultTable) = c("Allele","Weight")
resultTable = data.frame(resultTable)
resultTable$Allele = as.character(resultTable$Allele)
resultTable$Weight = as.character(resultTable$Weight)

resultTable[nrow(resultTable)+1,1] = "Satty"
resultTable[nrow(resultTable),2] = Results$Saaty
write.table(resultTable,"Weight", sep = "\t", quote = F , row.names = F)