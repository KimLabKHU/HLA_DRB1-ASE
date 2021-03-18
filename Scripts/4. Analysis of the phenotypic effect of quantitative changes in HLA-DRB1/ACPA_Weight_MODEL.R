Start = Sys.time()
#--------------------------------------------------------------------------------------------------------
#argument prepare
#--------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
  input = args[1]
  weight = args[2]
  output = args[3]
}

#--------------------------------------------------------------------------------------------------------
#Samples prepare
#--------------------------------------------------------------------------------------------------------
SampleTable = read.table(input, header = T)
HLATable = read.table(input, header = T, colClasses = "character")[,c(2,3)]
SampleTable[,c(2,3)] = HLATable

weight = read.table(weight, header = T)
weight[,1] = gsub("DRB1\\*","",weight[,1])
weight[,1] = gsub("\\:","",weight[,1])
weight[,1] = substr(weight[,1],1,4)

TableCol = ncol(SampleTable)
for (i in (TableCol+1) : (TableCol+nrow(weight))){
  SampleTable[,i] = 0
}
colnames(SampleTable)[c((TableCol+1):(TableCol+nrow(weight)))] = weight[,1]

for (i in 1:nrow(SampleTable)){
  tmpHLA1position = which(colnames(SampleTable) == as.character(SampleTable[i,2]))
  tmpHLA2position = which(colnames(SampleTable) == as.character(SampleTable[i,3]))
  if(length(tmpHLA1position) == 1 && length(tmpHLA2position) == 1){ 
  SampleTable[i,tmpHLA1position] = SampleTable[i,tmpHLA1position] + 1 
  SampleTable[i,tmpHLA2position] = SampleTable[i,tmpHLA2position] + 1
  }
}


removeUndefinedSamples = c()

for (i in 1:nrow(SampleTable)){
  if(sum(SampleTable[i,c((TableCol+1):ncol(SampleTable))]) != 2){
    removeUndefinedSamples = c(removeUndefinedSamples,i)
  }
}
if(length(removeUndefinedSamples) != 0){
SampleTable = SampleTable[-removeUndefinedSamples,]



write.table(removeUndefinedSamples,"./removeUndefinedSamples", sep = "\t", quote= F, row.names = F, col.names = F)

}

TotalTableCol = ncol(SampleTable)
for (i in 1:nrow(SampleTable)){
SampleTable[i,(TotalTableCol+1)] = sum(SampleTable[i,c((TableCol+1):TotalTableCol)] * weight[,3])
}
colnames(SampleTable)[TotalTableCol+1] = "A74"

AltSampleTable = SampleTable
for (i in 1:nrow(AltSampleTable)){
    tmpValue = c()
    for (j in ((TableCol+1):(ncol(AltSampleTable)-1))){
      if(AltSampleTable[i,j] == 1){
        tmpValue = c(tmpValue,j)
      }
    }
    if(length(tmpValue) == 2){
      tmpPosition1 = tmpValue[1]
      tmpPosition2 = tmpValue[2]
      
      tmpPosition1Weight = weight[which(as.character(weight$Allele) == colnames(AltSampleTable)[tmpPosition1]),2]
      tmpPosition2Weight = weight[which(as.character(weight$Allele) == colnames(AltSampleTable)[tmpPosition2]),2]
      
      AltSampleTable[i,tmpPosition1] = AltSampleTable[i,tmpPosition1] * 2 * tmpPosition1Weight / (tmpPosition1Weight + tmpPosition2Weight)
      AltSampleTable[i,tmpPosition2] = AltSampleTable[i,tmpPosition2] * 2 * tmpPosition2Weight / (tmpPosition1Weight + tmpPosition2Weight)
    }
}

for (i in 1:nrow(AltSampleTable)){
  AltSampleTable[i,(TotalTableCol+1)] = sum(AltSampleTable[i,c((TableCol+1):TotalTableCol)] * weight[,3])
}

write.table(SampleTable,"SampleTable", quote = F, sep = "\t", row.names = F)
write.table(AltSampleTable,"AltSampleTable", quote = F, sep = "\t", row.names = F)


#--------------------------------------------------------------------------------------------------------
#Modeling
#--------------------------------------------------------------------------------------------------------

xnam = paste("SampleTable[,", c(6:TableCol,(TotalTableCol+1)),"]", sep="")
fmla = as.formula(paste("SampleTable$INTCCP ~ ", paste(xnam, collapse= "+")))

ContModel = lm(fmla)

altxnam = paste("AltSampleTable[,", c(6:TableCol,(TotalTableCol+1)),"]", sep="")
altfmla = as.formula(paste("AltSampleTable$INTCCP ~ ", paste(altxnam, collapse= "+")))

AltModel = lm(altfmla)

Heteronam = paste(colnames(AltSampleTable)[c(6:TableCol,(TotalTableCol+1))], sep="")
Heterofmla = as.formula(paste("INTCCP ~ ", paste(Heteronam, collapse= "+")))

HeteroA74Model = lm(Heterofmla, data = subset(AltSampleTable, SampleTable$A74 == 1))

print(nrow(SampleTable))
print(nrow(subset(SampleTable, SampleTable$A74 == 1)))

print(summary(ContModel))
print(summary(ContModel)$coefficients)
print(confint(ContModel))


print(summary(AltModel))
print(summary(AltModel)$coefficients)
print(confint(AltModel))

print(summary(HeteroA74Model))
print(summary(HeteroA74Model)$coefficients)
print(confint(HeteroA74Model))

results = data.frame()
results[1,1] = AIC(ContModel)
results[1,2] = AIC(AltModel)
results[2,1] = BIC(ContModel)
results[2,2] = BIC(AltModel)
results[3,1] = summary(ContModel)$adj.r.squared
results[3,2] = summary(AltModel)$adj.r.squared
results[4,1] = "estimate"
results[4,2] = summary(HeteroA74Model)$coefficients[which(rownames(summary(HeteroA74Model)$coefficients) =="A74"),1]
results[5,1] = "stdError"
results[5,2] = summary(HeteroA74Model)$coefficients[which(rownames(summary(HeteroA74Model)$coefficients) =="A74"),2]
results[6,1] = "t-value"
results[6,2] = summary(HeteroA74Model)$coefficients[which(rownames(summary(HeteroA74Model)$coefficients) =="A74"),3]
results[7,1] = "p-value"
results[7,2] = summary(HeteroA74Model)$coefficients[which(rownames(summary(HeteroA74Model)$coefficients) =="A74"),4]

rownames(results) = c("AIC","BIC","Adjusted R2", "Hetero A74_1", "Hetero A74_2", "Hetero A74_3", "Hetero A74_4")
colnames(results) = c("Conv","Alt")
write.table(results[6,2],"OBS.t", quote=F)
write.table(results, output, quote= F)

End = Sys.time()

print(End-Start)
print("FINISH!")
