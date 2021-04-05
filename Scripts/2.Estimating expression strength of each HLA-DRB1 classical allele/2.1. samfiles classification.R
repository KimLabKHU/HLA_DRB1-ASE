#!/usr/bin/Rscript

# Lab Research - HLA-DRB1 ASE
# 20210328 Sehwan Chun at Corestem, Inc.
# 2.1. samfiles classification 

#### 0.argument setup ####
args = commandArgs(trailingOnly = TRUE)
if(length(args) == 1){ARG1 = args[1]}

test = read.csv(paste0('./SamFiles/',ARG1), header = FALSE, sep = "") # used in samfiles folders 

colnames(test) = c("QName","Flag","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","ASI","XSI","XNI","XMI","XOI","XGI","NMI","MDZ","YSI","YTZ")

#### 1. Singleton-unread filtered ####

#First Filter,  Singleton-unread
test_singleton_unread = subset(test, ASI == "YT:Z:UP")
test_singleton_unread = test_singleton_unread[,c(1:13)]
colnames(test_singleton_unread) = c("QName","Flag","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","YTZ","YFZ")
test_SU_Removed = subset(test, ASI != "YT:Z:UP")

# status correction
test_SU_Removed$ASI = as.character(test_SU_Removed$ASI)
test_SU_Removed$XSI = as.character(test_SU_Removed$XSI)
test_SU_Removed$XNI = as.character(test_SU_Removed$XNI)
test_SU_Removed$XMI = as.character(test_SU_Removed$XMI)
test_SU_Removed$XOI = as.character(test_SU_Removed$XOI)
test_SU_Removed$XGI = as.character(test_SU_Removed$XGI)
test_SU_Removed$NMI = as.character(test_SU_Removed$NMI)
test_SU_Removed$MDZ = as.character(test_SU_Removed$MDZ)
test_SU_Removed$YSI = as.character(test_SU_Removed$YSI)
test_SU_Removed$YTZ = as.character(test_SU_Removed$YTZ)

#XNI Correction
for (i in 1:nrow(test_SU_Removed)){
  if(test_SU_Removed[i,13] == "XN:i:0") {
      if(test_SU_Removed[i,19] == "YT:Z:UP"){
        test_SU_Removed[i,21] = test_SU_Removed[i,19]
        test_SU_Removed[i,20] = "YS:i:NULL"
        test_SU_Removed[i,19] = test_SU_Removed[i,18]
        test_SU_Removed[i,18] = test_SU_Removed[i,17]
        test_SU_Removed[i,17] = test_SU_Removed[i,16]
        test_SU_Removed[i,16] = test_SU_Removed[i,15]
        test_SU_Removed[i,15] = test_SU_Removed[i,14]
        test_SU_Removed[i,14] = test_SU_Removed[i,13]
        test_SU_Removed[i,13] = "XS:i:NULL"
      }
      else if(test_SU_Removed[i,19] != "YT:Z:UP"){
        test_SU_Removed[i,21] = test_SU_Removed[i,20]
        test_SU_Removed[i,20] = test_SU_Removed[i,19]
        test_SU_Removed[i,19] = test_SU_Removed[i,18]
        test_SU_Removed[i,18] = test_SU_Removed[i,17]
        test_SU_Removed[i,17] = test_SU_Removed[i,16]
        test_SU_Removed[i,16] = test_SU_Removed[i,15]
        test_SU_Removed[i,15] = test_SU_Removed[i,14]
        test_SU_Removed[i,14] = test_SU_Removed[i,13]
        test_SU_Removed[i,13] = "XS:i:NULL"
      }  
  } 
}

#YSI
for (i in 1:nrow(test_SU_Removed)){
  if (test_SU_Removed[i,20] == "YT:Z:UP"){
    test_SU_Removed[i,21] = test_SU_Removed[i,20]
    test_SU_Removed[i,20] = "YS:i:NULL"
  }
}

#### 2.classfied each samfiles into 3 subgroups ####
#Second Filter

test_SU_Removed_Corrected = subset(test_SU_Removed, YTZ != "YT:Z:DP")
test_singleton_read = subset(test_SU_Removed_Corrected, YSI == "YS:i:NULL")
test_remain_read = subset(test_SU_Removed_Corrected, YSI != "YS:i:NULL")

write.table(test_singleton_unread, paste0("./SamFiles/",ARG1,"_singleton_unread"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(test_singleton_read, paste0("./SamFiles/",ARG1,"_singleton_read"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(test_remain_read, paste0("./SamFiles/",ARG1,"_remain"), quote = FALSE, row.names = FALSE, sep = "\t")

