rm(list = ls())
library(dplyr)
library(tibble)
library(ggVennDiagram)
####################################
#       set your working directory
####################################
setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5'seq_061922/Quinnoutput/")

####################################
#Read your TSSpreditor file output.
####################################

Ksg <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5'seq_061922/5'seqdata_TSSpreditor/KsgData_Cleaned.csv", header = T)
top <- Ksg %>% filter(SuperStrand == "+")
comp <- Ksg %>% filter(SuperStrand == "-")

####################################
#Read your annotation file
#################################### 

genes <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/Nikhil-T4pipeline-master/NC_003028.v3.17.ncrna.genes", header = F)
genes <- genes %>% add_column(UTR="", .before = "V2") %>% 
  add_column(termUTR = "", .after = "V3")
genes <- genes[complete.cases(genes),]
for(i in 1:nrow(genes)){
  genes$UTR[i] <- genes$V2[i] - 150
  genes$termUTR[i] <- genes$V3[i] + 150
}
genes.top <- genes %>% filter(V4 == "+")
colnames(genes.top) <- c("genome", "UTR", "from", "to","termUTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )
genes.comp <- genes %>% filter(V4 == "-")
colnames(genes.comp) <- c("genome", "termUTR", "to", "from","UTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )

genes.top$UTR <- as.numeric(genes.top$UTR)
genes.top$from <- as.numeric(genes.top$from)
genes.top$to <- as.numeric(genes.top$to)
genes.top$termUTR <- as.numeric(genes.top$termUTR)
genes.top <- genes.top[complete.cases(genes.top),]
genes.comp$UTR <- as.numeric(genes.comp$UTR)
genes.comp$from <- as.numeric(genes.comp$from)
genes.comp$to <- as.numeric(genes.comp$to)
genes.comp$termUTR <- as.numeric(genes.comp$termUTR)
genes.comp <- genes.comp[complete.cases(genes.comp),]
comp$SuperPos <- as.numeric(comp$SuperPos)
top$SuperPos <- as.numeric(top$SuperPos)

for(i in 2:nrow(genes.top)){
  if(genes.top$UTR[i] < genes.top$to[i-1]){
    genes.top$UTR[i] <- genes.top$to[i-1] + 50
  }
  if(genes.top$termUTR[i-1] > genes.top$from[i]){
    genes.top$termUTR[i-1] <- genes.top$from[i] - 50
  }
}

for(i in 2:nrow(genes.comp)){
  if(genes.comp$UTR[i-1] > genes.comp$to[i]){
    genes.comp$UTR[i-1] <- genes.comp$to[i] - 50
  }
  if(genes.comp$termUTR[i] < genes.comp$from[i-1]){
    genes.comp$termUTR[i] <- genes.comp$from[i-1] + 50
  }
}

rownames(genes.top) <- 1:nrow(genes.top)
rownames(genes.comp) <- 1:nrow(genes.comp)

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")

####################################
#Annotate your file
####################################

for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg <- rbind(comp, top)
rm(comp,top)

Ksg <- Ksg %>% arrange(SuperPos) %>% 
  filter(!Genome == "NDCt0")
Ksggenes <- Ksg[,-c(15:33)]

Ksggenes <- Ksggenes %>% filter(!Locus == "" | !TSSUTR == "" | !termUTR == "")
Ksggenes$Locus <- paste(Ksggenes$TSSUTR, Ksggenes$termUTR, Ksggenes$Locus) 
Ksggenes <- Ksggenes[,-c(2,4)]
KsgNDC <- Ksggenes %>% filter(Genome == "NDCt60")
Ksg1Q <-  Ksggenes %>% filter(Genome == "Ksg1q_t60")
Ksg3Q <-  Ksggenes %>% filter(Genome == "Ksg3q_t60")

KsgNDC <- KsgNDC[complete.cases(KsgNDC),]
Ksg1Q <- Ksg1Q[complete.cases(Ksg1Q),]
Ksg3Q <- Ksg3Q[complete.cases(Ksg3Q),]
####################################
#       Ksg NDC primary TSS
####################################

listOfCats <- unique(KsgNDC$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(KsgNDC)))
colnames(outputframe) <- colnames(KsgNDC)
for(x in 1:length(listOfCats)){
  
  tempFrame <- KsgNDC[KsgNDC$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempFrame <- tempFrame[tempFrame$stepHeight == tempMaximum,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
KsgNDCprimary <- outputframe
####################################
#           Ksg 1Q primary TSS
####################################
listOfCats <- unique(Ksg1Q$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(Ksg1Q)))
colnames(outputframe) <- colnames(Ksg1Q)
for(x in 1:length(listOfCats)){
  
  tempFrame <- Ksg1Q[Ksg1Q$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempFrame <- tempFrame[tempFrame$stepHeight == tempMaximum,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
Ksg1Qprimary <- outputframe
####################################
#       Ksg 3Q primary TSS
####################################       

listOfCats <- unique(Ksg3Q$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(Ksg3Q)))
colnames(outputframe) <- colnames(Ksg3Q)
for(x in 1:length(listOfCats)){
  
  tempFrame <- Ksg3Q[Ksg3Q$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempFrame <- tempFrame[tempFrame$stepHeight == tempMaximum,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
Ksg3Qprimary <- outputframe
####################################
#       Annotate Ksg NDC primary
####################################
KsgNDCprimary <- KsgNDCprimary[,-c(2)]
Ksg1Qprimary <- Ksg1Qprimary[,-c(2)]
Ksg3Qprimary <- Ksg3Qprimary[,-c(2)]
top <- KsgNDCprimary %>% filter(SuperStrand == "+")
comp <- KsgNDCprimary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

KsgNDCprimary <- rbind(comp, top)

####################################
#    Annotate Ksg 1Q Primary TSS
####################################
top <- Ksg1Qprimary %>% filter(SuperStrand == "+")
comp <- Ksg1Qprimary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg1Qprimary <- rbind(comp, top)
####################################
#    Annotate Ksg 3Q Primary TSS
####################################
top <- Ksg3Qprimary %>% filter(SuperStrand == "+")
comp <- Ksg3Qprimary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")

for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg3Qprimary <- rbind(comp, top)
####################################
#    Ksg NDC secondary TSS
####################################

listOfCats <- unique(KsgNDC$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(KsgNDC)))
colnames(outputframe) <- colnames(KsgNDC)
for(x in 1:length(listOfCats)){
  
  tempFrame <- KsgNDC[KsgNDC$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempsecondary <- max(tempFrame$stepHeight[tempFrame$stepHeight != max(tempFrame$stepHeight)])
  tempFrame <- tempFrame[tempFrame$stepHeight == tempsecondary,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
KsgNDCsecondary <- outputframe

####################################
#    Ksg 1Q secondary TSS 
####################################

listOfCats <- unique(Ksg1Q$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(Ksg1Q)))
colnames(outputframe) <- colnames(Ksg1Q)
for(x in 1:length(listOfCats)){
  
  tempFrame <- Ksg1Q[Ksg1Q$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempsecondary <- max(tempFrame$stepHeight[tempFrame$stepHeight != max(tempFrame$stepHeight)])
  tempFrame <- tempFrame[tempFrame$stepHeight == tempsecondary,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
Ksg1Qsecondary <- outputframe

####################################
#    Ksg 3Q secondary TSS 
####################################

listOfCats <- unique(Ksg3Q$Locus)
tempFrame <- c()
outputframe <- as.data.frame(matrix(nrow = 0, ncol = ncol(Ksg3Q)))
colnames(outputframe) <- colnames(Ksg3Q)
for(x in 1:length(listOfCats)){
  
  tempFrame <- Ksg3Q[Ksg3Q$Locus %in% listOfCats[x],]
  tempMaximum <- max(tempFrame$stepHeight)
  tempsecondary <- max(tempFrame$stepHeight[tempFrame$stepHeight != max(tempFrame$stepHeight)])
  tempFrame <- tempFrame[tempFrame$stepHeight == tempsecondary,]
  outputframe <- rbind(tempFrame, outputframe)
}

rownames(outputframe) <- 1:nrow(outputframe)
outputframe <- outputframe %>% arrange(SuperPos)
Ksg3Qsecondary <- outputframe

KsgNDCsecondary <- KsgNDCsecondary[,-c(2)]
Ksg1Qsecondary <- Ksg1Qsecondary[,-c(2)]
Ksg3Qsecondary <- Ksg3Qsecondary[,-c(2)]

####################################
#   Annotate Ksg NDC secondary TSS
####################################
top <- KsgNDCsecondary %>% filter(SuperStrand == "+")
comp <- KsgNDCsecondary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")

for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

KsgNDCsecondary <- rbind(comp, top)

####################################
#   Annotate Ksg 1Q secondary TSS
####################################
top <- Ksg1Qsecondary %>% filter(SuperStrand == "+")
comp <- Ksg1Qsecondary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg1Qsecondary <- rbind(comp, top)

####################################
#   Annotate Ksg 3Q secondary TSS
####################################

top <- Ksg3Qsecondary %>% filter(SuperStrand == "+")
comp <- Ksg3Qsecondary %>% filter(SuperStrand == "-")

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>%
  add_column(termUTR ="", .after = "Locus")
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg3Qsecondary <- rbind(comp, top)


Ksgprimary <- rbind(KsgNDCprimary, Ksg1Qprimary, Ksg3Qprimary)
primary.genic <- Ksgprimary[((!Ksgprimary$Locus == "") & (Ksgprimary$TSSUTR == "") & (Ksgprimary$termUTR == "")),]
primary.TSS <- Ksgprimary[((!Ksgprimary$TSSUTR == "") & (Ksgprimary$termUTR == "")),]
primary.termUTR <- Ksgprimary[((Ksgprimary$TSSUTR == "") & (!Ksgprimary$termUTR == "")),]

Ksgsecondary <- rbind(KsgNDCsecondary, Ksg1Qsecondary, Ksg3Qsecondary)
secondary.genic <- Ksgsecondary[((!Ksgsecondary$Locus == "") & (Ksgsecondary$TSSUTR == "") & (Ksgsecondary$termUTR == "")),]
secondary.TSS <- Ksgsecondary[((!Ksgsecondary$TSSUTR == "") & (Ksgsecondary$termUTR == "")),]
secondary.termUTR <- Ksgsecondary[((Ksgsecondary$TSSUTR == "") & (!Ksgsecondary$termUTR == "")),]
