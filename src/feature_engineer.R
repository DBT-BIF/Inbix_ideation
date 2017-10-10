# Ideation challenge code for InBIx-2017
# R script to create features from raw data and feature selection using Boruta
# DBT-BIF team
# Author: Rakesh Saraswat, Deepak Sharma
# Date: 16 September, 2017

# Install packages if not installed ---------------------------------------

# install.packages("foreign")
install.packages("Boruta")
install.packages("Peptides")
source("https://bioconductor.org/biocLite.R")
install.packages("/tmp/RtmpCglJm8/downloaded_packages/BiocInstaller_1.24.0.tar.gz", repos = NULL, type = "source")

# Import required libraries -----------------------------------------------

library(foreign)
library(Boruta)
library(Peptides)
library(bioMedR)

# Feature engineering -----------------------------------------------------

swapped = readFASTA("/home/depthgr8/nnonpeptides/Inbix_ideation/data/3dswap-pred_negative_dataset.fasta",legacy.mode = FALSE,seqonly=TRUE)
length(swapped)
index1 <- c()
for(i in 1:length(swapped)-count){
  c <- checkProt(swapped[[i]])
  if(c==FALSE){
    index1 <- cbind(index1,i)
    }
}
length(swapped)
swapped <- swapped[-index1]
length(swapped)
index2 <- c()
for(i in 1:length(swapped)){
  c <- lengthpep(swapped[i]) > 10
    if(c==FALSE){
      index2 <- cbind(index2,i)
    }
}
length(swapped)
swapped <- swapped[-index2]
length(swapped)
aacomp <- aaComp(swapped)
aacomp
str(aacomp)
x1 <- data.frame()
for(i in 1:length(aacomp)){
  x <- as.data.frame(aacomp[i])
  x1 <-rbind(x1,t(x$Mole.))
}
colnames(x1) <- c("Tiny","Small","Aliphatic","Aromatic","NonPolar","Polar","Charged","Basic","Acidic")
head(x1)
x2 <-cbind(x1,aIndex(swapped))
x3 <- cbind(x2,boman(swapped))
x3 <- cbind(x3,charge(swapped))
fasgai <- fasgaiVectors(swapped)
x5 <- data.frame()
for(i in 1:length(fasgai)){
  x4 <- unlist(fasgai[i])
  x5 <-as.data.frame(rbind(x5,t(x4)))
}
colnames(x5) <- c("Hydrophobicity index","Alpha and turn propensities","Bulky properties","Compositional characteristic index","Local flexibility","Electronic properties")
x6 <- cbind(x3,x5)
x7 <- cbind(x6,hmoment(swapped))
xx <- hydrophobicity(swapped, scale = "KyteDoolittle")
x8 <- cbind(x7,hydrophobicity(swapped, scale = "KyteDoolittle"))
kiderafactor <- kideraFactors(swapped)
x9 <- data.frame()
for(i in 1:length(kiderafactor)){
  x10 <- unlist(kiderafactor[i])
  x9 <-as.data.frame(rbind(x9,t(x10)))
}
colnames(x9) <- c("Helix/bend preference","Side-chain size","Extended structure preference","Hydrophobicity","Double-bend preference","Partial specific volume","Flat extended preference","Occurrence in alpha region","pK-C","Surrounding hydrophobicity")
x11 <- cbind(x8,x9)
head(x11)
x12 <- cbind(x11,lengthpep(swapped))
mswhim <- mswhimScores(swapped)
x14 <- data.frame()
for(i in 1:length(mswhim)){
  x13 <- unlist(mswhim[i])
  x14 <-as.data.frame(rbind(x14,t(x13)))
}
x15 <- cbind(x12,x14)
x16 <- cbind(x15,pI(swapped, pKscale = "EMBOSS"))
protfp <- protFP(swapped)
x18 <- data.frame()
for(i in 1:length(protfp)){
  x17 <- unlist(protfp[i])
  x18 <-as.data.frame(rbind(x18,t(x17)))
}
x19 <- cbind(x16,x18)
head(x19)
x21 <- data.frame()
vhse_scales <- vhseScales(swapped)
for(i in 1:length(vhse_scales)){
  x20 <- unlist(vhse_scales[i])
  x21 <-as.data.frame(rbind(x21,t(x20)))
}
x22 <- cbind(x19,x21)
x24 <- data.frame()
z_scales <- zScales(swapped)
for(i in 1:length(z_scales)){
  x23 <- unlist(z_scales[i])
  x24 <-as.data.frame(rbind(x24,t(x23)))
}
x25 <- cbind(x22,x24)
write.arff(data_8477,file="8477.arff",eol = "\n", relation = deparse(substitute(data_inbix)))
names(data_inbix)[1:56] <- paste(seq(1:56))
factors <- c(colnames(data_inbix))
formula_data <- as.formula(paste("output~", paste(factors, collapse="+")))
positive_fasta = readFASTA("~/Downloads/0.fasta")
x = readFASTA("~/Downloads/0.fasta"[[1]],legacy.mode = FALSE,seqonly=TRUE)
x <- swapped
x27 <- data.frame()
for(i in 1:length(x)){
  x26 <- extrProtAAC(x[[i]])
  x27 <-as.data.frame(rbind(x27,t(x26)))
}
x28 <- cbind(x25,x27)
x30 <- data.frame()
for(i in 1:length(x)){
  x29 <- extrProtDC(x[[i]])
  x30 <-as.data.frame(rbind(x30,t(x29)))
}
x31 <- cbind(x28,x30)
x33 <- data.frame()
for(i in 1:length(x)){
  x32 <- extrProtTC(x[[i]])
  x33 <-as.data.frame(rbind(x33,t(x32)))
}
x34 <- cbind(x31,x33)
negative <- x34
negative$output <- 0
data_filtered_8477 <- rbind(positive,negative)
data_filtered_8477$output <- as.factor(data_filtered_8477$output)
write.arff(data_filtered_8477,file="data_filtered_8477.arff",eol = "\n", relation = deparse(substitute(data_inbix)))
checkProt(positive_fasta)

# Feature selection using boruta ------------------------------------------

train <- data_filtered_8477
str(train$output)
set.seed(123)
boruta.train <- Boruta(output ~., data = train, doTrace = 2)
lz <- lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)
final.boruta <- TentativeRoughFix(boruta.train)
boruta.df <- attStats(final.boruta)
levels(boruta.df$decision)
confirmed <- subset(boruta.df,decision=="Confirmed")
confirmed_variable <- c(row.names(confirmed))
c <- as.vector(confirmed_variable)
idx <- c()
for(i in 1:8476){
  if(boruta.df[i,]$decision=="Confirmed"){
    idx <- rbind(idx,i)
  }
}

selected_boruta <-  data.frame(row.names = seq(1:1185))
id=0
for(i in 1:length(idx)){
  id <- idx[i]
  selected_boruta <- cbind(selected_boruta,data_filtered_8477[id])
  }

selected_boruta <- cbind(selected_boruta,data_filtered_8477[8477])
write.arff(selected_boruta,file="selecte_boruta.arff",eol = "\n", relation = deparse(substitute(selected_boruta)))




