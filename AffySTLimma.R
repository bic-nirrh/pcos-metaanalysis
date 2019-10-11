#Read affy cel files and analyze differentially expressed genes to generate .rnk files
#please go to the directory with cel files first
source("https://bioconductor.org/biocLite.R")
biocLite("oligo")
library(oligo)
require(oligo)
require(limma)

NamMat<-function(x) {
  return(UseGenes[which(rownames(UseGenes)==x[1]),1])
}

#read sample attributes in Group.txt file
g1<-list.files(pattern="Group.txt")
group<-read.table(g1,header=TRUE,sep="\t")
gc1<-table(group[,2])

#read, rma background correction, quantile normalization, and pm correction
celfiles <- list.files(pattern = "CEL", full = TRUE)
dataRaw<- read.celfiles(celfiles)
datamQN<- rma(dataRaw)
featureData(datamQN) <- getNetAffx(datamQN, 'transcript')
RawGenes<-exprs(datamQN)
RefSeq<-apply(datamQN@featureData@data,1,function(x) unlist(strsplit(x[8],"//"))[1])
Symb1<-apply(datamQN@featureData@data,1,function(x) unlist(strsplit(x[8],"//"))[2])
RawGenes<-cbind(Symb1,RefSeq,RawGenes)
colnames(RawGenes)[1:2]<-c("GeneName","RefSeq")
UseGenes<-RawGenes
UseGenes<-UseGenes[!is.na(UseGenes[,1]),]

#Collapse duplicated genes by selecting the probe with highest SD
DupGenes<-unique(UseGenes[duplicated(UseGenes[,1]),1])
for (i in 1:length(DupGenes)){
  dat<-NULL
  sdv<-NULL
  selected<-NULL
  indx<-NULL
  indx <- which(UseGenes[,1] == DupGenes[i])
  dat<-UseGenes[indx,]
  sdv<-apply(dat[,3:ncol(dat)], 1, sd)
  selected<-names(sdv)[which(sdv==max(sdv))]
  dat<-dat[which(rownames(dat)==selected)[1],]
  for (j in 1:length(indx)){
  UseGenes[indx[j],]<-dat #replace all duplicated genes with values of max SD
}}
UseGenes <- UseGenes[!duplicated(UseGenes[,1]), ]
UseGenes[,3:ncol(UseGenes)]<-as.numeric(UseGenes[,3:ncol(UseGenes)])
group2<-group[order(group[,2]),]
ExpressOut<-UseGenes[,c(1,1)]
name3<-NULL
for (i in 1:nrow(group2)){
  ExpressOut<-cbind(ExpressOut,UseGenes[,which(colnames(UseGenes)==group2[i,1])])
  name3<-c(name3,colnames(UseGenes)[which(colnames(UseGenes)==group2[i,1])])
}
colnames(ExpressOut)<-c("Name","geneid",as.character(name3))
write.table(ExpressOut,"Expression_All_Normalized.txt",sep="\t",row.names=FALSE,quote=FALSE)
med<-apply(ExpressOut[,3:ncol(ExpressOut)],1,median)
ExpressOutMed<-matrix(as.numeric(ExpressOut[,3:ncol(ExpressOut)]),ncol=(ncol(ExpressOut)-2))-as.numeric(med)
ExpressOutMed<-cbind(ExpressOut[,1:2],ExpressOutMed)
colnames(ExpressOutMed)<-colnames(ExpressOut)
write.table(ExpressOutMed,"MedianCentered_All.txt",sep="\t",row.names=FALSE,quote=FALSE)

#Determine the number of comparisons, choose the groups, and process the data by limma
  exp1<-group[which(group[,2]==names(gc1)[2]),1]
  data1<-UseGenes[,which(colnames(UseGenes) %in% exp1)]
  con1<-group[which(group[,2]==names(gc1)[1]),1]
  data2<-UseGenes[,which(colnames(UseGenes) %in% con1)]
  #This is a simpler way to code ModT for two group comparison by generating design matrix with a contrast column (Exp vs Con) instead individual group label 
  design<-matrix(ncol=2,nrow=(gc1[2]+gc1[1]),data=1)
  rownames(design)<-c(as.character(exp1),as.character(con1))
  colnames(design)<-c("All","EvC")
  design[(length(exp1)+1):nrow(design),2]<-0
  compa<-cbind(data1,data2)
  compa2<-matrix(as.numeric(compa),ncol=ncol(compa))
  colnames(compa2)<-colnames(compa)
  rownames(compa2)<-rownames(compa)
  
  #ModT Test and output results, .rnk (T-values),Express, and cls files
  fit<-lmFit(compa2,design)
  fit<-eBayes(fit)
  result1<-topTable(fit,coef="EvC",adjust="BH",number=nrow(compa2))
  result1<-cbind(rownames(result1),result1)
  names<-apply(result1,1,NamMat)
  result1<-cbind(names,result1)
  write.table(result1,paste0("ModT-Results-",names(gc1)[2],"-vs-",names(gc1)[1],".txt"),sep="\t",row.names=FALSE,col.names=TRUE)
  rnkfile<-cbind(as.character(result1[,1]),as.numeric(result1[,5]))
  rnkfile<-rnkfile[!is.na(rnkfile[,2]),]
  write.table(rnkfile,paste0("ModT_Results_",names(gc1)[2],"_vs_",names(gc1)[1],".rnk"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  compa2<-cbind(UseGenes[,1],UseGenes[,1],compa2)
  colnames(compa2)[1:2]<-c("Name","geneid")
  write.table(compa2,paste0("Expression_",names(gc1)[2],"_vs_",names(gc1)[1],".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  clsfile<-c(rep(names(gc1)[2],gc1[2]),rep(names(gc1[1]),gc1[1]))
  fileConn<-file(paste0(names(gc1)[2],"_vs_",names(gc1)[1],"_classes.cls"))
  writeLines(c(paste(length(clsfile),"2 1"),paste("# ", unique(clsfile)[1], " ",unique(clsfile)[2])), fileConn)
  write.table(t(clsfile),paste(names(gc1)[2],"_vs_",names(gc1)[1],"_classes.cls",sep=""),col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)
  close(fileConn)