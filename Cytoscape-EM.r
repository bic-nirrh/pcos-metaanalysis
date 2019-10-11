if(!"RCy3" %in% installed.packages()){
install.packages("BiocManager")
BiocManager::install("RCy3")
}
library(RCy3)
cytoscapePing ()
cytoscapeVersionInfo ()

# Make sure to modify the file namesand threshold values below
pvalue_gsea_threshold <- 1.0
qvalue_gsea_threshold <- 0.1
similarity_threshold <- "0.375"
similarity_metric = "COMBINED"
cur_model_name <-"GSE6798_updatedGMT"
current_network_name <- paste(cur_model_name,pvalue_gsea_threshold,qvalue_gsea_threshold,sep="_")
current_network_name
setwd("C:/Users/ICMR/Desktop/GSEA_updatedGMT_16112018/GSE6798_updatedGMT.GseaPreranked.1542634917274/edb")

# Transfer the GMT file (used for GSEA, Expression file obtained from CELfile analysis and rank file to this directory)
list.files()
#Assigning paths to read the file for Cytoscape
gmt_gsea_file <- file.path(getwd(),"New_Human_GOBP_AllPathways_no_GO_iea_November_01_2018_symbol_updated_ids.gmt")
gsea_ranks_file <- file.path(getwd(),"ModT_Results_2_vs_1_GSE5850_new.rnk")
gsea_results_filename <- file.path(getwd(),"results.edb")

em_command = paste('enrichmentmap build analysisType="gsea" gmtFile=',gmt_gsea_file,
'pvalue=',pvalue_gsea_threshold, 'qvalue=',qvalue_gsea_threshold,
'similaritycutoff=',similarity_threshold,
'coefficients=',similarity_metric,'ranksDataset1=',
gsea_ranks_file,'enrichmentsDataset1=',gsea_results_filename)
response <- commandsGET(em_command)

#With expression data
gsea_expression_filename <- file.path(getwd(),"Expression_2_vs_1.txt")
em_command = paste('enrichmentmap build analysisType="gsea" gmtFile=',gmt_gsea_file,
'pvalue=',pvalue_gsea_threshold, 'qvalue=',qvalue_gsea_threshold,
'similaritycutoff=',similarity_threshold,
'coefficients=',similarity_metric,'ranksDataset1=',
gsea_ranks_file,'enrichmentsDataset1=',gsea_results_filename, 'expressionDataset1=',gsea_expression_filename)

#To export the network table
install.packages("dplyr")
library("dplyr")
check <-data.frame(getTableColumns())
#Remove the Gene column as it is a sublist in the table
check2<-select(check,-(EM2_Genes))
write.csv(check2, sep = "", col.names = TRUE, file=paste(cur_model_name, "_node_table.csv", sep=""))#file = "foo.csv")

#To export the network table - selected entries and unlist the gene names
symbols <- strsplit(as.character(check$EM2_Genes), ",")
res <- data.frame(EM2_GS_DESCR=rep(check$EM2_GS_DESCR, sapply(symbols, length)),EM2_gs_size=rep(check$EM2_gs_size, sapply(symbols, length)), EM2_NES=rep(check$EM2_NES..Dataset.1., sapply(symbols, length)),EM2_Genes=unlist(symbols))
res$EM2_Genes <- gsub("c", "",res$EM2_Genes)
res$EM2_Genes <- gsub("\\(", "",res$EM2_Genes)
res$EM2_Genes <- gsub("\"", "",res$EM2_Genes)
res$EM2_Genes <- gsub("\\)", "",res$EM2_Genes)
View(res)
write.csv(res, sep = "", col.names = TRUE, file=paste(cur_model_name, "_gene_list.csv", sep=""))#file = "foo2.csv")

#To Save and export the network visualisation

full.path=paste(getwd(),'vignette_session',sep='/')
saveSession(filename=full.path) #.cys

full.path=paste(getwd(),'vignette_image',sep='/')
exportImage(filename=full.path, type = 'PDF') #.pdf




