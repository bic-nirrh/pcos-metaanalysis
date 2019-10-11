
library(GSVA)
library(ComplexHeatmap)

## load in data
# G:/AllGSE_GPL570.txt for all genes of GPL570 without any p value or logFC criteria
# G:/AllGSE_absolute_GPL570.txt for all genes with absolute values of GPL570 without any p value or logFC criteria

# expr.data <- read.table(file = "G:/PCOS/cumulus_metaphase_tissues.txt", header = TRUE, sep = "\t", row.names = 1)
expr.data <- read.table(file = "G:/PCOS/Allgenes_3tissues_gpl570.txt", header = TRUE, sep = "\t", row.names = 1)
expr.data.mat <- as.matrix(expr.data)
genesets <- read.table(file = "G:/PCOS/GMT_only5_diseases.txt", header = TRUE, sep = "\t")

geneset.data=data.frame((genesets))



#GSVA and ssgsea Scoring

GSC.gsva <- gsva(expr = expr.data.mat,
                 gset.idx.list = genesets,
                 method = 'gsva',
                 kcdf = 'Gaussian',
                 min.sz = 1,
                 parallel.sz = 0,
                 verbose = F)


#Heatmap

#gsva

hm.gsva <- Heatmap(GSC.gsva,
                   row_names_gp = gpar(fontsize = 10),
                   show_column_names = T)
draw(hm.gsva)


