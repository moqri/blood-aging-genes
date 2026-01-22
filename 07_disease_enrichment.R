library(DOSE)
library(readxl)
library(clusterProfiler)
library(enrichplot)
maindb = read.csv("~/Projects/Collaboration/Mahdi/genes_and_cpgs_for_network.csv", header=TRUE)

id <- as.data.frame(maindb[1:1])
rna <- as.data.frame(maindb[2:106])
dnam <- as.data.frame(maindb[,108:ncol(maindb)])
dnam <- dnam[,-which(colnames(dnam) == "age")]
age_column <- as.data.frame(maindb["age"])  

#Input and conversion of names
list_rna <- colnames(rna)
eg = bitr(list_rna, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- eg[,-1]
head(eg)
#data(geneList)
#gene <- names(geneList)[abs(geneList) > 1.5]
#head(gene)

teg <- eg[1:10]
x <- enrichDO(gene          = eg,
              ont           = "HDO",
              pvalueCutoff  = 0.05,
              readable      = TRUE)
head(x)

dgn <- enrichDGN(eg) 
head(dgn)

barplot(dgn, showCategory=15) 
