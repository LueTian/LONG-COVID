library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GOSemSim)
library(mygene)
library(AnnotationDbi)
library(WGCNA)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mygene")
BiocManager::install("WGCNA")
BiocManager::install("org.Mm.eg.db")
n### Part Source
###### 路径设置 
setwd("C:\\Users\\31435\\Desktop\\ACESO Paper\\Appendix")
data <- read.csv("Supplementary table 6 Protein_Frequency.csv")

# KEGG Path
kegg_dir <- paste('KEGG',file_list,sep = '_')
kegg_dir <- paste('C:\\Users\\31435\\Desktop\\ACESO 数据备份\\Research Project\\ADHD\\GO_KEGG\\Result\\Disease\\KEGG\\',kegg_dir,sep = '')
go_dir <- paste('GO',file_list,sep = '_')
go_dir <- paste('C:\\Users\\31435\\Desktop\\ACESO 数据备份\\Research Project\\ADHD\\GO_KEGG\\Result\\Disease\\GO\\',go_dir,sep = '')


genes <- as.character(data$Protein) #转换成字符格式
go <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'ALL', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")
write.table(as.data.frame(go),"GO_HFDROSATNC.csv",sep=",",row.names =F)

### WGCNA 
rm(list = ls())
expressiondata = read.csv("gene_expression_filtered_df_HFDROSATNC.csv", sep = ",")
sampleTree = hclust(dist(expressiondata), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 31000, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 30, minSize = 10)
table(clust)

MEList = moduleEigengenes(expression)
MEs = MEList$eigengenes

#clust

#Cluster 1 contains the samples we want to keep.
keepSamples = (clust==1)
expression = expression[keepSamples, ]
#dim(expression0)
nGenes = ncol(expression)
nSamples = nrow(expression)



