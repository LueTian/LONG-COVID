###### 包调用
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)

### Part Source

###### 路径设置 
setwd("Data\\Code_07_Outcome\\Source Protein")
file_list <- list.files()
# KEGG Path
kegg_dir <- paste('KEGG',file_list,sep = '_')
kegg_dir <- paste('Data\\Code_08_Outcome\\Source Protein\\GO\\',kegg_dir,sep = '')
go_dir <- paste('GO',file_list,sep = '_')
go_dir <- paste('Data\\Code_08_Outcome\\Source Protein\\GO\\',go_dir,sep = '')
for (i in 1:length(file_list)) {
  genelist <- read.csv(file_list[i],sep = ',')
  genes <- as.character(genelist$Entrez) #转换成字符格式
  kegg <- enrichKEGG(gene = genes,keyType = "kegg",organism = 'hsa',pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "fdr")
  go <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'ALL', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")
  write.table(as.data.frame(kegg),kegg_dir[i],sep=",",row.names =F)
  write.table(as.data.frame(go),go_dir[i],sep=",",row.names =F)
}


### Part Drug
###### 路径设置 
setwd("Data\\Code_07_Outcome\\Source Protein Add")
file_list <- list.files()
# KEGG Path
kegg_dir <- paste('KEGG',file_list,sep = '_')
kegg_dir <- paste('Data\\Code_08_Outcome\\Source Protein Add\\KEGG\\',kegg_dir,sep = '')
go_dir <- paste('GO',file_list,sep = '_')
go_dir <- paste('Data\\Code_08_Outcome\\Source Protein Add\\GO\\',go_dir,sep = '')
for (i in 1:length(file_list)) {
  genelist <- read.csv(file_list[i],sep = ',')
  genes <- as.character(genelist$Entrez) #转换成字符格式
  kegg <- enrichKEGG(gene = genes,keyType = "kegg",organism = 'hsa',pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "fdr")
  go <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'ALL', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")
  write.table(as.data.frame(kegg),kegg_dir[i],sep=",",row.names =F)
  write.table(as.data.frame(go),go_dir[i],sep=",",row.names =F)
}

### Part Target
###### 路径设置 
setwd("Data\\Code_07_Outcome\\Target Protein")
file_list <- list.files()
# KEGG Path
kegg_dir <- paste('KEGG',file_list,sep = '_')
kegg_dir <- paste('\Data\\Code_08_Outcome\\Target Protein\\KEGG\\',kegg_dir,sep = '')
go_dir <- paste('GO',file_list,sep = '_')
go_dir <- paste('Data\\Code_08_Outcome\\Target Protein\\GO\\',go_dir,sep = '')
for (i in 1:length(file_list)) {
  genelist <- read.csv(file_list[i],sep = ',')
  genes <- as.character(genelist$protein) #转换成字符格式
  kegg <- enrichKEGG(gene = genes,keyType = "kegg",organism = 'hsa',pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "fdr")
  go <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'ALL', pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")
  write.table(as.data.frame(kegg),kegg_dir[i],sep=",",row.names =F)
  write.table(as.data.frame(go),go_dir[i],sep=",",row.names =F)
}