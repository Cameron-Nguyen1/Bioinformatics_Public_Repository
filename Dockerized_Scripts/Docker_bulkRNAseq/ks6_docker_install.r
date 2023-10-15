#!/usr/bin/env Rscript
install.packages("BiocManager",INSTALL_opts = '--no-lock')
dep1=c("GSEABase","Biobase","GSVA","gprofiler2","clusterProfiler","msigdbr","enrichplot","IsoformSwitchAnalyzeR")
lapply(X=dep1,FUN=BiocManager::install,character.only=TRUE,INSTALL_opts = '--no-lock')
dep2= c("methods","tidyverse","tximport","ensembldb","EnsDb.Hsapiens.v86","EnsDb.Mmusculus.v79","edgeR","matrixStats","cowplot","DT","plotly","gt","limma","pheatmap")
#for (i in dep2){install.packages(i,INSTALL_opts = '--no-lock')}
lapply(X=dep2,FUN=install.packages,character.only=TRUE,INSTALL_opts = '--no-lock')