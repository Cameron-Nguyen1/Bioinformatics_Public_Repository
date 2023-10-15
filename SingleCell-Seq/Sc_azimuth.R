#!/bin/env Rscript
dep = c("plyr","data.table","devtools","Azimuth","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","ggrepel","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid")
lapply(dep,library,character.only=TRUE)
rem_cmo = function(x){
    x = x[!grepl("CMO",rownames(x)),]
    return(x)
}
aziPredict = function(dataset,reference){
    dataset = RunAzimuth(query=dataset, reference)
    mzl = list("Map_Score.50"=table(dataset$mapping.score > .50),"Map_Score.75"=table(dataset$mapping.score > .75))
    newl_s = list("L1_Annot"=dataset$predicted.ann_level_1.score, "L2_Annot"=dataset$predicted.ann_level_2.score, "L3_Annot"=dataset$predicted.ann_level_3.score, "L4_Annot"=dataset$predicted.ann_level_4.score, "L5_Annot"=dataset$predicted.ann_level_5.score, "L6_Annot"=dataset$predicted.ann_finest_level.score)
    for (i in 1:length(names(newl_s))){
        mz_50 = list(name50=table(newl_s[[i]] >= .50))
        mz_75 = list(name75=table(newl_s[[i]] >= .75))
        mzl = c(mzl, mz_50, mz_75)
    }
    lnum = 1
    for (i in 3:length(names(mzl))){
        if (i %% 2 != 0 && i != 3){
            lnum = lnum + 1
        }
        if (i %% 2 == 0){
            names(mzl)[i] = paste0("L",lnum,"_Pred.Score_75")
        }else{
            names(mzl)[i] = paste0("L",lnum,"_Pred.Score_50")
        }
    }
    df_pred = data.frame(mzl)
    df_pred = data.frame(t(df_pred[,!grepl("Var1",colnames(df_pred))]))
    df_pred = cbind(df_pred,list(pct.True=df_pred[,2]/(df_pred[,2]+df_pred[,1])))
    colnames(df_pred) = c("Below.Scr","Above.Scr","Pct.Above")
    write.csv(x=df_pred,file="Azimuth_Prediction_Scores.csv")
    return(dataset)
}
make_UMAP_Azimuth = function(x,pdf_name){
    ref = colnames(x@meta.data)
    plst = list()
    for (i in ref[grepl("^predicted",ref) & !grepl("score",ref)]){
        plst[[i]] = DimPlot(x,reduction="umap", label=T, repel=T, label.size=2.5, group.by=i) + ggtitle(paste0("Azimuth "),i)
    }
    g=marrangeGrob(grobs=plst,ncol=2,nrow=1)
    ggsave(paste0(pdf_name,".pdf"),plot=g,width=23,height=10,device="pdf", limitsize = FALSE)
}
hm_prop = function(x){
    strang = deparse(substitute(x))
    if (strang %like% "harm"){
        c_title = "log-Norm + Harmony"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$RNA_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$RNA_snn_res.0.2)))}
    else if(strang %like% "SCT"){
        c_title = "SCT + CCA"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$integrated_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$integrated_snn_res.0.2)))}
    else{
        c_title = "log-Norm + Merged"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$RNA_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$RNA_snn_res.0.2)))}
    nr <- length(levels(as.factor(x@meta.data$orig.ident)))
    hm1 <- ComplexHeatmap::Heatmap(tab, cluster_columns = FALSE,
                            show_column_dend = TRUE, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = 12), column_names_gp = grid::gpar(fontsize=12),
                            cluster_rows = TRUE, show_row_dend = FALSE, col = heatmap.cols, border_gp = gpar(col = "black", lty = 1),
                            rect_gp = gpar(col = "white", lwd = 2), column_title = c_title, height = unit(5, "mm")*nr, width = unit(5, "mm")*nc, name="Proportion of cells")
    return(hm1)
}
bdim = function(x){
  if (x[["pca"]]@assay.used == "integrated"){
    DefaultAssay(x) = "integrated"
    pct=x[['pca']]@stdev/sum(x[['pca']]@stdev)*100 #Make PCT var for PCA dim selection
    csum = cumsum(pct)
    c1 = which(pct < 5 & csum > 90)[1] #PCA dim must contribute <5 variance and > 90 cumulative variance
    c2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1 #Find dim with less than .1 difference in step
    if (is.na(c1) | is.na(c2)){
      z = c(c1,c2)
      bestdims = z[!is.na(z)]
    }else{
    bestdims = min(c1,c2) #We take the lowest dimension of the two as the range of dims that capture a majority of variance
    }
    return(bestdims)
  }
  pct=x[['pca']]@stdev/sum(x[['pca']]@stdev)*100 #Make PCT var for PCA dim selection
  csum = cumsum(pct)
  c1 = which(pct < 5 & csum > 90)[1] #PCA dim must contribute <5 variance and > 90 cumulative variance
  c2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1 #Find dim with less than .1 difference in step
  if (is.na(c1) | is.na(c2)){
    z = c(c1,c2)
    bestdims = z[!is.na(z)]
  }else{
  bestdims = min(c1,c2) #We take the lowest dimension of the two as the range of dims that capture a majority of variance
  }
  return(bestdims)
}

annotate = function(x,assay){
  #ImmGen (Mouse Bulk Expression)
  immgen.ref <- celldex::ImmGenData(ensembl=FALSE)
  #table(immgen.ref$label.main)
  # perform predictions
  immgen.pred <- SingleR::SingleR(test = x@assays[[assay]]@data, ref = immgen.ref, labels = immgen.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  #table(immgen.pred$labels)
  #table(immgen.pred$pruned.labels)
  x@meta.data$immgen.pred <- immgen.pred$labels
  x@meta.data$immgen.pred_pruned <- immgen.pred$pruned.labels # contains NA

  #MRSD (Mouse RNA-seq Data)
  mrsd.ref <- celldex::MouseRNAseqData(ensembl=FALSE)
  mrsd.pred <- SingleR::SingleR(test = x@assays[[assay]]@data, ref = mrsd.ref, labels = mrsd.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  x@meta.data$mrsd.pred <- mrsd.pred$labels
  x@meta.data$mrsd.pred_pruned <- mrsd.pred$pruned.labels # contains NA

  return(c(x,immgen.ref,mrsd.ref))
}

clean = function(x){ ##Optimize this later, use  a for loop to hit all *.pred columns.
  x@meta.data$immgen.pred.clean <- forcats::fct_collapse(x@meta.data$immgen.pred,`B cells` = "B cells, pro")
  freq <- as.data.frame(table(x@meta.data$immgen.pred))
  freq <- freq %>% mutate("prop" = prop.table(Freq))
  colnames(freq) <- c("cells", "n", "proportion")
  freq$percent <- freq$proportion * 100
  freq <- freq[order(freq$n, decreasing = TRUE), ]
  others <- as.vector(freq$cells[which(freq$n < 50)])
  x@meta.data$immgen.pred.clean <- forcats::fct_collapse(x@meta.data$immgen.pred.clean, Others = others)

  x@meta.data$mrsd.pred.clean <- forcats::fct_collapse(x@meta.data$mrsd.pred,`B cells` = "B cells, pro")
  freq <- as.data.frame(table(x@meta.data$mrsd.pred))
  freq <- freq %>% mutate("prop" = prop.table(Freq))
  colnames(freq) <- c("cells", "n", "proportion")
  freq$percent <- freq$proportion * 100
  freq <- freq[order(freq$n, decreasing = TRUE), ]
  others <- as.vector(freq$cells[which(freq$n < 50)])
  x@meta.data$mrsd.pred.clean <- forcats::fct_collapse(x@meta.data$mrsd.pred.clean, Others = others)

  return(x)
}

marker_predict = function(x,reference){
#Marker Gene Detection
#Classic approach / Cell Level / Bulk RNAseq Reference
classic.pred <- SingleR::SingleR(test = x@assays[["RNA"]]@data,
                ref = reference, labels = reference$label.main,
                fine.tune = TRUE, de.method = "classic",
                de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
#Pairwise Wilcoxon ranked sum test / Cell Level / Single-cell Reference
wilcox.pred <- SingleR::SingleR(test = x@assays[["RNA"]]@data, ref = reference,
                labels = reference$label.main, fine.tune = TRUE, de.method = "wilcox",
                de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
return(list(classic.pred,wilcox.pred))
}

Louvain_Annotation_HM = function(x,assay_full,file_name){
  mrsd_1 <- prop.table(table(x@meta.data$mrsd.pred.clean, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  immgen_2 <- prop.table(table(x@meta.data$immgen.pred.clean, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  azimuth_3 <- prop.table(table(x@meta.data$predicted.ann_level_3, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  pdf(file_name)
  for (obj in objects()[objects() %like% "_[1-9]$"]){
    mhm = ComplexHeatmap::Heatmap(get(obj), cluster_columns = FALSE, show_column_dend = TRUE, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = 12),
                        column_names_gp = grid::gpar(fontsize=12), cluster_rows = TRUE, show_row_dend = FALSE,  col = heatmap_ann.cols, border_gp = gpar(col = "black", lty = 1),
                        rect_gp = gpar(col = "white", lwd = 2), column_title = paste0("Louvian 0.2 clusters vs. ",substr(obj,1,nchar(obj)-2)),  name="Proportion of cells")
    plot(mhm)
  }
  dev.off()
}

make_UMAP_Annotation = function(x,full_assay,pdf_name){
p1 <- DimPlot(x, reduction = "umap", label=T, group.by = full_assay, cols = cluster.cols) + ggtitle(paste0("Louvain clusters - resolution ",substr(full_assay,nchar(full_assay)-3,nchar(full_assay))))
p3 <- DimPlot(x, reduction = "umap", label=T, group.by = "mrsd.pred.clean", cols = mrsd.cols) + ggtitle("MRSD Pred")
p4 <- DimPlot(x, reduction = "umap", label=T, group.by = "immgen.pred.clean", cols = immgen.cols) + ggtitle("ImmGen - cell level")
g1 = p1 + p3 + p4 + patchwork::plot_layout(ncol = 1, nrow = 3)
ggsave(pdf_name,width=15,height=32,plot=g1)
}

theme_90 <- function() {
  theme_bw(base_size=12)+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5),axis.text=element_text(color="black"),panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),strip.background = element_rect(fill="white"))
}
find_plot_markers = function(x,lfc_pdf,annot_pdf,dotplot){
    x = NormalizeData(x,normalization.method="LogNormalize",scale=10000,assay="RNA")
    iter = 0
    res = x@meta.data[,colnames(x@meta.data) %like% "*res.[0-9].[0-9]"]
    for (i in res){
        iter = iter+1
        res_str = str_match(colnames(res)[iter],".[0-9].[0-9]")[1]
        print(paste0("Working on resolution: ",res_str))
        Idents(x) = i
        markers = Seurat::FindAllMarkers(x, assay="RNA", only.pos = TRUE, slot="counts", min.pct = 0.3, logfc.threshold = .5, max.cells.per.ident=5000, random.seed = 888)
        #Organize results into table
        markers %>%
            group_by(cluster) %>%
            top_n(n = 5, wt = avg_log2FC) -> top5
        marker.list = list(AllMarkers=markers, Top5Markers=top5)
        marker.dt = as.data.table(markers)
        smry = marker.dt[, list(N=.N, mean_log2FC_cluster=mean(avg_log2FC), mean_pct_cluster=mean(pct.1)),  by=cluster]
        p1 = ggplot(data=marker.dt, aes(x=avg_log2FC, y=cluster, fill=cluster, color=cluster))+
                    ggridges::geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
                                                point_shape = "|", point_size = 2, size = 0.25, position = position_points_jitter(height = 0), alpha=0.2)+
                    geom_vline(xintercept=0.5, color="grey20", linetype="dashed", lwd=0.2)+
                    labs(y=paste0("Louvain clusters, resolution ",res_str), x="Average log2 fold change of cluster markers")
        p2 = ggplot(data=marker.dt, aes(x=pct.1, y=cluster, fill=cluster, color=cluster))+
            ggridges::geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
                                        point_shape = "|", point_size = 2, size = 0.25, position = position_points_jitter(height = 0), alpha=0.2)+
            labs(y=paste0("Louvain clusters, resolution ",res_str), x="Percentage of cells expressing marker in cluster")+
            geom_vline(xintercept=0.3, color="grey20", linetype="dashed", lwd=.2)
        pw = p1 / p2 + patchwork::plot_layout(guides="collect") & scale_fill_manual(values=cluster.cols) & scale_color_manual(values=cluster.cols) & theme_linedraw()
        ggsave(paste0(lfc_pdf,res_str,".pdf"),height=10,width=10,plot=pw,device='pdf')

        heatmap_data.cols = circlize::colorRamp2(breaks=c(0:5),hcl_palette="Inferno")
        heatmap_scale.cols = circlize::colorRamp2(breaks=c(-3.5,0,3.5),hcl_palette = "Lisbon")
        heatmap_marker.cols = circlize::colorRamp2(breaks=c(0:5), hcl_palette="Inferno")

        marker.genes = top5$gene
        x.hm = subset(x, downsample=100)

        mat = x.hm[["RNA"]]@data[marker.genes, ] %>% as.matrix()

        ra = ComplexHeatmap::rowAnnotation(`Cluster`=as.character(top5$cluster), col=list(`Cluster`=cluster.cols), show_legend=FALSE)

        ta_immgen = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[colnames(res)[iter]]], `ImmGen`=x.hm@meta.data$immgen.pred.clean, col=list(`ImmGen`=immgen.cols,`Cluster`=cluster.cols))
        ta_mrsd = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[colnames(res)[iter]]], `MRSD`=x.hm@meta.data$mrsd.pred.clean, col=list(`MRSD`=mrsd.cols,`Cluster`=cluster.cols))
        ta_azi = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[colnames(res)[iter]]], `L3_Azi`=x.hm@meta.data$predicted.ann_level_3, col=list(`L3_Azi`=L3_azi.cols,`Cluster`=cluster.cols))

        hm_list = list()
        for (ta in objects()[grepl("^ta_",objects())]){
            hm_list[[ta]] = ComplexHeatmap::Heatmap(mat, name = "Expression", cluster_columns = FALSE, cluster_rows=FALSE, column_split=x.hm@meta.data[[colnames(res)[iter]]],
                        row_split=top5$cluster, row_names_gp = grid::gpar(fontsize = 11), column_title=character(0), column_gap = unit(0.5, "mm"),
                        col = heatmap_data.cols, top_annotation = get(ta), show_column_names = FALSE, left_annotation = ra)
        }
        nc = length(unique(markers$cluster))
        pdf(paste0(annot_pdf,res_str,".pdf"),height=round_any(nc,8,f=ceiling),width=round_any(nc,8,f=ceiling))
        lapply(hm_list,ComplexHeatmap::draw,merge_legends = TRUE, padding = unit(c(2, 2, 2, 8), "mm"))
        dev.off()

        p = Seurat::DotPlot(x, assay="RNA", features=unique(marker.genes), group.by=colnames(res)[iter], scale=TRUE)+
            theme(axis.text.x=element_text(size=6), axis.title.x=element_blank())+
            labs(y=paste0("Louvain2 clusters - resolution ",res_str))+ theme_90()
        ggsave(paste0(dotplot,res_str,".pdf"),plot=p,height=8,width=14,device="pdf")
    }
}
feature_plots = function(x, barplot, ridgeplot,assay){ #Doesn't consider integration type i.e. Harm vs SCT
    dat <- x@meta.data
    if (assay == "RNA"){
      p1 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=predicted.ann_level_3))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("L3_Azimuth", values=L3_azi.cols)
      p2 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=mrsd.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("MRSD", values=mrsd.cols)
      p3 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=immgen.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("IMMGEN", values=immgen.cols)
      p4 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=predicted.ann_level_3))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("L3_Azimuth", values=L3_azi.cols)
      p5 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=mrsd.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("MRSD", values=mrsd.cols)
      p6 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=immgen.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("IMMGEN", values=immgen.cols)
      g = cowplot::plot_grid(p1, p2, p3, ncol = 1)
      g2 = cowplot::plot_grid(p4, p5, p6, ncol = 1)
      ggsave(paste0(barplot,"_resolution_0.2.pdf"),plot=g,width=10,height=12,device="pdf")
      ggsave(paste0(barplot,"_resolution_0.5.pdf"),plot=g2,width=10,height=12,device="pdf")
      p1 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="RNA_snn_res.0.2", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.2")
      p2 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="RNA_snn_res.0.5", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.5")
      p = cowplot::plot_grid(p1, p2, ncol = 1, nrow= 2)
      ggsave(paste0(ridgeplot,".pdf"),device='pdf',height=12,width=8,plot=p)
    }
    if (assay=="integrated"){
      p1 = ggplot(data=dat, aes(x=integrated_snn_res.0.2, fill=predicted.ann_level_3))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("L3_Azimuth", values=L3_azi.cols)
      p2 = ggplot(data=dat, aes(x=integrated_snn_res.0.2, fill=mrsd.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("MRSD", values=mrsd.cols)
      p3 = ggplot(data=dat, aes(x=integrated_snn_res.0.2, fill=immgen.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("IMMGEN", values=immgen.cols)
      p4 = ggplot(data=dat, aes(x=integrated_snn_res.0.5, fill=predicted.ann_level_3))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("L3_Azimuth", values=L3_azi.cols)
      p5 = ggplot(data=dat, aes(x=integrated_snn_res.0.5, fill=mrsd.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("MRSD", values=mrsd.cols)
      p6 = ggplot(data=dat, aes(x=integrated_snn_res.0.5, fill=immgen.pred.clean))+
          geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("IMMGEN", values=immgen.cols)
      g = cowplot::plot_grid(p1, p2, p3, ncol = 1)
      g2 = cowplot::plot_grid(p4, p5, p6, ncol = 1)
      ggsave(paste0(barplot,"_resolution_0.2.pdf"),plot=g,width=10,height=12,device="pdf")
      ggsave(paste0(barplot,"_resolution_0.5.pdf"),plot=g2,width=10,height=12,device="pdf")
      p1 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="integrated_snn_res.0.2", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.2")
      p2 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="integrated_snn_res.0.5", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.5")
      p = cowplot::plot_grid(p1, p2, ncol = 1, nrow= 2)
      ggsave(paste0(ridgeplot,".pdf"),device='pdf',height=12,width=8,plot=p)
    }
}
subcluster = function(x,tx,integ,parent){ #tx should be a list() with names of CMO groups attached to CMO labels, i.e. list('model'=c("CMO_1","CMO_2"))
    q = NULL ### Assign every sample to a group name
    for (ele in tx){
        q = c(q,ele)
    }
    iter = 0
    sample2.cols = NULL ###Set sample colors
    for (ele in q){
        iter = iter + 1
        sample2.cols = c(sample2.cols,ele = sample.cols[iter])
    }
    names(sample2.cols) = q ###Ugly way to set names

    for (resolution in colnames(x@meta.data)[colnames(x@meta.data) %like% "*res.[0-9].[0-9]"]){
        for (i in unique((x@meta.data[[resolution]]))){
            eval(parse(text=paste0("sub",i,"=subset(x,subset=`",resolution,"`=='",i,"')"))) #Programatically assign subcluster objects to subclusters
        }
        for (i in objects()[objects() %like% "sub\\d{1}"]){
            iter = as.character(str_match(i,"\\d{1,}"))
            wod  = paste0(parent,"/",integ,"/",resolution,"/subcluster",iter,"/")
            print(wod)
            dir.create(file.path(wod),recursive = TRUE)
            setwd(wod)
            sub_temp = get(i)
            DefaultAssay(sub_temp) = "RNA"
            print(paste0("Starting ",i," dimensional reduction"))
            sub_temp <- Seurat::NormalizeData(sub_temp, normalization.method = "LogNormalize", scale.factor = 10000, verbose=FALSE)
            sub_temp <- Seurat::FindVariableFeatures(sub_temp, selection.method="vst", nfeatures=2000)
            sub_temp <- Seurat::ScaleData(sub_temp, verbose=FALSE)
            sub_temp <- Seurat::RunPCA(sub_temp, assay="RNA",npcs=bdim(sub_temp))
            g1=DimHeatmap(sub_temp, dims=1:9, nfeatures=10, fast=TRUE, slot="scale.data", assays="RNA")
            ggsave(paste0("Cluster_",iter,"_PCA_Loadings_harm_HM.pdf"),width=8,height=8,plot=g1)
            print(paste0("Dim Reduction Complete"))
            q = DimPlot(sub_temp, group.by="orig.ident", reduction="pca") + scale_color_manual(values=sample.cols)
            ggsave(paste0("PCA_by_CMO_harm_Cluster_",iter,".pdf"),plot=q,device='pdf')

            sub_temp_meta <- sub_temp@meta.data
            sub_temp_meta$barcode <- unlist(tstrsplit(rownames(sub_temp_meta),"_",keep=1))
            #sub_temp_meta$cell_id <- ifelse(sub_temp_meta$orig.ident=="untreated",
            #                          paste0("Untreated_", sub_temp_meta$barcode),
            #                          paste0("Treated_", sub_temp_meta$barcode))
            sub_temp_meta$rn <- rownames(sub_temp_meta)
            p1 <- DimPlot(sub_temp, group.by="orig.ident", reduction="pca") +
                    labs(title="Treatment") +
                    scale_color_manual("Treatment", values=sample.cols)

            p2 <- DimPlot(sub_temp, group.by="mrsd.pred", reduction="pca") +
                    labs(title="MRSD annotation") +
                    scale_color_manual("Annotation", values=mrsd.cols)

            p3 <- DimPlot(sub_temp, group.by="immgen.pred", reduction="pca") +
                    labs(title="ImmGen annotation") +
                    scale_color_manual("Annotation", values=immgen.cols)

            p4 <- DimPlot(sub_temp, group.by="predicted.ann_level_3", reduction="pca") +
                    labs(title="L3_Azi Annot") +
                    scale_color_manual("L3_Azi", values=L3_azi.cols)

            p5 <- DimPlot(sub_temp, group.by="Phase", reduction="pca") +
                    labs(title="Cell cycle") +
                    scale_color_manual("Cell cycle", values=cell_cycle.cols)

            pw <- p1 + p2 + p3 + p4 + p5 + patchwork::plot_layout(ncol=2,nrow=3)
            ggsave(paste0("Sub_",iter,"_PCA_harm.pdf"),width=15,height=10,plot=pw)
            sub_temp@meta.data$orig.ident = factor(sub_temp@meta.data$orig.ident)
            if (nrow(sub_temp@meta.data) < 30){
              neighbors = round(nrow(sub_temp@meta.data)*.10,0)
            }else{
              neighbors = 30
            }
            sub_temp <- harmony::RunHarmony(sub_temp, dims=1:bdim(sub_temp), group.by.vars = "orig.ident")
            sub_temp <- Seurat::RunUMAP(sub_temp, assay.use="RNA", dims=1:bdim(sub_temp), reduction="harmony",n.neighbors=neighbors)
            sub_temp <- Seurat::FindNeighbors(sub_temp, reduction="harmony", dims=1:bdim(sub_temp))
            sub_temp <- Seurat::FindClusters(sub_temp, algorithm=2, resolution=c(0.1, 0.2))

            p1 <- DimPlot(sub_temp, group.by="RNA_snn_res.0.1")+
            rcartocolor::scale_color_carto_d(palette = "Bold")

            p2 <- DimPlot(sub_temp, group.by="predicted.ann_level_3")+
            labs(title="L3_Azimuth Prediction") +
            scale_color_manual("L3_Azi", values=L3_azi.cols)

            p3 <- DimPlot(sub_temp, group.by="RNA_snn_res.0.2")+
            rcartocolor::scale_color_carto_d(palette = "Bold")

            p4 <- DimPlot(sub_temp, group.by="orig.ident")+
            labs(title="Treatment")+
            scale_color_manual("Treatment", values=sample.cols)

            p5 <- FeaturePlot(sub_temp, features="percent.mt")+
            viridis::scale_color_viridis(option="mako")

            p6 <- FeaturePlot(sub_temp, features="Phase")+
            labs(title="Cell cycle")+
            scale_color_manual("Cell cycle", values=cell_cycle.cols)

            pw <- p1 + p2 + p3 + p4 + p5 + p6 + patchwork::plot_layout(ncol=2, nrow=3)
            ggsave(paste0("Sub_",iter,"_subclustering_All.pdf"),height=10,width=10,plot=pw)

            #Single-cell differential GEX approach (Not as good as "pseudobulk")
            Idents(sub_temp) <- sub_temp@meta.data$orig.ident
            sub_temp@meta.data$orig.ident=factor(sub_temp@meta.data$orig.ident)
            c(v1,v2) %<-% c(1,2)
            flag = FALSE
            scde_list = list()
            scde_flag = TRUE
            for (i in names(tx)){
              tx[[i]] = tx[[i]][tx[[i]] %in% levels(sub_temp@meta.data$orig.ident)]
            }
            scde_flag = tryCatch({
            for (i in 1:(length(names(tx))*(length(names(tx))-1)/2)){
                if (flag == FALSE){
                    v2_og = v2
                    flag = TRUE
                }
                scde_results <- Seurat::FindMarkers(sub_temp, slot="data",
                ident.1=tx[[v1]],
                ident.2=tx[[v2]],
                logfc.threshold=0.25,test.use="wilcox", min.pct=0.1, min.diff.pct=-Inf, only.pos=FALSE)
                scde_results <- cbind(gene=rownames(scde_results),scde_results)
                scde_results$Sig <- data.table::fcase(scde_results$p_val_adj <= 0.05, "FDR",
                                                    scde_results$p_val <= 0.05, "Nominal",
                                                    default="NS")
                scde_list[[paste0(names(tx)[v1],".",names(tx)[v2])]] = scde_results
                if (v2 < length(tx)){
                    v2 = v2+1
                }
                else{
                    v1 = v1+1
                    v2_og = v2_og+1
                    v2 = v2_og}
                }
                scde_flag = TRUE},
                error=function(cond){
                  scde_flag = FALSE
                  message(paste0("SCDE failure: Cluster ", iter, " // GROUPS: ",names(tx)[v1]," ",names(tx)[v2]))
                  message(paste0(cond))
                  return(scde_flag)
                },
                warning=function(cond){
                  message(cond)
                },
                finally={}
            )

            #Pseudobulk GEX approach
            ##Checking data metrics
            dat <- as.data.table(sub_temp@meta.data)
            dat[, .N, by=orig.ident]

            counts <- Seurat::AggregateExpression(sub_temp,assays="RNA",group.by="orig.ident",slot="counts")$RNA
            filt_counts <- counts[rowSums(counts>=10)>1,]

            meta <- sub_temp@meta.data
            meta <- meta[!duplicated(meta$orig.ident),]
            rownames(meta) <- meta$orig.ident
            colnames(counts) = meta$orig.ident
            colnames(filt_counts) = meta$orig.ident
            meta <- meta[match(colnames(counts), rownames(meta)),]
            # set factor levels for DE analysis
            meta$group <- factor(meta$orig.ident)#, levels=levels(sub_temp$orig.ident))
            dds <- DESeq2::DESeqDataSetFromMatrix(countData=filt_counts, colData=meta, design=~group)
            dds <- DESeq2::estimateSizeFactors(dds, type="poscounts")
            lognorm_counts <- log2(DESeq2::counts(dds, normalized=TRUE)+1)
            mlognorm_counts <- reshape2::melt(lognorm_counts)

            p = ggplot(data=mlognorm_counts, aes(x=value, y=Var2, fill=Var2, color=Var2)) + ggridges::geom_density_ridges(alpha=0.8)+
                theme_linedraw() + labs(x="log2 (normalized counts) + 1", y="") + scale_fill_manual("Condition", values=sample.cols) + scale_color_manual("Condition", values=sample.cols)
            ggsave(paste0("Sub_",iter,"_Lognorm_Counts.pdf"),plot=p)

            pdf(paste0("Sub_",iter,"_Expression_PCA.pdf"))
            pca <- FactoMineR::PCA(t(lognorm_counts), scale.unit=TRUE, ncp=10, graph=FALSE)
            plot(factoextra::fviz_screeplot(pca))
            plot(factoextra::fviz_pca_ind(pca))
            dev.off()

            if (scde_flag == TRUE){
              for (comp in names(scde_list)){
                  hm_results = scde_list[[comp]][!is.na(scde_list[[comp]]$gene),]
                  hm_results = hm_results[order(hm_results$p_val, decreasing=FALSE),]
                  if (nrow(hm_results)>=75){
                    hm_results = hm_results[1:75,]
                  }else{
                    hm_results = hm_results[1:nrow(hm_results),]
                  }
                  de_genes = hm_results$gene

                  mat = sub_temp[["RNA"]]@data[de_genes, ] %>% as.matrix()
                  mat = t(scale(t(mat)))
                  top_ann = ComplexHeatmap::HeatmapAnnotation(`Sample`=sub_temp@meta.data$orig.ident, col=list(`Sample`=sample2.cols), annotation_name_gp = grid::gpar(fontsize=9, fontface="bold"))

                  # row colors
                  perm <- match(rownames(mat), hm_results$gene)
                  hm_results <- hm_results[perm,]
                  row_ann <- ComplexHeatmap::rowAnnotation(`Significance`=hm_results$Sig, col=list(`Significance`=sig.cols), annotation_name_gp = grid::gpar(fontsize=9, fontface="bold"))

                  pdf(paste0("Sub_",iter,"_",comp,"_Significance_Expression_HM.pdf"),height=15,width=15)
                  ch <- ComplexHeatmap::Heatmap(mat, name = "Expression", cluster_columns = TRUE, column_split=sub_temp@meta.data$orig.ident,
                                          cluster_rows=TRUE, row_split=hm_results$Sig, col=heatmap_scale.cols, row_names_gp = grid::gpar(fontsize = 11),
                                          column_title=character(0), top_annotation = top_ann, show_column_names = FALSE, left_annotation = row_ann)
                  ComplexHeatmap::draw(ch, merge_legends = TRUE, padding = unit(c(2, 2, 2, 8), "mm"))
                  dev.off()

                  pdf(paste0("Sub_",iter,"_",comp,"_Vln_Top5_Sigdiff.pdf"),width=15)
                  top5_genes <- hm_results$gene[1:6]
                  q = VlnPlot(sub_temp, features=top5_genes, assay="RNA", group.by="orig.ident") &
                      scale_fill_manual(values=sample.cols) & theme_linedraw() & theme(axis.text.x=element_blank(), axis.title.x=element_blank())
                  plot(q)
                  dev.off()

                  x_axis_lim <- abs(max(scde_list[[comp]]$avg_log2FC))+0.75
                  p <- ggplot(data=scde_list[[comp]],
                              aes(x=avg_log2FC, y=-log10(p_val), color=Sig))+
                              geom_point(size=2.4, pch=19, alpha=0.7)+
                              labs(x="Average Log2 Fold Change", y="-log10(pvalue)")+
                              xlim(-x_axis_lim, x_axis_lim)+
                              scale_color_manual("Significance", values=sig.cols)


                  #GSEA geneset expression analysis
                  pathways <- fgsea::gmtPathways(paste0(parent,"/GMT.gmt"))
                  ranks <- scde_list[[comp]]$avg_log2FC
                  names(ranks) <- scde_list[[comp]]$gene
                  ranks <- ranks[order(ranks, decreasing=TRUE)]
                  gsea_results <- fgsea::fgsea(pathways=pathways,
                                              stats=ranks,
                                              minSize=1)
                  setorder(gsea_results, pval, na.last=TRUE)
                  fwrite(gsea_results,paste0("Sub_",iter,"_",comp,"_GSEA_HALLMARK.csv"))
                  gsea_plot <- gsea_results
                  scale_lim <- max(abs(gsea_plot$NES))
                  gsea_plot$Pathway <- gsub("HALLMARK_", "", gsea_plot$pathway)
                  gsea_plot$pval <- ifelse(gsea_plot$padj > 0.05, NA, gsea_plot$padj)
                  p <- ggplot(data=gsea_plot, aes(x=NES, y=reorder(Pathway, NES), fill=-log10(pval)))+
                  geom_bar(stat="identity", aes(color=""))+
                  scale_fill_viridis_c(option="mako", na.value="#DBE2E9")+
                  scale_color_manual(values="#DBE2E9")+
                  guides(color=guide_legend("Not significant", override.aes=list(fill="#DBE2E9")))+
                  theme_linedraw()
                  pdf(paste0("Sub_",iter,"_",comp,"_Barplot_NES_Paths.pdf"))
                  plot(p)
                  dev.off()
                }
              }else{
                message("Skipping GSEA due to SCDE failure.")
              }
            setwd(parent)
        }
    remove(list=objects()[objects() %like% "^sub\\d{1,}"])
    }
}

base="/path/to/outs/per_sample_outs/"
parent=getwd()
for (sample in list.files(base)){ #Generate Seurat objects and assign percent.mt
  print(paste0("S",sample,"=ReadMtx(cells='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz',features='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/features.tsv.gz',mtx='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz'"))
  eval(parse(text=paste0("S",sample,"=ReadMtx(cells='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz',features='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/features.tsv.gz',mtx='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz')")))
  eval(parse(text=paste0("S",sample,"=rem_cmo(S",sample,")")))
  eval(parse(text=paste0("S",sample,"_Seurat=CreateSeuratObject(S",sample,",project='CMO_",sample,"',assay='RNA')")))
  eval(parse(text=paste0("S",sample,"_Seurat=PercentageFeatureSet(S",sample,"_Seurat,pattern = '^mt-', col.name='percent.mt')")))
}


### Add a function here to find variable features, good to know before we normalize everything ###
sample.cols = paletteer_d("ggsci::category20_d3")[1:length(list.files(base))]

store = objects()[objects() %like% "Seurat$"]
iter = 0
for (obj in store){
  iter = iter+1
  if (iter > (length(store))){
    break
  }
  if (iter >= 3){
    S_Combined = merge(S_Combined,y=get(store[iter]),project="REPC")
  }
  if (iter == 1){
    S_Combined = merge(get(store[1]),y=get(store[2]),project="REPC")
  }
}

S_Combined = PercentageFeatureSet(S_Combined, pattern = "^mt-|^MT-", col.name = "percent.mt")  # Mitochondria % for merged set
pdf("Merged_preQC_Violin.pdf",width=15)
VlnPlot(object = S_Combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols=sample.cols, pt.size=0, ncol = 3) & geom_jitter(alpha = 0.25, size = 0.1) #preQC stats
dev.off()

nUMI.low = 500                   #QC SECTION
nUMI.high = 40000
percent_mito.low = -Inf
percent_mito.high = 20
nGene.low = 250
nGene.high = 5000

S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$nFeature_RNA > nGene.high | S_Combined@meta.data$nFeature_RNA < nGene.low, "Remove", "Keep")
S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$percent.mt > percent_mito.high, "Remove", S_Combined@meta.data$Filter)
S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$nCount_RNA > nUMI.high | S_Combined@meta.data$nCount_RNA < nUMI.low, "Remove", S_Combined@meta.data$Filter)

filter.cols = c(`Keep`="grey25", `Remove`="firebrick")
S_Combined@meta.data$orig.ident = factor(S_Combined@meta.data$orig.ident)
S_Combined@meta.data$orig.ident = factor(S_Combined@meta.data$orig.ident)
p1 = ggplot(data = S_Combined@meta.data,
       aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
p2 = ggplot(data = S_Combined@meta.data,
       aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
p3 = ggplot(data = S_Combined@meta.data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
g1 = cowplot::plot_grid(p1, p2, p3, ncol = 3)
p1 = ggplot(data = S_Combined@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = Filter)) + geom_point() +
  geom_hline(yintercept = c(nGene.high, nGene.low)) + geom_vline(xintercept = c(nUMI.high, nUMI.low)) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position = "none")
p2 = ggplot(data = S_Combined@meta.data, aes(x = percent.mt, y = log10(nFeature_RNA), color = Filter)) + geom_point() +
  geom_vline(xintercept = c(percent_mito.high)) + geom_hline(yintercept = c(log10(nGene.high), log10(nGene.low))) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position ="none")
p3 = ggplot(data = S_Combined@meta.data, aes(x = percent.mt, y = log10(nCount_RNA), color = Filter)) + geom_point() +
  geom_vline(xintercept = c(percent_mito.high)) + geom_hline(yintercept = c(log10(nUMI.high), log10(nUMI.low))) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position = "none")
g2 = cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(filename="Violin_Features.pdf",plot=g1,height=10,width=18)
ggsave(filename="Scatter_Features.pdf",plot=g2,height=10,width=18)

S_Combined = subset(x = S_Combined, subset = nFeature_RNA > nGene.low & nFeature_RNA < nGene.high & nCount_RNA > nUMI.low & nCount_RNA < nUMI.high & percent.mt < percent_mito.high)
pdf("Merged_postQC_Violin.pdf",width=15) #POSTQC VIOLIN
VlnPlot(object = S_Combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols=sample.cols, pt.size=0, ncol = 3) & geom_jitter(alpha = 0.25, size = 0.1) #preQC stats
dev.off()

S_Combined = aziPredict(S_Combined,"lungref")

S_Combined = Seurat::NormalizeData(S_Combined, normalization.method = "LogNormalize", scale.factor = 10000)
S_Combined = Seurat::FindVariableFeatures(S_Combined, selection.method = "vst", nfeatures = 4000)
S_Combined = Seurat::ScaleData(S_Combined)
S_Combined=RunPCA(S_Combined,npcs=50) #Test 50 dimensions
bestdims_c = bdim(S_Combined)

h2m.df = read.table(paste0(getwd(),"/HMD_HumanPhenotype.rpt.txt"), sep = "\t") #Cell Cycle Scoring
h2m.df = h2m.df[, c("V1","V3")]
colnames(h2m.df) = c("Gene.name", "ortholog_name")
h2m.df = h2m.df[!duplicated(h2m.df$Gene.name), ]
h2m.df = h2m.df[!duplicated(h2m.df$ortholog_name), ]
mmus_s = h2m.df[which(h2m.df$Gene.name%in%cc.genes$s.genes), ]$ortholog_name
mmus_g2m = h2m.df[which(h2m.df$Gene.name%in%cc.genes$g2m.genes), ]$ortholog_name
S_Combined = CellCycleScoring(S_Combined, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = FALSE) #Species specific

myl = SplitObject(S_Combined,split.by="orig.ident")
myl = lapply(myl,SCTransform,method='glmGamPoi',vst.flavor='v2',vars.to.regress='percent.mt',return.only.var.genes=FALSE,verbose=TRUE)
features = SelectIntegrationFeatures(object.list = myl, nfeatures = 4000)
myl = PrepSCTIntegration(object.list = myl, anchor.features = features)
anchors = FindIntegrationAnchors(object.list = myl, normalization.method = "SCT", anchor.features = features)
S_Combined.SCT = IntegrateData(anchorset = anchors, normalization.method = "SCT")  ###This can bug out if cells per sample is < 100. Apply k.weight = min(sample)
S_Combined.SCT=RunPCA(S_Combined.SCT,npcs=50)
bestdims_sct = bdim(S_Combined.SCT)
S_Combined.SCT <- RunUMAP(S_Combined.SCT, dims = 1:bestdims_sct, verbose = FALSE)
S_Combined.SCT <- FindNeighbors(S_Combined.SCT, dims = 1:bestdims_sct, verbose = FALSE)
S_Combined.SCT <- FindClusters(S_Combined.SCT, resolution=0.2,verbose = FALSE)

S_Combined@meta.data$CC.Diff = S_Combined@meta.data$S.Score - S_Combined@meta.data$G2M.Score ##PCA Cell Cycle Scores
S_Combined = RunPCA(S_Combined, npcs = bestdims_c, verbose = FALSE)
S_Combined@meta.data$Full = "All CMOs"
pc.df = cbind(S_Combined@reductions$pca@cell.embeddings, S_Combined@meta.data)
p1 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p2 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
p3 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2M.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p4 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2m.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
pdf("Cell_Cycle_Scores_PCA.pdf",width=21,height=12)
p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CELL CYCLE SCORES PCA')}

S_Combined = RunUMAP(S_Combined, reduction = "pca", dims = 1:bestdims_c, verbose = F) ##UMAP Cell Cycle Scores
umap.df = cbind(S_Combined@reductions$umap@cell.embeddings, S_Combined@meta.data)
p1 = ggplot(data = umap.df, aes(x = umap_1, y = umap_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p2 = ggplot(data = umap.df, aes(x = umap_1, y = umap_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
p3 = ggplot(data = umap.df, aes(x = umap_1, y = umap_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2M.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p4 = ggplot(data = umap.df, aes(x = umap_1, y = umap_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2m.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
pdf("Cell_Cycle_Scores_UMAP.pdf",width=21,height=12)
p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CELL CYCLE SCORES UMAP')}

pdf("Cells_By_Phase_Barplot.pdf",width=8)
cell_cycle.cols <- c(`G1`="#003049", `G2M`="#d62828", `S`="#f77f00")
ggplot(data = S_Combined@meta.data, aes(x = orig.ident, fill = Phase)) +
  geom_bar(position = "fill", color = "black") + scale_fill_manual("", values = cell_cycle.cols) + xlab("") + ylab("Proportion of cells")
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CELLS BY PHASE BARPLOT')}

pdf("Cells_By_Phase_PCA_UMAP.pdf",width=16)
p1 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, reduction = "pca", split.by = "Full") + ggtitle("")
p2 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, reduction = "pca", split.by = "orig.ident") + ggtitle("")
a1 = p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides = "collect")
p3 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, split.by = "Full") + ggtitle("")
p4 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, split.by = "orig.ident") + ggtitle("")
a2 = p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides = "collect")
ggsave(filename="Cells_By_Phase_PCA.pdf",plot=a1,width=16)
ggsave(filename="Cells_By_Phase_UMAP.pdf",plot=a2,width=16)

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CELLS BY PHASE UMAP')}

S_Combined = RunPCA(S_Combined, npcs = bestdims_c, verbose = FALSE)
S_Combined = RunUMAP(S_Combined, reduction = "pca", dims = 1:bestdims_c, verbose = F)
S_Combined = FindNeighbors(S_Combined, reduction = "pca", dims = 1:bestdims_c)
S_Combined = FindClusters(S_Combined, resolution = 0.2)

S_Combined.harm = S_Combined %>%
    RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE)
S_Combined.harm = RunPCA(S_Combined.harm, npcs = 50, verbose = FALSE)
S_Combined.harm = RunUMAP(S_Combined.harm, reduction = "harmony", dims = 1:bdim(S_Combined.harm))
S_Combined.harm = FindNeighbors(S_Combined.harm, reduction = "harmony", dims = 1:bdim(S_Combined.harm))
S_Combined.harm = FindClusters(S_Combined.harm, resolution = 0.2)

objs = list(S_Combined=S_Combined,S_Combined.harm=S_Combined.harm,S_Combined.SCT=S_Combined.SCT)

p1 = DimPlot(S_Combined, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("Merged + logNorm")
p2 = DimPlot(S_Combined.harm, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("logNorm + PCA + Harmony")
p3 = DimPlot(S_Combined.SCT, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("SCT + CCA")
g1 = p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave(filename="UMAP_CMOs.pdf",plot=g1,width=15,height=11)

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: UMAP CMO')}

p1.1 = DimPlot(S_Combined, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("Merged + logNorm")
p2.1 = DimPlot(S_Combined.harm, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("logNorm + PCA + Harmony")
p2.2 = DimPlot(S_Combined.SCT, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("SCT + CCA")
g1 = p1.1 + p2.1 +  p2.2 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave(filename="PCA_CMOs.pdf",plot=g1,width=15,height=11)

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: PCA CMO')}

tab <- as.data.frame(table(S_Combined@meta.data$RNA_snn_res.0.2, S_Combined@meta.data$orig.ident))
p1.2 <- ggplot(tab, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("Merged + log-Norm")
tab2 <- as.data.frame(table(S_Combined.harm@meta.data$RNA_snn_res.0.2, S_Combined.harm@meta.data$orig.ident))
p2.2 <- ggplot(tab2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("log-Norm + Harmony")
tab3 <- as.data.frame(table(S_Combined.SCT@meta.data$integrated_snn_res.0.2, S_Combined.SCT@meta.data$orig.ident))
p3.2 <- ggplot(tab3, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("SCT + CCA")
g1 = p1.2 + p2.2 + p3.2 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave("CMO_Barplot.pdf",width=15,height=11,plot=g1)

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CMO BARPLOT')}

heatmap.cols <- circlize::colorRamp2(c(0, 1), c("white", "black"))
heatmap_ann.cols <- circlize::colorRamp2(c(0, 0.1, 1), c("#EDD9A3", "#EA4F88", "#4B2991"))

pdf("CMO_Cluster_Proportions.pdf",width=10)
hm1 = hm_prop(S_Combined)
hm2 = hm_prop(S_Combined.harm)
hm3 = hm_prop(S_Combined.SCT)
hm1 + hm2 + hm3+ patchwork::plot_layout(ncol = 3, guides="collect")
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: CMOCLUSTERPROP')}

pdf("By_PC_Axis.pdf",width=12,height=12)
lapply(objs,function(x){
    qr = names(x@commands)
    if(length(qr[grepl("harmony",qr)])>=1){title_wow="Harmony"}
    if(length(qr[grepl("integrated",qr)])>=1){title_wow="SCT"}
    if(length(qr[grepl("harmony|integrated",qr)])==0){title_wow="Non-Integrated"}
    pc <- as.data.frame(x[["pca"]]@cell.embeddings)
    pc$orig.ident <- x@meta.data$orig.ident
    q = GGally::ggpairs(pc, columns = 1:4, upper = "blank",
                            mapping = ggplot2::aes(color = orig.ident, alpha = 0.5),
                            progress = FALSE) +
        scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = sample.cols) + ggtitle(title_wow)
    q
})
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: BY PC AXIS')}

# Begin LOUVAIN clustering #
S_Combined.harm = Seurat::FindNeighbors(S_Combined.harm, reduction = "harmony", dims = 1:bdim(S_Combined.harm))
S_Combined.SCT = Seurat::FindNeighbors(S_Combined.SCT,reduction="pca",dims=1:bestdims_sct)
S_Combined.harm = Seurat::FindClusters(S_Combined.harm, algorithm = 1, resolution = c(0.2, 0.5))
S_Combined.SCT = Seurat::FindClusters(S_Combined.SCT, algorithm = 1, resolution = c(0.2, 0.5))

pdf("Cluster_Louvain.pdf",width=15)
cowplot::plot_grid(ncol = 2,
    DimPlot(S_Combined.harm, reduction = "umap", label = T, group.by = "RNA_snn_res.0.2") + ggtitle("Louvain Harm - res 0.2"),
    DimPlot(S_Combined.harm, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5") + ggtitle("Louvain Harm - res 0.5")
)
cowplot::plot_grid(ncol = 2,
    DimPlot(S_Combined.SCT, reduction = "umap", label = T, group.by = "integrated_snn_res.0.2") + ggtitle("Louvain SCT - res 0.2"),
    DimPlot(S_Combined.SCT, reduction = "umap", label = T, group.by = "integrated_snn_res.0.5") + ggtitle("Louvain SCT - res 0.5")
)
dev.off()

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: ClusterLouvain')}

## Annotation + Clean Annotations ##
c(S_Combined.SCT,immgen.ref,mrsd.ref) %<-% annotate(S_Combined.SCT,"integrated")[c(1,2,3)]
c(S_Combined.harm,immgen.ref,mrsd.ref) %<-% annotate(S_Combined.harm,"RNA")[c(1,2,3)]
S_Combined.SCT = clean(S_Combined.SCT)
S_Combined.harm =  clean(S_Combined.harm)

#Marker Gene Detection - [1] if using bulk-RNAseq reference, [2] if using sc-RNAseq reference
c(preds.harm.immgen,preds.sct.immgen) %<-% c(marker_predict(S_Combined.harm,immgen.ref)[1],marker_predict(S_Combined.SCT,immgen.ref)[1])
c(preds.harm.mrsd,preds.sct.mrsd) %<-% c(marker_predict(S_Combined.harm,mrsd.ref)[1],marker_predict(S_Combined.SCT,mrsd.ref)[1])

Louvain_Annotation_HM(S_Combined.SCT,"integrated_snn_res.0.2","Annot_Cluster_Prop_SCT.pdf")
Louvain_Annotation_HM(S_Combined.harm,"RNA_snn_res.0.2","Annot_Cluster_Prop_Harm.pdf")

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: LouvainAnnotHM')}

cluster.cols <- pals::polychrome(n = 30) #Supports up to 30 cluster colors
names(cluster.cols) = as.character(0:29)
c(cluster.cols["1"],cluster.cols["4"],cluster.cols["8"],cluster.cols["9"]) %<-% c("#CDC1C5","#43CD80","#97FFFF","#9AFF9A")
# immgen annotation colors
immgen.cols = c(`B cells` = "#F0A0FF", `Basophils` = "#9370DB",
                `Eosinophils`= "#B452CD", `DC` = "#9DCC00",  `Epithelial cells` = "#87CEFA",   `Endothelial cells` = "#6495ED", `Fibroblasts` = "navy",
                `ILC` = "#191919", `Macrophages` = "#1C8356", `Mast cells` = "#E9Debb", `Monocytes` = "#81C57A",`Microglia` = "#EEB4B4",
                `Neutrophils` = "#b10da1",`NK cells` = "#8F7C00", `NKT` = "#FEAF16", `Others` = "gray40", `Stem cells` = "red", `Stromal cells` = "#F08080",`T cells` = "#FE902EFF", `Tgd` = "#FFEE33")

mrsd.cols = c('Adipocytes'=paletteer_d("ggsci::category20_d3")[1],'Astrocytes'=paletteer_d("ggsci::category20_d3")[2],'B cells'=paletteer_d("ggsci::category20_d3")[3],'Cardiomyocytes'=paletteer_d("ggsci::category20_d3")[4],
'Dendritic cells'=paletteer_d("ggsci::category20_d3")[5],'Endothelial cells'=paletteer_d("ggsci::category20_d3")[6],'Epithelial cells'=paletteer_d("ggsci::category20_d3")[7],'Erythrocytes'=paletteer_d("ggsci::category20_d3")[8],'Microglia'=paletteer_d("ggsci::category20_d3")[9],
'Fibroblasts'=paletteer_d("ggsci::category20_d3")[10],'Granulocytes'=paletteer_d("ggsci::category20_d3")[11],'Hepatocytes'=paletteer_d("ggsci::category20_d3")[12],'Macrophages'=paletteer_d("ggsci::category20_d3")[13],'Monocytes'=paletteer_d("ggsci::category20_d3")[14],
'Neurons'=paletteer_d("ggsci::category20_d3")[15],'NK cells'=paletteer_d("ggsci::category20_d3")[16],'Oligodendrocytes'=paletteer_d("ggsci::category20_d3")[17],'T cells'=paletteer_d("ggsci::category20_d3")[18],'Others'=paletteer_d("ggsci::category20_d3")[19])

L3_azi.cols = ""
namae2 = unique(S_Combined@meta.data$predicted.ann_level_3)
for (i in 1:length(namae2)){L3_azi.cols[namae2[i]]=paletteer_d("ggsci::category20_d3")[i]}
L3_azi.cols = L3_azi.cols[2:length(L3_azi.cols)]


make_UMAP_Annotation(S_Combined.harm,"RNA_snn_res.0.2","UMAP_ANNOT_HARM.pdf")
make_UMAP_Annotation(S_Combined.SCT,"integrated_snn_res.0.2","UMAP_ANNOT_SCT.pdf")
if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: UMAP_ANNOT')}

make_UMAP_Azimuth(S_Combined.SCT,"UMAP_AZIMUTH_SCT")
make_UMAP_Azimuth(S_Combined.harm,"UMAP_AZIMUTH_HARM")

if (TRUE %in% grepl("Rplots.pdf",list.files())){print('DEBUG ::: UMAP_AZI')}

# Annotation Confidence Heatmap + Delta Scores #
myl2 = list()
myl3 = list()
for (obj in objects()[grepl("preds.\\D{1,}\\.",objects())]){
  namae = str_split_1(obj,"\\.")
  myl2[[obj]] = SingleR::plotScoreHeatmap(get(obj),show.pruned=TRUE,main=paste(namae[2],namae[3]))[[4]]
  myl3[[obj]] = SingleR::plotDeltaDistribution(get(obj)) + ggtitle(paste(namae[2],namae[3]))
}
matr = matrix(c(1,2,3,4),nrow=2,ncol=2)
g=grid.arrange(arrangeGrob(grobs= myl2,ncol=2,layout_matrix=matr))
g2=grid.arrange(arrangeGrob(grobs= myl3,ncol=2,layout_matrix=matr))
ggsave("HM_CONF.pdf",plot=g,width=16,height=15,device='pdf')
ggsave("Delta_CONF.pdf",plot=g2,width=16,height=15,device='pdf')

sig.cols <- c(`FDR`="#FF00CC", `Nominal`="#ffb3fd", `NS`="#C8C8CD")
heatmap_data.cols = circlize::colorRamp2(breaks=c(0:5),hcl_palette="Inferno")
heatmap_scale.cols = circlize::colorRamp2(breaks=c(-3.5,0,3.5),hcl_palette = "Lisbon")
heatmap_marker.cols = circlize::colorRamp2(breaks=c(0:5), hcl_palette="Inferno")

find_plot_markers(S_Combined.harm,"HARM_LFC","HARM_ANNOT_HM","Dotplot_Features_Harm")
find_plot_markers(S_Combined.SCT,"SCT_LFC","SCT_ANNOT_HM","Dotplot_Features_SCT")

feature_plots(S_Combined.harm,"Annotation_Cluster_Barplot_Harm","Ridgeplot_Features_Harm","RNA")
feature_plots(S_Combined.SCT,"Annotation_Cluster_Barplot_SCT","Ridgeplot_Features_SCT","SCT")

tx = list("Mock"=c("CMO_1","CMO_2","CMO_3","CMO_4"),"Infected"=c("CMO_5","CMO_6","CMO_7","CMO_8"))
subcluster(S_Combined.harm,tx=tx,integ="Harmony",parent)
subcluster(S_Combined.SCT,tx=tx,integ="SCT",parent)

module_genes <- list(Monocyte=c("Ly6c+","Ly6c-","Lyz2"),AT1=c("Akap5"),AT1_2=c("Akap5","Sftpd","Slc34a2","Lamp3"),Secretory=c("Sftpd","Cyp2f2","Scgb3a2"),
Basal=c("Krt5","Krt13"),Ciliated=c("Tmem212","Sntn"),AF1=c("Tcf21","Adh1","Dpt"),AF2=c("Mfap5","Dpt"),Pericyte=c("Higd1b","Cox4i2"),AEC=c("Dkk2","Gja5"),
AT2=c("Sftpd","Lamp3","Slc34a2","Lyz2"),CD4_T=c("Ms4a4b","Ccr7"),CD8_T=c("Ms4a4b","Cd8b1","Cd8a"),ILC=c("Cxcr6","Arg1"),NK=c("Gzma","Klrk1","Ms4a4b"),
B_Cell=c("Ms4a1","Igkc"),DCs=c("H2-Aa","ifi205","Ccl22","Ccr7"),Neutrophil=c("Stfa2","Lyz2","Gm5483"),Basophil=c("Mcpt8","Cd200r3"),Mast=c("Mcpt4","Kit","Tpsb2"),
Macrophage=c("F7","Atp6v0d2","Lyz2","C1qa","C1qc"),Megakaryocyte=c("Gp5","Gp6"))
#problem genes LY6C,LYZ2,LAMP3,CYP2F2,KRT13,ADH1,MS4A4B,CD8B1,H2-AA,IFI205,STFA2,GM5483,MCPT8,CD200R3,MCPT4,TPSB2
#problem cell types Monocyte, Basal, Neutrophil, Basophil
#module_genes = lapply(module_genes,toupper)
S_Combined.harm <- UCell::AddModuleScore_UCell(S_Combined.harm,features=module_genes, assay="RNA",name=NULL)
S_Combined.SCT <- UCell::AddModuleScore_UCell(S_Combined.SCT,features=module_genes, assay="integrated",name=NULL)
for (name in names(module_genes)){
  p1 <- FeaturePlot(S_Combined.harm, features=name, order=TRUE) + labs(title=paste0(name," Markers"))
  p2 <- DimPlot(S_Combined.harm, group.by="predicted.ann_level_3") +
    labs(title="L3_Azi annotation") +
    scale_color_manual("L3_Azi", values=L3_azi.cols)
  g1 = p1 + p2 + patchwork::plot_layout(ncol=2)
  ggsave(filename=paste0(name,"_MARKER_HARM_UMAP.pdf"),plot=g1,width=18)
  p3 <- FeaturePlot(S_Combined.SCT, features=name, order=TRUE) + labs(title=paste0(name," Markers"))
  p4 <- DimPlot(S_Combined.SCT, group.by="predicted.ann_level_3") +
    labs(title="L3_Azi annotation") +
    scale_color_manual("L3_Azi", values=L3_azi.cols)
  g2 = p3 + p4 + patchwork::plot_layout(ncol=2)
  ggsave(filename=paste0(name,"_MARKER_SCT_UMAP.pdf"),plot=g2,width=18)
}