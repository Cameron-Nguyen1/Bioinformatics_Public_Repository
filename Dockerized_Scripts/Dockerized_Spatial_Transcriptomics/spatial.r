library(optparse)
option_list = list(
  #make_option(c("--inputdir"), type="character", default=NULL,help="Input directory, contains spaceranger output for >1 slides.", metavar="character"),
  make_option(c("--outputdir"), type="character", default=NULL,help="Output directory, will contain products of pipeline", metavar="character"),
  make_option(c("--slides"), type="character", default=NULL,help="Comma seperated input designating slide folder i.e.: A1,D2 ", metavar="character")
)
p = parse_args(OptionParser(option_list=option_list))
dep = c("Seurat","SeuratData","BPCells","ggplot2","patchwork","dplyr","SeuratDisk","stringr")
lapply(X=dep,library,character.only=TRUE)
print(as.character(p[1]))
dir.create(as.character(p[1]))

load_merge = function(slide_vector){
    for (slide in slide_vector){
        eval(parse(text=paste0(slide,"=Load10X_Spatial(\"",slide,"\",slice=\"",slide,"\")")))
        eval(parse(text=paste0(slide,"=",slide,"[,unname(which(colSums(GetAssayData(",slide,"))!=0))]")))
        eval(parse(text=paste0(slide,"$orig.ident=\"",slide,"\"")))
        eval(parse(text=paste0(slide,"=SCTransform(",slide,",assay=\"Spatial\",verbose=FALSE)")))
    }
    for (i in 1:length(slide_vector)){
        if (i == 1){
            ret = merge(get(slide_vector[1]),get(slide_vector[2]))
        }
        else if (i == length(slide_vector)){
            break
        }
        else{
            ret = merge(ret,get(slide_vector[i+1]))
        }
    }
    return(list(ret,get(slide_vector[1]),get(slide_vector[2])))
}
Spatial_Prediction_Plots = function(merged_obj,name_vector,ncols,name,width,height){
    mFP = SpatialFeaturePlot(merged_obj, features = name_vector, pt.size.factor = 2.5, ncol = ncols, crop = TRUE)
    for (i in 1:length(mFP)){
        mFP[[i]]$theme$legend.key.width = unit(1.5,"cm")
    }
    ggsave(plot=mFP,filename=paste0("/usr/local/work/",p[1],"/",name),width=width,height=height)
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

slide_vector = unlist(strsplit(as.character(p[2]),","))

x = load_merge(slide_vector)
#saveRDS(x,"FirstMerge.rds")
#x = readRDS("FirstMerge.rds")
Merged = x[1][[1]]

DefaultAssay(Merged) = "SCT"
VariableFeatures(Merged) = c(VariableFeatures(x[2][[1]]), VariableFeatures(x[3][[1]]))
Merged = RunPCA(Merged, verbose = TRUE)
m_bdim = bdim(Merged)
Merged = FindNeighbors(Merged, dims = 1:m_bdim)
Merged = FindClusters(Merged, verbose = TRUE)
Merged = RunUMAP(Merged, dims = 1:m_bdim)

umap_c = DimPlot(Merged, reduction = "umap", group.by = c("ident", "orig.ident"))
ggsave(plot=umap_c,filename=paste0("/usr/local/work/",p[1],"/Umap_Comparison.pdf"),width=20,height=10)
spatial_cluster_c = SpatialDimPlot(Merged,label=TRUE,repel=TRUE,label.size=4)
ggsave(plot=spatial_cluster_c,filename=paste0("/usr/local/work/",p[1],"/Spatial_Cluster_Comparison.pdf"),width=20,height=12)
Fplot = SpatialFeaturePlot(Merged, features = c("Sftpc","Cxcr2"))
ggsave(plot=Fplot,filename=paste0("/usr/local/work/",p[1],"/FeaturePlot.pdf"),width=20,height=12)

#Convert("LungMAP_MouseLung_CellRef.v1.1.h5ad", "LungMAP_Mouse.h5seurat")
MouseRef = LoadH5Seurat("LungMAP_Mouse.h5seurat")
MouseRef = SCTransform(MouseRef, ncells = 3000, verbose = TRUE)
MouseRef = RunPCA(MouseRef,verbose = FALSE)
MouseRef = RunUMAP(MouseRef, dims = 1:bdim(MouseRef))

anchors <- FindTransferAnchors(reference = MouseRef, query = Merged, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = MouseRef$celltype_level2, prediction.assay = TRUE,
    weight.reduction = Merged[["pca"]], dims = 1:m_bdim)
Merged[["predictions"]] <- predictions.assay
DefaultAssay(Merged) <- "predictions"
#saveRDS(Merged,"Merged_annot.rds")
#Merged = readRDS("Merged_annot.rds")

EndoC = c("AEC","CAP1/EPC","CAP2","LEC","VEC")
EpiC = c("AT1","AT1/AT2","AT2","Basal","Ciliated","Deuterosomal","PNEC","Secretory","Sox9_Epi")
MesenC = c("AF1","AF2","ASMC","Chondrocyte","Mesothelial","Pericyte","PMP","SCMF","VSMC")
ImmuneC = c("AM","B","Basophils","CD4 T","CD8 T","cDC1","cDC2","ILC","IM","iMON","maDC","MAST","Neutrophil","NK","Megakaryocyte/Platelet","Treg")

Spatial_Prediction_Plots(Merged,EndoC,2,"EndothelialCells_pred.pdf",20,25)
Spatial_Prediction_Plots(Merged,EpiC,6,"EpithelialCells_pred.pdf",35,25)
Spatial_Prediction_Plots(Merged,MesenC,6,"MesenchymalCells_pred.pdf",35,25)
Spatial_Prediction_Plots(Merged,ImmuneC,4,"ImmuneCells_pred.pdf",40,35)