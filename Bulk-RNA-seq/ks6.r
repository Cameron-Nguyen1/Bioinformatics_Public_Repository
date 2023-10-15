#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
option_list = list(
  make_option(c("--studydesign"), type="character", default=NULL, dest="studydesign",
              help="Study Design file name.", metavar="character"),
    make_option(c("--database"), type="character", default=NULL, dest="database",
    help="Database like EnsDb.Hsapiens.v86 or EnsDb.Mmusculus.v79", metavar="character"),
    make_option(c("--expression"), type="character", default=NULL, help="Expression that defines a comparison using sample group names. e.g. infection=Group1-Group2", metavar="character"),
   make_option(c("--ex_eh"), type="character", default=NULL,
              help="Expression that defines: Average of Group 1 samples // e.g.
              (Group1_97+Group1_98+Group1_99)/3", metavar="character"),
  make_option(c("--ex_uh"), type="character", default=NULL,
              help="Expression that defines: Average of Group 2 samples // e.g.
              (Group2_85+Group2_86+Group2_87)/3", metavar="character"),
  make_option(c("--species"), type="character", default=NULL,
              help="Genus Species from which samples are taken from. Format: \"Mus Musculus\"", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="Directs output to designated directory.", metavar="character")
)
p <- parse_args(OptionParser(option_list=option_list))
print(paste("Starting ----- ",p$expression,sep=""))
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(ensembldb))
#suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))


dir.create(p$out,showWarnings=FALSE)
kalli_flow = setRefClass("kalli_flow",fields = list(), methods = list(
  read_data = function(wdir, study_file, db)
  {
    setwd(as.character(wdir))
    targets <- read_tsv(as.character(study_file))# read in your study design
    path <- file.path(targets$Sample, "abundance.tsv") # set file paths to your mapped data
    Tx <- as_tibble(transcripts(db, columns=c("tx_id", "gene_name")))
    Tx <- dplyr::rename(Tx, target_id = tx_id)
    Tx <- dplyr::select(Tx, "target_id", "gene_name")
    Txi_gene <- tximport(path,
                         type = "kallisto",
                         tx2gene = Tx,
                         txOut = FALSE, #determines whether your data represented at transcript or gene level
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)
    return(list(Txi_gene,targets)) #Return Gene list with Targets tibble
  },
  visualize = function(read_res)
  {
    sampleLabels <- as.character(read_res[[2]]$Sample)
    myDGEList <- DGEList(read_res[[1]]$counts)

    cpm <- cpm(myDGEList)
    keepers <- rowSums(cpm>1)>=5 #user defined
    myDGEList.filtered <- myDGEList[keepers,]
    myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
    log2.cpm <- cpm(myDGEList.filtered.norm, log=TRUE)
    log2.cpm.filtered.norm.df <- as_tibble(log2.cpm, rownames = "geneID")
    colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

    log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
    colnames(log2.cpm.df) <- c("geneID", sampleLabels)
    log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                      cols = as.character(sampleLabels[1]):as.character(sampleLabels[length(sampleLabels)]), # column names to be stored as a SINGLE variable // Auto
                                      names_to = "samples", # name of that new variable (column)
                                      values_to = "expression") # name of new variable (column) storing all the values (data)
    pdf(paste(p$out,"/Log2_CPM_Chart.pdf",sep=""))
   pl0t = ggplot(log2.cpm.df.pivot) +
      aes(x=samples, y=expression, fill=samples) +
      geom_violin(trim = FALSE, show.legend = FALSE) +
      stat_summary(fun = "median",
                   geom = "point",
                   shape = 95,
                   size = 10,
                   color = "black",
                   show.legend = FALSE) +
      labs(y="log2 expression", x = "sample",
           title="Log2 Counts per Million (CPM)",
           #subtitle="unfiltered, non-normalized",
           caption=paste0("produced on ", Sys.time())) +
      theme_bw() +
      coord_flip()
plot(pl0t)
    dev.off()
    return(list(log2.cpm,myDGEList.filtered.norm,sampleLabels,log2.cpm.filtered.norm.df))

  },
  pca = function(data_from_read, data_from_vis)
  {
    Group <- factor(data_from_read[[2]]$Group)
    pca.res <- prcomp(t(data_from_vis[[1]]), scale.=F, retx=T)
    pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
    pc.per <- round(pc.var/sum(pc.var)*100, 1)
    pca.res.df <- as_tibble(pca.res$x)
    pdf(paste(p$out,"/PCA_Chart.pdf",sep=""))
    pca.plot = ggplot(pca.res.df) +
      aes(x=PC1, y=PC2, label=data_from_read[[2]]$Sample, color = Group) +
      geom_point(size=3)+#, shape = Shape) +
      #stat_ellipse() +
      xlab(paste0("PC1 (",pc.per[1],"%",")")) +
      ylab(paste0("PC2 (",pc.per[2],"%",")")) +
      labs(title="PCA plot",
           caption=paste0("produced on ", Sys.time())) +
      coord_fixed() +
      theme_bw()
    plot(pca.plot)
    dev.off()
  },
  vp = function(data_from_read, data_from_vis, labels)
  {
    Group <- data_from_read[[2]]$Group
    Group <- factor(Group)
    design <- model.matrix(~0 + Group)
    colnames(design) <- levels(Group)

    v.DEGList.filtered.norm <- voom(data_from_vis[[2]], design, plot = FALSE)
    fit <- lmFit(v.DEGList.filtered.norm, design)
    contrast.matrix <- makeContrasts(ex,
                                     levels=design)

    fits <- contrasts.fit(fit, contrast.matrix)
    ebFit <- eBayes(fits)
    myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
    myTopHits.df <- myTopHits %>%
      as_tibble(rownames = "geneID")
    results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1.5) #Pulls DEGs from VP

    num_neg = length(myTopHits.df$logFC[myTopHits.df$logFC<=(-1.5) & myTopHits.df$adj.P.Val <= .05])
    num_pos = length(myTopHits.df$logFC[myTopHits.df$logFC>=1.5 & myTopHits.df$adj.P.Val <= .05])
    vplot <- ggplot(myTopHits.df) +
      aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
      geom_point(size=2) +
      geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
      geom_vline(xintercept = 1.5, linetype="longdash", colour="#BE684D", size=1) +
      geom_vline(xintercept = -1.5, linetype="longdash", colour="#2C467A", size=1) +
      #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
      #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
      labs(title=paste("Volcano plot",substr(ex,11,nchar(ex))),
           caption=paste("Num. Genes Positive:",num_pos,"/// Num. Genes Negative:",num_neg,sep=" ")) +
      theme_bw()

    r4 = myTopHits.df
    q = r4[r4[,6] <= .05 & abs(r4[,2]) >= 1.5,]
    write.csv(r4,paste(p$out,"/",substr(ex,11,nchar(ex)),"_All_Hits.csv",sep=""))
    write.csv(q,paste(p$out,"/",substr(ex,11,nchar(ex)),"_myTopHits.csv",sep=""))
    colnames(v.DEGList.filtered.norm$E) <- labels
    write.csv(v.DEGList.filtered.norm$E,file=paste(p$out,"/",substr(ex,11,nchar(ex)),"_Enrichment_Scores.csv",sep=""))
    pdf(paste(p$out,"/",substr(ex,11,nchar(ex)),"_Volcano_Plot.pdf",sep=""))
    plot(vplot)
    dev.off()
    return(list(v.DEGList.filtered.norm,vplot, results, q))
  },
  hm = function(labels, vp_results)
  {
    colnames(vp_results[[1]]$E) <- labels
    diffGenes <- vp_results[[1]]$E[vp_results[[3]][,1] !=0,] #Rows where expression data is DEG, has to be -1 or 1 from $E
    write.csv(diffGenes,file=paste(p$out,"/",substr(ex,11,nchar(ex)),"_Counts_output.csv",sep="")) ##diffGenes = counts for analysis, lfc=1.5, p=.05
    clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
    module.assign <- cutree(clustRows, k=2)

    pheatmap(diffGenes, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename=paste(p$out,"/",substr(ex,11,nchar(ex)),"_Allreg_heatmap.pdf",sep=""), width=10, height=16, show_rownames = FALSE,
             main="Gene Enrichment Heatmap", scale="row")
    mu_n = names(module.assign[module.assign==2])
    mu_n = mu_n[!mu_n %in% ""]
    myModule_up <- diffGenes[mu_n,]

    pheatmap(myModule_up, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename=paste(p$out,"/",substr(ex,11,nchar(ex)),"_Upreg_heatmap.pdf",sep=""), width=10, height=16, show_rownames = FALSE,
             main="Upregulated Gene Heatmap", scale="row")

    md_n =  names(module.assign[module.assign==1])
    md_n = md_n[!md_n %in% ""]
    myModule_down <- diffGenes[md_n,]

    pheatmap(myModule_down, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename=paste(p$out,"/",substr(ex,11,nchar(ex)),"_Downreg_heatmap.pdf",sep=""), width=10, height=16, show_rownames = FALSE,
             main="Downregulated Gene Heatmap", scale="row")

    return(list(myModule_up, myModule_down))
  },
  gsea = function(spec, cat,m_up, m_down, mydata.df)
  {
    gost.res_up <- gost(rownames(m_up), organism = "mmusculus", correction_method = "fdr")
    gp_up = gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
    gost.res_down <- gost(rownames(m_down), organism = "mmusculus", correction_method = "fdr")
    gp_down = gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
    hs_gsea_c2 <- msigdbr(species = spec, # change depending on species your data came from
                          category = cat) %>% # choose your msigdb collection of interest
      dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature
    mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
    mydata.gsea <- mydata.df.sub$LogFC
    names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
    mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
    myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
    myGSEA.df <- as_tibble(myGSEA.res@result)
    tableu = datatable(myGSEA.df,
                       extensions = c('KeyTable', "FixedHeader","Buttons"),
                       #caption = 'Signatures enriched in leishmaniasis',
                       options = list(dom="Bfrtip",buttons=c("excel","csv"),keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
      formatRound(columns=c(3:10), digits=2)
    sea = gseaplot2(myGSEA.res,
                    geneSetID = 1:5, #Lets you choose what you plot, comes from the data table "tableu"
                    pvalue_table = FALSE) #can set this to FALSE for a cleaner plot
                    #title = myGSEA.res$Description[]) #can also turn off this title
    myGSEA.df <- myGSEA.df %>%
      mutate(phenotype = case_when(
        NES < 0 ~ "healthy",
        NES > 0 ~ "disease"))
    bubble = ggplot(myGSEA.df, aes(x=phenotype, y=ID)) +
      geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
      scale_color_gradient(low="blue", high="red") +
      ggtitle(paste(substr(ex,11,nchar(ex)),"_Enrichment_Bubble_Chart",sep="")) +
      ylab("Pathway")+
      theme_bw()
    pdf(paste(p$out,"/",substr(ex,11,nchar(ex)),"_Bubble_Enrichment.pdf",sep=""))
    plot(bubble)
    dev.off()
    return(list(gp_up,gp_down,myGSEA.df,sea,bubble,myGSEA.res))
  },
  DTU_Analysis = function(expression_dir, design_file, annotation_file, cdna_genome, res_num, qval_bool)
  {
    xrd = importIsoformExpression(expression_dir)
    xrddesign = read.csv(design_file,sep="\t")
    mySwitchList <- importRdata(
      isoformCountMatrix = xrd$counts,
      isoformRepExpression = xrd$abundance,
      designMatrix = xrddesign,
      removeNonConvensionalChr = TRUE,
      addAnnotatedORFs=TRUE,
      ignoreAfterPeriod=TRUE,
      isoformExonAnnoation = annotation_file,
      isoformNtFasta = cdna_genome,
      showProgress = FALSE)

    mySwitchList <- isoformSwitchAnalysisCombined(
      n=50,
      switchAnalyzeRlist = mySwitchList,
      genomeObject=cdna_genome,
      pathToOutput = p$out) # directory must already exist
    extractSwitchSummary(mySwitchList)

    x = extractTopSwitches(
      mySwitchList,
      filterForConsequences = TRUE, # these 'consequences' related to the annotations I reference above.
      n = res_num,
      sortByQvals = qval_bool) #change to TRUE if you want this list sorted by FDR-adusted Pval (a.k.a., q value)
    write.csv(x,file=paste0(p$out,"/Switch_top_50.csv"))

    return(mySwitchList)
  },

  plot_switch = function(switch_list, gene_name, plot_name, cond1, cond2)
  {
    x = switchPlot(
      switch_list,
      gene=gene_name,
      condition1 = cond1,
      condition2 = cond2,
      localTheme = theme_bw())
    ggsave(plot=x,filename=plot_name,width=10,height=10)
  }
))
k = kalli_flow()
print("Initialized...")
res1 = k$read_data(getwd(),p$studydesign,eval(parse(text=p$database)))
print("Read complete.")
res2 = k$visualize(res1)
res3 = k$pca(res1,res2)
print("Visualization complete.")
ex = p$expression
print("Starting VP.")
res4 = k$vp(res1,res2,res2[[3]])
print("VP - OK / HM Begin")
res5 = k$hm(res2[[3]],res4)
print("HM Done, frame -> GSEA")
mydata.df <- mutate(res2[[4]],
                    disease.AVG = eval(parse(text=p$ex_uh)),
                    healthy.AVG = eval(parse(text=p$ex_eh)),
                    LogFC = (disease.AVG - healthy.AVG)) %>%
  mutate_if(is.numeric, round, 2)
res6 = k$gsea(p$species,"H",res5[[1]],res5[[2]],mydata.df)
print("GSEA complete.")
write.csv(res6[[3]],paste(p$out,"/",substr(ex,11,nchar(ex)),"_PATHS.csv",sep=""))
pdf(paste(p$out,"/",substr(ex,11,nchar(ex)),"_Running_Enrichment_Score.pdf",sep=""))
res6[[4]]
dev.off()

iso_pre = read_tsv(p$studydesign)
colnames(iso_pre)=c("sampleID","condition")
write.table(iso_pre,file="isoform.txt",sep="\t",row.names=FALSE)
iso = read_tsv("isoform.txt")
conds = unique(iso$condition)
myswitch = k$DTU_Analysis(expression_dir = getwd(),
                          design_file = "isoform.txt",
                          annotation_file = "/path/to/annotationfile/gencode.vM32.chr_patch_hapl_scaff.annotation.gtf",
                          cdna_genome = "/path/to/annotationfile/gencode.vM32.transcripts.fa",
                          res_num = 50,
                          qval_bool = FALSE)
#Specific gene switchplots:
#k$plot_switch(switchlist from DTU analysis, gene of interest, Outputfile, Condition1, Condition2)
k$plot_switch(myswitch, 'Tab1', "Tab1_SwitchPlot.pdf", conds[1],conds[2])
print("--- DONE ---")