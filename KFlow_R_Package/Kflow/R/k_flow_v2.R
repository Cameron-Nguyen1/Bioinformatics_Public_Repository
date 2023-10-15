library(methods) #----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
#library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
library(EnsDb.Mmusculus.v79)
library(edgeR)
library(matrixStats)
library(cowplot)
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables
library(limma)
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tool s for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(pheatmap)
library(enrichplot) # great for making the standard GSEA enrichment plots #Imports
library(IsoformSwitchAnalyzeR)
#library(BSgenome.Hsapiens.UCSC.hg38)
#-----
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
    sampleLabels <- read_res[[2]]$Sample
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
                                      cols = as.character(colnames(log2.cpm.df)[1]):as.character(colnames(log2.cpm.df)[length(log2.cpm.df)]), # column names to be stored as a SINGLE variable // Auto
                                      names_to = "samples", # name of that new variable (column)
                                      values_to = "expression") # name of new variable (column) storing all the values (data)
    png(filename="something.png",width=1200,height=800)
    pl0t <- ggplot(log2.cpm.df.pivot) +
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
    pca.plot <- ggplot(pca.res.df) +
      aes(x=PC1, y=PC2, label=data_from_read[[2]]$Sample, color = Group) +
      geom_point(size=3)+#, shape = Shape) +
      #stat_ellipse() +
      xlab(paste0("PC1 (",pc.per[1],"%",")")) +
      ylab(paste0("PC2 (",pc.per[2],"%",")")) +
      labs(title="PCA plot",
           caption=paste0("produced on ", Sys.time())) +
      coord_fixed() +
      theme_bw()
    png("PCA_Plot.png")
    ggplotly(pca.plot)
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

    vplot <- ggplot(myTopHits.df) +
      aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
      geom_point(size=2) +
      geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
      geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
      geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
      #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
      #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
      labs(title=paste("Volcano plot",substr(ex,11,nchar(ex))),
           caption=paste0("produced on ", Sys.time())) +
      theme_bw()

    r4 = myTopHits.df
    q = r4[r4[,6] <= .05 & abs(r4[,2]) >= 1.5,]
    write.csv(r4,paste("All_Hits_",substr(ex,11,nchar(ex)),".csv",sep=""))
    write.csv(q,paste("myTopHits_",substr(ex,11,nchar(ex)),".csv",sep=""))
    colnames(v.DEGList.filtered.norm$E) <- labels
    write.csv(v.DEGList.filtered.norm$E,file=paste("Enrichment_Scores_",substr(ex,11,nchar(ex)),".csv",sep=""))
    png(paste("Volcano Plot ",substr(ex,11,nchar(ex)),".png",sep=""),500,500)
    ggplotly(vplot)
    dev.off()
    return(list(v.DEGList.filtered.norm,vplot, results, q))
  },
  hm = function(labels, vp_results)
  {
    colnames(vp_results[[1]]$E) <- labels
    diffGenes <- vp_results[[1]]$E[vp_results[[3]][,1] !=0,] #Rows where expression data is DEG, has to be -1 or 1 from $E
    clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
    module.assign <- cutree(clustRows, k=2)

    pheatmap(diffGenes, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename="Allreg_heatmap.pdf", width=10, height=15, show_rownames = FALSE,
             main="Gene Enrichment Heatmap", scale="row")
    mu_n = names(module.assign[module.assign==2])
    mu_n = mu_n[!mu_n %in% ""]
    myModule_up <- diffGenes[mu_n,]

    pheatmap(myModule_up, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename="Upreg_heatmap.pdf", width=10, height=15, show_rownames = FALSE,
             main="Upregulated Gene Heatmap", scale="row")

    md_n =  names(module.assign[module.assign==1])
    md_n = md_n[!md_n %in% ""]
    myModule_down <- diffGenes[md_n,]

    pheatmap(myModule_down, color = colorRampPalette(c("#253494", "#2C7FB8", "#dcfaf3", "#fcfc92", "#fcfc05"))(100),
             cluster_cols=TRUE, cluster_rows=FALSE, clustering_distance_cols = "maximum",
             filename="Downreg_heatmap.pdf", width=10, height=15, show_rownames = FALSE,
             main="Downregulated Gene Heatmap", scale="row")

    write.csv(diffGenes,file="Counts_output.csv") ##diffGenes = counts for analysis, lfc=1.5, p=.05
    return(list(myModule_up, myModule_down))
  },
  gsea = function(spec, cat,m_up, m_down, mydata.df)
  {
    gost.res_up <- gost(rownames(m_up), organism = "hsapiens", correction_method = "fdr")
    gp_up = gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
    gost.res_down <- gost(rownames(m_down), organism = "hsapiens", correction_method = "fdr")
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
        NES > 0 ~ "disease",
        NES < 0 ~ "healthy"))
    bubble = ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) +
      geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
      scale_color_gradient(low="blue", high="red") +
      theme_bw()
    return(list(gp_up,gp_down,tableu,sea,bubble,myGSEA.res))
  },
  DTU_Analysis = function(expression_dir, design_file, annotation_file, cdna_genome, res_num, qval_bool)
  {
    xrd = importIsoformExpression(expression_dir)
    xrddesign = read.csv(design_file,sep="\t")
    mySwitchList <- importRdata(
      isoformCountMatrix   = xrd$counts,
      isoformRepExpression = xrd$abundance,
      designMatrix         = xrddesign,
      removeNonConvensionalChr = TRUE,
      addAnnotatedORFs=TRUE,
      ignoreAfterPeriod=TRUE,
      isoformExonAnnoation = annotation_file,
      isoformNtFasta       = cdna_genome,
      showProgress = TRUE)

    mySwitchList <- isoformSwitchAnalysisCombined(
      switchAnalyzeRlist   = mySwitchList,
      pathToOutput = getwd()) # directory must already exist
    extractSwitchSummary(mySwitchList)

    x = extractTopSwitches(
      mySwitchList,
      filterForConsequences = TRUE, # these 'consequences' related to the annotations I reference above.
      n = res_num,
      sortByQvals = qval_bool) #change to TRUE if you want this list sorted by FDR-adusted Pval (a.k.a., q value)
    write.csv(x,file="Switch_top_50.csv")

    return(mySwitchList)
  },

  plot_switch = function(switch_list, gene_name, plot_name, cond1, cond2)
  {
    png(plot_name,width=900,height=800)
    switchPlot(
      switch_list,
      gene=gene_name,
      condition1 = cond1,
      condition2 = cond2,
      localTheme = theme_bw())
    dev.off()
  }
))
