library(targets)
source("conf/params.R")
source("R/functions.R")
source("R/visualization.R")
options(tidyverse.quiet = TRUE)
tar_option_set(
  packages = 
    c("knitr", "tidyverse", "Hmisc", "gdata", 
      # Data processing. tidyverse core: "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats".
      "RColorBrewer", "ggplot2", "ggpubr", "cowplot", "pheatmap", "scales", # Visualization.
      "Rsamtools", "ATACseqQC", "ChIPseeker", "clusterProfiler", # Analysis tools.
      "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "DESeq2", "gprofiler2",
      "biglm", "readr", "targets", "AnnotationHub", "ggrepel", "dplyr")) # Databases.

list(
  
  #####################
  ## Peak annotation ##
  #####################
  tar_target(
    peaksFile,
    peaks.file.path,
    format = "file"
  ),
  tar_target(
    peaksDF,
    read.table(peaksFile, sep = "\t", header=T) %>% as.data.frame()
  ),
  tar_target(
    ensDB_AHO,
    AnnotationHub() %>% query(db.query.pattern) 
  ),
  tar_target(
    ensDB,
    ensDB_AHO %>% `[[`(`$`(ensDB_AHO,"ah_id"))
  ),
  tar_target(
    peakGRanges,
    peaksDF %>% makeGRangesFromDataFrame(keep.extra.columns = T) 
  ),
  tar_target(
    annotatedPeaks,
      ChIPseeker::annotatePeak(peak = peakGRanges, tssRegion = c(-3000, 3000),
                   TxDb = ensDB,
                   level = "transcript",
                   assignGenomicAnnotation = TRUE,
                   genomicAnnotationPriority = 
                     c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                   annoDb = "org.Hs.eg.db",
                   addFlankGeneInfo = FALSE,
                   flankDistance = 5000,
                   sameStrand = FALSE,
                   ignoreOverlap = FALSE,
                   ignoreUpstream = FALSE,
                   ignoreDownstream = FALSE,
                   overlap = "all",
                   verbose = TRUE)
  ),
  tar_target(
    processedAnnotatedPeaks,
    annotatedPeaks %>% 
      `@`(anno) %>% 
      as.data.frame() %>% 
      dplyr::select(c((ncol(.)-2):ncol(.),1:(ncol(.)-3))) %>% 
      filter(!is.na(.$geneId)) %>% 
      dplyr::rename(Symbol=SYMBOL, Geneid=interval_id)
    
  ),
  
  ####################################
  ## Preprocessing Counts in Peaks. ##
  ####################################
  tar_target(
    peaksCountsFile,
    peaks.counts.file.path,
    format = "file"
  ),
  tar_target(
    peaksCountsDF,
    read.table(peaksCountsFile, sep="\t", header=T) %>% 
      dplyr::filter(Geneid %in% processedAnnotatedPeaks$Geneid) %>%
      merge(processedAnnotatedPeaks, by="Geneid")
  ),
  tar_target(
    peaksCountsDF.reduced,
    peaksCountsDF %>% 
      dplyr::select(c(1,7:12)) %>% 
      tibble::column_to_rownames(var="Geneid") %>% 
      rename_with(~str_remove(., '.mLb.clN.bam'))
  ),
  # Metadata
  tar_target(
    colData,
    data.frame(
      samples=colnames(peaksCountsDF.reduced),
      condition=as.factor(gsub(pattern="_R[1-9]+", 
                               replacement = "", 
                               x = colnames(peaksCountsDF.reduced)))) %>%
      tibble::column_to_rownames(var="samples") %>%
      dplyr::arrange(condition) %>%
      add_column(Original.Name=c("SCC21", "SCC22", "SCC9", "wt", "scr_ori", "scr_1")) %>%
      add_column(ID=c("156437", "156439", "156435", "156429", "156433", "156431"))
  ), 
  
  #########################################
  ## Differential accessibility analysis ##
  #########################################
  # DESeq workflow
  tar_target(
    DDS,
    DAASeq(data=peaksCountsDF.reduced[, rownames(colData)], 
           colData = colData, 
           design = ~condition, 
           reference = "Lyn.WT", 
           method = "DESeq2")
  ),
  # Signal normalization using variance stabilization transformation.
  tar_target(
    VSD, 
    vst(DDS)
  ),
  # Processing of results and adding annotation.
  tar_target(
    results,
    results(DDS) %>% 
      as.data.frame() %>% 
      tibble::add_column(condition=if_else(.$log2FoldChange >= 0, "Lyn.KO", "Lyn.WT")) %>%
      na.omit() %>%
      tibble::add_column(Geneid=rownames(.)) %>%
      merge(processedAnnotatedPeaks, by="Geneid") %>% 
      dplyr::select(c(8, 10, 1, 3, 7, 74, 80, 81, 82, 83)) 
    
  ),
  # Selecting significant results based on predefined tresholds.
  tar_target(
    sign.results,
    results %>% filter(abs(log2FoldChange) >= 0.56 & padj <= 0.05) %>%
    tibble::add_column(reduced.annotation = gsub(pattern=" (\\(.*\\))", replacement="", x=.$annotation)) %>%
    arrange(condition)
  ),
  
  ####################################
  ## Functional enrichment analysis ##
  ####################################
  
  tar_target(
    fe.results.up,
    getFEResults(
      de_results = sign.results %>% 
                    filter(log2FoldChange > 0) %>% 
                      na.omit() %>% filter(condition=="Lyn.KO"), 
      direction = "UP", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  tar_target(
    fe.results.dn,
    getFEResults(
      de_results = sign.results %>% 
                    filter(log2FoldChange < 0) %>% 
                      na.omit() %>% filter(condition=="Lyn.WT"), 
      direction = "DOWN", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  tar_target(
    customDB.c2,
    upload_GMT_file(gmtfile = "data/c2.all.v7.5.1.symbols.gmt")
    ## Please be aware if the target up to date as the returned ID could
    ## be obsolete after a time. If you like to run this target always, please
    ## add the following parameter to the target.
    #, cue = tar_cue(mode="always")
  ),
  tar_target(
    c2.fe.results.up,
    getFEResults(
      de_results = sign.results %>% 
        filter(log2FoldChange > 0) %>% 
        na.omit() %>% filter(condition=="Lyn.KO"), 
      DB=customDB.c2, 
      direction = "UP", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  tar_target(
    c2.fe.results.dn,
    getFEResults(
      de_results = sign.results %>% 
        filter(log2FoldChange < 0) %>% 
        na.omit() %>% filter(condition=="Lyn.WT"), 
      DB=customDB.c2, 
      direction = "DOWN", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  tar_target(
    customDB.h,
    upload_GMT_file(gmtfile = "data/h.all.v7.5.1.symbols.gmt")
    ## Please be aware if the target up to date as the returned ID could
    ## be obsolete after a time. If you like to run
    #,cue = tar_cue(mode="always")
  ),
  tar_target(
    h.fe.results.up,
    getFEResults(
      de_results = sign.results %>% 
        filter(log2FoldChange > 0) %>% 
        na.omit() %>% filter(condition=="Lyn.KO"), 
      DB=customDB.h, 
      direction = "UP", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  tar_target(
    h.fe.results.dn,
    getFEResults(
      de_results = sign.results %>% 
        filter(log2FoldChange < 0) %>% 
        na.omit() %>% filter(condition=="Lyn.WT"), 
      DB=customDB.h, 
      direction = "DOWN", 
      ordered=F, 
      organism = "hsapiens", 
      gene.column.name = "Symbol")
  ),
  
  #####################
  ## Quality Control ##
  #####################
  
  # PCA
  tar_target(
    PCA_plot,
    pca_plot(VSD)
  ),
  # Similarity heatmap
  tar_target(
    SIM_hm,
    similarity_heatmap(VSD)
  ),
  
  #####################
  ##  Volcano Plots  ##
  #####################
  tar_target(
    VPI,
    plot_Volcano(DEG.table = results, selected.genes = vp.selected.genes, 
                 colors.selected.genes = c("black"), 
                 xlab="log2(FC) - Regulation of accessibility",
                 pt.size = 2,repel.text.size = 3, point.colors = vp.point.colors, 
                 pt.border.color = "lightgrey", repel.max.overlap = 200000, 
                 repel.force = 0.25, repel.force_pull = 3, 
                 repel.min.segment.length = 0.1, repel.segment.linetype = "dotted", 
                 repel.segment.alpha = 1, repel.segment.ncp = 3, 
                 repel.segment.color = "white", ylim=c(0, 20), xlim=c(-5,5), 
                 title="Interval/Peak Level", title.color = "darkred")
  ),
  tar_target(
    VPT,
    plot_Volcano(DEG.table = results %>%
                   group_by(condition, Symbol, transcriptId) %>%
                   summarise_at(vars(log2FoldChange, padj), mean), 
                 selected.genes = vp.selected.genes, colors.selected.genes = "black", 
                 xlab="log2(FC) - Regulation of accessibility",pt.size = 2,
                 repel.text.size = 3, point.colors = vp.point.colors, 
                 pt.border.color = "lightgrey", repel.max.overlap = 200000, 
                 repel.force = 0.25, repel.force_pull = 3, repel.min.segment.length = 0.1, 
                 repel.segment.linetype = "solid", repel.segment.alpha = 1, 
                 repel.segment.ncp = 3, repel.segment.color = "white", 
                 ylim=c(0, 20), xlim=c(-5,5), title="Transcript Level", 
                 title.color = "darkred")
  ),
  tar_target(
    VPG,
    plot_Volcano(DEG.table = results %>%
                   group_by(condition, Symbol) %>%
                   summarise_at(vars(log2FoldChange, padj), mean), 
                 selected.genes = vp.selected.genes, colors.selected.genes = "black", 
                 xlab="log2(FC) - Regulation of accessibility",pt.size = 2,
                 repel.text.size = 3, point.colors = vp.point.colors, 
                 pt.border.color = "lightgrey", repel.max.overlap = 200000, 
                 repel.force = 0.25, repel.force_pull = 3, repel.min.segment.length = 0.1, 
                 repel.segment.linetype = "solid", repel.segment.alpha = 1, 
                 repel.segment.ncp = 3, repel.segment.color = "white", 
                 ylim=c(0, 20), xlim=c(-5,5), title="Transcript Level", 
                 title.color = "darkred")
    ),
  
  #################
  ##    Heatmap  ##
  #################
  tar_target(
    HM_padj,
    plot_hm(counts=peaksCountsDF.reduced, 
            daa_results = results,
            colData = colData,
            annotated.peaks.df = processedAnnotatedPeaks,
            order.by = "padj")
  ),
  tar_target(
    HM_l2fc,
    plot_hm(counts=peaksCountsDF.reduced, 
            daa_results = results,
            colData = colData,
            annotated.peaks.df = processedAnnotatedPeaks,
            order.by = "log2FoldChange")
  ),
  
  tar_target(
    TowerPlot_25,
    TowerPlot(deg_data = DEA,
              dap_data = results, 
              move.gene.name.down = 1.3,
              only.significant.peaks = F, 
              top_n=25,
              point.shape=23, ylim=c(-6,8), 
              min.point.size = .1, 
              max.point.size = 12,
              color.gene.name = "black", 
              roi=c("Promoter"),
              color.region.of.interest = "yellow", 
              stack.inc=.3,
              color.border.sign.peaks = "black", 
              stroke.high = .5, 
              stroke.std = 0.5,
              color.border.unsign.peaks = "lightgrey", 
              #title="Differential accessibility of Peaks in top 25 UP- and top 25 DOWN-regulated genes",
              theme=theme_cowplot(), 
              return="plot")
  ),
  
  
  ### Integration
  
  ## Transcriptomics
  tar_target(DEA_file, transcriptomics_dea, format = "file"),
  
  tar_target(
    DEA,
    read.table(file=DEA_file, sep="\t", header=T)  %>% 
    rename(Symbol=SYMBOL, 
           Description=GENENAME, 
           log2FoldChange=logFC, 
           pvalue=P.Value, 
           padj=qval)
  ),
  
  ## Proteomics
  tar_target(PDEA_file, proteomics_dea, format = "file"),
  
  tar_target(
    PDEA,
    read.table(PDEA_file, header=T) %>% 
      rename(Symbol=SYMBOL, log2FoldChange=logFC, pvalue=p.val, padj=p.adj)
    
  )
  
)


















































