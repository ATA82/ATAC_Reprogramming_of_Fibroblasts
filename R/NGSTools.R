peaks2saf <- function(peaks){
  annotated.merged.peaks.saf           <- peaks[, c("SYMBOL", "seqnames", "start", "end")]
  annotated.merged.peaks.saf$Strand    <- "-"
  colnames(annotated.merged.peaks.saf) <- c("GeneID","Chr","Start","End","Strand")
  return(annotated.merged.peaks.saf)
}

dist <- function(peak.start, peak.end, tss){
  start = peak.start # min(peak.start, peak.end)
  end = peak.end # max(peak.start, peak.end)
  if(tss > end){
    dist <- tss - end
  } else {
    dist <- start-end
  }
  return(dist)
}

overlap <- function(start_1, end_1, start_2, end_2){
  return(max(start_1, start_2) <= min(end_1, end_2))
}


counts <- function(peaks.list, bam_dir, bam_labels, modus, sample_labels) {
  ## Reduce peak file to saf format.
  #annotated.merged.peaks.saf           <- peaks.list[, c("SYMBOL", "seqnames", "start", "end")]
  #annotated.merged.peaks.saf$Strand    <- "-"
  #colnames(annotated.merged.peaks.saf) <- c("GeneID","Chr","Start","End","Strand")
  
  annotated.merged.peaks.saf <- peaks.list
  
  readsDir <- bam_dir
  labels   <- bam_labels
  bams     <- paste(readsDir, labels, ".sorted.",modus,"bam",sep="")

  fc_PE    <- Rsubread::featureCounts(files=bams, annot.ext = annotated.merged.peaks.saf, isPairedEnd=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE), nthreads=2)  
  counts   <- as.data.frame(fc_PE$counts)
  colnames(counts)  <- sample_labels
  counts$Gene <- rownames(counts)
  return(counts)
}

pca.plot <- function(deseq.object, title, scale, label) {
  ntop = 500
  Pvars <- rowVars(assay(deseq.object))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
  PCA <- prcomp(t(assay(deseq.object)[select, ]), scale = scale)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  #dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4],  Group  = colData(deseq.object)$Group)
  dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],  Group  = colData(deseq.object)$Group)
  if(label=="label"){
    group.plot <- # ggplot(dataGG,aes(x=PC1,y=PC2,color=Group )) +
      
       (qplot(PC1, PC2, data = dataGG, color =  Group) + 
                labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)), y = paste0("PC2, VarExp:", round(percentVar[2],4))) +
      scale_colour_brewer(type="qual", palette=2)) + 
      ggrepel::geom_label_repel(aes(label=rownames(dataGG)))  +
      ggtitle(paste(title,": PC1 vs PC2, top 500 variable genes",sep=""))+ 
      theme(plot.title=element_text(size=7, colour="red")) + theme_minimal()
  } else {
    group.plot <- (qplot(PC1, PC2, data = dataGG, color =  Group, main = paste(title,": PC1 vs PC2, top 500 variable genes",sep=""), size = I(3)) + 
                     labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)), y = paste0("PC2, VarExp:", round(percentVar[2],4))) +
                     scale_colour_brewer(type="qual", palette=2)) + ggrepel::geom_text_repel(aes(label=rownames(dataGG)))  + theme(plot.title=element_text(colour="red")) + theme_minimal()
  }
  return(group.plot)
}

generate.heatmaps <- function(raw.deseq, vsd.deseq, rlog.deseq, top=100, selected.genes, cluster_rows=T, cluster_cols=F, main="", order=NULL, gene.mask=NULL) {
  # data
  norm_strategies     <- c("raw", "vst", "rlog")
  types               <- c("top100", "top100DEG", "stem")
  filtering_strategy  <- c("rmdup", "keepdup")
  # data description
  norm_strategies_text    <- c("No Norm. (Scaled)", "No Norm. (Scaled)", "No Norm. (Scaled)", "Variance stabilization", "Variance stabilization","Variance stabilization", "Regularized logarithm",  "Regularized logarithm", "Regularized logarithm")
  types_text              <- c("Top 100 variable genes", "Top 100 DEG's", "Stem cell genes")
  filtering_strategy_text <- c("Removing duplicates", "Keeping duplicates")
  
  heatmaps <- vector("list", 9)
  
  i <- 1
  for(norm_strategy in norm_strategies){
    
    if(norm_strategy=="raw") { 
      data <- assay(raw.deseq)
      ntd <- raw.deseq
    }
    if(norm_strategy=="vst") { 
      data <- assay(vsd.deseq)
      ntd <- vsd.deseq
    }
    if(norm_strategy=="rlog"){
      data <- assay(rlog.deseq)
      ntd <- rlog.deseq
    }
    
    # presented genes
    
    
    counts.df <- as.data.frame(data)
    # Scaling results.
    center=T
    scale=T
    if(norm_strategy != "raw") scale=F
    scaled.data <- as.data.frame(t(as.data.frame(scale(t(counts.df), center = center, scale = scale))))
    
    for(type in types){
      
      
      if(type=="top100"){
        data <- data[!startsWith(rownames(data), "MT"),]
        data <- data[!startsWith(rownames(data), "LOC"),]
        data <- data[!startsWith(rownames(data), "MIR"),]
        data <- data[!startsWith(rownames(data), "LINC"),]
        ntop=top
        Pvars <- rowVars(as.matrix(data))
        select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
        scaled.data.1 <-  scaled.data[select, ]
      }
      
       else if(type=="top100DEG"){
        de.data <- as.data.frame(results(raw.deseq))
        de.data$Gene <- rownames(de.data)
        if(!is.null(gene.mask)){
          de.data <- de.data[de.data$Gene %in% gene.mask,]
        }
        de.data <- de.data[!startsWith(de.data$Gene, "MT"),]
        de.data <- de.data[!startsWith(rownames(de.data), "LOC"),]
        de.data <- de.data[!startsWith(de.data$Gene, "MIR"),]
        de.data <- de.data[!startsWith(de.data$Gene, "LINC"),]
        
        de.data <- na.omit(de.data[abs(de.data$log2FoldChange) >=0.56 & de.data$padj <=0.05,])
        de.data$absFC <- abs(de.data$log2FoldChange)
        de.data.pos <- de.data %>% top_n(50, log2FoldChange)
        de.data.neg <- de.data %>% top_n(-50, log2FoldChange)
        
        #  de.data <- de.data[order(abs(de.data$absFC)),]
        # de.data <- de.data[1:100,]
        
        select <- unique(c(de.data.pos$Gene,de.data.neg$Gene))
        scaled.data.1 <- scaled.data[rownames(scaled.data) %in% select, ]
        print(paste(type, ": Number of selected genes", length(select), sep=""))
      }
      
      else if(type=="stem"){
        select <- selected.genes
        scaled.data.1 <- scaled.data[rownames(scaled.data) %in% select, ]
        if(!is.null(order)){
          scaled.data.1$Gene <- rownames(scaled.data.1)
          scaled.data.1$Gene  <- factor(scaled.data.1$Gene , levels = order)
          scaled.data.1 <- scaled.data.1[order(scaled.data.1$Gene),]
        }
      }
      
      colnames(scaled.data) <- colnames(counts.df)
      scaled.data <- na.omit(scaled.data)
      library(pheatmap)
      df <- as.data.frame(colData(ntd)[,c("Group")])
      
      rownames(df) <- colnames(scaled.data)
      df$Sample <- rownames(df)
      colnames(df) <- c("Group","Sample")
      show=T
      font_size_row=12
      if(type=="top100") font_size_row=8
      if(type=="top100" & ntop>500) show=F
      
      library(ggplotify)
      heatmaps[[i]] <- as.grob(pheatmap(scaled.data.1[, c(1:(ncol(scaled.data.1)-1))], 
                                        main = ifelse(main=="", main, paste("Normalization: ", norm_strategies_text[[i]], sep="")), 
                                        border_color = "#FFFFFF",
                                        gaps_col = c(3),
                                        cluster_rows=F, 
                                        show_rownames=show,
                                        cluster_cols=F, 
                                        color= colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                 "RdYlBu")))(100),
                                        fontsize_row = font_size_row, 
                                        annotation_col=df, 
                                        scale = "row"))
      if(nrow(scaled.data[rownames(scaled.data) %in% selected.genes,]) > 2)
        pheatmap(scaled.data[rownames(scaled.data) %in% selected.genes,], cluster_rows=T, show_rownames=T,cluster_cols=F, fontsize_row = font_size_row, annotation_col=df)
      i <- i + 1
    }
  }
  
  return(heatmaps)
}