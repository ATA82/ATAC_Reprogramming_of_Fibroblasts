source("R/correlation_analysis.R")
source("R/footprinting.R")

### IO ###
read_peaks <- function(){
  
}

select_peaks_by_region <- function(){
  
}

### Quantification ###

count_reads <- function(){
  
}

### Downstream Analysis ###

create_design <- function(){
  
}

# differential accessibility analysis
DAASeq <- function(data, colData, design, reference=NULL, method="DESeq2"){
  if(method=="DESeq2"){
    DDS <- DESeqDataSetFromMatrix( countData = data[, rownames(colData)],
                                   colData = colData, 
                                   design = ~condition)
    if(!is.null(reference)){
      DDS$condition <- relevel(DDS$condition, ref = reference)
    }
    DDS <- DESeq(DDS)
    return(DDS)
  } else {
    stop("Method is not supported!")
  }
}

top_genes <- function(){
  
}

# function enrichment analysis
FESeq <- function(){
  
}

top_FETerms <- function(){
  
}

top_FEGenes <- function(){
  
}

### Visualization ###

# QC
pca_plot <- function(DDS, ntop=400){
  
  ### Principal component analysis.
  Pvars <- rowVars(assay(DDS))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
  PCA <- prcomp(t(assay(DDS)[select, ]), scale = T)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                      PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                      condition = colData(DDS)$condition,
                      Original.Name = colData(DDS)$Original.Name,
                      sample = rownames(colData(DDS)))
  
  PC_ids <- combn(4, 2, simplify = F)
  
  pca_plot <- ggplot(data =dataGG, aes(x = PC1, y = PC2, color =  condition, label=paste(sample," (",Original.Name,")",sep=""))) + geom_point(size=2)+geom_text_repel() + 
    labs(x = paste0(paste("PC",PC_ids[[1]][1],", VarExp:"), round(percentVar[PC_ids[[1]][1]],4)), 
         y = paste0(paste("PC",PC_ids[[1]][2],", VarExp:"), round(percentVar[PC_ids[[1]][2]],4))) + 
    scale_colour_brewer(type="qual", palette=2) + theme_pubr()
  
  return(pca_plot)
  
}

similarity_heatmap <- function(DDS){
  sampleDists <- stats::dist(t(assay(DDS)), method = "euclidean")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(DDS$condition, DDS$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( (brewer.pal(11, "RdYlBu")) )(255)
  sim_hm <- pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
                               labels_row = colnames(DDS), labels_col = colData(DDS)$Original.Name, col=colors)
  sim_hm <- ggplotify::as.ggplot(sim_hm)  + 
    theme(plot.title = element_text(hjust = 0.5, face="bold")) 
  return(sim_hm)
}

ma_plot <- function(){
  
}

# Gene level plots
volcano_plot <- function(){
  
}

topGenes_heatmap <- function(){
  
}

selection_heatmap <- function(){
  
}

# Functional level plots
bar_plot <- function(){
  
}

topFEGenes_heatmap <- function(){
  
}


compare <- function(a, b, c=">"){
  if(c == ">"){
    return(a>b)
  } else if (c == ">="){
    return(a>=b) 
  } else if (c == "<"){
    return(a < b)
  } else if (c == "<="){
    return(a <= b)
  } else if (c == "=="){
    return(a == b)
  } else {
    stop("comparison not supported")
  }
}

DoHeatmap <- function(data, meta.data, selection,
                      do.selection=F, do.selection.top=15, do.selection.group.by="Condition", 
                      do.selection.order.by="logFC", do.selection.order.direction ="descending", do.selection.feature.name="Symbol",
                      do.selection.significance.variables = c("logFC", "adj.P.Val"), do.selection.significance.cutoffs = c(.56, .05),
                      do.selection.cutoff.relation = c(">=", "<="), remove.batch.variable = NULL,
                      do.selection.cutoff.l2fc = .56, do.selection.cutoff.padj = .05, 
                      do.selection.regulation = "up", do.selection.regulation.variable = "logFC",
                      selected.meta.data.columns, sample_name_column, gaps_col=NULL, heatmap.plot.package = "pheatmap",
                      cell.colors=colorRampPalette(c("darkblue", "white", "red2"))(50), fontsize_row = 10, fontsize_col=10, 
                      cluster_rows=F, cluster_cols=F, scale="row", show_rowname=T, show_colnames=T, 
                      annotation_names_col=T, space.between.title.and.heatmap = c(0,100), return.data = F,
                      title = "", subtitle = "", plot.title.face = "bold", plot.title.color = "darkred",
                      move.title.vertically = 0, move.title.horizontally = .2, plot.title.size = 14,
                      plot.subtitle.face = "bold", plot.subtitle.color = "#000080", move.subtitle.vertically = -1,
                      move.subtitle.horizontally = .3,plot.subtitle.size = 11, border_color="lightgrey",
                      col.annotation.colors=NA){
  ### Order columns based on order in meta.data
  data <- data[, meta.data[,sample_name_column]]
  ### Store original selection
  selected.DE.table <- selection
  if(!is.null(remove.batch.variable)){
    ### Remove batch effect (The step before is important, otherwise the order of the batch information would be wrong)
    data <- limma::removeBatchEffect(data, batch = meta.data[,remove.batch.variable])
  }
  ### Copy data to heatmap matrix.
  hm.matrix <- as.data.frame(data)
  order_indicator <- ifelse(do.selection.order.direction=="descending", -1, 1)
  
  ### Define significance values of significance variables.
  if(do.selection){
    
    cutoff.filter <- T
    for(v in 1:length(do.selection.significance.variables)){
      cutoff.filter <- cutoff.filter & 
        (compare(selection[,do.selection.significance.variables[v]], 
                 do.selection.significance.cutoffs[v],
                 do.selection.cutoff.relation[v]))
    }
    selection <- selection[cutoff.filter, ]
    
    ### In case that genes are not from multipe groups of samples (e.g. condition), then column order is not important.
    if(do.selection.group.by == "none"){
      
      ### Create a dummy group for all genes.
      selection[, do.selection.group.by] <- 1
      ### Select top genes based on order-variable.
      regulation.filter <- ifelse(do.selection.regulation == "up", selection[,do.selection.regulation.variable] > 0, 
                                  ifelse(do.selection.regulation == "down",selection[,do.selection.regulation.variable] < 0, 
                                         selection[,do.selection.regulation.variable] != 0))
      selection <- selection[regulation.filter,] %>% 
        group_by(eval(parse(text = do.selection.group.by))) %>% 
        slice_max(order_by = -order_indicator*eval(parse(text = do.selection.order.by)), n = do.selection.top)
      
      selected.DE.table <- selection
      
      ### Reduce heatmap matrix to selected genes.
      hm.matrix <- hm.matrix[rownames(hm.matrix) %in% unlist(selection[, do.selection.feature.name]),]
    } else {
      
      ### Order table by group in custom group order and then by the desired factor (L2FC or Padj)
      selection <- selection[order(selection[,do.selection.group.by], order_indicator*selection[, do.selection.order.by]),]
      selection <- selection %>% 
        group_by(eval(parse(text = do.selection.group.by))) %>% 
        slice_max(order_by = -order_indicator*eval(parse(text = do.selection.order.by)), n = do.selection.top)
      selected.DE.table <- selection
      
      # define custom order based on meta.data.
      meta.data[, do.selection.group.by] <- as.character(meta.data[, do.selection.group.by])
      order <- data.frame(unique(unlist(meta.data[, do.selection.group.by])))
      order$order_id <- 1:nrow(order)
      colnames(order) <- c(do.selection.group.by, "order_id")
      selection <- merge(selection, order, by=do.selection.group.by)
      selection <- selection[order(selection[,"order_id"]),]
      
      groups <- unique(unlist(selection[, do.selection.group.by]))
      matrix.chunks <- vector("list", length(groups))
      for(g in 1:length(groups)){
        selected.group <- selection[selection[,do.selection.group.by]==groups[g], ]
        matrix.chunks[[g]] <-  hm.matrix[rownames(hm.matrix) %in% unlist(selected.group[, do.selection.feature.name]),]
      }
      
      hm.matrix <- do.call("rbind", matrix.chunks)
      
      
    }
    selection <- unlist(selection[,do.selection.feature.name])
  } else {
    
    hm.matrix <- hm.matrix[rownames(hm.matrix) %in% selection,]
  }
  ### Generate annotation of heatmap columns bei selecting the desired column annotation
  annotation_col = data.frame(a1=factor(meta.data[,selected.meta.data.columns[1]]))
  if(length(selected.meta.data.columns) >= 2){
    for(a in 2:length(selected.meta.data.columns)){
      tmp <- selected.meta.data.columns[a]
      annotation_col[tmp] <- meta.data[,selected.meta.data.columns[a]]
    }
  }
  colnames(annotation_col) <- selected.meta.data.columns
  rownames(annotation_col) <- rownames(meta.data)
  
  ### Order columns based on row order in meta.data table (any desired order should prepared before passing the meta.data
  ### parameter)
  hm.matrix <- hm.matrix[, meta.data[,sample_name_column]]
  
  ## change row order based on order of selected genes
  hm.matrix <- as.data.frame(hm.matrix)
  hm.matrix$genes <- rownames(hm.matrix)
  hm.matrix <- (na.omit(hm.matrix[match(selection, hm.matrix$genes),]))
  hm.matrix$genes <- NULL
  hm.matrix <- as.matrix(hm.matrix)
  main=""
  ### Generate raw heatmap.
  
  #-gaps-# automatically compute gaps if gaps between columns is not defined by user.
  if(is.null(gaps_col)){
    counted <-  ((meta.data %>% add_count((eval(parse(text = do.selection.group.by))), sort=F)))
    gaps <- unique(counted[, c(do.selection.group.by, "n")])$n
    for(gap in 1:length(gaps)){
      gaps[gap] <- gaps[gap] + sum(gaps[gap-1])
    }
    gaps_col <- gaps
  }
  
  if(heatmap.plot.package == "pheatmap"){
    ann_colors <- col.annotation.colors
    if(!is.na(col.annotation.colors)){
      m <- length(selected.meta.data.columns)
      if(m <= 8 && col.annotation.colors=="viridis"){
        library("viridis")
        types <- c("B", "D", "F", "H", "A", "C", "E", "G")
        ann_colors <- vector("list", m)
        for(s in 1:m){
          ann_colors[[s]] <- viridis_pal(option = types[s])(length(unique(unlist(meta.data[, selected.meta.data.columns[s]]) )))
          ann_colors[[s]] <- rev(ann_colors[[s]])
          names(ann_colors[[s]]) <- unlist(unique(meta.data[, selected.meta.data.columns[s]]))
        }
        names(ann_colors) <- selected.meta.data.columns
      }
    }
    
    hm <- pheatmap::pheatmap(mat = hm.matrix,
                             # color of cell borders
                             border_color = border_color,
                             main = main,  
                             ### Scale configuration
                             color = cell.colors, scale = scale, 
                             ### Hierarchical clustering
                             cluster_rows=cluster_rows,   cluster_cols=cluster_cols, 
                             ### Columns configuration
                             annotation_names_col = annotation_names_col, labels_col = rownames(meta.data), 
                             show_colnames = show_colnames,annotation_col = annotation_col, gaps_col = gaps_col,
                             fontsize_col = fontsize_col, annotation_colors = ann_colors,
                             ### Rows configuration
                             labels_row = selection, show_rownames=show_rowname, fontsize_row = fontsize_row
    )
  } else {
    stop(paste("The provided heatmap plotting package ", heatmap.plot.package , " is not supported."))
  }
  
  ### Postprocessing of heatmap
  hm.pp <-cowplot::plot_grid(ggplot()+theme_void(), ggplotify::as.ggplot(hm), rel_heights=space.between.title.and.heatmap, ncol=1) + 
    ggtitle(title) + labs(subtitle=subtitle) + 
    theme(plot.title = element_text(face=plot.title.face, color=plot.title.color, 
                                    size=plot.title.size, vjust=move.title.vertically,
                                    hjust=move.title.horizontally), 
          plot.subtitle = element_text(face=plot.subtitle.face,color=plot.subtitle.color,
                                       size = plot.subtitle.size, vjust=move.subtitle.vertically, 
                                       hjust=move.subtitle.horizontally))
  
  if(return.data){
    return(list(hm.plot=ggplotify::as.ggplot(hm.pp), hm.selection.list=as.data.frame(selected.DE.table), hm.matrix=as.data.frame(hm.matrix)))
  } else {
    return(ggplotify::as.ggplot(hm.pp))
  }
}


plot_hm <- function(counts, daa_results, colData, annotated.peaks.df, order.by="padj"){
  
  counts <- counts %>% 
            tibble::add_column(Geneid=rownames(.)) %>% 
            merge(tar_read("processedAnnotatedPeaks"), by="Geneid") %>% 
            filter(!is.na(.$Symbol)) %>%
            tibble::column_to_rownames("Geneid") %>%
            group_by(Symbol) %>%
            summarise_at(colnames(counts[,1:6]), mean) %>%
            as.data.frame() %>%
            tibble::add_column(Rowname=.$Symbol) %>% 
            tibble::column_to_rownames(var="Rowname")
  
  results <- na.omit(daa_results) 
  results[results$condition=="Lyn.WT",]$log2FoldChange <- -results[results$condition=="Lyn.WT",]$log2FoldChange
  colData$samples <- rownames(colData)
  
  if(order.by=="padj"){
    results <- results %>% group_by(condition) %>% top_n(50, -padj)
    hm <- DoHeatmap(data=counts[, rownames(colData)], meta.data = colData, selection = results, 
                    do.selection = T, do.selection.top = 25, do.selection.group.by = "condition", 
                    do.selection.order.by="padj", do.selection.order.direction ="ascending", 
                    do.selection.feature.name="Symbol", do.selection.significance.variables = c("log2FoldChange","padj"), 
                    do.selection.significance.cutoffs = c(1, .05), do.selection.cutoff.relation = c(">=", "<="),
                    remove.batch.variable = NULL, do.selection.cutoff.l2fc = 1, do.selection.cutoff.padj = .05, 
                    do.selection.regulation = "up", do.selection.regulation.variable = "log2FoldChange",
                    selected.meta.data.columns=c("condition"), sample_name_column="samples", gaps_col=NULL, heatmap.plot.package = "pheatmap",
                    cell.colors=colorRampPalette(c("darkblue", "white", "red2"))(50), fontsize_row = 10, fontsize_col=10, 
                    cluster_rows=F, cluster_cols=F, scale="row", show_rowname=T, show_colnames=T, 
                    annotation_names_col=T, space.between.title.and.heatmap = c(0,100), return.data = F,
                    title ="Top 25 peaks per condition", subtitle = "Ordered by adjusted p-value", plot.title.face = "bold", plot.title.color = "darkred",
                    move.title.vertically = 0, move.title.horizontally = .2, plot.title.size = 14,
                    plot.subtitle.face = "bold", plot.subtitle.color = "#000080", move.subtitle.vertically = -1,
                    move.subtitle.horizontally = .3,plot.subtitle.size = 11, border_color="lightgrey",
                    col.annotation.colors=NA)
  } else if (order.by %in% c("l2fc","log2FoldChange")){
    hm <- DoHeatmap(data=counts[, rownames(colData)], meta.data = colData, selection = results, 
              do.selection = T, do.selection.top = 25, do.selection.group.by = "condition", 
              do.selection.order.by="log2FoldChange", do.selection.order.direction ="descending", 
              do.selection.feature.name="Symbol", do.selection.significance.variables = c("log2FoldChange","padj"), 
              do.selection.significance.cutoffs = c(1, .05), do.selection.cutoff.relation = c(">=", "<="),
              remove.batch.variable = NULL, do.selection.cutoff.l2fc = 1, do.selection.cutoff.padj = .05, 
              do.selection.regulation = "up", do.selection.regulation.variable = "log2FoldChange",
              selected.meta.data.columns=c("condition"), sample_name_column="samples", gaps_col=NULL, heatmap.plot.package = "pheatmap",
              cell.colors=colorRampPalette(c("darkblue", "white", "red2"))(50), fontsize_row = 10, fontsize_col=10, 
              cluster_rows=F, cluster_cols=F, scale="row", show_rowname=T, show_colnames=T, 
              annotation_names_col=T, space.between.title.and.heatmap = c(0,100), return.data = F,
              title = "Top 25 peaks per condition", subtitle = "Ordered by log2FoldChange", plot.title.face = "bold", plot.title.color = "darkred",
              move.title.vertically = 0, move.title.horizontally = .2, plot.title.size = 14,
              plot.subtitle.face = "bold", plot.subtitle.color = "#000080", move.subtitle.vertically = -1,
              move.subtitle.horizontally = .3,plot.subtitle.size = 11, border_color="lightgrey",
              col.annotation.colors=NA)
  } else {
    Stop("Ordering option not supported!")
  }
  
  return(hm)
  
}


plot_Volcano <- function(  DEG.table, l2fc =1,  padj = 0.05, pt.trans = 0.5, pt.size = 2, repel.max.overlap = 200, xlim = NULL, ylim=NULL, pt.border.color="black",
                           repel.pt.padding = 0, xlab.color = "red", ylab.color = "blue", vertical.linetype = "dashed",
                           horizontal.linetype = "dashed", label.only.significant =T, show.only.significant=F, repel.nudge_x = .15, repel.segment.alpha	= 0.5,
                           repel.box.padding = 0.1, repel.nudge_y = 0.25, repel.segment.curvature = -0.1, segment_color="#000080", repel.segment.color = "lightgrey",
                           repel.segment.ncp = 1, repel.segment.angle = 20, xlab = NULL, ylab = NULL, coord_flip=F, repel.min.segment.length=0.5,
                           selected.genes, colors.selected.genes, theme=theme_classic(), add.lines=T, point.colors,  repel.force_pull=1, repel.force=1,
                           title="Volcano plot", title.face="bold", title.color="#000080", title.hjust=0.5,   max.iter	=50000, max.time=5, repel.segment.linetype="dashed",
                           repel.text.size=4, scale_x_countinuous=NULL, subtitle=NULL, subtitle.face="bold", subtitle.color="steelblue", subtitle.hjust=0.5, repel=T ){ 
  
  if(label.only.significant){
    selected.genes <- selected.genes[selected.genes %in% DEG.table[abs(DEG.table$log2FoldChange) >= l2fc & DEG.table$padj <= padj,]$Symbol]
  }
  if(show.only.significant){
    DEG.table <- DEG.table[DEG.table$padj <= padj & abs(DEG.table$log2FoldChange >=l2fc),]
  }
  
  lownegsign  <- DEG.table$padj <= padj & DEG.table$log2FoldChange > -l2fc & DEG.table$log2FoldChange < 0
  highnegsign <- DEG.table$padj <= padj & DEG.table$log2FoldChange <= -l2fc
  lowpossign  <- DEG.table$padj <= padj & DEG.table$log2FoldChange < l2fc & DEG.table$log2FoldChange >=0
  highpossign <- DEG.table$padj <= padj & DEG.table$log2FoldChange >= l2fc
  highnonsign <- DEG.table$padj > padj & abs(DEG.table$log2FoldChange) >= l2fc

  are.selected <- DEG.table$Symbol %in% selected.genes
  DEG.table$category <-ifelse(are.selected, "sel", ifelse(lownegsign, "lownegsign",
                              ifelse(highnegsign, "highnegsign",
                                     ifelse(lowpossign, "lowpossign",
                                            ifelse (highpossign,"highpossign",
                                                    ifelse(highnonsign,"highnonsign","nonsign"))))))
  
  DEG.table$color <- "black"
  
  vp <- ggplot(DEG.table, aes(x=log2FoldChange, y=-log10(padj))) + geom_point(size=pt.size, aes( fill=factor(category)), color=pt.border.color, alpha=pt.trans, shape=21)+ 
    scale_fill_manual(values = point.colors) +theme
  
  if(is.null(ylab)){ vp <- vp + ylab("-log10(adjusted P-value) - Significance") } else { vp <- vp + ylab(ylab) }
  if(is.null(xlab)){ vp <- vp + xlab("log2(FC) - Regulation") } else { vp <- vp + xlab(xlab) }
  
  if(add.lines){
    vp <- vp +  geom_vline(xintercept = -l2fc, linetype="dashed"  )+
      geom_vline(xintercept = l2fc, linetype="dashed" )  +
      geom_hline(yintercept = -log10(padj), linetype="dashed" )
  }
  if(class(selected.genes) != "list"){
    selected.genes <- list(selected.genes)
    colors.selected.genes <- list(colors.selected.genes)
  }
  
  pa <- padj
  
  if(repel){
  for(t in 1:length(selected.genes)){
    
    vp <- vp + geom_text_repel(aes(label=ifelse(Symbol %in% selected.genes[[t]] & ((abs(log2FoldChange) >= l2fc) & DEG.table$padj <= pa), as.character(Symbol),'')), force_pull = repel.force_pull, force = repel.force,
                               nudge_x = repel.nudge_x, box.padding = repel.box.padding, nudge_y = repel.nudge_y, min.segment.length = repel.min.segment.length,
                               segment.curvature = repel.segment.curvature, segment.ncp = repel.segment.ncp, segment.linetype=repel.segment.linetype, 
                               segment.angle = repel.segment.angle, max.overlaps = repel.max.overlap, segment.color = repel.segment.color, segment.alpha=repel.segment.alpha,
                               color= colors.selected.genes[[t]], size=repel.text.size, segment.inflect	=T, lineheight=0.2)
  }
  }
  vp <- vp   +
    ggtitle(title)+
    theme(plot.title = element_text(face = title.face, hjust = title.hjust, colour = title.color)) +
    theme(legend.position = "none")
  
  
  if(!is.null(subtitle)){
    vp <- vp + labs(subtitle = subtitle)+theme(plot.subtitle = element_text(face=subtitle.face, hjust = subtitle.hjust, color = subtitle.color))
  }
  
  if(!is.null(xlim)){
    vp <- vp + xlim(xlim)
  }
  
  if(!is.null(ylim)){
    vp <- vp + ylim(ylim)
  }
  
  if(coord_flip){
    vp <- vp + coord_flip()
  }
  
  if(!is.null(scale_x_countinuous)){
    vp <- vp + scale_x_countinuous
  }
  
  return(vp)
}


getFEResults <- function(de_results, DB="standard", direction="UP", ordered=F, organism="hsapiens", gene.column.name="gene", simple.table=T){
  tmp <- de_results
  tmp.go <- NULL
  if(nrow(tmp) > 0){
    if(DB != "standard"){
      tmp.go <- gost( query = unique(tmp[, gene.column.name]),
                      organism = DB, ordered_query = ordered,
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = T,
                      user_threshold = 0.05, correction_method = "g_SCS",
                      domain_scope = "annotated", custom_bg = NULL,
                      numeric_ns = "", as_short_link = FALSE)
    } else {
      tmp.go <- gost( query = unique(tmp[, gene.column.name]),
                      organism = organism, ordered_query = ordered,
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = T,
                      user_threshold = 0.05, correction_method = "g_SCS",
                      domain_scope = "annotated",
                      sources = c("GO:MF", "GO:BP", "GO:CC", "REAC", "TF", "CORUM", "HP"),
                      custom_bg = NULL,
                      numeric_ns = "", as_short_link = FALSE)
    }
  }
  if(is.null(tmp.go)) return(NULL)
  if(simple.table){return(tmp.go$result)}else{return(tmp.go)}
}

