TowerPlot <- function(deg_data, dap_data, annotation.name="annotation", region.name="", gene.name="Symbol",
                      padj.trans=0.05, padj.atac=0.05, only.significant.peaks=F, color.border.sign.peaks="black",
                      l2fc.trans=1, color.border.unsign.peaks="darkgrey", color.region.of.interest="yellow", 
                      roi=c("Promoter","3' UTR"), gene.label.angle=90, point.shape=23, min.point.size=3, 
                      max.point.size=10, atac.regulation.color.low="steelblue",  roi.rev=F, ylim=c(-10,10), 
                      stack.inc=0.5, atac.regulation.color.high="darkred", x.lab="Top 50 Trans DEG", 
                      y.lab="Log2FC Trans", top_n=25, stroke.high=1.5, stroke.std=1, move.gene.name.down=1.2,
                      color.gene.name="black", order.by="padj", order.dir="desc", theme=theme_pubr(),
                      title="", title.face="bold", title.color="black", title.hjust=0.5,
                      legend.title.color="black", legend.title.face="bold", legend.title.background.fc="Log2FC ATAC",
                      legend.title.background.padj="-Log10 P-Value", legend.title.point.types="Peak type",
                      legend.background.padj.point.color="black", return=c("plot", "plot+data","data"),
                      legend.padj.point.size.percentiles=c(.05, .25, .5, .75, 1)){
  
  ### Combine and reshape data for downstream usage.
  if(is.null(deg_data$RegulationDirection)){
    deg_data <- deg_data %>% add_column(RegulationDirection=if_else(.$log2FoldChange>=0, 1, -1), .before = gene.name)
  }
  deg_data$log2FoldChange.gw <- abs(deg_data$log2FoldChange)
  deg_data$RegulationDirection <- as.factor(deg_data$RegulationDirection)
  deg_data <- deg_data %>%
    filter((log2FoldChange.gw) >= l2fc.trans & padj <= padj.trans) %>%
    group_by(RegulationDirection) %>%
    dplyr::top_n(wt = (log2FoldChange.gw), n= top_n) %>% 
    arrange(RegulationDirection, log2FoldChange)
  if(is.null(deg_data$ID)){
    deg_data <- deg_data %>%
      add_column(ID=1:nrow(deg_data), .before="RegulationDirection")
  }
  
  dap_data <- dap_data[dap_data$Symbol %in% deg_data$Symbol,]
  
  if(roi.rev){corr_trans_epi <- deg_data %>% merge(dap_data[!grepl(region.name,dap_data[,annotation.name]),], by=gene.name)}
  else{corr_trans_epi <- deg_data %>% merge(dap_data[grepl(region.name,dap_data[,annotation.name]),], by=gene.name)}
  
  
  print(corr_trans_epi)
  
  if(only.significant.peaks){ corr_trans_epi <- corr_trans_epi %>% filter(padj.y <= padj.atac)} 
  order.var <- ifelse(order.by=="l2fc", "log2FoldChange.y", ifelse(order.by=="padj", "padj.y", ""))
  corr_trans_epi$s.annotation <- gsub(x=corr_trans_epi$annotation, pattern = " \\(.*\\)", "")
  
  ### Compute position of peaks on each stack across y axis. Variable stac.inc define the density of the points on the
  ### the stack
  corr_trans_epi <- corr_trans_epi  %>% arrange(get(gene.name), abs(get(order.var))) %>% group_by(get(gene.name)) %>%
    mutate(groupid = (cur_group_id())) %>%  mutate(gid = cumsum(groupid)/cur_group_id()) %>% ungroup() %>% 
    add_column(position=(ifelse(.$gid==1,0,stack.inc*(.$gid-1))+.$log2FoldChange.x), .before="log2FoldChange.x")
  
  ### Only annotate first peak with gene name for each gene.
  corr_trans_epi <- corr_trans_epi %>% group_by(groupid) %>% mutate(groupid=cur_group_id())
  corr_trans_epi$label <- ifelse(corr_trans_epi$gid==1, corr_trans_epi$Symbol, "")
  ### Define relative size of points based on si
  corr_trans_epi$size <- 10+(-log10(corr_trans_epi$padj.y)/max(-log10(corr_trans_epi$padj.y)))
  
  
  labels <- ifelse(corr_trans_epi$s.annotation %in% roi & corr_trans_epi$padj.y <= padj.atac, paste(roi,collapse="|"),
                   ifelse(corr_trans_epi$padj.y <= padj.atac,
                          "Significant",
                          "NS")
  )
  colors <- c(color.region.of.interest, color.border.sign.peaks, color.border.unsign.peaks)
  names(colors) <- c(paste(roi,collapse="|"), "Significant", "NS")
  corr_trans_epi$color <- ifelse(corr_trans_epi$s.annotation %in% roi & corr_trans_epi$padj.y <= padj.atac,
                                 color.region.of.interest,
                                 ifelse(corr_trans_epi$padj.y <= padj.atac,
                                        color.border.sign.peaks,
                                        color.border.unsign.peaks))
  
  corr_trans_epi$stroke <- ifelse(corr_trans_epi$s.annotation %in% roi & corr_trans_epi$padj.y <= padj.atac, stroke.high, stroke.std)
  
  
  q <- (quantile(corr_trans_epi$padj.y, legend.padj.point.size.percentiles))
  q <- q[order(q)]
  
  .convert <- function(corr_trans_epi, pv){
    return(10+(-log10(pv)/max(-log10(corr_trans_epi$padj.y))))
  }
  
  ### Return only data.
  if(return[1]=="data"){ return(corr_trans_epi)}
  ### Or either plot+data or data.
  else {
    ### Generate plot.
    p <- ggplot(corr_trans_epi, aes(x=ID, y=position, label=str_pad(label,12, side = "left"), fill=log2FoldChange.y, color=labels)) + 
      
      geom_point(shape=point.shape, aes(size=size,  stroke=stroke)) + ylim(ylim[1],ylim[2]) + 
      
      ### Define multiple scales of geometric points.
      scale_color_manual(values=colors) + 
      scale_size_continuous(corr_trans_epi$size, range = c(min.point.size, max.point.size), 
                            breaks=.convert(corr_trans_epi, q), 
                            labels=paste(str_pad(as.character(as.double(round(-log10(q),1))), width = 5, side="right"), 
                                         str_pad(paste("(",names(q),")",sep=""), width = 8, side = "right"),sep="")) + 
      scale_fill_gradient2(low=atac.regulation.color.low, mid="white", high=atac.regulation.color.high) +
      
      # Format legends
      guides(size=guide_legend(override.aes = aes(fill="black"), title.theme = element_text(face=legend.title.face, color=legend.title.color), title=legend.title.background.padj, order=2)) +
      guides(fill=guide_colorbar(title.theme = element_text(face=legend.title.face, color=legend.title.color), title=legend.title.background.fc, order = 1)) +
      guides(color=guide_legend(override.aes= aes(size=3),title.theme = element_text(face=legend.title.face, color=legend.title.color), title=legend.title.point.types, order=3))
    
    ### Theme adjusted
    p <- p + theme
    
    ### Remove x-axis ticks, text and line of both axis)
    p <- p +  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+ theme(axis.line = element_blank()) +
      ### Set axis titles
      xlab(x.lab) + ylab(y.lab) +
      ### Gene name format
      geom_text(angle=90, hjust=move.gene.name.down, color=color.gene.name) +
      ### Add horizontal line at l2fc=0
      geom_hline(yintercept=0, linetype="dashed", color="black", alpha=0.5)
    
    ### Set and format plot title.
    p <- p + ggtitle(title) + theme(plot.title=element_text(face=title.face, color=title.color, hjust = title.hjust))
    
    ### Return plot and data.
    if(return[1]=="plot+data"){return(list(plot=p,data=corr_trans_epi))}
    ### Or only plot.
    else {return(p)}
  }
}


correlation_plot <- function(DAA, DEA, DAA.name="ATAC-Seq", DEA.name="Transcriptome"){
  ## Data
  Table1 <- na.omit(DAA)
  Table1 <- Table1 %>%
    group_by(condition, Symbol,  transcriptId) %>%
    summarise_at(vars(log2FoldChange, padj), mean)
  Table1 <- as.data.frame(Table1)
  Table2 <- DEA
  table1_name <- DAA.name
  table2_name <- DEA.name
  
  colnames(Table1) <- paste("first.",  colnames(Table1) , sep="") 
  colnames(Table1)[2] <- gsub(pattern = "first.", replacement = "", x = colnames(Table1)[2])
  colnames(Table1)[3] <- gsub(pattern = "first.", replacement = "", x = colnames(Table1)[3])
  colnames(Table2) <- paste("second.",  colnames(Table2) , sep="") 
  colnames(Table2)[1] <- gsub(pattern = "second.", replacement = "", x = colnames(Table2)[1])
  
  ### Configuration
  ES1 <- 4
  PV1 <- 5
  ES2 <- 3
  PV2 <- 7
  SV <- c("Symbol")
  
  ## Logic
  cut_lfc_1 <- 0
  cut_padj_1 <- 1
  cut_lfc_2 <- 0
  cut_padj_2 <- 1
  cut <- .56
  include.background <- F
  
  ## UI
  
  # Main correlation plot
  pt.size=2
  pt.shape=21
  pt.border.color="black"
  theme <- theme_classic()
  both.up.color <- "red"
  both.dn.color <- "steelblue"
  up.dn.color <- "yellow"
  dn.up.color <- "lightblue"
  nol.sign.color <- "lightgrey"
  xlab <- paste("Log2FoldChange (",table1_name,")",sep="")
  ylab <- paste("Log2FoldChange (",table2_name,")",sep="")
  
  # Marginal distribution
  bins = 100
  hist.border.color <- "black"
  hist.color <- "white"
  density.lwd <- 1
  density.lt <- "solid"
  density.lc <- "steelblue"
  
  # Title
  title.color = "black"
  title.face = "bold"
  title.hjust <- 0.4
  title.font.size <- 20
  title.text <- paste("Correlation plot - ",DAA.name," vs. ",DEA.name,sep="")
  free.label.color <- "black"
  free.label.fill <- "white"
  label.pad.width <- 5
  
  
  Table1 <- Table1[abs(Table1[,ES1]) >= cut_lfc_1 & Table1[,PV1] <= cut_padj_1,]
  Table2 <- Table2[abs(Table2[,ES2]) >= cut_lfc_2 & Table2[,PV2] <= cut_padj_2,]
  Table1 <- Table1[, c(SV, colnames(Table1)[ES1], colnames(Table1)[PV1])]
  Table2 <- Table2[, c(SV, colnames(Table2)[ES2], colnames(Table2)[PV2])]
  
  if(include.background){
    Background.Table1 <- DEG.ext[[2]]
    Background.Table2 <- DEG.ox.ext[[2]]
    colnames(Background.Table1)[1:8] <- paste("first.bg.",  colnames(Background.Table1)[1:8] , sep="") 
    colnames(Background.Table1)[1] <- gsub(pattern = "first.bg.", replacement = "", x = colnames(Background.Table1)[1])
    colnames(Background.Table1)[8] <- gsub(pattern = "first.bg.", replacement = "", x = colnames(Background.Table1)[8])
    colnames(Background.Table2)[1:8] <- paste("second.bg.",  colnames(Background.Table2)[1:8] , sep="") 
    colnames(Background.Table2)[1] <- gsub(pattern = "second.bg.", replacement = "", x = colnames(Background.Table2)[1])
    colnames(Background.Table2)[8] <- gsub(pattern = "second.bg.", replacement = "", x = colnames(Background.Table2)[8])
    BES1 <- 3
    BES2 <- 3
    BPV1 <- 7
    BPV2 <- 7
    Background.Table1 <- Background.Table1[abs(Background.Table1[, BES1]) >= cut & Background.Table1[, BPV1] <= cut_padj_1, ]
    Background.Table2 <- Background.Table2[abs(Background.Table2[, BES2]) >= cut & Background.Table2[, BPV2] <= cut_padj_2, ]
    Background.Table1 <- Background.Table1[, c(SV, colnames(Background.Table1)[BES1], colnames(Background.Table1)[BPV1])]
    Background.Table2 <- Background.Table2[, c(SV, colnames(Background.Table2)[BES2], colnames(Background.Table2)[BPV2])]
  }
  ES1 <- 3
  PV1 <- 4
  ES2 <- 3
  PV2 <- 4
  combined <- merge(Table1, Table2, by=SV)
  
  if(include.background){
    combined <- merge(combined, Background.Table1, by=SV)
    combined <- merge(combined, Background.Table2, by=SV)
  }
  
  #combined[combined$first.log2FoldChange  0 & combined$first.bg.log2FoldChange > 0 & (combined$second.log2FoldChange < 0 & combined$second.bg.log2FoldChange > 0), ]
  
  
  
  ES1 <- 2
  PV1 <- 3
  ES2 <- 4
  PV2 <- 5
  
  combined$up.up <- ifelse((combined[,ES1]) >= cut & (combined[,ES2]) >= cut, T, F)
  combined$dn.dn <- ifelse((combined[,ES1]) <= -cut & (combined[,ES2]) <= -cut, T, F)
  combined$up.dn <- ifelse((combined[,ES1]) >= cut & (combined[,ES2]) <= -cut, T, F)
  combined$dn.up <- ifelse((combined[,ES1]) <= -cut & (combined[,ES2]) >= cut, T, F)
  
  combined$up.ns <- ifelse((combined[,ES1]) >= cut & ((combined[,ES2]) >= -cut &(combined[,ES2]) <= cut), T, F)
  combined$ns.up <- ifelse(((combined[,ES1]) >= -cut &(combined[,ES1]) <= cut) & (combined[,ES2]) >= cut, T, F)
  combined$dn.ns <- ifelse((combined[,ES1]) <= -cut & ((combined[,ES2]) >= -cut &(combined[,ES2]) <= cut), T, F)
  combined$ns.dn <- ifelse(((combined[,ES1]) >= -cut &(combined[,ES1]) <= cut) & (combined[,ES2]) <= -cut, T, F)
  
  nr.up.up <- nrow(combined[combined$up.up,])
  nr.dn.dn <- nrow(combined[combined$dn.dn,])
  nr.up.dn <- nrow(combined[combined$up.dn,])
  nr.dn.up <- nrow(combined[combined$dn.up,])
  nr.up.ns <- nrow(combined[combined$up.ns,])
  nr.ns.up <- nrow(combined[combined$ns.up,])
  nr.dn.ns <- nrow(combined[combined$dn.ns,])
  nr.ns.dn <- nrow(combined[combined$ns.dn,])
  nr.else <- nrow(combined)-(nr.up.up+nr.dn.dn+nr.up.dn+nr.dn.up)-(nr.up.ns+nr.ns.up+nr.dn.ns+nr.ns.dn)
  
  
  df <- data.frame(Upregulated=c(nr.up.up, nr.dn.up, nr.ns.up), Downregulated=c(nr.up.dn, nr.dn.dn, nr.ns.dn), NS_L2FC=c(nr.up.ns, nr.dn.ns, nr.else))
  colnames(df) <- paste(colnames(df), " (",table1_name,")",sep="")
  rownames(df) <- c("Upregulated", "Downregulated", "NS_L2FC")
  rownames(df) <- paste(rownames(df), " (",table2_name,")",sep="")
  df <- gridExtra::tableGrob(df, theme=gridExtra::ttheme_minimal())
  
  combined$Relation <- ifelse(combined$up.up, "Up.Up", 
                              ifelse(combined$dn.dn, "Dn.Dn", 
                                     ifelse(combined$up.dn, "Up.Dn", 
                                            ifelse(combined$dn.up, "Dn.Up", 
                                                   ifelse(combined$up.ns, "Up.Ns", 
                                                          ifelse(combined$ns.up, "Ns.Up", 
                                                                 ifelse(combined$dn.ns, "Dn.Ns", 
                                                                        ifelse(combined$ns.dn, "Ns.Dn", "Ns.Ns"))))))))
  combined$Relation <- as.factor(combined$Relation)
  named.color.vector <- c(both.up.color, both.dn.color, up.dn.color, dn.up.color, nol.sign.color, nol.sign.color, nol.sign.color, nol.sign.color, nol.sign.color)
  names(named.color.vector) <- c("Up.Up", "Dn.Dn", "Up.Dn", "Dn.Up", "Up.Ns", "Ns.Up", "Dn.Ns", "Ns.Dn", "Ns.Ns")
  
  
  linetype= "dashed"
  
  merge(combined, tar_read(PDEA), by="Symbol")
  
  ## Main plot
  p <- ggplot(combined, aes(x=first.log2FoldChange, y=second.log2FoldChange, fill=Relation)) + 
    geom_point(size=pt.size, shape = pt.shape, color = pt.border.color) + theme_classic() + 
    scale_fill_manual(values = named.color.vector) + 
    geom_hline(yintercept = cut, linetype=linetype) + geom_hline(yintercept = -cut, linetype=linetype) + 
    geom_vline(xintercept = cut, linetype=linetype) + geom_vline(xintercept = -cut, linetype=linetype) +
    xlab(xlab) + ylab(ylab)
  
  
  
  tp <- ggplot(combined, aes(x=first.log2FoldChange)) + 
    geom_histogram(aes(y = ..density..), bins = bins, colour = hist.border.color, fill = hist.color) + 
    geom_density(lwd = density.lwd, linetype = density.lt,  colour = density.lc)+ 
    theme_classic() + theme(axis.title.x = element_blank())
  
  rp <- ggplot(combined, aes(x=second.log2FoldChange)) + 
    geom_histogram(aes(y = ..density..), bins = bins, colour = hist.border.color, fill = hist.color) + 
    geom_density(lwd = density.lwd, linetype = density.lt,  colour = density.lc)+ 
    theme_classic()  + theme(axis.title.y = element_blank()) + coord_flip()
  
  return(cowplot::plot_grid(
    ggplotify::as.ggplot(cowplot::plot_grid(
      cowplot::plot_grid(tp, ggplot()+theme_void(), rel_widths = c(12,2)),
      cowplot::plot_grid(  p + theme(legend.position = "none"), rp, rel_widths = c(12,2)),
      rel_heights = c(2,8), ncol=1)) + 
      ggtitle(title.text) + theme(plot.title = element_text(face=title.face, 
                                                            size = title.font.size, 
                                                            hjust=title.hjust, 
                                                            colour = title.color)), df, rel_heights = c(4,1), ncol=1))
}
