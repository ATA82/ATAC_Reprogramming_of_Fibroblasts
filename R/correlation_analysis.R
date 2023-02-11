library(ggplot2)
library(ggrepel)
library(ggvenn)
library(dplyr)
library(tidyr)
# ali 
core.filter <- function(combined, l2fc, padj, which, regulation){
  
  if(regulation=="up"){
    combined <- combined %>%
      filter(!!sym(paste(which,".log2FoldChange",sep="")) >=l2fc & !!sym(paste(which,".padj",sep="")) <=padj)
  }
  else if(regulation=="dn"){
    combined <- combined %>%
      filter(!!sym(paste(which,".log2FoldChange",sep="")) <= -l2fc & !!sym(paste(which,".padj",sep="")) <=padj)
  }
  else if(regulation=="de"){
    combined <- combined %>%
      filter(abs(!!sym(paste(which,".log2FoldChange",sep=""))) >= l2fc & !!sym(paste(which,".padj",sep="")) <=padj)
  }
  else if(regulation=="ns"){
    combined <- combined %>%
      filter((!!sym(paste(which,".log2FoldChange",sep="")) > -l2fc & !!sym(paste(which,".log2FoldChange",sep="")) < l2fc) |
               !!sym(paste(which,".padj",sep="")) > padj)
  }
  else{
    stop("regulation type not supported")
  }
}

combined.filter <- function(combined, l2fc, padj, bg.l2fc, bg.padj, which, regulation, bg.which, bg.regulation, bg=T){
  if(bg){
    return(combined %>% core.filter(l2fc, padj, which, regulation) %>% core.filter(bg.l2fc, bg.padj, bg.which, bg.regulation))
  } else{
    return(combined %>% core.filter(l2fc, padj, which, regulation))
    
  }
}

core.category <- function(combined, l2fc, padj, which){
  selection <- combined %>% select(!!sym(paste(which,".log2FoldChange",sep="")), !!sym(paste(which,".padj",sep="")))
  selection$core.category <- ifelse(selection[,1] >= l2fc & selection[,2] <= padj, "up", ifelse (selection[,1] <= -l2fc & selection[,2] <= padj, "dn", "ns"))
  return(selection$core.category)
}

combined.category <- function(combined, l2fcs, padjs, whichs){
  
  list <- vector("list", length(l2fcs))
  
  for(l in 1:length(l2fcs)){
    list[[l]] <- (core.category(combined, l2fcs[l], padjs[l], whichs[l]))
  }
  t <- mapply(list, FUN = "paste")
  t <- as.data.frame(t)
  regulation <- t %>% unite("Regulation", sep = ".")
  
  return(regulation$Regulation)
}

overlap.table <- function(combined, dea_1_name, dea_2_name){
  df <- combined$CorrelationVariables %>% table()  %>% as.data.frame()  %>% t() %>% as.data.frame()
  colnames(df) <- df[1,]
  df <- df[-1,]
  CorrelationVariablesCategories <- c("dn.dn", "dn.up", "dn.ns", "up.dn", "up.up", "up.ns", "ns.dn", "ns.up", "ns.ns")
  sde <- setdiff(1:9, which(CorrelationVariablesCategories %in% colnames(df)))
  for(sd in sde){
    df[,CorrelationVariablesCategories[sd]] <- 0
  }
  df <- df[, CorrelationVariablesCategories]
  df <- data.frame(Upregulated=c(df$up.up, df$dn.up, df$ns.up),
                   Downregulated=c(df$up.dn, df$dn.dn, df$ns.dn),
                   NS_L2FC=c(df$up.ns, df$dn.ns, df$ns.ns))
  colnames(df) <- paste(colnames(df), " (",dea_1_name,")",sep="")
  rownames(df) <- c("Upregulated", "Downregulated", "NS_L2FC")
  rownames(df) <- paste(rownames(df), " (",dea_2_name,")",sep="")
  return(df)
}


correlation_analysis.2 <- function(
  
  ### Analysis data
  # dea results
  dea_1, dea_2, bg.dea_1=NULL, bg.dea_2=NULL,
  # positions of effect sizes and adjusted p-value.
  l2fc_pos_1=3, l2fc_pos_2=3, padj_pos_1=7, padj_pos_2=7,
  Bl2fc_pos_1=3, Bl2fc_pos_2=3, Bpadj_pos_1=7, Bpadj_pos_2=7,
  gene_fields = c("ensembl_gene_id", "symbol"),
  
  # cutoffs
  l2fc_1 = 1, padj_1 = 0.05,
  l2fc_2 = 1, padj_2 = 0.05,
  rl2fc_1 = 0, rpadj_1 = .05,
  rl2fc_2 = 0, rpadj_2 = .05,
  
  
  fg_regulation_criteria  = c("de", "up", "dn", "ns"),
  fg_regulation_criteria_1= c("de", "up", "dn", "ns"),
  fg_regulation_criteria_2= c("de", "up", "dn", "ns"),
  
  ### Involving Background Analysis
  include.background=T,
  bg_l2fc_1 = 1, bg_padj_1 = 0.05,
  bg_l2fc_2 = 1, bg_padj_2 = 0.05,
  bg_regulation_criteria  = c("de", "up", "dn", "ns"),
  bg_regulation_criteria_1= c("de", "up", "dn", "ns"),
  bg_regulation_criteria_2= c("de", "up", "dn", "ns"),
  
  ### Include only genes, which are regulated the same direction in both conditions of interest.
  directed = T,
  ### Include only genes, which are regulated the same way with respect to background in both
  ### conditions of interest.
  mode="both",
  
  ### GUI
  # Main correlation plot
  dea_1_name="DEA_1", dea_2_name="DEA_2", pt.size=2, pt.shape=21, pt.border.color="black",
  theme = theme_classic(), linetype="dotted", both.up.color = "red", xlab=NULL, ylab=NULL,
  dn.up.color = "lightblue", ns.color = "lightgrey", both.dn.color = "steelblue",
  up.dn.color = "yellow",
  # Marginal distribution
  bins = 100, hist.border.color = "black", hist.color = "white", density.lwd = 1,
  density.lt = "solid", density.lc = "steelblue",
  # Venn Diagram
  inlcude.venn.diagram=T, venn.diagram.colors=c("#0073C2FF", "#EFC000FF"), venn.set.name.size=4,
  # Title
  title.color = "black", title.face = "bold", title.hjust = 0.4, title.font.size = 20,
  title.text = "Correlation plot - OsxNTR vs. ColNTR", free.label.color = "black",
  free.label.fill = "white", label.pad.width = 5){
  
  ### Convert table to unified format.
  
  ## Convert table of condition of interest from first experiment
  dea_1 <- dea_1[, c(which(colnames(dea_1) %in% gene_fields), l2fc_pos_1, padj_pos_1)]
  colnames(dea_1) <- c(gene_fields, "first.log2FoldChange", "first.padj")
  ## Convert table of condition of interest from second experiment
  dea_2 <- dea_2[, c(which(colnames(dea_2) %in% gene_fields), l2fc_pos_2, padj_pos_2)]
  colnames(dea_2) <- c(gene_fields, "second.log2FoldChange", "second.padj")
  combined <- merge(dea_1, dea_2, by=gene_fields)
  
  if(include.background){
  ## Convert table of background condition from first experiment
  bg.dea_1 <- bg.dea_1[, c(which(colnames(bg.dea_1) %in% gene_fields), Bl2fc_pos_1, Bpadj_pos_1)]
  colnames(bg.dea_1) <- c(gene_fields, "first.bg.log2FoldChange", "first.bg.padj")
  ## Convert table of background condition from secon experiment
  bg.dea_2 <- bg.dea_2[, c(which(colnames(bg.dea_2) %in% gene_fields), Bl2fc_pos_2, Bpadj_pos_2)]
  colnames(bg.dea_2) <- c(gene_fields, "second.bg.log2FoldChange", "second.bg.padj")
  ### Combine all conditions.
  combined.bg <- merge(bg.dea_1, bg.dea_2, by=gene_fields)
  combined <- merge(combined, combined.bg, by=gene_fields)
  }
  bg = include.background
  
  
  ## Filter table by specified cutoffs (relaxed and strict)
  if(mode=="one"){
    combined.relaxed.1 <- combined %>%
      combined.filter(l2fc = rl2fc_1, padj = rpadj_1, bg.l2fc = bg_l2fc_1, bg.padj = bg_padj_1,
                      which = "first", bg.which = "first.bg", regulation = fg_regulation_criteria_1,
                      bg.regulation = bg_regulation_criteria_1, bg=bg)
    combined.relaxed.2 <- combined %>%
      combined.filter(l2fc = rl2fc_2, padj = rpadj_2, bg.l2fc = bg_l2fc_2, bg.padj = bg_padj_2,
                      which = "second", bg.which = "second.bg", regulation = fg_regulation_criteria_2,
                      bg.regulation = bg_regulation_criteria_2, bg=bg)
    
    combined.relaxed   <- combined %>% filter(combined[,gene_fields[1]] %in% c(combined.relaxed.1[,gene_fields[1]], combined.relaxed.2[,gene_fields[1]]))
    combined.strict.1  <- combined %>%
      combined.filter(l2fc = l2fc_1, padj = padj_1, bg.l2fc = bg_l2fc_1, bg.padj = bg_padj_1,
                      which = "first", bg.which = "first.bg", regulation = fg_regulation_criteria_1,
                      bg.regulation = bg_regulation_criteria_1, bg=bg)
    combined.strict.2  <- combined %>%
      combined.filter(l2fc = l2fc_2, padj = padj_2, bg.l2fc = bg_l2fc_2, bg.padj = bg_padj_2,
                      which = "second", bg.which = "second.bg", regulation = fg_regulation_criteria_2,
                      bg.regulation = bg_regulation_criteria_2, bg=bg)
    
    combined.strict  <- combined %>% filter(combined[,gene_fields[1]] %in% c(combined.strict.1[,gene_fields[1]], combined.strict.2[,gene_fields[1]]))
    combined.relaxed$RegulationCriteria <- combined.relaxed[,gene_fields[1]] %in% combined.strict[,gene_fields[1]]
    combined.strict$RegulationCriteria <- T
    
    
  #  print(combined.strict)
  #  print(combined.relaxed)
    
  } else if(mode=="both"){
    combined.relaxed <- combined %>%
      combined.filter(l2fc = rl2fc_1, padj = rpadj_1, bg.l2fc = bg_l2fc_1, bg.padj = bg_padj_1,
                      which = "first", bg.which = "first.bg",
                      regulation = fg_regulation_criteria_1, bg.regulation = bg_regulation_criteria_1, bg=bg) %>%
      combined.filter(rl2fc_2, rpadj_2, bg_l2fc_2,
                      bg_padj_2, which = "second", bg.which = "second.bg", regulation = fg_regulation_criteria_2,
                      bg.regulation = bg_regulation_criteria_2, bg=bg)
    combined.strict <- combined %>%
      combined.filter(l2fc_1, padj_1, bg_l2fc_1,
                      bg_padj_1, which = "first", bg.which = "first.bg", regulation = fg_regulation_criteria_1,
                      bg.regulation = bg_regulation_criteria_1, bg=bg) %>%
      combined.filter(l2fc_2, padj_2, bg_l2fc_2,
                      bg_padj_2, which = "second", bg.which = "second.bg", regulation = fg_regulation_criteria_2,
                      bg.regulation = bg_regulation_criteria_2, bg=bg)
    

    combined.relaxed$RegulationCriteria <- combined.relaxed[,gene_fields[1]] %in% combined.strict[,gene_fields[1]]
    combined.strict$RegulationCriteria <- T
  }
  combined.strict$CorrelationVariables <- combined.category(combined.strict,
                                                            l2fcs=c(l2fc_1,l2fc_2), padjs=c(padj_1,padj_2),
                                                            whichs =  c("first", "second"))
  
  combined.strict <- combined.strict[, c(ncol(combined.strict), ncol(combined.strict)-1, 1:(ncol(combined.strict)-2))]
  combined.relaxed$CorrelationVariables <- combined.category(combined.relaxed,
                                                             l2fcs=c(l2fc_1,l2fc_2), padjs=c(padj_1,padj_2),
                                                             whichs = c("first", "second"))

  combined.relaxed <- combined.relaxed[, c(ncol(combined.relaxed), ncol(combined.relaxed)-1, 1:(ncol(combined.relaxed)-2))]

  
  ### Compute size of each category of combined regulation.
  df.relaxed <- overlap.table(combined.relaxed, dea_1_name = dea_1_name, dea_2_name = dea_2_name)
  rownames(df.relaxed) <- paste(rownames(df.relaxed), " (",dea_2_name,")",sep="")
  df.relaxed <- gridExtra::tableGrob(df.relaxed, theme=gridExtra::ttheme_minimal())
  
  df <- overlap.table(combined.relaxed, dea_1_name = dea_1_name, dea_2_name = dea_2_name)
  rownames(df) <- paste(rownames(df), " (",dea_2_name,")",sep="")
  df <- gridExtra::tableGrob(df, theme=gridExtra::ttheme_minimal())
  
  combined.relaxed$CorrelationVariables <- as.factor(combined.relaxed$CorrelationVariables)
  combined.strict$CorrelationVariables <- as.factor(combined.strict$CorrelationVariables)
  
  named.color.vector <- c(both.up.color, both.dn.color, up.dn.color, dn.up.color,
                          ns.color, ns.color, ns.color, ns.color, ns.color)
  names(named.color.vector) <- c("up.up", "dn.dn", "up.dn", "dn.up", "up.ns", "ns.up", "dn.ns", "ns.dn", "ns.ns")
  
  ## Main plot
  plot.relaxed <- generate_correlation_plot(combined.relaxed, dea_1_name=dea_1_name, dea_2_name=dea_2_name, xlab=xlab, ylab=ylab,
                                            pt.size=pt.size, pt.shape=pt.shape, pt.border.color=pt.border.color, named.color.vector=named.color.vector,
                                            selected_l2fc_cut_1=rl2fc_1, selected_l2fc_cut_2=rl2fc_2, line_l2fc_cut_1=l2fc_1, line_l2fc_cut_2=l2fc_2,  linetype=linetype, hist.border.color=hist.border.color, hist.color=hist.color,
                                            density.lwd=density.lwd, density.lt=density.lt, density.lc=density.lc, gene_fields = gene_fields,
                                            title.text=title.text, title.hjust=title.hjust, title.font.size=title.font.size, title.face=title.face, title.color=title.color)
  
  plot.strict <- generate_correlation_plot(combined.strict, dea_1_name=dea_1_name, dea_2_name=dea_2_name, xlab=xlab, ylab=ylab,
                                           pt.size=pt.size, pt.shape=pt.shape, pt.border.color=pt.border.color, named.color.vector=named.color.vector,
                                           selected_l2fc_cut_1=l2fc_1, selected_l2fc_cut_2=l2fc_2, line_l2fc_cut_1=l2fc_1, line_l2fc_cut_2=l2fc_2, linetype=linetype, hist.border.color=hist.border.color, hist.color=hist.color,
                                           density.lwd=density.lwd, density.lt=density.lt, density.lc=density.lc, gene_fields=gene_fields,
                                           title.text=title.text, title.hjust=title.hjust, title.font.size=title.font.size, title.face=title.face, title.color=title.color)
  
  ### Venn plots
  vp.strict <- generate_venn_plot(combined.strict, venn.diagram.colors = venn.diagram.colors, venn.set.name.size = venn.set.name.size, l2fc_1 = l2fc_1, l2fc_2 = l2fc_2, padj_1 = padj_1, padj_2 = padj_2)
  
  vp.relaxed <- generate_venn_plot(combined.relaxed, venn.diagram.colors = venn.diagram.colors, venn.set.name.size = venn.set.name.size, l2fc_1 = l2fc_1, l2fc_2 = l2fc_2, padj_1 = rpadj_1, padj_2 = rpadj_2)
  
  
  results.relaxed <- list(plot.relaxed, combined.relaxed, df.relaxed, vp.relaxed)
  results.strict <- list(plot.strict, combined.strict, df, vp.strict)
  names(results.relaxed) <- c("correlation.plot", "combined.table", "overlap.table", "venn")
  names(results.strict) <- c("correlation.plot", "combined.table", "overlap.table", "venn")
  
  results <- list(results.relaxed, results.strict)
  names(results) <- c("relaxed", "strict")
  
  return(results)
}


generate_correlation_plot <- function(combined, dea_1_name, dea_2_name, xlab, ylab, named.color.vector,
                                      selected_l2fc_cut_1, selected_l2fc_cut_2, line_l2fc_cut_1, line_l2fc_cut_2, linetype, hist.border.color, hist.color,
                                      density.lwd, density.lt, density.lc, pt.size, pt.shape, pt.border.color,
                                      title.text, title.hjust, title.font.size, title.face, title.color, gene_fields){
  ## Main plot
  if(!is.null(ylab)){ ylab = paste("Log2FoldChange (",dea_2_name,")",sep="") }
  if(!is.null(xlab)){ xlab = paste("Log2FoldChange (",dea_1_name,")",sep="") }
  
  
  
  p <- ggplot(combined, aes(x=first.log2FoldChange, y=second.log2FoldChange, fill=CorrelationVariables, label=!!sym(colnames(combined)[grep(gene_fields[1], x=colnames(combined))]))) +
    geom_point(size=pt.size, shape = pt.shape, color = pt.border.color)  +
    geom_text_repel(max.overlaps = 20, text.size = 3, pt.border.color = "lightgrey", max.overlap = 200000,
                    force = 0.25, force_pull = 3, min.segment.length = 0.1, segment.linetype = "dotted",
                    segment.alpha = 1, segment.ncp = 3, segment.color = "white") +  theme_classic() +
    scale_fill_manual(values = named.color.vector) +
    geom_hline(yintercept =  line_l2fc_cut_1, linetype=linetype) +
    geom_hline(yintercept = -line_l2fc_cut_2, linetype=linetype) +
    geom_vline(xintercept =  line_l2fc_cut_1, linetype=linetype) +
    geom_vline(xintercept = -line_l2fc_cut_2, linetype=linetype) +
    xlab(xlab) + ylab(ylab) + theme(legend.position = "none")
  
  ## Histogram and density plot of expression of selected genes for first de analysis.
  top.plot <- ggplot(combined, aes(x=first.log2FoldChange)) +
    geom_histogram(aes(y = ..density..), bins = bins, colour = hist.border.color, fill = hist.color) +
    geom_density(lwd = density.lwd, linetype = density.lt,  colour = density.lc)+
    theme_classic() + theme(axis.title.x = element_blank())
  
  ## Histogram and density plot of expression of selected genes for second de analysis.
  rp <- ggplot(combined, aes(x=second.log2FoldChange)) +
    geom_histogram(aes(y = ..density..), bins = bins, colour = hist.border.color, fill = hist.color) +
    geom_density(lwd = density.lwd, linetype = density.lt,  colour = density.lc)+
    theme_classic()  + theme(axis.title.y = element_blank()) + coord_flip()
  
  vp <-  ggplot()+theme_void()
  
  ## Combined plot
  cpl <- ggplotify::as.ggplot(cowplot::plot_grid(
    cowplot::plot_grid(top.plot, vp, rel_widths = c(12,2)),
    cowplot::plot_grid(p , rp, rel_widths = c(12,2)), rel_heights = c(2,8), ncol=1)) +
    ggtitle(title.text) + theme(plot.title = element_text( face=title.face,
                                                           size = title.font.size,
                                                           hjust=title.hjust,
                                                           colour = title.color))
  
  return(cpl)
  
}

generate_venn_plot <- function(combined, venn.diagram.colors, venn.set.name.size, l2fc_1, l2fc_2, padj_1, padj_2){
  
  print(combined)
  upregulated_genes_first    <- combined[combined$first.log2FoldChange  >=  l2fc_1 & combined$first.padj  <= padj_1,]$ensembl_gene_id
  upregulated_genes_second   <- combined[combined$second.log2FoldChange >=  l2fc_2 & combined$second.padj <= padj_2,]$ensembl_gene_id
  downregulated_genes_first  <- combined[combined$first.log2FoldChange  <= -l2fc_1 & combined$first.padj  <= padj_1,]$ensembl_gene_id
  downregulated_genes_second <- combined[combined$second.log2FoldChange <= -l2fc_2 & combined$second.padj <= padj_2,]$ensembl_gene_id
  
  venn_up=list(first.Upregulated=upregulated_genes_first, second.Upregulated=upregulated_genes_second)
  venn_up_corr <- ggvenn(venn_up, fill_color = venn.diagram.colors, stroke_size = 0.5, set_name_size = venn.set.name.size)
  venn_dn=list(first.Downregulated=downregulated_genes_first, second.Downregulated=downregulated_genes_second)
  venn_dn_corr <- ggvenn(venn_dn, fill_color = venn.diagram.colors, stroke_size = 0.5, set_name_size = venn.set.name.size)
  
  venn_up.dn=list(first.Upregulated=upregulated_genes_first, second.Downregulated=downregulated_genes_second)
  venn_up.dn_corr= ggvenn(venn_up.dn, fill_color = venn.diagram.colors, stroke_size = 0.5, set_name_size = venn.set.name.size)
  venn_dn.up=list(first.Downregulated=downregulated_genes_first, second.Upregulated=upregulated_genes_second)
  venn_dn.up_corr= ggvenn(venn_dn.up, fill_color = venn.diagram.colors, stroke_size = 0.5, set_name_size = venn.set.name.size)
  
  venn <- list(venn_up_corr, venn_dn_corr, venn_up.dn_corr, venn_dn.up_corr)
  names(venn) <- c("up", "dn", "up.dn", "dn.up")
  
  return(venn)
}
