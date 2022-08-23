matching <- function(symbols, selected.genes, dimers=T, dim.sep="::"){
  matches <- vector("list", length(symbols))
  for(s in 1:length(symbols)){
    symbols_list <- unlist(str_split(symbols[s], pattern="::"))
    for(sg in selected.genes){
      if(sg %in% symbols_list){
        matches[[s]] <- T
        break
      }
    }
    if(is.null(matches[[s]])){
      matches[[s]]=F
    }
  }
  return(unlist(matches))
}

footprinting_plot <- function(activity_df, padj=0.05, z.score = 0, 
                              conditions=c("WT", "KO"), activity="Z_score", 
                              second.variable=NULL, log10.second.variable=F, 
                              coord_flip=F, selected.genes=NULL,
                              pt.border.color="black", pt.sizes=c(2,1,2), 
                              pt.colors=c("red", "grey", "orange", "steelblue"), 
                              xlab="Points randomly (normally) distributed", 
                              ylab="KO <-- Activity z-Score --> WT", 
                              title.color="#800000", title.hjust=0.5, 
                              title.face="bold", source="hint", ylim=c(-5,5), 
                              title="Footprinting analysis - Overview", 
                              repel.max.overlaps=20000, repel.text.size=3, 
                              repel.font.face="bold", repel.text.color="black",
                              repel="text", sign.hline="both", sign.hcut=1, 
                              sign.vline="none", sign.vcut=1){
  l <- -log10(padj)
  if(is.null(selected.genes)){
    selected.genes <- ""
  }
  if(source=="hint"){
    activity_df <- activity_df %>% separate(Motif, c("MotifID", "MotifGenes"), "\\.(?=[A-Z])")
    activity_df <- activity_df %>% separate(MotifGenes, c("Symbol"), "\\(", remove = F)
    activity_df$Enrichment <- 
      ifelse(activity_df$Z_score > z.score & activity_df$P_values < padj, 
             conditions[1],
             ifelse(activity_df$Z_score < z.score & activity_df$P_values < padj,
             conditions[2], 
             ifelse(matching(symbol = activity_df$Symbol, selected.genes = selected.genes),
             "NS.ROI",
             "NS")))
    activity_df$label <-
      ifelse((activity_df$P_values < padj | 
                matching(symbol = activity_df$Symbol, 
                         selected.genes = selected.genes)) , 
             activity_df$MotifGenes, "")
    activity_df$extended.label <- 
      ifelse(activity_df$label != "", 
             paste(activity_df$MotifGenes," (", activity_df$MotifID,")", sep=""), 
             "")
    
    if(is.null(second.variable) & is.null(activity_df$random.x)){
      set.seed(123)
      activity_df$random.x <- sample(1:10000, nrow(activity_df), replace=TRUE)
    } else {
      activity_df$random.x <- activity_df[,(second.variable)]
    }
    if(log10.second.variable){
      activity_df$random.x <- -log10(activity_df$random.x)
    }
    fp <- ggplot(activity_df, 
                 aes(x=random.x, y=(get(activity)), label=label, fill=Enrichment)) + 
      theme_pubr() + 
      geom_point(color=pt.border.color, aes(size=Enrichment), pch=21) + 
      scale_fill_manual(values=pt.colors) 
    
    if(repel=="text"){
      fp <- fp+ geom_text_repel(size=repel.text.size,
                                max.overlaps = repel.max.overlaps,
                                color=repel.text.color,
                                fontface=repel.font.face)
    } else {
      fp <- fp + geom_label_repel(size=repel.text.size,
                                  max.overlaps = repel.max.overlaps,
                                  alpha=0.5,
                                  color=repel.text.color,
                                  fontface=repel.font.face) 
    }
    fp <- fp + scale_size_manual(values=pt.sizes) + xlab(xlab) + ylab(ylab) + 
      ylim(ylim[1],ylim[2]) + ggtitle(title) + 
      theme(plot.title=element_text(hjust=title.hjust, 
                                    color=title.color, 
                                    face=title.face))
    
    if(is.null(second.variable)){
      fp <- fp + theme(axis.text.x=element_blank(), 
                       axis.ticks.x=element_blank())
    }
    if(coord_flip){
      fp <- fp + coord_flip()
    }
    if(sign.hline=="both"){
      fp <- fp + geom_hline(yintercept = sign.hcut, linetype="dotted")
      fp <- fp + geom_hline(yintercept = -sign.hcut, linetype="dotted")
    } else if(sign.hline=="one"){
      fp <- fp + geom_hline(yintercept = sign.hcut, linetype="dotted")
    }
    if(sign.vline=="both"){
      fp <- fp + geom_vline(xintercept = sign.vcut, linetype="dotted")
      fp <- fp + geom_vline(xintercept = -sign.vcut, linetype="dotted")
    } else if(sign.vline=="one"){
      fp <- fp + geom_vline(xintercept = sign.vcut, linetype="dotted")
    }
    return(list(fp, activity_df))
  } else {
    stop(paste("Provided footprinting analysis software '",source,"' is not supported!"))
  }
}