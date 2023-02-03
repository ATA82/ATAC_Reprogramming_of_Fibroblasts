### Configuration

## Peaks
distance.to.tss <- Inf  
db.query.pattern <- c("EnsDb", "sapiens", "104")
mount.point <- "~/Desktop/CHEOPS"
server.dir <- ""
project.dir <- "AG_Hallek/xxxxxx_vomstein_atacseq/results/bwa/mergedLibrary/macs/narrowPeak/consensus/"
peaks.file.name <-  "consensus_peaks.mLb.clN.boolean.annotatePeaks.txt"
peaks.dir <- paste(mount.point, "/", server.dir, "/", project.dir, sep="")
peaks.file.path <- paste(peaks.dir, "/", peaks.file.name, sep="")

## Counts
peaks.counts.file.name <- "consensus_peaks.mLb.clN.featureCounts.txt"
peaks.counts.file.path <- paste(peaks.dir, "/", peaks.counts.file.name, sep="")

### Volcano Plots
#vp.point.colors <- list("darkgrey", "dodgerblue3", "dodgerblue4", "coral2", "coral3", "lightgrey")
vp.point.colors <- list("#A80000", "lightgrey", "dodgerblue3", "dodgerblue4", "lightgrey", "black", "lightgrey")

names(vp.point.colors) <- c("sel", "nonsign", "lownegsign", "highnegsign", "lowpossign", "highpossign", "highnonsign")
vp.selected.genes <- c("THBS1", "JUN", "FOS", "STAT1", "NFKB1")

### Other OMICS
transcriptomics_dea <- "results/transcriptomics/transcriptomics.csv"
proteomics_dea <- "results/proteomics/proteomics.csv"

### pdf plots
vpm_p <- "volcano.plot.motifs.pdf"
pca_p <- "pca.and.simhm.pdf"
vpi_p <- "volcano.plot.intervall.pdf"
vpt_p <- "volcano.plot.transcript.pdf"
tp_p <- "towerplot_top25_on_both_sides.pdf"
cp_all_p <- "correlation.plot.pdf"
cp_das_p <- "correlation.plot.DA.significant.pdf"
cp_des_p <- "correlation.plot.DE.significant.pdf"
cp_bos_p <- "correlation.plot.both.significant.pdf"
fap_p <- "footprinting.activity.plot.pdf"
fvp_p <- "footprinting.plot.volcano.pdf"
