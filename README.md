# ATAC-Seq analysis related to the paper: "Lyn in fibroblasts promotes leukemia"

## Introduction
This is a repository of the code related to the analysis of the ATAC-Seq dataset from the publication: 
"Alexander F. vom Stein, Rocio Rebollido-Rios, Anna Lukas, Maximilian Koch, Anton von Lom, Sebastian Reinartz, Daniel Bachurski, France Rose, Katarzyna Bozek, Ali T. Abdallah, Viktoria Kohlhas, Julia Saggau, Rebekka Zölzer, Yue Zhao, Christiane Bruns, Paul J. Bröckelmann, Philipp Lohneis, Reinhard Büttner, Björn Häupl, Thomas Oellerich, Phuong-Hien Nguyen, and Michael Hallek. **LYN kinase programs stromal fibroblasts to facilitate leukemic survival via regulation of c-JUN and THBS1**. _Nature Communications_. 2023".

Reproducibility of data analysis in scientific research is important because it helps ensure the validity and reliability of the results, as well as promoting transparency and trust in the scientific community. When others are able to reproduce the same findings using the same methods and data, it provides evidence that the results are robust and not due to chance or error. This also allows for further investigation and improvement of the study, leading to the advancement of knowledge and the discovery of new insights. Therefore, this README is a description of data used in the analysis and the method used to analyze the data as well as the steps needed to reproduce the same results.

The outline of the description is as follows:
- Introduction
- Description of used data
- Method description
- Running the analysis
- References
- Citation

## Description of used data

### Raw data
6 ATAC-Seq libraries generated from samples extracted from HS-5 LYNWT and LYNKO cells cultured under normal growth conditions:
- **KO_R1:** LYN KO cells from human HS-5 cell line.
- **KO_R2:** LYN KO cells from human HS-5 cell line.
- **KO_R3:** LYN KO cells from human HS-5 cell line.
- **WT_R1:** LYN WT cells from human HS-5 cell line.
- **WT_R2:** LYN WT cells from human HS-5 cell line.
- **WT_R3:** LYN WT cells from human HS-5 cell line.

### Processed data


## Method description

**ATAC-sequencing Experimental procedures:** HS-5 LYNWT and LYNKO cells were cultured under normal
growth conditions. Cells were used to perform tagmentation and preparations for
ATAC-sequencing according to the manufacturer’s instructions of the used ATAC-kit
(Active Motif). DNA was isolated respectively and processed and sequenced with 50M
reads, PE100 per sample in the Cologne Center for Genomics following in-house
standards.

**Primary analysis:** The ATAC-seq dataset was analysed on the CHEOPS HPC cluster
of the University of Cologne using nf-core’s (v1.2.1)97 atac-seq pipeline 
(nf-core/atacseq pipeline) in a Singularity environment and the corresponding workflow
management software Nextflow (v21.04.1). Peaks were called using MACS2
(v2.2.7.1) with standard parameters (narrowPeak mode). Consensus peaks
generated by the nf-core/atac-seq pipeline were used for downstream analysis. Counts 
file of consensus peaks was generated using featureCounts (v2.0.1) with standard
parameters using nf-core/atac-seq.

**Differential accessibility analysis:** Peak annotation was performed using the
annotatePeak function of ChIPseeker (v1.30.3) with overlap = ”all”, Homo Sapiens
ensembl database version 104, and otherwise standard parameters. Differential
accessibility analysis was performed using DESeq2 (v1.34.0) with standard
parameters. Normalization of counts for visualization was performed using DESeq2
built-in variance stabilization transformation (vst method) with standard parameters.
Intervals with log2-FC >= 0.56 (Fold Change >= 1.5) and adjusted p-value <= 0.05 were
considered as significant for downstream analyses.

The **motif enrichment** was performed using the motif analysis tool of the python
package Regulatory Analysis Toolbox (RGT) (v0.13.2). In short, the tool implements a
simple Fischer’s exact test for each known human transcription factor from known
databases like Hocomoco104 and JASPAR105. Input region provided are all intervals
generated by ATAC-Seq primary analysis with log2-FC >= 0.56 (Fold Change >= 1.5)
and adjusted p-value <= 0.05 as computed by the differential accessibility analysis.
Background regions used are the universe of all called peaks. The fold change is
computed as ratio of foreground and background relative frequency, where a relative
frequency is defined as number of peaks with a motif match relative to number of peaks
without motif match. Motifs with 15% positive fold change and adjusted p-value <= 0.05
were considered as top significantly enriched motifs in LYNKO versus LYNWT.

The **footprinting analysis** was performed using the footprinting tool from RGT (v0.13.2)
in the atac-seq mode (See HINT-ATAC106). As input, we have used the sorted and
filtered bam files of merged replicates and corresponding narrowPeak files generated
by the nf-core/atac-seq pipeline (internally using macs2). We followed the standard
footprinting procedure to generate tracks as well as lists of motif-predicted binding sites
(mpbs) in bed format for both LYNKO and LYNWT conditions. Finally, we ran the
differential footprinting command of RGT to generate differential activity scores and
corresponding adjusted p-values for each motif. 

### QC


### Peak Annotation


### Differential accessbility analysis


### Motif enrichment analysis


### Footprinting analysis


## Running the analysis

### Installing the environment.
The following was tested on Ubuntu 20.04. We will need an R version 4.1.2, Bioconductor package 3.16 and RStudio 
1. Install conda (For more information see: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html): 
    
    `bash Miniconda3-py38_22.11.1-1-Linux-x86_64.sh`
    
3. An updated R version (>=4.0) is available through the conda-forge channel, so first add the channel using

    `conda config --add channels conda-forge`

3. We will install R through conda-forge and not through the default channel, so we need to set its priority over the default channel.

     `conda config --set channel_priority strict` 
   NOTE: You can undo this change by setting strict priority to the default channel as described here.

4. Check whether an updated R version is added in the conda search space or not.

    `conda search r-base`

5. Now, it is always a good practice (recommoned here) to create a new conda environment, which will help to debug the package-specific and compatibility issues without disrupting the base environment.

    `conda create -n atac_lyn_in_fibroblasts_promotes_leukemia python=3.8`

6. Let's activate the newly create conda environment.

    `conda activate atac_lyn_in_fibroblasts_promotes_leukemia`

7. And finally install the R package.

    `conda install -c conda-forge r-base=4.1.2`

### Installing the github repository


### Running the analysis



## Cite the paper
