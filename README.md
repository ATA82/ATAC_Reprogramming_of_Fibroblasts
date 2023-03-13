# ATAC-Seq data analysis related to the paper: "Lyn in fibroblasts promotes leukemia" (running title)

**Authors:**

-   **Ali T. Abdallah**: Implementation and data analysis.
-   **Alexander F. vom Stein**: Experiment performer.

## Introduction

This is a repository of the code and processed data related to the analysis of the ATAC-Seq dataset from the publication: *"Alexander F. vom Stein, Rocio Rebollido-Rios, Anna Lukas, Maximilian Koch, Anton von Lom, Sebastian Reinartz, Daniel Bachurski, France Rose, Katarzyna Bozek, **Ali T. Abdallah**, Viktoria Kohlhas, Julia Saggau, Rebekka Zölzer, Yue Zhao, Christiane Bruns, Paul J. Bröckelmann, Philipp Lohneis, Reinhard Büttner, Björn Häupl, Thomas Oellerich, Phuong-Hien Nguyen, and Michael Hallek. **LYN kinase programs stromal fibroblasts to facilitate leukemic survival via regulation of c-JUN and THBS1**. Nature Communications*. 2023".

### About reproducibility

*Reproducibility* of data analysis is an important quality metric of scientific research because it helps ensure the validity and reliability of the results, as well as promoting transparency and trust in the scientific community. When others are able to reproduce the same findings using the same methods and data, it provides evidence that the results are robust and not due to chance or error. This also allows for further investigation and improvement of the study, leading to the advancement of knowledge and the discovery of new insights (See the wiki article on [reproducibility](https://en.wikipedia.org/wiki/Reproducibility) for more details).

If you follow the steps described in this README, you should be able to reproduce the results in the Paper.

### Outline

This README is a description of (1) the data used in the analysis, (2) the method used to analyze the data, as well as (3) the steps needed to reproduce the same results.

The following is the outline of the description:

-   [Introduction](#introduction)
-   [Description of used data](#description-of-used-data)
-   [Method description](#method-description)
-   [Running the analysis](#running-the-analysis)
-   [References](#references)
-   [Citation](#citation)

## Description of used data

### Raw data

6 ATAC-Seq libraries generated from samples extracted from HS-5 LYNWT and LYNKO cells cultured under normal growth conditions:

-   **KO_R1:** LYN KO cells from human HS-5 cell line.
-   **KO_R2:** LYN KO cells from human HS-5 cell line.
-   **KO_R3:** LYN KO cells from human HS-5 cell line.
-   **WT_R1:** LYN WT cells from human HS-5 cell line.
-   **WT_R2:** LYN WT cells from human HS-5 cell line.
-   **WT_R3:** LYN WT cells from human HS-5 cell line.

The raw data are stored in the functional genomics data collection (ArrayExpress) with the following accession number [E-MTAB-12531](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12531#) and could be directly downloaded from the [samples and data section](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12531/sdrf)

### Processed data

-   **Consensus Peaks**: These are the consensus peaks called by the nf-core pipeline.
-   **Annotated peaks**: These are the annotated peak using the targets pipeline.
-   **Differential accessibility table**: This is the list of differentially accessible peaks.
-   **Motif table**: This is the list of motifs generated using the RGT toolbox.
-   **Footprinting table**: This is the list of footprints detected by the ATAC-HINT method in the RGT toolbox.

## Reproduce the results

### Getting started

If you do not have docker installed on your host system, install it first see: <https://docs.docker.com/get-docker/>, if you are going to remote into the host running the RStudio server instance in docker, only the remote host and not your ssh client needs to have docker installed.

-   Download raw data (fastq files) from the above-mentioned repository ([E-MTAB-12531](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12531/sdrf) )
-   Clone this repo to your desired directory: `git clone git@github.com:ATA82/ATAC_Reprogramming_of_Fibroblasts.git`
-   Run the nf-core pipeline using the [start_nfcore_atacseq.sh](https://github.com/ATA82/ATAC_Reprogramming_of_Fibroblasts/blob/main/basic_analysis/start_nfcore_atacseq.sh) script on our raw data, using the fastq files downloaded in the previous step, the nf.conf file, and the samplesheet.csv file (these files could be found in the [here](https://github.com/ATA82/ATAC_Reprogramming_of_Fibroblasts/tree/main/basic_analysis))
    - The start_nfcore_atacseq.sh, the nf.conf file, the samplesheet file, and the raw data should be in the same folder.
    - If you want to use different folders you should adjust the code or the samplesheet file accordingly.
    - There are comments within the start_nfcore_atacseq.sh script regarding the requirements.
    - **This step is only needed if you want to start from scratch otherwise you could work with the results of the nf-core pipeline.**
-   Create a file in this directory called `.env` which contains `USERNAME=<your username on your host system>`
-   Customise the docker-compose.yml file to set the:
    -   project/container name - set this to a name for the project to make a new container for the project and to make it easier to identify your containers in e.g. `docker ps`
    -   project directory - set this to the path on your host system where you have just cloned this repo It is important that you create this directory prior to starting the container as if the project path does not exist it will be owned by the root user by default which will cause permissions issues
    -   port - setting a unique port for the container makes it easier to have multiple projects up on your host at the same time as they will not clash
    -   map the cache of renv on the host system to the container. If you don't have a cache you could just select a directory on the host system where the user of rstudio have access to, such as the one used in the docker-compose.yml file. The advantage of selecting an existing cache is given in case you already have most of the package versions on your system. There are additional comments in `docker-compose.yml` to aid in setting it up

### Running the analysis

-   Start your container with `docker compose up -d`, `-d` is for detach so the container will now run in the background and will even resume following a reboot. to stop session use `docker compose down`, add the `--build` flag when starting if you want to trigger a re-build of the container image.
-   Connect to your RStudio server session at: localhost:port or 0.0.0.0:port, where port is the port you set in `docker-compose.yml`, in this example it is 4444 The default user and password are 'rstudio' & '1rstudio' respectively, these are set in `rstudio_docker/rstudio.env` note that you can connect to an RStudio server instance running on a remote system to which you have access over ssh by using ssh port forwarding. You can forward a port over ssh with a command of the form: `ssh -nNT -L <local port>:<host>:<remote port> <host>` e.g. `ssh -nNT -L 4444:abdallah.cecad.uni-koeln.de:4444 aabdallah@aabdallah.cecad.uni-koeln.de` you should now be able to access the RStudio server instance at localhost of the ssh client in a web browser. (Tip: If working with multiple RStudio containers we may sometimes have issues logging into different sessions at the same time opening different sessions using Firefox container tabs fixes this issue)
-   Switch to the "project" directory using: setwd("project")
-   Use `renv::restore()` to load the appropriate versions of the R package dependencies from the renv lockfile
-   If your cache does not contain these dependencies, it will be installed via renv.
-   After finishing the restore process, restart the R session.
-   If you are notified about the synchronisation status of the project, run renv::status().
-   You will be notified to run either renv::snapshot() or renv::restore() to get the project synchronized. The first one works.
-   load the {targets} R package `library(targets)` and run the targets pipeline with `tar_make()`
-   After all target objects are created you could switch to R Notebook within the project and Run all chunks.

Customisations can be made within the `_targets.R` file and the R functions files located in `R/`. In case you have used our code and analysis to get some insights, I would likely to kindly ask for a citation of the related paper. The full citation is mentioned below.

## Method description

**ATAC-sequencing Experimental procedures:** HS-5 LYNWT and LYNKO cells were cultured under normal growth conditions. Cells were used to perform tagmentation and preparations for ATAC-sequencing according to the manufacturer’s instructions of the used ATAC-kit (Active Motif). DNA was isolated respectively and processed and sequenced with 50M reads, PE100 per sample in the Cologne Center for Genomics following in-house standards.

**Primary analysis:** The ATAC-seq dataset was analysed on the CHEOPS HPC cluster of the University of Cologne using nf-core’s (v1.2.1)97 atac-seq pipeline (nf-core/atacseq pipeline) in a Singularity environment and the corresponding workflow management software Nextflow (v21.04.1). Peaks were called using MACS2 (v2.2.7.1) with standard parameters (narrowPeak mode). Consensus peaks generated by the nf-core/atac-seq pipeline were used for downstream analysis. Counts file of consensus peaks was generated using featureCounts (v2.0.1) with standard parameters using nf-core/atac-seq.

**Differential accessibility analysis:** Peak annotation was performed using the annotatePeak function of ChIPseeker (v1.30.3) with overlap = ”all”, Homo Sapiens ensembl database version 104, and otherwise standard parameters. Differential accessibility analysis was performed using DESeq2 (v1.34.0) with standard parameters. Normalization of counts for visualization was performed using DESeq2 built-in variance stabilization transformation (vst method) with standard parameters. Intervals with log2-FC \>= 0.56 (Fold Change \>= 1.5) and adjusted p-value \<= 0.05 were considered as significant for downstream analyses.

The **motif enrichment** was performed using the motif analysis tool of the python package Regulatory Analysis Toolbox (RGT) (v0.13.2). In short, the tool implements a simple Fischer’s exact test for each known human transcription factor from known databases like Hocomoco104 and JASPAR105. Input region provided are all intervals generated by ATAC-Seq primary analysis with log2-FC \>= 0.56 (Fold Change \>= 1.5) and adjusted p-value \<= 0.05 as computed by the differential accessibility analysis. Background regions used are the universe of all called peaks. The fold change is computed as the ratio of foreground and background relative frequency, where a relative frequency is defined as a number of peaks with a motif match relative to a number of peaks without a motif match. Motifs with 15% positive fold change and adjusted p-value \<= 0.05 were considered as top significantly enriched motifs in LYNKO versus LYNWT.

The **footprinting analysis** was performed using the footprinting tool from RGT (v0.13.2) in the atac-seq mode (See HINT-ATAC106). As input, we have used the sorted and filtered bam files of merged replicates and corresponding narrowPeak files generated by the nf-core/atac-seq pipeline (internally using macs2). We followed the standard footprinting procedure to generate tracks as well as lists of motif-predicted binding sites (mpbs) in bed format for both LYNKO and LYNWT conditions. Finally, we ran the differential footprinting command of RGT to generate differential activity scores and corresponding adjusted p-values for each motif.

## References

1.  Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PMID: 32055031.
2.  Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
3.  Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PMID: 28398311.
4.  Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PubMed PMID: 18798982; PubMed Central PMCID: PMC2592715.
5.  Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PubMed PMID: 24227677.
6.  Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics. 2015 Jul 15;31(14):2382-3. doi: 10.1093/bioinformatics/btv145. Epub 2015 Mar 11. PMID: 25765347.
7.  Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. PubMed PMID: 25516281; PubMed Central PMCID: PMC4302049.
8.  Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H. gprofiler2 -- an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler. F1000Res. 2020 Jul 15;9:ELIXIR-709. doi: 10.12688/f1000research.24956.2. PMID: 33564394; PMCID: PMC7859841.
9.  Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015 Dec 23;1(6):417-425. doi: 10.1016/j.cels.2015.12.004. PMID: 26771021; PMCID: PMC4707969.
10. Kulakovskiy IV, Medvedeva YA, Schaefer U, Kasianov AS, Vorontsov IE, Bajic VB, Makeev VJ. HOCOMOCO: a comprehensive collection of human transcription factor binding sites models. Nucleic Acids Res. 2013 Jan;41(Database issue):D195-202. doi: 10.1093/nar/gks1089. Epub 2012 Nov 21. PMID: 23175603; PMCID: PMC3531053.
11. Fornes O, Castro-Mondragon JA, Khan A, van der Lee R, Zhang X, Richmond PA, Modi BP, Correard S, Gheorghe M, Baranašić D, Santana-Garcia W, Tan G, Chèneby J, Ballester B, Parcy F, Sandelin A, Lenhard B, Wasserman WW, Mathelier A. JASPAR 2020: update of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2020 Jan 8;48(D1):D87-D92. doi: 10.1093/nar/gkz1001. PMID: 31701148; PMCID: PMC7145627.
12. Li Z, Schulz MH, Look T, Begemann M, Zenke M, Costa IG. Identification of transcription factor binding sites using ATAC-seq. Genome Biol. 2019 Feb 26;20(1):45. doi: 10.1186/s13059-019-1642-2. PMID: 30808370; PMCID: PMC6391789.


## Citation

Vom Stein AF, Rebollido-Rios R, Lukas A, Koch M, von Lom A, Reinartz S, Bachurski D, Rose F, Bozek K, **Abdallah AT**, Kohlhas V, Saggau J, Zölzer R, Zhao Y, Bruns C, Bröckelmann PJ, Lohneis P, Büttner R, Häupl B, Oellerich T, Nguyen PH, Hallek M. LYN kinase programs stromal fibroblasts to facilitate leukemic survival via regulation of c-JUN and THBS1. Nat Commun. 2023 Mar 10;14(1):1330. doi: 10.1038/s41467-023-36824-2. PMID: 36899005.
