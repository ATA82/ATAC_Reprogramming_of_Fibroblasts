#!/bin/bash

export PATH="/path/to/nexflow/bin/:$PATH" # add nextflow to the path
export NXF_WORK="/scratch/$USER/nf_work" # NB these may not exist by default and may be reset when scratch's data are deleted.
export NXF_TEMP="/scratch/$USER/nf_tmp"
export NXF_EXECUTOR=slurm # select appropriate executor
export NXF_SINGULARITY_CACHEDIR="/path/to/caches/nextflow_singularity_cache" # Here you should set the path to the singularity cache.

#"Nextflow Options" ##########################################################################################################################################
# -bg: run nextflow in background.                                                                                                                          #
# -r 3.0: use version 3.0                                                                                                                                   #
# -profile singularity: run pipeline in singularity                                                                                                         #
#############################################################################################################################################################

#"Used nf-core RNA-Seq Pipeline Options" #####################################################################################################################
# --input: a csv listing all samples, replication information, file names and strandedness information.                                                     #
# --outdir: directory where to store results.                                                                                                               #
# --multiqc_title: the title use in the multiqc report.                                                                                                     #
# --fasta: the URL/path of the genome sequence file (fasta format)                                                                                          #
# --gff: the URL/path of the annotation of the genome sequence file (ggf3 or gtf)                                                                           #
# --star_index: the location of the pre-built star index (usefule in case of resource bottlenecks or if using old versions or the like is explicity needed  #
# --gtf_extra_attributes: define which attributes should be used as additional identifier (default: gene_name)                                              #
# --gtf_group_features: define the name of the feature in the gtf/gff file to group upon in salmon (default: gene_id)                                       #
# --gtf_count_type: define the feature used to make the quantification (default: exon) | in newer versions replaced by  "--featurecounts_feature_type"      #
# --save_reference: saves reference sequences.                                                                                                              #
# --seq_center: define sequencing center to store this information in the bam files, which is useful, for tracing back.                                     #
# --deseq2_vst: use variance stabilization transformation of deseq2 rather than the rlog transformation.                                                    #
# --skip_biotype_qc: skips the biotype qc step - here due to inconsistency as biotype field is not present.                                                 #
#############################################################################################################################################################

fasta='http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
gtf='http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz'
fasta='/projects/cecad/bioinformatics-core-facility/public/ensembl/main/release-104/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa'
fasta='/projects/cecad/bioinformatics-core-facility/public/ensembl/main/release-104/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
gtf='/projects/cecad/bioinformatics-core-facility/public/ensembl/main/release-104/homo_sapiens/Homo_sapiens.GRCh38.104.gtf'
blacklist='/projects/cecad/bioinformatics-core-facility/public/blacklists/hg38-blacklist.v2.bed/'
gsize=2913022398
mito_name='MT'
run_name='ATACSeq_211123_MH'
email='your-e-mail@domain.xx'
mqc_title='\\"ATAC-Seq\\ -\\ AG\\ Hallek\\ -\\ Dr.\\ Alexander\\ vom\\ Stein\\"'
seq_center='\\"Cologne\\ Center\\ For\\ Genomics\\"'

nextflow run nf-core/atacseq -bg -r 1.2.1 -c nf.conf -profile singularity --input samplesheet.csv --narrow_peak --fasta $fasta --gtf $gtf --blacklist $blacklist --macs_gsize $gsize --mito_name $mito_name --email $email
