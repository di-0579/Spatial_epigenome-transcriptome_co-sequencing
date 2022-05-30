# Spatial_epigenome-transcriptome_co-sequencing
### Introduction

This repository aims to share the raw data processing and visualization codes used in Spatial_epigenome-transcriptome_co-sequencing project.

![](https://github.com/di-0579/Spatial_epigenome-transcriptome_co-sequencing/blob/main/workflow/workflow.jpg?raw=true)

### Data processing
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq 6000 sequencer (pair-end 150 bp mode).
#### Spatial_ATAC-seq
##### 1.Raw Fastq data
Read 1: contains the spatial Barcode A and Barcode B
Read 2: contains the genome sequences
##### 2. Reformat raw Fastq file to Cell Ranger ATAC format (10x Genomics)
**Raw read 1 -> New Read 1 + New Read 2**
- New Read 1: contains the genome sequences
- New Read 2: contains the spatial Barcode A and Barcode B

**Raw read 2 -> New Read 3**

Reformatting raw data was implemented by BC_process.py in the Data_preprocessing folder.


##### 3. Sequence alignment and generation of fragments file
The reformated data was processed using Cell Ranger ATAC v1.2 with following references:
Mouse reference (mm10):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz

Human reference (GRCh38):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz

**A preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder. To run the pipeline, use the command:**

    sbatch Snakemake.sh

#### Spatial_RNA-seq

**1. Raw Fastq data processing using ST pipeline and generate expression matrix**

The DBiT-seq Raw fastq file
Read 1: Contains the cDNA sequence
Read 2: Contains the spatial Barcode A, Barcode B and UMIs

**2.Reformat Fastq Read 2 file**
To reformat the Raw data, run the fastq_process.py in Rawdata_processing folder and gzip the resulted fastq file to save space:

    python fastq_process.py
    gzip sample_R2_processed.fastq

The reformated data was processed following ST pipeline.  

**3、Run ST pipeline**
Run st_pipeline.sh to start the ST pipeline: The input is processed_R2.fastq.gz and Raw R1.fastq.gz. It also requires a "spatial_barcodes_index.txt" to decode the spatial location information. Genome references and annotatation files were aslo needed.

    #!/bin/bash
    # FASTQ reads
    FW=$tmp/${sample}_R2_processed.fastq
    RV=$tmp/${sample}_raw_qc_R1.fastq.gz
    # References for mapping, annotation and nonRNA-filtering
    MAP=/Dropseq_Alignment_References/mm10/
    ANN=/Dropseq_Alignment_References/mm10/mm10.gtf 
    CONT=/Spatial_omics_references/mouse/GRCm38_86/ncRNA/StarIndex/

    # Barcodes settings
    ID=/useful/spatial_barcodes.txt

    # Output folder and experiment name
    OUTPUT=../output/
    mkdir -p $OUTPUT

    TMP_ST=$OUTPUT/tmp
    mkdir -p $TMP_ST

    # Running the pipeline
    st_pipeline_run.py \
      --output-folder $OUTPUT \
      --ids $ID \
      --ref-map $MAP \
      --ref-annotation $ANN \
      --expName $sample \
      --htseq-no-ambiguous \
      --verbose \
      --log-file $OUTPUT/${sample}_log.txt \
      --demultiplexing-kmer 5 \
      --threads 20 \
      --temp-folder $TMP_ST \
      --no-clean-up \
      --umi-start-position 16 \
      --umi-end-position 26 \
      --demultiplexing-overhang 0 \
      --min-length-qual-trimming 20 \
      $FW $RV

**4.Convert Ensemble to Gene Names**
Then, Run converttoname.sh to annotate the resulting sample_stdata.tsv.
    
    #!/bin/bash
    tsv_E=$OUTPUT/${sample}_stdata.tsv
    path_to_annotation_file=/Dropseq_Alignment_References/mm10/mm10.gtf

    convertEnsemblToNames.py --annotation $path_to_annotation_file --output $OUTPUT/${sample}_stdata_names.tsv $tsv_E

####  Identify useful pixels (pixel on tissue) from microscope image using Matlab
link：
https://github.com/edicliuyang/DBiT-seq_FFPE/tree/master/Figure_Processing



### Data visualization
The data visualization were completed with R language. The package used extensively the functions in Signac V1.6， Seurat V4.0 and ArchR v1.0.2.














