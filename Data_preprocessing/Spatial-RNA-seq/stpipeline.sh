#!/bin/bash

# FASTQ reads
FW=PATH_TO_PROCESSED_R2/sample_R2_processed.fastq.gz
RV=PATH_TO_R1/R1.fastq.gz

# References for mapping, annotation and nonRNA-filtering
MAP=PATH_TO_ALIGNMENT_REF/Dropseq_Alignment_References/mm10/
ANN=PATH_TO_ALIGNMENT_REF_GTF/Dropseq_Alignment_References/mm10/mm10.gtf
#CONT=/mouse/GRCm38_86/ncRNA/StarIndex/

# Barcodes settings
ID=PATH_TO_BARCODE_INDEX/spatial_barcodes_index.txt 

# Output folder and experiment name
OUTPUT=PATH_TO_OUTPUT/st_pipeline_new/
mkdir -p PATH_TO_OUTPUT/st_pipeline_new/

TMP=PATH_TO_TEMP/st_pipeline_new/tmp
mkdir -p PATH_TO_TEMP/st_pipeline_new/tmp

# Do not add / or \ to the experiment name
EXP=P21

# Running the pipeline
st_pipeline_run.py \
  --output-folder $OUTPUT \
  --ids $ID \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $EXP \
  --htseq-no-ambiguous \
  --verbose \
  --log-file $OUTPUT/${EXP}_log.txt \
  --allowed-kmer 5 \
  --mapping-threads 20 \
  --temp-folder $TMP \
  --no-clean-up \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --overhang 0 \
  --min-length-qual-trimming 10 \
  $FW $RV