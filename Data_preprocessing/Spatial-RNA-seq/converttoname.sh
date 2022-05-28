tsv_E=3t3_stdata.tsv
path_to_annotation_file=/Dropseq_Alignment_References/mm10/mm10.gtf

convertEnsemblToNames.py $tsv_E --annotation $path_to_annotation_file --output 3t3_exp_matrix.tsv