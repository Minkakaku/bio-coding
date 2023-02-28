# DESeq2_params.sh: to define RepExpress env variables
# for running DESeq2 expression comparison on featureCounts
# unique counts files, or TPM files from that, generated
# by map_count_reads.sh runs.
#
# Variables that are consistent across a series of runs
# are already defined in basic_params.sh
# 
# Edit this to reflect the parameters required for each DESeq2 run.

# the output directory for DESeq2: leave blank for
# output into the current directory.

deseq2_output_dir="../05.deseq2_out/";

# We need to define the count or tpm files by having a list of each set in
# a text file, one per line.

deseq2_file_list_1="../05.deseq2_out/control_tpm_list.txt";

sample_name_1="control";

deseq2_file_list_2="../05.deseq2_out/treatment_tpm_list.txt";

sample_name_2="treatment";

# the sample_name values can be left blank and will default
# to "control" and "treatment"

# name for the final output file: can be left empty for a fairly awkward default

DESeq2_output_file="../05.deseq2_out/EWAT_deseq2_express.txt";

delete_temp_files="";
