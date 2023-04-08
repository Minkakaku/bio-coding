#!/bin/bash

# check for parameter:

if [[ -z $1 || -z $2 ]]; then

printf "This script needs two parameters:\n";
printf "  1. Name of basic run parameter file\n";
printf "  2. Name of run-specific info file\n";

exit 1;

fi

# Check that we can read this file

if [[ -f $1 && -f $2 ]]; then

# pick up definitions for the run from parameter files
. "$1";
. "$2";

# Do some sanity checks for necessary values

if [[ -z "${read1_fastq}" ]]; then

printf "\nError: read1_fastq file not defined in run parameter file '$2'\n";

exit 1;

fi

# check that mapping output dir exists, create if not

if [[ -z "${mapping_output_dir}" && "${mapping_output_dir}" != "." ]]; then

mkdir -p "${mapping_output_dir}""TE/";
mkdir -p "${mapping_output_dir}""GEN/";

fi

# generate a file name header for STAR output

if [[ "${read1_fastq}" == *".gz" ]]; then

r1_fq_base=$(basename "${read1_fastq}" ".fq.gz");

else

r1_fq_base=$(basename "${read1_fastq}" ".fq");

fi

# generate output file name

if [[ -n "${mapping_output_dir}" ]]; then

star_te_out_prefix="${mapping_output_dir}""TE/""${r1_fq_base}""_";
star_gen_out_prefix="${mapping_output_dir}""GEN/""${r1_fq_base}""_";
else

star_te_out_prefix="TE/""${r1_fq_base}""_";
star_gen_out_prefix="GEN/""${r1_fq_base}""_"
fi

star_out_gene_name="${star_gen_out_prefix}"".gene.""Aligned.sortedByCoord.out.bam";
star_out_te_name="${star_te_out_prefix}"".te.""Aligned.sortedByCoord.out.bam";
#star_out_bam_name="${star_out_prefix}""Aligned.out.bam";

# check if STAR has already been run, to avoid repetition
if [[ -f "${star_out_gene_name}" ]]; then

printf "STAR output file '${star_out_gene_name}'already exists, using this file\n";

else
# # run STAR on reads:

STAR --runThreadN "${starthreads}" \
--genomeDir "${star_genome_dir}" \
--outSAMtype BAM SortedByCoordinate \
--sjdbGTFfile "${ucsc_gene_gtf}" \
--outFileNamePrefix "${star_gen_out_prefix}"".gene." \
--readFilesIn "${read1_fastq}" "${read2_fastq}" \
--readFilesCommand gunzip -c;

fi

if [[ -f "${star_out_te_name}" ]]; then

printf "STAR output file '${star_out_te_name}'already exists, using this file\n";

else

STAR --runThreadN "${starthreads}" \
--genomeDir "${star_genome_dir}" \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix "${star_te_out_prefix}"".te." \
--readFilesIn "${read1_fastq}" "${read2_fastq}" \
--outFilterMultimapNmax "${starfiltermax}" \
--winAnchorMultimapNmax "${staranchormax}" \
--readFilesCommand gunzip -c;

fi
# run FeatureCounts on results, Multicount, then Unique count

if [[ -n ${verbose} ]]; then
printf "Running ${path_to_featurecounts}featureCounts on\n  '${star_out_te_name}', multi-run\n"
fi

featureCounts -a "${ucsc_repeat_gtf}" \
  -o "${star_te_out_prefix}""M_FC.txt" -t "exon" -f -p -M -O --fraction \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  ${featurecounts_overlap_partial} \
  -T "${featurecounts_threads}" "${star_out_te_name}"

if [[ -n ${verbose} ]]; then
printf "Running ${path_to_featurecounts}featureCounts on\n  '${star_out_te_name}', unique-run\n"
fi

featureCounts -a "${ucsc_repeat_gtf}" \
  -o "${star_te_out_prefix}""U_FC.txt" -t "exon" -f -p -O \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  ${featurecounts_overlap_partial} \
  -T "${featurecounts_threads}" "${star_out_te_name}"


featureCounts -a "${ucsc_repeat_gtf}" \
  -o "${star_te_out_prefix}""meta.txt" -t "exon" -p -O \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  ${featurecounts_overlap_partial} \
  -T "${featurecounts_threads}" "${star_out_te_name}"
# GENECODE runs

if [[ -n ${verbose} ]]; then

printf "Running ${path_to_featurecounts}featureCounts on\n  '${star_out_gene_name}'\n"

fi

featureCounts -a "${ucsc_gene_gtf}" \
  -o "${star_gen_out_prefix}""gene_counts.txt" -t "exon" -p -O \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  ${featurecounts_overlap_partial} \
  -T "${featurecounts_threads}" "${star_out_gene_name}"

# put unique RE element ID at end of line,
# then run identgeneloc on this, to locate proximal genes

cat << 'SCRIPT2' > reorder_cols.awk
# put col1 to end of line

{for (i = 2; i<= NF; i++)
  printf("%s\t",$i);
printf("%s\n",$1);
}
SCRIPT2


# generate a header line for output
printf "#Chromosome\tstart\tend\tstrand\tlength\thits\tRE_uniq_ID\tDistToGene\tOccupy\tLocationWRTgene\tGeneSense\tGeneCoord\tGeneName\tgenelocation\n" > "${star_te_out_prefix}""U_FC.genloc";

# the track for combining expression values

cat "${star_te_out_prefix}""U_FC.txt" | awk -f reorder_cols.awk | \
"${path_to_dmap}"identgeneloc -T -f "${ucsc_gene_gtf}" -i -C 9 -a "transcript" -A gene_id  -r - | \
cut -f 1-11,14- >> "${star_te_out_prefix}""U_FC.genloc";
# tidy up scripts if required

if [[ -n ${delete_temp_files} ]]; then

rm append_tpm_FC.awk
rm get_tpm_parameter.awk
rm reorder_cols.awk
rm add_gene_names_to_genloc.awk
rm appnd_stringtie_ga.awk

fi

exit 0;

else

exit 1;

fi