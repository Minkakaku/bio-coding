#!/bin/bash
# get params from path file
. "$1";

# check if ref_te is gz

ucsc_repeat_gtf = ''
ucsc_repeat=$(basename "${ucsc_repeat_gtf}" ".gtf")

ucsc_gene_gtf = ''
ucsc_gene=$(basename "${ucsc_gene_gtf}" ".gtf")

if [[ "${ucsc_repeat_gtf}" == *".gz" ]]; then

  gzip -dc "${ucsc_repeat_gtf}" > "${ucsc_gene}"".gtf";

fi

# check for uncompressed ucsc gene gtf file - if not generate it.

if [[ ! -f "${ucsc_gene_gtf}" ]]; then

  gzip -dc "${ucsc_gene_gtf_src}" > "${repeat_gene_gtf_dir}""${ucsc_gene_gtf}";

fi

# append modified repeat gtf to ucsc genomic gtf

cat "${ucsc_gene_gtf}" "${ucsc_repeats_uniq_gtf}" > "${gene_repeat_gtf}";

if [[ -n ${verbose} ]]; then

printf "Files '${ucsc_gene_gtf}' and\n";
printf " '${ucsc_repeat_src}' combined to\n";
printf " '${gene_repeat_gtf}'\n";

fi

# check if we really need to build the STAR index

if [[ -f "${star_genome_dir}""/Genome" && -f "${star_genome_dir}""/geneInfo.tab" ]]; then

if [[ -n ${verbose} ]]; then

printf "STAR index found in '${star_genome_dir}' - using this index\n";

fi

else

# generate STAR index

# first check if genome file is .gz compressed, if so, we must uncompress it

if [[ -n ${verbose} ]]; then

printf "Beginning STAR run on '${genome_fasta_file}'\n";
printf "with '${gene_repeat_gtf}'\n";

fi

  if [[ "${genome_fasta_file}" == *".gz" ]]; then

  uncomp_genome_fasta=$(basename "${genome_fasta_file}" ".gz");
  gzip -dc ${genome_fasta_file} > ${uncomp_genome_fasta};

  "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${uncomp_genome_fasta}" --sjdbGTFfile "${gene_repeat_gtf}";

  else

  mkdir -p "${star_genome_dir}";

  "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${genome_fasta_file}" --sjdbGTFfile "${gene_repeat_gtf}";

  fi
fi

# create Ensembl ID to gene name list file

cat << 'ENS_SCRIPT' > ensembl_ID_to_gname.awk
# scan the ucsc_gene_gtf file, extracting Ensemble IDs
# and gene names from the attributes, write a tab separated
# list of these

BEGIN{FS="\t";}

$1!~/#/&&$3=="gene"{ns = split($9,s9," ");
for (i = 1; i < ns; i+=2)
  {
  if (index(s9[i],"gene_id") > 0)
    ensid = substr(s9[i+1],2,length(s9[i+1])-3);
  if (index(s9[i],"gene_name") > 0)
    {
    printf("%s\t%s\n",ensid,substr(s9[i+1],2,length(s9[i+1])-3));
    break;
    }
  }
}
ENS_SCRIPT

if [[ -n ${verbose} ]]; then

printf "Creating Ensembl ID vs gene name file '${ensid_vs_gname}'\n";

fi

awk -f ensembl_ID_to_gname.awk "${ucsc_gene_gtf}" > "${ensid_vs_gname}";

printf "genome setup complete\n";

if [[ -n delete_temp_files ]]; then

echo "Deleting scripts";

rm make_uniq_gtf_ex_allfields.awk ensembl_ID_to_gname.awk

fi

if [[ -n ${verbose} ]]; then

printf "\nRepExpress build completed\n";

fi

exit 0;

else

printf "Can't open parameter file '$1'\n";

exit 1;

fi
