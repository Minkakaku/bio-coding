#!/bin/bash
# get params from path file
. "$1";

# check if we really need to build the STAR index

if [[ -f "${star_genome_dir}""/Genome" && -f "${star_genome_dir}""/geneInfo.tab" ]]; then

  if [[ -n ${verbose} ]]; then

  printf "STAR index found in '${star_genome_dir}' - using this index\n";

  fi

else

# generate STAR index

# first check if genome file is .gz compressed, if so, we must uncompress it

  if [[ "${genome_fasta_file}" == *".gz" ]]; then

  uncomp_genome_fasta=$(basename "${genome_fasta_file}" ".gz");
  gzip -dc ${genome_fasta_file} > ${uncomp_genome_fasta};

 "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${uncomp_genome_fasta}"

  else

  mkdir -p "${star_genome_dir}";

  "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${genome_fasta_file}";

  fi
fi
