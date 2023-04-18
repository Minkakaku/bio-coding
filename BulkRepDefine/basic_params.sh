# basic_params.sh: define basic RepExpress env variables to show:
# 
# 1. Control verbosity as script runs
# 2. Locations of required executables (e.g. STAR, Stringtie, featureCounts)
# 3. Locations of raw genomic sequence file(s) and generated index for mapping:
# 4. Locations of gtf annotation files
#
# Edit this file to reflect the parameters required on your target system.
# This information is expected to suit a whole series of runs which would
# be made against the same genome and annotations.  Definitions for each
# individual run should be made in a run parameter script (run_params.sh,
# for instance).

# 1. Control verbosity: set to empty string to reduce feedback

verbose="yes";

# Control retention of working scripts, set this string to empty to retain them

delete_temp_files="";

# 2.  Locations of required executables and invariant run parameters:

# path to STAR executable: can be left empty if STAR is already on your exec PATH
#  The command 'which STAR' will indicate this for you.

# STAR mapping run parameters:
starfiltermax="150";
staranchormax="150";
starthreads="16";

# path to featureCounts executable: empty if on your path
# 'which featureCounts' will indicate this

path_to_featurecounts="";

# featureCounts parameters

featurecounts_threads="16";
featurecounts_overlap="--minOverlap 25";
featurecounts_overlap_partial="--fracOverlap 1"


# path to DMAP executables, particularly identgeneloc.
#  leave empty if these are already on your path
# 'which identgeneloc' will indicate this

path_to_dmap="/data/gpfs02/wmo/software/DMAP/src/";

# 3. Locations of raw genomic sequence file(s) and generated index for mapping:

# genome fasta files - this is for all sequences in one large fasta file.
#  It is possible to have separate files for each chromosome, but it
#  tends to get a bit messy.

genome_fasta_file="/data/gpfs02/wmo/genome/GENCODE/GRCm39/ref/GRCm39.genome.fa";

# the location where the index files are written and read from

star_genome_dir="/data/gpfs02/wmo/genome/UCSC/index/star";

# 4. Locations of gtf annotation files

# name of the UCSC repeat element source file downloaded from
# https://genome.ucsc.edu/cgi-bin/hgTables
# with web page settings:
#
# 'assembly': must be set to the required build
#
# 'group': set to Repeats
#
# 'output format': set to 'all fields from selected table' in order to retrieve
#   details that are omitted from the GTF - gene transfer format (limited) format.
#
# 'output file': should be set to something to avoid having the data sent directly
#  for display by the browser (e.g. hg19_ucsc_repeats.txt).
#
# The file can be gzip compressed or not.  The gzip compression at the web interface
#  didn't seem to work.
ucsc_repeat_gtf="/data/gpfs02/wmo/genome/UCSC/rep_gtf/rep_RepeatMasker_GRCm39.gtf";

# name of gencode gene annotation gtf file, available from
# https://www.gencodegenes.org/human/release_32lift37.html or
# by ftp from ftp.ebi.ac.uk at /pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping
#
# The desired file for GRCh37/hg19 is gencode.v32lift37.annotation.gtf.gz
#
# The GRCh38 release is available from the related directory.
#
ucsc_gene_gtf="/data/gpfs02/wmo/genome/UCSC/gen_gtf/mm39.ncbiRefSeq.gtf"
#uncompress first
