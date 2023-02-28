# build_deseq2_files.awk: script to take two input files of
# featureCounts output files and generate join commands to
# build the two input files required for DESeq2 running
# under R.

BEGIN{file1_descr = "control";
file2_descr = "treatment";
meta_file = "DESeq2_metadata.txt";
printf("ColID\tdescr\n") > meta_file;
nullfld = "0";
fcount = 0;
id_list_file = "all_deseq2_TEs.txt";
ostring = "0,";
countcol = 7;
basecmd = sprintf("join --nocheck-order -a 1 -e \"%s\"",nullfld);
prvcmd = sprintf("%s -o '%s2.%d' %s ",basecmd,ostring,countcol,id_list_file);
outfile = "DESeq2_count_matrix.txt";
cmd = "";
prvfile = "";
hdr = "#TE_id";
}

FNR==NR{printf("%s\t%s\n",basename($0),file1_descr) > meta_file;
add_join($0);
}

FNR < NR{printf("%s\t%s\n",basename($0),file2_descr) > meta_file;
add_join($0);
}

function basename(fullpath)
{
ns = split(fullpath,fps,"/");
if (ns > 0)
  return(fps[ns]);
else
  return("");
}

function add_join(jfilename)
# to add the next join section to existing
# command string
{
hdr = hdr "\t" basename(jfilename);
ostring = "0,";
for (i = 1; i <= NR; i++)
  ostring = ostring sprintf("1.%d,",i+1);
cmd = cmd prvcmd;
prvcmd = sprintf(" '%s' | %s -o '%s2.%d' - ",jfilename,basecmd,ostring,countcol);
prvfile = jfilename;
}

END{printf("printf \"%s\\n\" > %s\n",hdr,outfile);
printf("%s '%s' | tr \" \" \"\\t\" | awk '$1!~/==>/&&$1!=\"Geneid\"' >> '%s'\n",cmd,prvfile,outfile);
}
