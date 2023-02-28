# generate a script to select required columns from
# geneloc output files
# Note: this version uses awk

# input file is a list of files containing
# column-selected geneloc values.
# reject line 1 of input files.

BEGIN{lcnt=0;
sorted_id_file="all_TE_ids.txt";}

{ns = split($0,s0,"/");
flist[NR] = sprintf("%s_c946.txt",s0[ns]);
printf("awk 'NR>1{printf(\"%%s\\t%%s\\t%%s\\n\",$9,$4,$6);}' %s > %s\n",
        $0,flist[NR]);
lcnt = NR;
}

END{printf("cut -f 1,2");
for (i=1; i <= lcnt; i++)
  printf(" \\\n %s",flist[i]);
printf(" | sort -u -k 1,1 > %s\n",sorted_id_file);
for (i = 1; i <= lcnt; i++)
  {
  printf("rm %s\n",flist[i]) > "remove_tmp_files.sh";
  printf("%s\n",flist[i]) > "genloc_colsel_files.txt";
  }
close("remove_tmp_files.sh");
close("genloc_colsel_files.txt");
}
