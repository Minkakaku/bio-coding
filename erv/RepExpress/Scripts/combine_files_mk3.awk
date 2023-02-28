# to combine tpm or feature files with join, filling absent
# fields.
# required input: file of sample tpm names, one
# per line.

BEGIN{nullfld = "-";
  fcount = 0;
  id_list_file = "all_TE_ids.txt";
  ostring = "0,1.2,";
  countcol = 3;
  basecmd = sprintf("join -a 1 -e \"%s\"",nullfld);
  prvcmd = sprintf("%s -o '%s2.%d' %s ",basecmd,ostring,countcol,id_list_file);
  outfile = "combined.txt";
  cmd = "";
  prvfile = "";
  hdr = "#TE_id\tSense";
}

{
ns = split($0,s0,"/");
hdr = hdr "\t" substr(s0[ns],1,length(s0[ns])-9);
ostring = "0,";
for (i = 1; i <= NR+1; i++)
  ostring = ostring sprintf("1.%d,",i+1);
cmd = cmd prvcmd;
prvcmd = sprintf(" %s | %s -o '%s2.%d' - ",$0,basecmd,ostring,countcol);
prvfile = $0;
}

END{printf("echo \"%s\" > %s\n",hdr,outfile);
printf("%s %s | tr \" \" \"\\t\" >> %s\n",cmd,prvfile,outfile);
}
