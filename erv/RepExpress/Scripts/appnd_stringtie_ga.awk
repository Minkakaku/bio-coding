# appnd_stringie_ga.awk: script to read counts and TPMs from Stringtie gene abundance file
# into an array, indexed by the Ensemble gene ID, then append matching values to a
# file from identgeneloc with (or without) gene names appended.  This script
# will also append the gene name from Stringtie
#
# Peter Stockwell: Sep-2020
#

BEGIN{if (length(stie_file_name) <= 0)
  {
  printf("script needs a Stringtie file name as a -v stie_file_name= command option\n");
  exit(1);
  }
else
  {
  while (getline ret < stie_file_name > 0)
    {
    ns = split(ret,rsplit);
    ens_gname = rsplit[1];
    if (!(ens_gname in st_tpms))
      {
      st_tpms[ens_gname] = rsplit[ns];
      gnames[ens_gname] = rsplit[2];
      }
    }
  close(stie_file_name);
  ensid_col = 16;
  }
}

$1~/#/{printf("%s\tgene\tst_tpm\n",$0);}

$1!~/#/{ens_name = substr($ensid_col,2,length($ensid_col) - 2);
if (ens_name in st_tpms)
  printf("%s\t%s\t%s\n",$0,gnames[ens_name],st_tpms[ens_name]);
else
  printf("%s\t-\t-\n",$0);
}
