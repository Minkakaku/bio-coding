# add_gene_names_to_genloc.awk: read an ENSEMBL gene_id
# from a column (16) in input and search for the corresponding
# gene name.  Append to line.
#
#

BEGIN{ensid_col = 15;
ensid_to_gene_file = "./ensid_vs_gname.txt";
while (getline ret < ensid_to_gene_file > 0)
  {
  ns = split(ret,rsplit);
  if (!(rsplit[1] in gnames))
    gnames[rsplit[1]] = rsplit[ns];
  }
close(ensid_to_gene_file);
}

$1~/#/{printf("%s\tgene\n",$0);}

$1!~/#/{ens_name = substr($ensid_col,2,length($ensid_col)-2);
if (ens_name in gnames)
  printf("%s\t%s\n",$0,gnames[ens_name]);
else
  printf("%s\t-\n",$0);
}
