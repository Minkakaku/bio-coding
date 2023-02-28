# append_tpm_FC.awk: take a FeatureCount gene abundance file
# and calculate the TPM for each line, appending the value
# to the output.  This process requires deriving a parameter
# for TPM calculations, requiring the file to be pre-scanned
# 
# usage: awk -f append_tpm_FC.awk <FeatureCount_abundance_file>
#
# command line options:
#   colhdr: header for appended column
#   nonzerotpms=1 : only write non-zero tpm values
#
# Peter Stockwell: Aug-2020

BEGIN{atotal = 0.0;
nonzerotpms = 0;
}

NR==2{printf("%s\t%s\n",$0,colhdr);
atotal = get_tpm_param(FILENAME);
if (atotal <= 0.0)
  {
  printf("TPM calculation failed, check setting of featurecounts_strandedness (-p) in mapping parameters file\n");
  exit(1);
  }
}

NR>2{if($(NF-1)+0!=0)
  {aparam=$NF/$(NF-1);
  tpm = aparam * 1.0e6 / atotal;
  }
  else
    tpm = 0.0;
if ((tpm > 0.0) || (!nonzerotpms))
  printf("%s\t%.4f\n",$0,tpm);
}

function get_tpm_param(filenam)
# scan filenam to derive tpm parameter
{
cmd = sprintf("awk -f get_tpm_parameter.awk %s",filenam);
if (cmd | getline ret > 0)
  {
  close(cmd);
  ns = split(ret,rsplit);
# printf("scanned '%s' getting param=%s\n",filenam,ret);
  return(rsplit[ns] + 0.0);
  }
else
  {
  printf("TPM parameter scan for '%s' failed\n",filenam);
  exit(1);
  }
}
