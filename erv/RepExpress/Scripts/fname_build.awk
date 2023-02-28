# take a file of sample fastq names and
# generate a file of corresponding tpm output
# files for each

BEGIN{tpm_file_location="./";}

/.fastq.gz$/{printf("%s\n",rejigname($0,"_at_r1.fastq.gz","_at_r1_M_U_FC_tpm.genloc"));}

/.fastq$/{printf("%s\n",rejigname($0,"_at_r1.fastq","_at_r1_M_U_FC_tpm.genloc"));}

/.fq$/{printf("%s\n",rejigname($0,"_at_r1.fq","_at_r1_M_U_FC_tpm.genloc"));}

/.fq.gz$/{printf("%s\n",rejigname($0,"_at_r1.fq.gz","_at_r1_M_U_FC_tpm.genloc"));}

function rejigname(oldname,oldextension,newextension)
# return the name with the new extension appended in place of old
{
ns = split(oldname,namesplit,"/");
newname=sprintf("%s/%s%s",tpm_file_location,
                 substr(namesplit[ns],1,length(namesplit[ns])-length(oldextension)),
                 newextension);
return(newname);
}
