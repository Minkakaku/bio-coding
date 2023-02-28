# return the 'A' value for TPM from FeatureCounts output files.

BEGIN{totscaled=0.0;}

NR>2{if($(NF-1)!=0)totscaled+=$NF/$(NF-1);}

END{printf("%.2f\n",totscaled*1.0e3);}
