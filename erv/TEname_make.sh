wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -O hg38.te.txt.gz
zcat hg38.te.txt.gz | grep -E 'LINE|SINE|LTR|DNA|Retroposon' | cut -f 11 | sort | uniq > hg38.TEname.txt