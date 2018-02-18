#!bin/tcsh -ef

set run=$1
set chrom=$2

if ($chrom == TRANSLOC) then

   gunzip -c tmp/PHITS_genome/$run/*.introns.gz  tmp/PHITS_mito/$run/*.introns.gz | gawk -F '\t' -f scripts/d1.de_uno.awk run=$run | sort -k 1,1  -k 2,2  -k 4,4n | gzip > tmp/Transloc/$run/d1.$run.de_uno.txt.gz

else

  gunzip -c tmp/PHITS_genome/$run/*.introns.gz  tmp/PHITS_mito/$run/*.introns.gz  | gzip > tmp/OR/$run/d1.$run.de_uno.txt.gz
  ls -ls  tmp/OR/$run/d1.$run.de_uno.txt.gz

endif


exit 0

