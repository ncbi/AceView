#!bin/tcsh -ef

set run=$1

  gunzip -c tmp/PHITS_genome/$run/*.introns.gz  tmp/PHITS_mito/$run/*.introns.gz  | gzip > tmp/OR/$run/d1.$run.de_uno.txt.gz
  ls -ls  tmp/OR/$run/d1.$run.de_uno.txt.gz
  touch tmp/OR/$run/d1.de_uno.done

echo done


