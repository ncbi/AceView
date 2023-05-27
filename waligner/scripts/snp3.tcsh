#!/bin/tcsh
# may 2023   analyse the tsmps and the BRS snps using snpsummary for the SEQC2_capture project

set phase=$1
set  zone=$2

if ($zone == "") then
  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    scripts/snp3.tcsh $phase $zone
  end
  goto done
endif


if ($phase == 1) then 
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/$MAGIC.snp3a.$zone -snpType 0  -minSnpFrequency 20 -minSnpCover 20
   goto done
endif

if ($phase == 999) then 
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/BRS.$zone -snpType 0 -minSnpFrequency 20 -minSnpCover 20
   goto done
endif

if ($phase == p) then 
  set toto=toto
  echo "pparse MetaDB/$MAGIC/runs.ace" > $toto
  echo "pparse MetaDB/$MAGIC/groups.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/samples.ace" >> $toto
  echo "save" >> $toto
  echo "quit" >> $toto

  bin/tace tmp/TSNP_DB/$zone -no_prompt < $toto
  goto done
endif



done:
  echo done
  exit 0


