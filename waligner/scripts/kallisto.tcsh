#!bin/tcsh -f

set phase=$1
set species=$2
set target=$3
set run=$4

if ($phase == klst1) then
      set toto=`cat tmp/Kallisto/$run/fastc.list | gawk '{printf("Fastc/%s.fastc.gz ",$1);}'`
      echo "kallisto quant -i $MagicRootDir/tmp/Kallisto_index/$species.$target.kallisto_index -o $MagicRootDir/tmp/Kallisto/$run/$target --plaintext --single $toto -l 300 -s 50"
            kallisto quant -i $MagicRootDir/tmp/Kallisto_index/$species.$target.kallisto_index -o $MagicRootDir/tmp/Kallisto/$run/$target --plaintext --single $toto -l 300 -s 50
endif