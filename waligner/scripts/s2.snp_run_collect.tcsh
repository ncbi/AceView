#!bin/tcsh -f

set run=$1
set minF=$2
set minCover=$3
set maxCover=$4

foreach zone (mito SpikeIn rrna `ls tmp/SNP_ZONE/zone?.*.txt | sed -e 's/tmp\/SNP_ZONE\///' -e 's/\.txt//' | sort -k 1n`)
  echo $zone

    if (-e tmp/SNP/$run/$zone.detect.u.snp.gz) then
      gunzip -c tmp/SNP/$run/$zone.detect.u.snp.gz | bin/snp -merge -run $run -minFrequency $minF -minCover $minCover | gawk -F '\t' '{if($9 > mx)next;print;}' mx=$maxCover | gzip > tmp/SNP/$run/$zone.select.u.txt.gz
    endif

end

touch tmp/SNPGROUP/$run/s2.snp_run_collect.done


