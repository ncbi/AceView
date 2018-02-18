#!bin/tcsh -f

set group=$1

foreach zone (mito SpikeIn rrna `ls tmp/SNP_ZONE/zone?.*.txt | sed -e 's/tmp\/SNP_ZONE\///' -e 's/\.txt//' | sort -k 1n`)
  echo $zone
  date >  tmp/SNPGROUP/$group/toto.$zone
  foreach run (`cat MetaDB/$MAGIC/_q2 | gawk -F '\t' '{if($4 == g || $1 == g) print $1;}' g=$group | sort -u`)
    if (-e tmp/SNP/$run/$zone.count.u.txt.gz) then
      gunzip -c tmp/SNP/$run/$zone.count.u.txt.gz >> tmp/SNPGROUP/$group/toto.$zone
    endif
  end
  cat tmp/SNPGROUP/$group/toto.$zone | bin/snp -merge -run $group -minCover 5 | gzip > tmp/SNPGROUP/$group/$zone.count.u.txt.gz
  \rm tmp/SNPGROUP/$group/toto.$zone 

end

touch tmp/SNPGROUP/$group/s2.snp_group_collect.done

