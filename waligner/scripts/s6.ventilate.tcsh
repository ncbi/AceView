#!bin/tcsh -f
set MAGIC=$1
set group=$2
set minSnpCover=$3
set ercc=$4

setenv ici `pwd`
set mytmp=$TMPDIR/aceview.s1.snp_accumulate.$$.
# set mytmp=tmp/KL.$MAGIC/$group/RUNS

mkdir $mytmp

# ventilate in append mode all the runs of the group
foreach run (`cat MetaDB/$MAGIC/_q2 | gawk -F '\t' '{if($1 == g || $4 == g){z=$1;n[z]++;if(n[z]==1)print z;}}' g=$group | sort -u `)
  if (-e  tmp/NEWHITS_snp/$run/s590.multi.countEdited) then
    echo "Found file  tmp/NEWHITS_snp/$run/s590.multi.countEdited"
    cat tmp/NEWHITS_snp/$run/s590.multi.countEdited | gawk -F '\t' '/^#/{next}/Chrom/{next}{if (ercc==1 && substr($1,1,5) != "ERCC-") next ; split($1,aa,"|"); print >> cali aa[1] ;}' cali="$mytmp/confirmed."  ercc=$ercc
  else
    echo "missing file  tmp/NEWHITS_snp/$run/s590.multi.countEdited"
  endif

end

# regroup and add a threshold
pushd $mytmp
foreach ff (confirmed.*)
  cat $ff | bin/snp -merge -minCover $minSnpCover -run $group -gzo > $ici/tmp/KL.$MAGIC/$group/$ff.gz
end 
popd

touch tmp/KL.$MAGIC/$group/s6.ventilate.done
#\rm -rf  $mytmp

exit 0

