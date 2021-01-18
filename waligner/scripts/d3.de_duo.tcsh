#!bin/tcsh -f

set dd=$1
set chrom=$2
set run=$3
set ii=$4
set intronMaxLength=$5
set minX=$6
set maxX=$7
set oL=$8
set section=$9
set toto=$TMPDIR/aceview.de_duo.$$

echo $toto
if (-e $toto) \rm $toto

set mask=""
if (-e TARGET/Targets/$species.genome.mask.txt) set mask="-genomeMask TARGET/Targets/$species.genome.mask.txt"

if ($dd == Transloc)  then

 
    # ls -ls  tmp/$dd/$chrom/$run.clean.gz
    gunzip -c  tmp/$dd/$run/*/$run.clean.gz >> $toto

  cat $toto | $bin/geneelements -strategy $Strategy -newIntrons -minimalSupport $ii -non_gt_ag 8 -transsplicing 8 -intronMaxLength $intronMaxLength -overhangLength $oL $mask -o  tmp/Transloc/$run/$run.txt

else
  set ok=0

    if (-e tmp/$dd/$run/$chrom/$run.clean.gz) then
      ls -ls  tmp/$dd/$run/$chrom/$run.clean.gz
      gunzip -c  tmp/$dd/$run/$chrom/$run.clean.gz >> $toto
      set ok=1
    endif

  # 2015_03_11, by now most intron are found in de_uno mode and often know their strand
  # so we remove the de_duo search for non classic intron by removibg the option ' -non_gt_ag 8 '
  if ($ok == 1) then
    echo "bin/geneelements -hitFile $toto -strategy $Strategy -newIntrons -minimalSupport $ii -intronMaxLength $intronMaxLength -minX $minX -maxX $maxX -overhangLength $oL $mask -o tmp/$dd/$run/$chrom/$run.$ii.$section.txt"
          bin/geneelements -hitFile $toto -strategy $Strategy -newIntrons -minimalSupport $ii -intronMaxLength $intronMaxLength -minX $minX -maxX $maxX -overhangLength $oL $mask -o tmp/$dd/$run/$chrom/$run.$ii.$section.txt
  else
    echo "ok = 0 in de.de_duo.tcsh, dd=$dd no .clean.gz file"
    ls -ls  tmp/$dd/$run/$chrom/$run.clean.gz
  endif
endif

echo $toto
ls -ls $toto
#if (-e $toto) \rm $toto
