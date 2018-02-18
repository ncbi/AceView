#!bin/tcsh -f

# if i use -ef, i block between hello6a and hello7 on Rhs436/f.1
set run=$1
set lane=$2

set mytmp=$TMPDIR/aceview.s1.snp_confirmation.$$
#set mytmp=./mytmp

mkdir $mytmp
mkdir $mytmp/COUNT
mkdir $mytmp/COUNT/$run

phases4:
 echo -n "Phase s4: $run $lane  snp confirmation $mytmp "
 date

set nTarget=0
# apply the previous score limit and find
# reads confirming the SNP-edited sequences

set target=snp
set target_class=t_snp   # transcriptome SNPs
set target_class=s_snp
echo hello

  echo -n "$run $lane target=$target  "
  date

  if (! -e tmp/NEWHITS_snp/all.edited_sequence.m.snp.fasta.gz) then
    echo "Missing file   tmp/NEWHITS_snp/all.edited_sequence.m.snp.fasta.gz"
    exit 1
  endif
  if (! -e tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz) then
    echo "Missing file   tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz"
    exit 1
  endif

  set minAli8=`cat Fastc/$run/Max_probe_length | gawk '{x=int(.90 * $1); if(x>65) x = 65; if (x<35) x=35 ; printf("%d", x);}'`

  set ln1=`gunzip -c tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz | head -2 | tail -1 | gawk '{ln=int(1+length($1)/2);print ln - 2}'`
  set ln2=`gunzip -c tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz | head -2 | tail -1 | gawk '{ln=int(1+length($1)/2);print ln + 2}'`
  set Xsnp="-targetStart $ln1 -targetStop $ln2"

  if (! -d $mytmp/NEWHITS_$target) mkdir $mytmp/NEWHITS_$target
  if (! -d $mytmp/NEWHITS_$target/$run) mkdir $mytmp/NEWHITS_$target/$run

  set isSplice=" "
  set slBonus=" "
  set bonus=0
  set targetBonus="-targetBonus $bonus"

  set clipPolyA=" "
  set clipPolyT=" "
  foreach run2 (`cat MetaDB/$MAGIC/RunRnaList`)
    if ($run == $run2) then
      if ($target == genome || $target == gdecoy || $target == mito || $target == chloro) then
        set clipPolyA="-clipPolyA"
      else
        set clipPolyA="-SclipPolyA"
      endif
      foreach run3 (`cat MetaDB/$MAGIC/RunNonStrandedList`)
        if ($run == $run3) then
          if ($target == genome  || $target =~ chrom*  || $target == gdecoy || $target == mito || $target == chloro) then
            set clipPolyT="-clipPolyT"
          else
            set clipPolyT="-SclipPolyT"
          endif
        endif
      end
    endif
  end

  set isSolid=" "
  set iiSolid=0
  foreach run2 (`cat MetaDB/$MAGIC/RunSolidList`)
    if ($run == $run2) set iiSolid=1
  end

  if ($iiSolid == 1) set isSolid="-solid"

  set isStranded=" "
  
  set isSplice=" "

  if ($species == worm) set slBonus="-slBonus 5"

  set isSplice=" "
  
  echo "## $Strategy #$target $run $isSolid $clipPolyT  $isSplice $slBonus $isStranded $targetBonus"

  set minEntropy4=`echo $minAli | gawk '{printf("%d",int((2.0*$1)/3));}'`
  set off3="1 -seedShift 5"

  set v2=""
  if ($?exitAdaptorRaw) then
    set v2=" -exitAdaptor $exitAdaptorRaw " 
  endif

echo  ZZZZZ | gzip > $mytmp/ZZZZZ.gz
echo "prepare previous score PLUS ZERO, for wild type, PLUS 3 for mutant, PLUS 8 for brk"
date
  set prevScore=""
  if (-e tmp/COUNT/$lane.hits.gz)  then
echo hello1
# Here we compensate the target bonus of the mito chloro and rrna
    gunzip -c  tmp/COUNT/$lane.hits.gz $mytmp/ZZZZZ.gz  tmp/NEWHITS_snp/$run/too_many_hits.list.gz  | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next;}{if(z<1){p=$1;s=$2;if($8=="A_mito" || $8=="C_chloro" || $8 == "B_rrna")s--;if(p==oldp)next;oldp=p;printf("%s\t%d\n",p,s+0);}else{printf("%s\t%d\n",$1,$3);}}' | sort -u | gzip -1 > $mytmp/COUNT/$lane.best_score.w.snp.gz
echo hello2
    gunzip -c  $mytmp/COUNT/$lane.best_score.w.snp.gz | gawk '{printf("%s\t%d\n",$1,3 + $2);}' | gzip -1 > $mytmp/COUNT/$lane.best_score.m.snp.gz
echo hello3
    gunzip -c  $mytmp/COUNT/$lane.best_score.w.snp.gz | gawk '{printf("%s\t%d\n",$1,8 + $2);}' | gzip -1 > $mytmp/COUNT/$lane.best_score.m.brk.gz
echo hello4
    touch  $mytmp/COUNT/$lane.multi.list
echo hello5
    if (-e tmp/NEWHITS_snp/$run/too_many_hits.list.gz) then
      gunzip -c tmp/NEWHITS_snp/$run/too_many_hits.list.gz | grep Multi >> $mytmp/COUNT/$lane.multi.list
    endif
echo hello6
    if (-e tmp/PHITS_genome/$lane.overRepresentedSeeds.gz) then
echo hello6a
      gunzip -c tmp/PHITS_genome/$lane.overRepresentedSeeds.gz | grep Multi >> $mytmp/COUNT/$lane.multi.list
    endif
echo hello7
  else
    echo "missing file tmp/COUNT/$lane.hits.gz"
    exit 1
  endif
  set nTarget=1

echo "align "
date

# increase maxHit to 100 in case we have several overlapping SNPs to test
set maxHit2=100

  set prevScore=" -selectPreviousScore  $mytmp/COUNT/$lane.best_score.w.snp.gz "
  set reject="-reject tmp/NEWHITS_snp/$run/too_many_hits.list.gz"
  set minAlix=$minAli

  echo "bin/clipalign -i Fastc/$lane.fastc.gz -gzo  $isSolid -t tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz  -maxHit $maxHit2 $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAlix  -clipN $clipN -minAli $minAlix $isSplice $slBonus  $targetBonus  $Xsnp $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -o $mytmp/NEWHITS_$target/$lane.w.any.pre  -target_class $target_class $v2 $prevScore  -strategy $Strategy $reject"
  bin/clipalign -i Fastc/$lane.fastc.gz -gzo  $isSolid -t tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz  -maxHit  $maxHit2 $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAlix  -clipN $clipN -minAli $minAlix $isSplice $slBonus  $targetBonus  $Xsnp $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -o $mytmp/NEWHITS_$target/$lane.w.any.pre  -target_class $target_class $v2 $prevScore  -strategy $Strategy $reject

  set prevScore=" -previousScore  $mytmp/COUNT/$lane.best_score.w.snp.gz "
  set select="-select $mytmp/COUNT/$lane.multi.list"

  bin/clipalign -i Fastc/$lane.fastc.gz -gzo  $isSolid -t tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz  -maxHit $maxHit2  $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAlix  -clipN $clipN -minAli $minAlix $isSplice $slBonus  $targetBonus  $Xsnp $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -o $mytmp/NEWHITS_$target/$lane.w.multi  -target_class $target_class $v2 $prevScore  -strategy $Strategy $select

foreach sb (snp brk)
  set prevScore=" -selectPreviousScore  $mytmp/COUNT/$lane.best_score.m.$sb.gz "
  if ($sb == brk) set reject=" -reject tmp/PHITS_genome/$lane.overRepresentedSeeds.gz "
  if ($sb == brk) set minAlix=$minAli8

  echo "bin/clipalign -i Fastc/$lane.fastc.gz -gzo  $isSolid -t tmp/NEWHITS_snp/all.edited_sequence.m.$sb.fasta.gz  -maxHit $maxHit2  $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAlix  -clipN $clipN -minAli $minAlix  $isSplice $slBonus  $targetBonus  $Xsnp $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -o $mytmp/NEWHITS_$target/$lane.m.$sb.pre  -target_class $target_class $v2  $prevScore  -strategy $Strategy  $reject"
  bin/clipalign -i Fastc/$lane.fastc.gz -gzo  $isSolid -t tmp/NEWHITS_snp/all.edited_sequence.m.$sb.fasta.gz  -maxHit $maxHit2 $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAlix  -clipN $clipN -minAli $minAlix  $isSplice $slBonus  $targetBonus  $Xsnp $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -o $mytmp/NEWHITS_$target/$lane.m.$sb.pre  -target_class $target_class $v2  $prevScore  -strategy $Strategy  $reject

end

foreach wm (w m)
  gunzip -c $mytmp/NEWHITS_$target/$lane.$wm.*.pre.hits.gz | bin/bestali  -filter $species -maxHit $maxHit2 -countBest -strategy $Strategy -exportBest   -gzo -o $mytmp/NEWHITS_$target/$lane.$wm.hits
end


  touch $mytmp/NEWHITS_$target/$lane.done


# -doubleSeed
# exit 1
#######################################################################################
#######################################################################################
## cleanup

  if (! -d tmp/NEWHITS_snp) mkdir tmp/NEWHITS_snp
  if (! -d tmp/NEWHITS_snp/$run) mkdir tmp/NEWHITS_snp/$run

  mv $mytmp/NEWHITS_snp/$run/*.pre.hits.gz tmp/NEWHITS_snp/$run/PRE
  mv $mytmp/NEWHITS_snp/$run/* tmp/NEWHITS_snp/$run

\rm -rf $mytmp

echo -n "done "
date

exit 0


##############################################################

# debug test
Fastc/Rhs493/f.4.fastc.gz 

clipalign -i totop.fasta    -t tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz  -maxHit 20 -SclipPolyA  -SclipPolyT -minEntropy 16 -seedLength 16 -probeMinLength 24  -clipN 2 -minAli 24      -targetBonus 0  -targetStart 119 -targetStop 123   -seedOffset 1 -seedShift 5 -intronMaxLength  50000  -target_class s_snp  -exitAdaptor AGGGGATAGG,TCGTATGCCGTCTTCTGCTTGAAAAAA,ttagacatatctccgtcgtagggatccc   -selectPreviousScore  /export/home/TMP/aceview.s1.snp_confirmation.10395/COUNT/Rhs493/f.4.best_score.w.snp.gz   -strategy RNA_seq -reject tmp/NEWHITS_snp/Rhs493/too_many_hits.list.gz


clipalign -i totop.fasta    -t tmp/NEWHITS_snp/all.edited_sequence.w.fasta.gz  -minAli 24  -seedOffset 1 -seedShift 5 
