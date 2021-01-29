#! bin/tcsh -f

set phase=$1
set run=$2

# use the overhangs to associate each read to an SL type and a target_class
# then scan the COUNT/hits file to know the hierarchy and keep only the best ali
# this give us SL type and oordinates in a genome or mRNA target

mkdir  $TMPDIR/SLpA.$run

if ($phase == "merge") goto phaseMerge

echo -n 'Collect :'
date


touch  $TMPDIR/SLpA.$run/select
\rm $TMPDIR/SLpA.$run/*
foreach lane (`cat Fastc/$run/LaneList`)
  if (-e $TMPDIR/$run.SLpA.collect.txt.gz.$$) \rm $TMPDIR/$run.SLpA.collect.txt.gz.$$
  touch $TMPDIR/$run.SLpA.collect.txt.gz.$$
  foreach target ($Etargets genome) 
    set ff=tmp/PHITS_$target/$lane.overhangs.gz
    if (-e $ff) then
      source scripts/target2target_class.txt
      gunzip -c $ff | gawk -F '\t' '/^#/{next;}{type=substr($6,1,2);if(type == "SL" || type == "pA" || type == "pT") print $1 "\t" $2  "\t" $6 "\t" tc "\t" $7 "\t" $8 "\t" $9 ;}' tc=$target_class | gzip >>  $TMPDIR/SLpA.$lane.collect.gz
    endif
  end
  gunzip -c  $TMPDIR/SLpA.$lane.collect.gz ZZZZZ.gz tmp/COUNT/$lane.hits.gz | gawk -F '\t' '/ZZZZZ/{zz++;next;}{if(zz+0<1){type[$1,$4,$5]=$0;next;}}{if(length(type[$1,$8,$11])>2){print type[$1,$8,$11] ; type[$1,$8,$11] = 0 ;} next;}' >> $TMPDIR/SLpA.$run/select
end

# cumulate the counts
cat $TMPDIR/SLpA.$run/select |  gawk -F '\t' '{n[$2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 ] += $2;}END {for (k in n) printf ("%09d\t%s\n",n[k],k) }' | sort -k 5 -k 6 -k 8n > $TMPDIR/SLpA.$run/cumul

# split by target
cat $TMPDIR/SLpA.$run/cumul | grep Z_genome | gawk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\n", $5,$7,$3,$6,$1);}' > $TMPDIR/SLpA.$run/cumul.genome
cat $TMPDIR/SLpA.$run/cumul | grep -v Z_genome > $TMPDIR/SLpA.$run/cumul.mrnas

# remap to the genome using the coordinates of the mrnas
# export in single file  
#   remapping        mrna x1 Exon x2 chrom a1 a2
#   plus types       mrna x  type count
# sort by (mrna x) if type is inside a known exon we can remap it 
# in a single pass without looping  
echo -n 'Remap :'
date

gunzip -c tmp/METADATA/mrnaRemap.gz >  $TMPDIR/SLpA.$run/mrnaRemap
foreach target ($Etargets)
  source scripts/target2target_class.txt
  cat $TMPDIR/SLpA.$run/cumul.mrnas | gawk -F '\t' '{if($4 == tc)print;}' tc=$target_class  | gawk -F '\t' '{printf("%s\t%09d\t%s\t%s\t%s\n", $5,$7,$3,$6,$1);}' > $TMPDIR/SLpA.$run/cumul.$target
  cat $TMPDIR/SLpA.$run/mrnaRemap | gawk -F '\t' "/^$target_class/"'{printf("%s\t%09d\tExon\t%s\t%s\t%s\t%s\n",$2,$3,$4,$5,$6,$7);}' > $TMPDIR/SLpA.$run/mrnaRemap.$target
  cat  $TMPDIR/SLpA.$run/cumul.$target $TMPDIR/SLpA.$run/mrnaRemap.$target | sort > $TMPDIR/SLpA.$run/mix.$target 
end

echo -n 'Cumul :'
date
# add up all methods
cat $TMPDIR/SLpA.$run/cumul.genome | gawk -F '\t' '{if ($3=="pT"){$3="pA";if($4=="Forward")$4="Reverse";else $4="Forward";}gsub(/^pT/,"pA",$3);printf("%s\t%s\t%09d\t%s\t%d\n", $3,$1,$2,$4,$5);}' >  $TMPDIR/SLpA.$run/counts.unsorted
foreach target ($Etargets)
  cat  $TMPDIR/SLpA.$run/mix.$target | gawk -F '\t' '{type = $3; if(type == "Exon"){mrna=$1;x1=0+$2;x2=0+$4;chrom=$5;a1=0+$6;a2=0+$7;strand="Forward";if(a1>a2)strand="Reverse";next;}x=0+$2;if (mrna == $1 &&x >= x1 &&x  <= x2){if(strand == "Forward")b=x-x1+a1;else b=a1-(x-x1);gsub(/^pT/,"pA",type);printf("%s\t%s\t%09d\t%s\t%d\n",type,chrom,b,strand,0+$5);}}' >> $TMPDIR/SLpA.$run/counts.unsorted
end
cat  $TMPDIR/SLpA.$run/counts.unsorted | sort | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 "\t" $4;if (z != old){if (n>0)printf("%s\t%d\n",old,n);old = z ; n = 0 ;} n+= $5;}END{if(n>0)if (n>0)printf("%s\t%d\n",old,n);}'  > $TMPDIR/SLpA.$run/final

# clean up
cat $TMPDIR/SLpA.$run/final | tags | sort -k 2nr
gzip $TMPDIR/SLpA.$run/final
mv $TMPDIR/SLpA.$run/final.gz tmp/SLpA/$run.SLpA.gz

goto done


############################

phaseMerge:

# cumul
foreach run2 (`cat MetaDB/$MAGIC/r2sublib MetaDB/$MAGIC/g2r | gawk '{if($1==run)print $2;}' run=$run | sort -u`)
  if (-e tmp/SLpA/$run2.SLpA.gz) then
    zcat tmp/SLpA/$run2.SLpA.gz >>  $TMPDIR/SLpA.$run/cumul
  endif
end
# merge
cat $TMPDIR/SLpA.$run/cumul | sort | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 "\t" $4;if (z != old){if (n>0)printf("%s\t%d\n",old,n);old = z ; n = 0 ;} n+= $5;}END{if (n>0)printf("%s\t%d\n",old,n);}' > $TMPDIR/SLpA.$run/final
cat $TMPDIR/SLpA.$run/final | tags | sort -k 2nr
gzip  $TMPDIR/SLpA.$run/final 
mv  $TMPDIR/SLpA.$run/final.gz  tmp/SLpA/$run.SLpA.gz

############################

done:
# \rm -rf  $TMPDIR/SLpA.$run
echo -n 'done :'
date
