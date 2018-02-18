#!bin/tcsh -f

set MAGIC=$1
set group=$2
set chrom=$3
set minX=$4
set maxX=$5
set run=any

set minIntronSupport=`cat MetaDB/$MAGIC/RunList | gawk '{n++}END{printf("%d",int(n/20)+1);}'`
set sMin="-minimalIntronSupport $minIntronSupport"

if (! -d  tmp/EHITS.$MAGIC/$chrom) mkdir  tmp/EHITS.$MAGIC/$chrom


set sDuo=" -newDoubleIntrons "
if (! -e tmp/EHITS.$MAGIC/$chrom/introns.de_duo) then
  echo "collating de_duo"
  if (! -e tmp/EHITS.$MAGIC/$chrom/introns.de_duo) then
    # touch  tmp/EHITS.$MAGIC/any.$chrom.limit1
    foreach run (`cat MetaDB/$MAGIC/RunList`)
        # was OR pA pT SL1 SL2 SL3 SL4 SL5 SL6 SL7 SL8 SL9 SL10 SL11 SL12, 2016_04_01 we now prefer to do the XA in phase f3
       foreach dd (OR)
         echo $run $chrom $dd
         set ok=0
         foreach n (1 2 5 10 20 50 100 200 1000)
           if ($ok == 1) continue
           if (-e tmp/$dd/$run/$chrom/$run.$n.0.txt) then
             set ok=1
             cat  tmp/$dd/$run/$chrom/$run.$n.*.txt | gawk -F '\t' '{if($3 != chrom)next;}/^pA/{print}/^pT/{print}/^SL/{print}/^INTRON/{print}'  chrom=$chrom |  sort   >>  tmp/EHITS.$MAGIC/$chrom/introns.de_duo
           endif
         end
       end
    end
  endif
endif
  
set  iDuo=tmp/EHITS.$MAGIC/$chrom/introns.de_duo
if (-e $iDuo) then
  set sDuo="-sxxNewIntronsFileName  $iDuo -newDoubleIntrons"
  ls -ls $iDuo
endif


if (! -e tmp/EHITS.$MAGIC/$chrom/f1.introns.de_uno.gz) then
  echo "collating de_uno"
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    gunzip -c tmp/OR/$run/d1.$run.de_uno.txt.gz | gawk -F '\t' '{if($3==chrom)print;}' chrom=$chrom >> tmp/EHITS.$MAGIC/$chrom/f1.introns.de_uno
    gunzip -c tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz | gawk '/^Intron/{ok=0;split($2,aa,"_");if (aa[1]==chrom)ok=1;a1=aa[3];a2=aa[4];next;}/^Run_/{if(ok!=1)next;n=$6;if(a1<a2){s="Forward";b1=a1-35;b2=a1-1;c1=a2+1;c2=a2+35;ln=a2-a1+1;}else{s="Reverse";b1=a1+35;b2=a1+1;c1=a2-1;c2=a2-35;ln=a1-a1+1;}printf("INTRON\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t\t%d\t%d\n",s,chrom,b1,b2,chrom,c1,c2,ln,int(n));}' chrom=$chrom >>  tmp/EHITS.$MAGIC/$chrom/f1.introns.de_uno
  end
  gzip tmp/EHITS.$MAGIC/$chrom/f1.introns.de_uno
endif

set iUno=tmp/EHITS.$MAGIC/$chrom/f1.introns.de_uno.gz

set sUno=""

if (-e $iUno) then
  set sUno="-sxxDeUnoIntronsFileName  $iUno"
  ls -ls $iUno
endif

set chr=""
# if ($species == worm)   set chr="CHROMOSOME_"

set mylanes=tmp/EHITS.$MAGIC/$chrom/_lane.$$
if (-e $mylanes) \rm $mylanes
touch $mylanes

if (! -e tmp/EHITS.$MAGIC/$chrom/double.$group.$chrom.$minX.$maxX.txts ) then

  echo "running bin/geneelements"
  if (0) then
    echo " gunzip -c tmp/PHITS_genome/*/*.$chrom.[frn]*.minerr0.hits.gz | bin/geneelements $sUno $sDuo $sMin -sxxChromosome $chr$chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -minX $minX -maxX $maxX   -minimalSupport $minExonCover  -o tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txt "
    if (-e $TMPDIR/f1.$$) \rm $TMPDIR/f1.$$
    foreach run (`cat MetaDB/$MAGIC/RunList`)
      gunzip -c tmp/PHITS_genome/$run/*.$chrom.[frn]*.minerr0.hits.gz >> $TMPDIR/f1.$$
    end
    cat $TMPDIR/f1.$$ | bin/geneelements  $sUno $sDuo $sMin -sxxChromosome $chr$chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -minX $minX -maxX $maxX  -minimalSupport $minExonCover  -o $TMPDIR/f1.$$.out

    echo "bin/geneelements done"

    cat $TMPDIR/f1.$$.out  | sort -u > tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txts
    ls -ls tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txt* 
    \rm $TMPDIR/f1.$$ $TMPDIR/f1.$$.out

  else
    echo " bin/geneelements  $sUno $sDuo $sMin -sxxChromosome $chr$chrom -minX $minX -maxX $maxX  -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz"
    bin/geneelements  $sUno $sDuo $sMin -sxxChromosome $chr$chrom -minX $minX -maxX $maxX   -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -stranded | sort -u  > tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txts
    ls -ls  tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txts
  endif
endif

echo "collate ace"
  cat   tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txts | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$6;u2=$11;dnaa=$12;dnab=$13;dnac=$14;support=$15;feet=$16;if(length(feet)==5){if(feet == "gt_ag" || feet == "gc_ag" || feet == "at_ac" || feet == "ct_ac"){feet=feet  "\n";}else {feet= "Other "  feet  "\n" ;}}else feet="Julot\n" ;typea=$3;typeb=$4;typec=$5; if(typeb == "Exon")next;xx="XI_" ; col="PALEYELLOW";nam=xx group "_" chrom "__" u1 "_" u2 ;names[nam]++; n = names[nam] ;nam = nam "." n ; printf("Sequence %s\ncDNA_clone %s\nIntMap %s %d %d\nIs_Read\nForward\nComposite %d\nColour %s\n%s",nam, nam, $1,u1,u2,support, col,feet) ;sl=0;if(substr(typea,1,2)=="SL"){x1=length(dnaa);sl=substr(typea,3,1)+0;}if(substr(typeb,1,2)=="SL"){x1=length(dnaa);sl=substr(typeb,3,1)+0;};if(sl>0)printf("Transpliced_to SL%d %d\n",sl,x1);pA=0;if(typeb=="pA"){x1=length(dnaa);}if(typec=="pA"){x1=length(dnaa)+length(dnab);pA=length(dnac);};if(pA>0)printf("PolyA_after_base %d\n",x1);if(typeb=="Intron")printf("Intron %s__%d_%d\n",chrom,$8,$9);if($2=="Reverse")dx=-1;else dx=1;if(typea=="Intron")printf("Intron %s__%d_%d\n",chrom,$7+dx,$8-dx);if(typec=="Intron")printf("Intron %s__%d_%d\n",chrom,$9+dx,$10-dx);printf("\n");}' group=$group  > tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.ace
 
echo "collate fasta"

cat   tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.txts  | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$6;u2=$11;typeb=$4;if(typeb == "Exon") { next;xx="XJ_" ;} else {xx="XI_" ; }nam=xx group "_" chrom "__" u1 "_" u2 ; names[nam]++; n = names[nam] ;nam = nam "." n ; printf(">%s\n%s%s%s\n",nam, $12,$13,$14) ;}' group=$group  > tmp/EHITS.$MAGIC/$chrom/f1.double.$group.$chrom.$minX.$maxX.fasta

\rm   $mylanes


touch tmp/EHITS.$MAGIC/$chrom/f1.double_collect.done

exit 0

# test
# we critically miss the introns in known transcripts as collected by phase ii2a
# they should be reformatted as
# gunzip -c tmp_old/EHITS/NB.any.9.de_uno_introns.gz | head
INTRON  Reverse 9       15151   15081   9       14940   14892   gt_ag   140     9
INTRON  Reverse 9       17420   17344   9       17166   17120   gt_ag   177     13
INTRON  Reverse 9       17777   17719   9       17479   17424   gt_ag   239     5

 
