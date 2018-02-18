#!bin/tcsh -ef

echo -n "f2 tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt start : "
date

if (! -e tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt) then
  if (-e tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt1) \rm tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt1
  foreach target ($Etargets)
    foreach run (`cat MetaDB/$MAGIC/RunList`)
      gunzip -c  tmp/INTRONRUNS/$run/*.u.doubleIntronSupport.ace.gz | gawk '{gsub(/\"/,"",$0);}/^DoubleIntron/{g=$2;}/^Run_U/{gg[g]+= $6;uu[g]+= $6;next;}/^Run_nU/{gg[g]+= $6;nu[g]+= $6;next;}END{for(g in gg)if(gg[g]>0)printf("%s\t%d\t%d\t%d\n",g,gg[g],uu[g],nu[g]);}' >>   tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt1
    end 
  end
 
  cat tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt1 | gawk '{n2[$1]+=$2;n3[$1]+=$3;n4[$1]+=$4;}END{for(k in n2)printf("%s\t%d\t%d\t%d\n",k,n2[k],n3[k],n4[k]);}' > tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt
  \rm   tmp/EHITS.$MAGIC/f2.allDoubleIntrons.txt1 
endif

exit 0

set chrom=$1
set minIntronSupport=`cat MetaDB/$MAGIC/RunList | gawk '{n++}END{printf("%d",int(n/20)+1);}'`
set sMin="-minimalIntronSupport $minIntronSupport"

set mini=$minExonCover

set Wstranded=`cat MetaDB/$MAGIC/GroupW_strandedList | head -1`
set WnewExon=`cat MetaDB/$MAGIC/GroupW_new_exonList| head -1`
set wiggleStranded=""

date
echo "f2.geneelements.tcsh chrom=$chrom sMin=$sMin minExon=$mini "

if (1) then
  if (-e  tmp/WIGGLERUN/$Wstranded/$chrom/R.chrom.u.f.BF.gz) then
    set wiggleStranded="-wiggle_f tmp/WIGGLERUN/$Wstranded/$chrom/R.chrom.u.f.BF.gz   -wiggle_r tmp/WIGGLERUN/$Wstranded/$chrom/R.chrom.u.r.BF.gz"
  else if (-e  tmp/WIGGLEGROUP/$Wstranded/$chrom/R.chrom.u.f.BF.gz) then
    set wiggleStranded="-wiggle_f tmp/WIGGLEGROUP/$Wstranded/$chrom/R.chrom.u.f.BF.gz   -wiggle_r tmp/WIGGLEGROUP/$Wstranded/$chrom/R.chrom.u.r.BF.gz"
  endif
endif

if ($WnewExon == "" && $Wstranded != "") set WnewExon=$Wstranded

echo " Wstranded=$Wstranded"
echo "wiggleStranded=$wiggleStranded"
echo "WnewExon=$WnewExon"

if (-e tmp/WIGGLERUN/$WnewExon/$chrom/R.chrom.frns.u.BF.gz) then
  set wiggle_ns="-wiggle_ns tmp/WIGGLERUN/$WnewExon/$chrom/R.chrom.frns.u.BF.gz"
else if (-e tmp/WIGGLEGROUP/$WnewExon/$chrom/R.chrom.frns.u.BF.gz) then
  set wiggle_ns="-wiggle_ns tmp/WIGGLEGROUP/$WnewExon/$chrom/R.chrom.frns.u.BF.gz"
endif
echo "wiggle_ns=$wiggle_ns"

# do we need the same file as constructed in f1 or a different one ?
# what threshold should we use in both cases

set sDuo=""
if (1) then
  set  iDuo=tmp/EHITS.$MAGIC/$chrom/introns.de_duo
  if (! -e $iDuo) then
    echo "missing file $iDuo, please run phase f1"
    exit 1
  endif
  if (-e $iDuo) set sDuo="-sxxNewIntronsFileName  $iDuo "
endif

set iUno=tmp/EHITS.$MAGIC/$chrom/introns.de_uno.gz
if (! -e $iUno) then
  echo "missing file  $iUno. please run phase f1"
  exit 1
endif

set sUno=""

if (-e $iUno) set sUno="-sxxDeUnoIntronsFileName  $iUno"

 
if (0) then

      echo "bin/geneelements -newExons -minimalSupport $mini $wiggle_ns $wiggleStranded  $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  ----  tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.nostrand.txt"

            bin/geneelements -newExons -minimalSupport $mini $wiggle_ns $wiggleStranded  $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  > tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.nostrand.txt
else
      echo "bin/geneelements -newExons -minimalSupport $mini $wiggleStranded -strand     $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  > tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.strand.txt"
            bin/geneelements -newExons -minimalSupport $mini $wiggleStranded -strand     $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  > tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.strand.txt
      echo "bin/geneelements -newExons -minimalSupport $mini $wiggleStranded -antistrand $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  > tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.antistrand.txt"
            bin/geneelements -newExons -minimalSupport $mini $wiggleStranded -antistrand $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  > tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.antistrand.txt
endif

if (0) then

  cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.*strand.txt | gawk -F '\t' '/^EXON/{nn++;a1=0+$4;a2=0+$5;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);nam="XE__" chrom "__" a1 "_" a2 ;printf("Sequence %s\ncDNA_clone %s\nIs_read\nColour PALEMAGENTA\nForward\nIntMap %s %d %d\n-D Composite\n", nam,nam,$3,a1,a2);n=split($6,aa,",");for(i=1;i<=n;i++){split(aa[i],bb,":");printf("Composite %d %d %d\n",bb[1],bb[2],bb[3]);}printf("\n");if(length($7)>1)printf("DNA %s\n%s\n\n", nam,$7);}' >  tmp/EHITS.$MAGIC/$chrom/f2.newexons.XE.ace

  cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.*strand.txt | gawk -F '\t' '/^INTRON/{nn++;a1=0+$4;a2=0+$5;a3=0+$6;a4=0+$7;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);nam="XD__" chrom "__" a1 "_" a2 "__" a3 "_" a4 ;dx=1;if(a3>a4)dx=-1;printf("Sequence %s\ncDNA_clone %s\nIs_read\nColour Green2\nComposite %d %d %d %d %d\nForward\nIntMap %s %d %d\n-D Intron\nIntron %s__%s_%s\n\n", nam,nam,$8,$9,$10,$11,$12,$3,a1,a4,chrom,a2+dx,a3-dx);if(length($13) + length($14)>1)printf("DNA %s\n%s%s\n\n", nam,$13,$14);}' >   tmp/EHITS.$MAGIC/$chrom/f2.newexons.XD.ace
endif

  cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.*strand.txt | gawk -F '\t' '/^pA/{nn++;a1=0+$4;a2=0+$5;n=0+$6;n2=$7+$8;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);nam="XA__" chrom "__" a1 "_" a2 ;printf("Sequence %s\ncDNA_clone %s\nIs_read\nColour LIGHTORANGE\nComposite %d %d\nReverse\nIntMap %s %d %d\nPolyA_after_base 9\n\n", nam,nam,n,n2,$3,a1,a2,nam);if(length($9))printf("DNA %s\n%s%s\n\n", nam,"TTTTTTTT",$9);}' >   tmp/EHITS.$MAGIC/$chrom/f2.newexons.XA.ace

  cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.*strand.txt | gawk -F '\t' '/^SL[1-9]/{nn++;a1=0+$4;a2=0+$5;n=0+$6;n2=$7+$8;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);nam="X" $1 "__" chrom "__" a1 "_" a2 ;printf("Sequence %s\ncDNA_clone %s\nIs_read\nColour GREEN\nComposite %d %d\nForward\nIntMap %s %d %d\nTranspliced_to %s 1\n\n", nam,nam,n,n2,$3,a1,a2,$1);if(length($7))printf("DNA %s\n%s\n\n", nam, $7);}' >   tmp/EHITS.$MAGIC/$chrom/f2.newexons.XSL.ace

endif

touch tmp/EHITS.$MAGIC/$chrom/f2.collect.done

exit 0

# this hack may be suboptimal
cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.txt | gawk -F '\t' '/^EXON/{nn++;a1=0+$4;a2=0+$5;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);if (a1<a2)printf("newexon_%s_%d\t%d\t%s\t%d\t%d\n",cali,nn,1,chrom,a1,a2);}' cali=$MAGIC chrom=$chrom > tmp/METADATA/$MAGIC.newx.$chrom.f.sponge

cat tmp/EHITS.$MAGIC/$chrom/f2.newexons.min$mini.txt | gawk -F '\t' '/^EXON/{nn++;a1=0+$4;a2=0+$5;chrom=$3;gsub(/CHROMOSOME_/,"",chrom);if (a1>a2)printf("newexon_%s_%d\t%d\t%s\t%d\t%d\n",cali,nn,1,chrom,a1,a2);}'  cali=$MAGIC chrom=$chrom > tmp/METADATA/$MAGIC.newx.$chrom.r.sponge

cat  tmp/METADATA/$MAGIC.newx.$chrom.[fr].sponge | sort -u > tmp/METADATA/$MAGIC.newx.$chrom.ns.sponge

exit 0


