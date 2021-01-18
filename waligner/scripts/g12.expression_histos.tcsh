#!bin/tcsh -f

set target=$1
set GM=$2

# Histogram of all measured non-NA indexes in the different targets
# and
# Histogram of the number of genes with given cumulated support: number of tags accross the whole project

# g12a: Prepare the data in parallel
if ($target != export) then

  if ($GM == GENE) set GM2=genes
  if ($GM == MRNAH) set GM2=transcripts
  if ($GM == MA) set GM2=microarray
  
    set target2=$target
    if ($target == av) set target2=AceView
    echo "$GM $target"
ls -ls  RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz
    if (-e  RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz && ! -e RESULTS/Expression/unique/$target/$MAGIC.index.$GM.distrib.txt2) then
      gunzip -c  RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz | gawk 'BEGIN{printf("0\n400\n");}/^Gene_group/{ok=0;next;}/^Gene/{ok=1;next;}/^$/{ok=0;}/^if(ok<1)next;}/^Run_U/{if(1){if($11=="NA")next;print 2*int(5*$3);}}' | bin/histo -w 200 -plain -o RESULTS/Expression/unique/$target/$MAGIC.index.$target2.$GM.distrib -title "$target $GM2 not NA"
      \rm  RESULTS/Expression/unique/$target/index.$target2.$GM.distrib.correl.txt
      if (-e  RESULTS/Expression/unique/$MAGIC.expression_index.histo.txt) \rm  RESULTS/Expression/unique/$MAGIC.expression_index.histo.txt
    endif

  exit 0
  # 2019_05_30 Danielle says this is no longer interesting
   if (-e RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz && ! -e RESULTS/Expression/unique/$target/$MAGIC.cumulated_support_notNA.$target2.$GM.txt) then 
      if ($GM == GENE) set tt1=Gene
      if ($GM == MRNAH) set tt1=Transcript
      if ($GM == MA) set tt1=Probe
      set tt2=tt1
      gunzip -c RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz  | gawk '{if($1 == tt1 || $1==tt2){g=$2;next;}}/zzRhs1074/{next;}/^Run_U/{if(1){if($11=="NA")next;nn[g]+=$6;}}END{for (g in nn)printf("%s\t%d\n",g,nn[g]);}' tt1=$tt1 tt2=$tt2 | sort -k 2nr >  RESULTS/Expression/unique/$target/$MAGIC.cumulated_support_notNA.$target2.$GM.txt 

      if (-e  RESULTS/Expression/unique/$MAGIC.cumulated_support_of_annotated_genes_and_transcripts.histo.txt) \rm RESULTS/Expression/unique/$MAGIC.cumulated_support_of_annotated_genes_and_transcripts.histo.txt
    endif 


  exit 0
endif

# g12b: Export the complete table 
set toto=RESULTS/Expression/unique/$MAGIC.expression_index.histo.txt
echo -n "# " > $toto
date >> $toto
echo "Histogram of all measured non-NA indexes in the different targets" >> $toto
echo 1 | gawk '{printf("\t\t");for(i=0;i<=200;i++)printf("\t%.1f",i/5);printf("\n");}' > $toto.1

foreach GM (GENE MRNAH MA)
  if ($GM == GENE) set GM2=genes
  if ($GM == MRNAH) set GM2=transcripts
  if ($GM == MA) set GM2=microarray
  foreach target ($Etargets $other_targets )
    set target2=$target
    if ($target == av) set target2=AceView

     if (-e RESULTS/Expression/unique/$target/$MAGIC.index.$target2.$GM.distrib.txt) then
      echo -n RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz >> $toto.1
      ls -ls RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.ace.gz | gawk '{printf("\t%s",$0);}' >> $toto.1
      echo -n "\t$target $GM2\t" >> $toto.1
      cat RESULTS/Expression/unique/$target/$MAGIC.index.$target2.$GM.distrib.txt | tail -n +4 | cut -f 2 | scripts/transpose >> $toto.1
    endif
  end
end
cat $toto.1 | scripts/transpose >> $toto
\rm $toto.1

######### 
# collect just cumulative read count among not NA runs

# setenv MAGIC SEQC_D
# histo in true log2 scale
set toto=RESULTS/Expression/unique/$MAGIC.cumulated_support_of_annotated_genes_and_transcripts.histo.txt

echo -n "### $toto : " > $toto
date >> $toto
echo "Cumulated support of genes and transcripts in project $MAGIC, excluding NA runs in number of reads" >> $toto

echo '#' > $toto.1
foreach GM (GENE MRNAH)
  foreach target (RefSeq EBI av)

    cat   TYTY/$target/$MAGIC.cumulated_support_notNA.$GM.txt  | gawk -F '\t'  '/ERCC-/{next;}/ERCC_vector/{next;}{n=$2 ; if(n==0)next;k=n;j=0;m=1;while(k>0){j++;k=int(k/2);m=2*m;}j--;nn[j,1]++;if(j>jmax)jmax=j;}END{m=1;for(j=0;j<=jmax;j++){printf("%d\t%d\n",m,nn[j,1]);m=2*m;}}' | sort -k 1,1nr | gawk '{if(zz<1)printf("# Level\t%s\t%s Cumul\n",tt,tt);zz=1;}{n+=$2;printf("%s\t%s\t%d\t%d\t%d\n",GM,tt,$1,$2,n);}' GM=$GM tt=$target >> $toto.1
  end
end

cat  $toto.1 | gawk -F '\t'  '/^#/{next;}{z=$3;x=$1 " " $2;zz[z]=1;xx[x]=1;nn1[z,x]=$4;nn2[z,x]=$5;}END{for(x in xx)printf("\tcount %s\tcumul %s",x,x);for(z in zz){printf("\n%s",z);for(x in xx)printf("\t%s\t%s",nn1[z,x],nn2[z,x]);}printf("\n");}' > $toto.2

cat $toto.2 | scripts/transpose | sort  | scripts/transpose | sort -k 1n >> $toto

\rm $toto.[12]

exit 0
setenv MAGIC SEQC_NB_Fatigue
foreach GM (GENE MRNAH)
  foreach target (RefSeq EBI av)
    set target2=$target
    if ($target == av) set target2=AceView
 
    \cp RESULTS/Expression/unique/$target/$MAGIC.cumulated_support_notNA.$target2.$GM.txt TYTY/$MAGIC.cumulated_support_notNA.$target2.$GM.$target.txt
    \cp ~/NB/RESULTS/Expression/unique/$target/NB.cumulated_support_notNA.$target2.$GM.txt TYTY/NB.cumulated_support_notNA.$target2.$GM.$target.txt
    \cp ~/Fatigue_non_stranded/RESULTS/Expression/unique/$target/Fatigue12.cumulated_support_notNA.$GM.txt TYTY/Fatigue.cumulated_support_notNA.$GM.$target.txt

  end
end

foreach GM (GENE MRNAH)
  foreach target (RefSeq EBI av)
    set target2=$target
    if ($target == av) set target2=AceView

    mkdir TYTY/$target
    cat TYTY/*.cumulated_support_notNA.$target2.$GM.$target.txt | gawk '{nn[$1]+=$2}END{for(g in nn)printf("%s\t%d\n",g,nn[g]);}' > TYTY/$target/SEQC_NB_Fatigue.cumulated_support_notNA.$target2.$GM.txt
  end
end

mkdir TYTY/introns
 cat  tmp/introns/SEQC_main.cumulated_support.txt ~/NB_2013/tmp/introns/NB.cumulated_support.txt ~/Fatigue_non_stranded/tmp/introns/Fatigue12.cumulated_support.txt | gawk -F '\t' '{gg[$1]+=$2;}END{for(g in gg)if(gg[g]>0)printf ("%s\t%d\n", g,gg[g]);}' | sort -k 2nr > TYTY/introns/$MAGIC.cumulated_support.txt

