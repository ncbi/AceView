#!bin/tcsh -f

set target=$1
set t2=$2
if ($target == "go") goto go
if ($target == "AUC") goto AUC
if ($target == "plusbas2") goto plusbas2
set out=$3
set ratio_bound=$4
set histo_shifting=$5
set maxGene=$6

set uu=nu
set myace=OTHER_PIPELINES/$target.ace

set deepType=`echo $target | gawk '{if(index($1,"junction")+index($1,"intron")>1)print "-deepIntron"; else { if(index($1,"transcript")>0)print "-deepTranscript"; else print "-deepGene";}}'`
set uu=`echo $target | gawk '{if(index($1,"AGLuK")>0)print "-u"; else print "-nu";}'`
set deepMix=`echo $target | gawk '{if(index($1,"_GJ_") > 1)print "-intronGeneMix"; else print " ";}'`

if (1 && ! -e $out/$MAGIC.newSampleClassificationBySignature5.txt) then
   echo "bin/geneindex $uu $deepMix -$deepType $myace -runList tmp/GENERUNS/$MAGIC.RunListSorted -runAce tmp/GENERUNS/$MAGIC.snp.info.ace  -o $out/$MAGIC -gzo -pA -correl -compare -ratio_bound $ratio_bound -histo_shifting $histo_shifting -maxGene $maxGene -noHisto $gx_opt "
         bin/geneindex $uu $deepMix -$deepType $myace -runList tmp/GENERUNS/$MAGIC.RunListSorted -runAce tmp/GENERUNS/$MAGIC.snp.info.ace -o $out/$MAGIC -gzo -pA -correl -compare -ratio_bound $ratio_bound -histo_shifting $histo_shifting -maxGene $maxGene  -iterate $gx_opt 
else
  echo "found $out/$MAGIC.newSampleClassificationBySignature.txt "
endif

set toto1=$out/AUCBIG5
echo "RES\t$Results" > $toto1
echo "TARGET\t$target"  >> $toto1
echo "T2\t$t2"  >> $toto1
if ($maxGene == mxgenes) set maxGene=5
echo "PARAM\tR$ratio_bound/$maxGene"g  >> $toto1
  # cat $out/*.beta.?.txt | grep MCC >> $toto1
foreach ff (`ls $out/*.beta.?.txt | grep -v MYC`)
  cat $ff  | gawk '/AUC/{print;next;}/phenotype/{print;next;}/Title/{print;next;}/^Plus/{print;next;}/^Minus/{print;next;}/^Beta/{print;next}/^Alpha/{print;next}' >> $toto1
end
\rm $out/$MAGIC.diff.*

exit 0


#############################


##### A2 T0 T2 R2 E2 S2 Mag0 Mag2 Mat0 Mat2
##### T2 ST2 Mag0 Mag2 Mat0 Mat2)  R0 ST0 ST2)

go:
foreach ii (S0 S2 AT0 AT2)
  foreach mxgenes (010 050 100 300 500)
echo "### $ii $mxgenes"

  setenv RR "1.0 1.5 2.0 3"

  if ($ii == A0) then
    # RNA seq best in digital
    setenv targets "AceView"
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif
  if ($ii == A1) then
    continue  # RNA seq best in digital
    setenv targets "AceView"
    setenv Results R_score
    setenv gx_opt "-digital 1 -normalize 1 -local 0"
  endif
  if ($ii == A2) then
    setenv targets "AceView"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == AT0) then
    # RNA seq best in digital
    setenv targets "AceView_transcripts"
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif
  if ($ii == AT1) then
    continue  # RNA seq best in digital
    setenv targets "AceView_transcripts"
    setenv Results R_score
    setenv gx_opt "-digital 1 -normalize 1 -local 0"
  endif
  if ($ii == AT2) then
    setenv targets "AceView_transcripts"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == S2) then
    setenv targets "seqc "
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == S0) then
    setenv targets "seqc "
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif
  if ($ii == ST0) then
    # continue  # RNA seq best in digital
    setenv targets "seqc_transcripts "
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif
  if ($ii == ST1) then
    continue  # RNA seq best in digital
    setenv targets "seqc_transcripts "
    setenv Results R_score
    setenv gx_opt "-digital 1 -normalize 1 -local 0"
  endif
  if ($ii == ST2) then
    setenv targets "seqc_transcripts "
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0 -threshold 12"
  endif
  if ($ii == R0) then
    # RNA seq best in digital
    setenv targets RefSeq
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif
  if ($ii == R1) then
    continue  # RNA seq best in digital
    setenv targets RefSeq
    setenv Results R_score
    setenv gx_opt "-digital 1 -normalize 1 -local 0"
  endif
  if ($ii == R2) then
    setenv targets "RefSeq"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == E0) then
    setenv targets "Encode"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == E2) then
    setenv targets "Encode"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0"
  endif
  if ($ii == Mag0) then
    setenv targets "AGLuK_rescaled80"
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0 -keepIndex -threshold 10"
  endif
  if ($ii == Mag2) then
    setenv targets "AGLuK_rescaled80"
    setenv Results R_digital 
    setenv gx_opt "-digital 2 -normalize 1 -local 0 -keepIndex -threshold 10"
  endif
  if ($ii == Mat0) then
    setenv targets "AGLuK_at"
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0 -MA 60 -threshold 10"
  endif
  if ($ii == Mat2) then
    setenv targets "AGLuK_at"
    setenv Results R_digital
    setenv gx_opt "-digital 2 -normalize 1 -local 0 -MA 60 -threshold 10" 
  endif
  if ($ii == Ia) then
    setenv targets "AceView_introns2_Magic"
    setenv Results R_plain
    setenv gx_opt "-digital 0 -normalize 1 -local 0"
  endif

foreach target ($targets)
  set t2=$target


  if ($target == EBI_37_67_Magic)        set t2=EBI_37_67_Magic
  if ($target == AceView_genes_Fudan)    set t2=Fudan_AceView_Gene
  if ($target == AceView_genes_Magic)    set t2=Magic_AceView_Gene
  if ($target == AceView_introns2_Magic) set t2=Magic_AceView_Junction
  if ($target == AceView_GJ_Magic)       set t2=Magic_AceView_GJ
  if ($target == AceView_mRNA_Fudan)     set t2=Fudan_AceView_Transcript
  if ($target == Agilent_microarray)     set t2=Koln_Agilent_Microarray
  if ($target == FBK_genes_FBK)          set t2=FBK_FBK_Gene
  if ($target == FBK2_genes_FBK)         set t2=FBK2_FBK_Gene
  if ($target == RefSeq_genes_BGI)       set t2=BGI_RefSeq_Gene
  if ($target == RefSeq_genes_Magic)     set t2=Magic_RefSeq_Gene
  if ($target == RefSeq_genes_Su)        set t2=Su_RefSeq_Gene
  if ($target == RefSeq_exons_Su)        set t2=Su_RefSeq_Exon
  if ($target == RefSeq_mRNA_Hong)       set t2=Hong_RefSeq_Transcript
  if ($target == RefSeq_mRNA_Su)         set t2=Su_RefSeq_Transcript
  if ($target == RefSeq_exons_Lilly)     set t2=Lilly_RefSeq_Exon
  if ($target == UCSC_exons_Lilly)       set t2=Lilly_UCSC_Exon
  if ($target == UCSC_genes_BGI)         set t2=BGI_UCSC_Gene
  if ($target == Vega_exons_Lilly)       set t2=Lilly_Vega_Exon



  if (! -e OTHER_PIPELINES/$target.ace) then
    echo "missing file : OTHER_PIPELINES/$target.ace"
    continue
  endif
  if (! -d OTHER_PIPELINES/$Results) mkdir  OTHER_PIPELINES/$Results
  if (! -d OTHER_PIPELINES/$Results/$target) mkdir OTHER_PIPELINES/$Results/$target
  foreach maxGene ($mxgenes)
    foreach ratio_bound ($RR)
      foreach histo_shifting (1)
          set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
          if (! -d $out) mkdir $out
          if (-d $out && ! -e $out/$MAGIC.outcome.1.txt) then
            set mem=`echo $target | gawk '{i=index($1,"intron") + index($1,"GJ") + index($1,"exon") + index($1,"transcripts");if(i>1)print "32G";}'`
            if (-e $out/AUC) then
              cat $out/AUC >> ./AUC
            endif
            if (1) then
              scripts/submit $out/$MAGIC "scripts/geneindex.optimizer.tcsh $target $t2 $out $ratio_bound $histo_shifting $maxGene" $mem
              # scripts/geneindex.optimizer.tcsh $target $t2 $out $ratio_bound $histo_shifting $maxGene 
            endif
          else
            # echo "found $out/$MAGIC.outcome.1.txt"
          endif
      end
    end
  end
end

end
end

exit 0

######


##########################################

plusbas2:

# grep results directly in the C code output

 #   echo Title_line | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk title=1 >> $toto

set toto1=OTHER_PIPELINES/$Results/toto1
echo "RES\t$Results" >> $toto1
foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_mRNA_Fudan Agilent_microarray FBK_genes_FBK  FBK2_genes_FBK RefSeq_genes_BGI RefSeq_genes_Magic  RefSeq_mRNA_Hong RefSeq_mRNA_Su RefSeq_genes_Su   UCSC_genes_BGI  AceView_introns2_Magic AceView_GJ_Magic RefSeq_exons_Su RefSeq_exons_Lilly UCSC_exons_Lilly Vega_exons_Lilly)
  if (! -e OTHER_PIPELINES/$target.ace) continue

echo "TARGET\t$target"  >> $toto1
  if (! -d OTHER_PIPELINES/$Results/$target) continue
  foreach maxGene (500) 
    foreach ratio_bound (1.0 1.2 1.5 2.0 2.5 3 4 6)
      foreach histo_shifting (1)
        set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
        if (-d $out) then
           echo "PARAM\tR$ratio_bound/$maxGene"g  >> $toto1
	   cat $out/*.beta.?.txt | grep MCC  >> $toto1
        endif
      end
    end
  end
end

exit 0

AUC:
cat OTHER_PIPELINES/Results_*/*/*/AUCBIG5  > AUC
set toto=AUC.txt
date > $toto
echo "ROC AUC for various iterations of the Magic classifier" >> $toto
cat AUC | gawk -F '\t'  -f scripts/geneindex_optimizer_C_scanner.awk >> $toto

\cp $toto RESULTS

exit 0

####################
set toto=ROC_optimizer.maxGene.$Results.txt
echo -n "# " > $toto
date > $toto
echo "ROC AUC for various iterations of the Magic classifier" >> $toto

    echo Title_line | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk title=1 >> $toto
foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_mRNA_Fudan Agilent_microarray FBK_genes_FBK  FBK2_genes_FBK RefSeq_genes_BGI RefSeq_genes_Magic  RefSeq_mRNA_Hong RefSeq_mRNA_Su RefSeq_genes_Su   UCSC_genes_BGI  AceView_introns2_Magic AceView_GJ_Magic RefSeq_exons_Su RefSeq_exons_Lilly UCSC_exons_Lilly Vega_exons_Lilly)
  if (! -e OTHER_PIPELINES/$target.ace) continue
ls  OTHER_PIPELINES/$target.ace
  if (! -d OTHER_PIPELINES/$Results/$target) continue
  foreach maxGene (10 20  40 100 300 500) 
    foreach ratio_bound (2.0 2.5 3.0 3.5 4.0 5.0)
      foreach histo_shifting (1)
           set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
          if (-e $out/_roc.1) then
            cat $out/_roc.1 ZZZZZ $out/_roc.2 | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk m=$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene  >> $toto
           endif
       end
    end
  end
  echo >> $toto
end

  \cp $toto RESULTS

exit 0

set toto=ROC_optimizer.Ratio.txt
echo -n "# " > $toto
date > $toto
echo "ROC AUC for various iterations of the Magic classifier" >> $toto

    echo Title_line | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk title=1 >> $toto
foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_mRNA_Fudan Agilent_microarray FBK_genes_FBK  FBK2_genes_FBK RefSeq_genes_BGI RefSeq_genes_Magic  RefSeq_mRNA_Hong RefSeq_mRNA_Su RefSeq_genes_Su   UCSC_genes_BGI  AceView_introns2_Magic AceView_GJ_Magic RefSeq_exons_Su RefSeq_exons_Lilly UCSC_exons_Lilly Vega_exons_Lilly)
  if (! -e OTHER_PIPELINES/$target.ace) continue
ls  OTHER_PIPELINES/$target.ace
  if (! -d OTHER_PIPELINES/$Results/$target) continue
  foreach maxGene (10 20  40 100 300 500) 
      foreach histo_shifting (1)
    foreach ratio_bound (2.0 2.5 3.0 3.5 4.0 5.0)
           set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
          if (-e $out/_roc.1) then
            cat $out/_roc.1 ZZZZZ $out/_roc.2 | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk m=$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene  >> $toto
           endif
       end
    end
  end
  echo >> $toto
end

  \cp $toto RESULTS


set toto=ROC_optimizer.shift.txt
echo -n "# " > $toto
date > $toto
echo "ROC AUC for various iterations of the Magic classifier" >> $toto

foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_mRNA_Fudan  AceView_introns2_Magic AceView_GJ_Magic Agilent_microarray FBK_genes_FBK UCSC_genes_BGI RefSeq_genes_BGI RefSeq_genes_Magic   RefSeq_mRNA_Hong RefSeq_genes_Su RefSeq_mRNA_Su RefSeq_exons_Su RefSeq_exons_Lilly UCSC_exons_Lilly  Vega_exons_Lilly) 
  if (! -e OTHER_PIPELINES/$target.ace)  continue
  if (! -d OTHER_PIPELINES/$Results/$target) mkdir OTHER_PIPELINES/$Results/$target

  if (! -e OTHER_PIPELINES/$target.ace) continue
  if (! -d OTHER_PIPELINES/$Results/$target) mkdir OTHER_PIPELINES/$Results/$target
      foreach histo_shifting (1)
  foreach maxGene (10 20  40 100 300) 
    foreach ratio_bound (2.0 2.5 3.0 3.5 4.0 5.0)
          set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
          if (-e $out/_roc.1) then
            cat $out/_roc.1 ZZZZZ $out/_roc.2 | gawk -F '\t' -f scripts/geneindex_optimizer_scanner.awk m=$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene  >> $toto
           endif
       end
    end
  end
  echo >> $toto
end

  \cp $toto RESULTS

pushd OTHER_PIPELINES/$Results

set toto=MCC_optimizer.txt
echo -n "# " > $toto
date > $toto
echo "MCC for various iterations of the Magic classifier" >> $toto
  echo $dd >> $toto
    echo Title_line | gawk -F '\t' -f ../../scripts/geneindex_optimizer_scanner.awk title=1 >> $toto
  foreach dd (`ls -d */*`)
    if (-e $dd/_roc.1) then
      cat $dd/_roc.1 ../../ZZZZZ $dd/_roc.2 | gawk -F '\t' -f ../../scripts/geneindex_optimizer_scanner.awk m=$dd mm=MCC >> $toto
    endif
  end
  \cp $toto ../../RESULTS

popd

# we need to understand the sp != sn in AUC=89 MCC=62
# cat AceView_genes_Fudan/R2.5.Shift2.MaxGene40/NB.outcome.1.txt | sort -k 8n | cut -f 16
cat OTHER_PIPELINES/$Results/AceView_genes_Fudan/R2.5.Shift2.MaxGene40/NB.outcome.1.txt | grep _HR |  sort -k 8n | gawk -F '\t' -f scripts/NB.rocplot.awk kk=17 ii=8

################### Histograms

foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_junctions_Magic AceView_mRNA_Fudan Agilent_microarray FBK_genes_FBK RefSeq_genes_BGI RefSeq_genes_Magic RefSeq_genes_Su RefSeq_mRNA_Hong RefSeq_mRNA_Su RefSeq_exons_Lilly UCSC_exons_Lilly UCSC_genes_BGI Vega_exons_Lilly )
   cat $target.ace | gawk '/Rhs469/{printf("%d\n",int( $6));}' | sort -k 1nr | tail -n +10 | gawk '{printf("%d\n",10*log(1+$1)/log(2));}END{print 500;}' | bin/histo -plain -w 100 -o $target.LogTagCountHisto.Rhs469 &
end


################### Exportation

echo "#" > toto_export
foreach target (AceView_genes_Fudan AceView_genes_Magic AceView_mRNA_Fudan Agilent_microarray FBK_genes_FBK  FBK2_genes_FBK RefSeq_genes_BGI RefSeq_genes_Magic  RefSeq_mRNA_Hong RefSeq_mRNA_Su RefSeq_genes_Su   UCSC_genes_BGI  AceView_introns2_Magic AceView_GJ_Magic RefSeq_exons_Su RefSeq_exons_Lilly UCSC_exons_Lilly Vega_exons_Lilly)


  if ($target == AceView_genes_Fudan)    set t2=Fudan_AceView_Gene:55
  if ($target == AceView_genes_Magic)    set t2=Magic_AceView_Gene:60
  if ($target == AceView_introns2_Magic) set t2=Magic_AceView_Junction:95
  if ($target == AceView_GJ_Magic)       set t2=Magic_AceView_GJ:90
  if ($target == AceView_mRNA_Fudan)     set t2=Fudan_AceView_Transcript:65
  if ($target == Agilent_microarray)     set t2=Koln_Agilent_Microarray:45
  if ($target == FBK_genes_FBK)          set t2=FBK_FBK_Gene:40
  if ($target == FBK2_genes_FBK)         set t2=FBK2_FBK_Gene:42
  if ($target == RefSeq_genes_BGI)       set t2=BGI_RefSeq_Gene:25
  if ($target == RefSeq_genes_Magic)     set t2=Magic_RefSeq_Gene:50
  if ($target == RefSeq_genes_Su)        set t2=Su_RefSeq_Gene:15
  if ($target == RefSeq_exons_Su)        set t2=Su_RefSeq_Exon:70
  if ($target == RefSeq_mRNA_Hong)       set t2=Hong_RefSeq_Transcript:35
  if ($target == RefSeq_mRNA_Su)         set t2=Su_RefSeq_Transcript:20
  if ($target == RefSeq_exons_Lilly)     set t2=Lilly_RefSeq_Exon:75
  if ($target == UCSC_exons_Lilly)       set t2=Lilly_UCSC_Exon:80
  if ($target == UCSC_genes_BGI)         set t2=BGI_UCSC_Gene:30
  if ($target == Vega_exons_Lilly)       set t2=Lilly_Vega_Exon:85



\rm  _mm_export
echo "ALL_SEX\tNB_female/NB_male\tsex_sex_sex" >> _mm_export
echo "ALL_FAV\tOS_Died_TrainingSet/OS_Survived_TrainingSet:1\tMCC sort on 38 favorable versus unfavorable (extremes) clf=Training Overall survival" >> _mm_export
echo "ALL_OS\tOS_Died_TrainingSet/OS_Survived_TrainingSet:3\tMCC sort on 40 overal survival clf=Training Overall survival" >> _mm_export
echo "ALL_EFS\tEventFree_TrainingSet/EventOccured_TrainingSet\tMCC sort on 39 event free survival clf=Training Event free" >> _mm_export
echo "HR_OS\tHR_Died_TrainingSet_43/HR_Survived_TrainingSet_42\tClassifier: High Risk" >> _mm_export
echo "HR_EFS\tHR_EF/HR_EO\tClassifier: HR EF" >> _mm_export
echo "MYCN_FAV\tMYCN2_high/MYCN2_low:1\tsort on 39 event free survival clf=MYCN april iterated" >> _mm_export
echo "MYCN_OS\tMYCN2_high/MYCN2_low:2\tsort on 40 overal survival clf=MYCN april iterated" >> _mm_export

# ALL_EFS/OS\tOS_Died_TrainingSet/OS_Survived_TrainingSet:2 
# echo "ALL_FAV/uf\tFavorable_TrainingSet_87/Unfavorable_TrainingSet_45\tMCC sort on 20 favorable versus unfavorable (extremes) clf=Training f/u" >> _mm_export

  cat NB_parameters.txt | gawk -F '\t' '{n++;if(n == 1) {for(i = 2 ; i<= NF ; i++)mm[i] = $i;next;}for (i=2 ; i<= NF ; i++)printf("%s\t%s\t%s\n",mm[i],$1,$i);}' >  NB_parameters.txt2

  foreach mm (`cat _mm_export | gawk -F '\t' '/^#/{next;}{if($1)print $2}'`)
    set t1=`cat  _mm_export | gawk -F '\t' '/^#/{next;}{if($2==mm)print $1;}' mm=$mm`
    set km=`cat  _mm_export | gawk -F '\t' '/^#/{next;}{if($2==mm)print $3;}' mm=$mm`
    set alpha=0

    set histo_shifting=1
    set maxGene=300
    set ratio_bound=2.0
    set maxGene=`cat NB_parameters.txt2 | grep $target | gawk -F '\t' '{split(t1,aa,"_");tt1=aa[2];if(index($1,tt1 "_N")>0)print $3}' t1=$t1`
    set ratio_bound=`cat NB_parameters.txt2 | grep $target | gawk -F '\t' '{split(t1,aa,"_");tt1=aa[2];if(index($1,tt1 "_R")>0)print $3}' t1=$t1`

    set out=OTHER_PIPELINES/$Results/$target/R$ratio_bound.Shift$histo_shifting.MaxGene$maxGene
    if (! -e $out/NB.newSampleClassificationBySignature.txt) then
       echo "EROR cannot find $out" 
       continue
    endif

    if (-e  $out/_roc.1) then
      set alpha=`cat $out/_roc.1 $out/_km.1 | grep Bestj | grep -v stratified | grep "$km" |  gawk '{i=index($0,"alpha=");split(substr($0,i+6),aa," ");x=aa[1];}END{print x}' `
    endif

   cat $out/$MAGIC.newSampleClassificationBySignature.txt | gawk -F '\t' '/^#Title/{for (i=1;i<=NF;i++){if((length($i)==5 || substr($i,6,1)=="b") && substr($i,1,2)=="NB")i2r[i]=$i;}next;}/^#/{next;}{i=index(t2,":");t2b=substr(t2,1,i-1);score=substr(t2,i+1);i=index(mm,":");if(i>1)mm=substr(mm,1,i-1);sgn=1;if(t1=="ALL_EFS" || t1=="HR_EFS")sgn=-1;if($2==mm)for (i=1;i<=NF;i++)if(i2r[i])printf("SEQC_NB_%s_:%02d:NCBI_MAG_%s\t%s\t%s\t%.2f\n",t1,score,t2b,i2r[i],sgn*$i,0+sgn*alpha);}' t1=$t1 t2=$t2 mm="$mm" alpha=$alpha >> toto_export

  end

end

foreach alpha (0 1 2 3)
  cat toto_export | gawk -F '\t' '/^#/{next;}{m=mm2m[$1];if(0+m<1){mx++;m=mx;mm2m[$1]=m;mm[m]=$1;split($1,aa,":");score[m]=aa[2];};r=r2rr[$2];if(0+r<1){rx++;r=rx;r2rr[$2]=r;rr[r]=$2;};nn[m,r]=$3;if(0+$3>=$4)x=1;else {x=-1;if(alpha==3)x=0;}if (alpha==0 || alpha==3)nn[m,r]=x;if (alpha==1)nn[m,r]=x*score[m];if (alpha==2) nn[m,r]=$3-$4;next;}END{printf("NB");for(m=1;m<=mx;m++)printf("\t%s",mm[m]);for(r=1;r<=rx;r++){printf("\n%s",rr[r]);for(m=1;m<=mx;m++)printf("\t%s",nn[m,r]);}printf("\n");}' alpha=$alpha  > toto_export2.$alpha.txt
end

cat  toto_export2.0.txt | scripts/transpose | sort -r | gawk '{split($1,aa,":");if(old && aa[1]!=old){for(i=2;i<=NF;i++){printf("\t%d",z[i]);z[i]=0;}printf("\n\n");}old=aa[1];  for(i=2;i<=NF;i++)z[i]+=$i; printf("%s%s",aa[1],aa[3]); for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >  RESULTS/NB_predictions_NCBI.binary.txt

cat  toto_export2.1.txt | scripts/transpose | sort -r  | gawk '{split($1,aa,":");if(old && aa[1]!=old)printf("\n\n");old=aa[1]; printf("%s%s",aa[1],aa[3]); for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >  RESULTS/NB_predictions_NCBI.scale.txt
cat  toto_export2.2.txt | scripts/transpose | sort -r  | gawk '{split($1,aa,":");if(old && aa[1]!=old)printf("\n\n");old=aa[1]; printf("%s%s",aa[1],aa[3]); for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >  RESULTS/NB_predictions_NCBI.alpha.txt


# verify in the binary export the av genes
cat RESULTS/SEQC_NB_NCBI_Predictions_sept23.txt | cut -f 1,75 > toto


exit 0

ls -lsd OTHER_PIPELINES/Res*/AceView_genes_Magic/*MaxGene*  OTHER_PIPELINES/Res*/Encode*Magic/*MaxGene* OTHER_PIPELINES/R*/Agilent*/*MaxGen* | grep 'Feb' | gawk '{print $10}' > toto 
\rm tutu
foreach f2 (NB.HR_Died_TrainingSet_44_HR_Survived_TrainingSet_43...beta.0.txt  NB.OS_Died_TrainingSet_OS_Survived_TrainingSet...beta.0.txt)
  foreach ff ( `cat toto` )
    echo -n "$ff\t$f2\t" >> tutu
    head -3 $ff/$f2 | gawk '/AUC2/{printf("\t%s",$0);}' >> tutu
    head -3 $ff/$f2 | gawk '/Parameters/{print;}' >> tutu
  end
end

ls -lsd OTHER_PIPELINES/R_digital/AceView_genes_Magic/R*gene200/*/*   | gawk '{print $10}' > toto 


foreach tt ($MAGIC.OS_Died_TrainingSet_OS_Survived_TrainingSet...beta $MAGIC.HR_Died_TrainingSet_44_HR_Survived_TrainingSet_43...beta $MAGIC.Favorable_TrainingSet_87_Unfavorable_TrainingSet_45...beta)
  \rm tutu
  foreach iter (0 1)
    foreach ff ( `ls OTHER_PIPELINES/R_*/*/*.MaxGene*/$tt.$iter.txt`)
      echo -n "$ff\t" >> tutu
      head -3 $ff | gawk '/AUC2/{printf("\t%s",$0);}' >> tutu
      head -3 $ff | gawk '/Parameters/{print;}' >> tutu
    end
  end
  cat tutu |  sed -e 's/OTHER_PIPELINES\/R_//g' -e 's/_genes_Magic\///g' -e 's/Shift1\.//g' -e 's/MaxGene/Mx/g' -e 's/Fatigue\.//g' -e 's/Random/Rdm/g'  -e 's/\.\.beta\./b/g' -e 's/\.txt//g' -e 's/ChronicFatigueGroup/All/g'| sort > RESULTS/AUC2.$tt.txt
  cat tutu | gawk -F '\t' '/AUC2/{if ($8>58)printf("%s\t%s\n",$1,$8)}' | sort > OTHER_PIPELINES/best.AUC2.$tt.txt
end


cat OTHER_PIPELINES/R_*/AGLuK/R1.0.Shift1.MaxGene200/*gue.FC*beta.1.txt | gawk -F '\t' '/^Plus/{for(i=3;i<=NF;i++)n[$i]++;}END{for(k in n)printf("\t%d\t%s\n",n[k],k);}' | sort -k 1n
cat OTHER_PIPELINES/R_digital/AceView_genes_Magic/R1.0.Shift1.MaxGene200/*gue.FC*beta.1.txt | gawk -F '\t' '/^Minus/{for(i=3;i<=NF;i++)n[$i]++;}END{for(k in n)printf("\t%d\t%s\n",n[k],k);}' | sort -k 1n

# LOC399942 == mec-12 in worm, is the worst minus gene = overexpressed in FC47.2

###########################################################################################
######## Utility to  try to name the S_seqc genes via the embedded RefSeq models
if (! -e OTHER_PIPELINES/seqc.clone2mrna2geneName.txt) then

setenv ici `pwd`
foreach chrom($chromSetAll)
  tbly ~/NB/tmp/XH$chrom <<EOF
    query find cdna_clone IS _* && from_gene
    bql -a -o $ici/OTHER_PIPELINES/tutu.$chrom.txt select c,g,m from c in @,g in c->from_gene,m in g->mrna
    quit
EOF
end

cat  OTHER_PIPELINES/tutu.*.txt | gawk -F '\t' '{gsub(/\"/,"",$0); n=split($1,aa,"_");if(n!=3)print "ERREUR " $0 ; printf("%s\t%s\t%s\t%s\n",aa[2],$3,$2,$1);}' | sort -u > OTHER_PIPELINES/seqc.clone2mrna2geneName.txt
cat OTHER_PIPELINES/seqc.clone2mrna2geneName.txt | cut -f 1,2 | sort -u > OTHER_PIPELINES/seqc.mrna2geneName.txt
\rm  $ici/OTHER_PIPELINES/tutu.$chrom.txt

endif
##########

#############################################################################################
########### FREQUENCY ANALYSIS OF GENES , removing baddies

####################################
### list the best predictions


foreach endp (B D F)
  if($endp == F) then
    set tt=$MAGIC.HR_Died_TrainingSet_44_HR_Survived_TrainingSet_43...beta
    cat  OTHER_PIPELINES/best.AUC2.$tt.txt | gawk -F '\t' '{if ($2>58.5)print $1;}' > totook
  endif
  if($endp == D) then
    set tt=$MAGIC.Favorable_TrainingSet_87_Unfavorable_TrainingSet_45...beta
    cat  OTHER_PIPELINES/best.AUC2.$tt.txt | gawk -F '\t' '{if ($2>98)print $1;}' > totook
  endif
  if($endp == B) then
    set tt=$MAGIC.OS_Died_TrainingSet_OS_Survived_TrainingSet...beta 
    cat  OTHER_PIPELINES/best.AUC2.$tt.txt | gawk -F '\t' '{if ($2>88.5)print $1;}' > totook
  endif

\rm totook2*
foreach ff (`cat totook`)
  cat $ff | cut -f 1,2 | gawk '/^Genes/{zz++;next;}{if(zz>0 && $1)print $1}' >> totook2
  cat $ff | gawk '/^Plus gene/{for(i=9;i<=NF;i++)print $i}' >> totook2p
  cat $ff | gawk '/^Minus gene/{for(i=9;i<=NF;i++)print $i}' >> totook2m

end

cat totook2p | gawk '{gene=$1;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}} print gene ;}' > totook3p
cat totook2m | gawk '{gene=$1;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}} print gene ;}' > totook3m

cat OTHER_PIPELINES/seqc.mrna2geneName.txt ZZZZZ totook3p | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$2]=gg[$2] "_" $1;next;}}{g=$1;if(substr(g,1,2)=="S_"){gsub(/^S_/,"",g);if(gg[g])g=substr(gg[g],2);else{split(g,aa,".");g=aa[1];}}nn[g]++;}END{for (g in nn)printf("%d\t%s\n",nn[g],g);}' | sort -k 1nr | gawk '{printf("%s\t%s\t%d\n",$2,endp,$1);}' endp=$endp  > toto1
cat  MicroArray/Hits/AGLuK.av.probe2gene.quasi_unique.*txt ZZZZZ toto1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[substr($1,1)]=$2;next;}}{g=$1;if(gg[g])g=gg[g];printf("%s\t%s\t%d\n",g,$2,$3);}' > totook4p.$endp

cat OTHER_PIPELINES/seqc.mrna2geneName.txt ZZZZZ totook3m | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$2]=gg[$2] "_" $1;next;}}{g=$1;if(substr(g,1,2)=="S_"){gsub(/^S_/,"",g);if(gg[g])g=substr(gg[g],2);else{split(g,aa,".");g=aa[1];}}nn[g]++;}END{for (g in nn)printf("%d\t%s\n",nn[g],g);}' | sort -k 1nr | gawk '{printf("%s\t%s\t%d\n",$2,endp,$1);}' endp=$endp   > toto1
cat MicroArray/Hits/AGLuK.av.probe2gene.quasi_unique.*txt ZZZZZ toto1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[substr($1,1)]=$2;next;}}{g=$1;if(gg[g])g=gg[g];printf("%s\t%s\t%d\n",g,$2,$3);}' > totook4m.$endp


# enpd B D F
end 
set toto=totook5pm
echo -n '# '> $toto
date >> $toto
echo "# Number of occurences of the genes in the signatures of the best MAGIC ene/transcript/agilent predictors" >> $toto
echo "# Agilent probes, RefSeq NM and transcripts are uniformly remapped to their gene indentifier" >> $toto
echo "# Gene\tEndpoint\tA\tA+\tA-\tB\tB+\tB-\tC\tC+\tC-\tD\tD+\tD-\tE\tE+\tE-\tF\tF+\tF-" >> $toto
cat  RESULTS/NB_all_Models/toto3 ZZZZZ totook4p.? ZZZZZ totook4m.?  | gawk -F '\t' '/^ZZZZZ/{zz++;}{g=$1;if(length(g)>1){nn[g]+=$3;nt[1+zz,g,$2]+=$3;}}END{for(g in nn){printf("%s\t%d",g,nn[g]);for(i=1;i<=6;i++)printf("\t%d\t%d\t%d",nt[1,g,substr("ABCDEF",i,1)],nt[2,g,substr("ABCDEF",i,1)],nt[3,g,substr("ABCDEF",i,1)]);printf("\n");}}' | sort -k 2,2nr -k 1,1 >> $toto
\cp $toto RESULTS/BestGenePredictors.magic.txt




########## 2013_06_29   analysis of most frequent genes in best predictors from all NB teams

pushd RESULTS/NB_all_Models

# coount the size of the signatures

set justIntrons=1
set toto3=toto3.$justIntrons
set totoC=totoC.$justIntrons
# set toto3=toto3.p
\rm $toto3
\rm $totoC
foreach endp (A B C D E F)
  foreach model (`cat GoodModels.txt | grep _$endp'_'`)
    if (-e toto2) \rm toto2
    set lab=`echo $model | gawk '{split($1,aa,"_");print aa[4];}'`
    set annot=`echo $model | gawk '{split($1,aa,"_");print aa[5];}'`
    set GM=`echo $model | gawk '{split($1,aa,"_");print aa[6];}'`
    set mn=`echo $model | gawk '{split($1,aa,"_");print aa[8];}'`

    if ($justIntrons >= 0 && $GM != J) continue
    if (0 && $annot != AG1) continue
    if (0 && ($lab != SAS || $annot != TUC || $endp != A)) continue
  if (0 && $lab != FBK) continue

    echo -n "\n#################$endp lab=$lab $annot GM=$GM"

    set sigFile=`ls  SEQC_NB_ListsOfSignatureGenes_$lab'_'*.txt`

    set nLine=`grep -n '@' $sigFile | grep '_'$endp'_' | grep _$annot'_' | grep _$GM'_' | gawk '{split($1,aa,":");printf("%d",0+aa[1]);}'` 
    echo " $sigFile $nLine"
    echo "NA" > toto1
    if ($nLine >= 1) then
      cat $sigFile | tail -n +$nLine | sed -e 's/\r//g' | sed -e 's/\r//g'  > toto1
    endif

    echo -n "$endp\t$lab\t$annot\t$GM" >> $totoC
    cat toto1  | gawk '/^NA$/{next;}/^@/{na++;if(na==2){exit;}next;}{if(length($1)>1)n++;}END{printf("\t%d",0+n);}' >> $totoC
    echo  "\t$model\t$mn" >> $totoC
    if ($nLine < 1) continue 

    set ok=1
    set GM1=$GM
    if ($annot == AG1) then
      set GM1=P
      set ok=0
      cat toto1 | gawk '/^@/{na++;if(na==2)exit;}{g=$1;split(g,aa,"|");g=aa[1];split(g,aa,":");g=aa[1];gsub(/^AG_/,"",g);gsub(":","",g);print g;}' > toto1a
      cat  ~/NB_2013/MicroArray/Hits/AGLuK.av.probe2gene.quasi_unique.*txt  ZZZZZ toto1a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[substr($1,6)]=$2;next;}}{g=$1;if(gg[g])g=gg[g];print g;}' >> toto2
    endif
    if ($ok == 1 && ($GM == G ||  $GM == T)) then 
      set ok=0
      cat toto1  | gawk '/^@/{na++;if(na==2)exit;}{g=$1;split(g,aa,".");g=aa[1]; gsub(/[:_]Gene_AceView/,"",g);gsub(/[:_]Gene_RefSeq/,"",g); print g;}' > toto1a
      cat  ~/37_5/TARGET/MRNAS/hs.NM2gene.txt ZZZZZ toto1a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[substr($1,1)]=$4;next;}}{g=$1;if(gg[g])g=gg[g];print g;}' >> toto2
    endif
    if ($ok == 1 && $GM == J && $justIntrons == 0) then 
      set ok=0
      cat toto1 | sed -e '/_random/random/g'  | gawk '{i=index($1,"|NT:-");if(i>1){split($1,aa,":");$1=aa[1] ":" aa[3] ":" aa[3];}print $1;}' | gawk '/^@/{na++;if(na==2)exit;}{g=$1;gsub(/_/,":",g);gsub(/-/,":",g);split(g,aa,":");c=aa[1];gsub(/chr/,"",c); a1=aa[2];a2=aa[3];if(a2>1)print c "__" a1 "_" a2 ;}' > toto1a
  cat  NB/intron2gene.txt ZZZZZ toto1a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gsub(/\"/,"",$0);g=$2;gg[$1]=$2;split($1,aa,"_");c=aa[1];a1=10*int((aa[3]+5)/10);a2=10*int((aa[4]+5)/10);u=c "__" a1 ; gg[u]=g;u=c "__" a2 ; gg[u]=g;next;}}{g=$1;if(gg[g]){g=gg[g];}else{split(g,aa,"_");c=aa[1];a1=10*int((aa[3]+5)/10);a2=10*int((aa[4]+5)/10);u=c "__" a1 ; if(gg[u])g=gg[u];else{u=c "__" a2 ;if(gg[u])g=gg[u];else{u=c "__" a1+10 ; if(gg[u])g=gg[u];else{u=c "__" a1-10 ; if(gg[u])g=gg[u];else{u=c "__" a2+10 ; if(gg[u])g=gg[u];else{u=c "__" a2-10 ; if(gg[u])g=gg[u]}}}}}}print g;}' >> toto2
       wc toto1 toto1a toto2
    endif
    if ($ok == 1 && $GM == J && $justIntrons == 1) then 
      # in this option we rationalize the introns names to correct 3__123_456, rather than to gene names
      set ok=0
      cat toto1 | sed -e '/_random/random/g' | gawk '{i=index($1,"|NT:-");if(i>1){split($1,aa,":");$1=aa[1] ":" aa[3] ":" aa[3];}print $1;}' | gawk '/^@/{na++;if(na==2)exit;}{g=$1;gsub(/_/,":",g);gsub(/-/,":",g);split(g,aa,":");c=aa[1];gsub(/chr/,"",c); a1=aa[2];a2=aa[3];if(a2>1)print c "__" a1 "_" a2 ;}' > toto1a
      cat  NB/intron2gene.txt ZZZZZ toto1a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gsub(/\"/,"",$0);g=$2;gg[$1]=$2;split($1,aa,"_");c=aa[1];a1=10*int((aa[3]+5)/10);a2=10*int((aa[4]+5)/10);u=c "__" a1 ; gg[u]=g;u=c "__" a2 ; gg[u]=g;next;}}{g=$1;if(gg[g]){g=gg[g];}else{split(g,aa,"_");c=aa[1];a1=10*int((aa[3]+5)/10);a2=10*int((aa[4]+5)/10);u=c "__" a1 ; if(gg[u])g=gg[u];else{u=c "__" a2 ;if(gg[u])g=gg[u];else{u=c "__" a1+10 ; if(gg[u])g=gg[u];else{u=c "__" a1-10 ; if(gg[u])g=gg[u];else{u=c "__" a2+10 ; if(gg[u])g=gg[u];else{u=c "__" a2-10 ; if(gg[u])g=gg[u]}}}}}}print g;}' >> toto2
      wc toto1 toto1a toto2
    endif

    cat toto2 | gawk '/^@/{next;}{nn[$1]++}END{for(g in nn) printf("%s\t%s%s\t%d\n",g,endp,GM,nn[g]);}' GM=$GM1 endp=$endp >> $toto3
  end

end

cat  $totoC | gawk -F '\t' '{ printf("%d\t",$7-$5);print;}' | sort -k 1n > NB/RESULTS/NB.all_labs.counts.txt

if ($justIntrons == 1) then 
  set toto=toto4i
  echo -n '# '> $toto
  date >> $toto 
  echo "# Number of occurences of the exon junctions in the signatures of the best exon-junctions predictors" >> $toto
  echo "# Each junction is attributed to a gene and lists its parent annotations" >> $toto
  echo '# Coordinates are oriented, they represent the first (usually the G of  Gt..ag) and last (gt...aG) base of the intron" >> $toto
  echo "# Junction\tGene\tAnnotated in\tAny endpoint\tA\tB\tC\tD\tE\tF" >> $toto
  cat $toto3 | gawk -F '\t' '{g=$1;if(length(g)>1){nn[g]+=$3;nt[g,$2]+=$3;}}END{for(g in nn){printf("%s\t%d",g,nn[g]);for(i=1;i<=6;i++)for(j=3;j<=3;j++)printf("\t%d",nt[g,substr("ABCDEF",i,1) substr("GTJP",j,1)]);printf("\n");}}' | sort -k 2,2nr -k 1,1 > toto4ia

  cat  NB/intron2gene.txt ZZZZZ NB/totox.RefSeq.txt ZZZZZ NB/totox.av.txt ZZZZZ  NB/totox.encode.txt ZZZZZ NB/totox.seqc.txt NB/totox.introns.txt ZZZZZ toto4ia | gawk -F'\t' 'BEGIN{nam[1]="RefSeq,";nam[2]="AceView,";nam[3]="Encode,";nam[4]="Discovery,";}/^ZZZZZ/{zz++;next;}{if(zz<1){gsub(/\"/,"",$0);gg[$1]=$2;next;}}{if(zz<1){gsub(/\"/,"",$0);av[$1]=1;next;}}{if(zz<5){annot[$1,zz]=1;av[$1]=1;next;}}{g=$1;split(g,aa,"_");printf("chr%s:%s-%s\t%s\t",aa[1],aa[3],aa[4],gg[g]);ok=1;for(i=1;i<5;i++)if(annot[g,i]==1 && (i<4 || ok==1)){if(i<3)ok=0;printf("%s",nam[i]);}for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto
   
  \cp $toto NB/RESULTS/BestIntronPredictors.all_teams.good_modelstxt
else
  set toto=toto4
  echo -n '# '> $toto
  date >> $toto 
  echo "# Number of occurences of the genes in the signatures of the best gene/transcript/exon-junctions predictors" >> $toto
  echo "# Agilent probes, RefSeq NM, transcripts and junctions are uniformly remapped to their gene identifier" >> $toto
  echo "# Gene\tGeneId\tLocus\tChromosome\t5'\t3'\tAny endpoint\tA\tB\tC\tD\tE\tF\t\tAG\tAT\tAJ\tAP\tBG\tBT\tBJ\tBP\tCG\tCT\tCJ\tCP\tDG\tDT\tDJ\tDP\tEG\tET\tEJ\tEP\tFG\tFT\tFJ\tFP" >> $toto
  cat NB/TARGET/GENES/*.gene2intMap.txt ZZZZZ NB/TARGET/GENES/*.gene2geneid.txt ZZZZZ $toto3 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gsub(/\"/,"",$0);g2intmap[$1]="chr" $2 ":" $3 "-" $4;g2c[$1]=$2;g2a1[$1]=$3;g2a2[$1]=$4;next;}}{if(zz<2){gsub(/\"/,"",$0);g2gid[$1]=g2gid[$1] "," $2;next;}}{g=$1;if(length(g)>1){nn[g]+=$3;nt[g,substr($2,1,1)]+=$3;nt[g,$2]+=$3;}}END{for(g in nn){printf("%s\t%s\t%s\t%s\t%d\t%d\t%d",g,substr(g2gid[g],2),g2intmap[g],g2c[g],g2a1[g],g2a2[g],nn[g]);for(i=1;i<=6;i++)printf("\t%d",nt[g,substr("ABCDEF",i,1)]);printf("\t");for(i=1;i<=6;i++)for(j=1;j<=4;j++)printf("\t%d",nt[g,substr("ABCDEF",i,1) substr("GTJP",j,1)]);printf("\n");}}' | sort -k 7,7nr -k 1,1 >> $toto
  \cp $toto NB/RESULTS/BestGenePredictors.all_teams.best_models.txt
endif

######################


grep Minus RESULTS/Expression/quasi_unique/RefSeq/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )next;for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 > RESULTS/Fatigue.RefSeq.baddy.txt
grep Plus RESULTS/Expression/quasi_unique/RefSeq/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )next;for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 >> RESULTS/Fatigue.RefSeq.baddy.txt


grep Plus RESULTS/Expression/quasi_unique/AceView/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )next;for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 > RESULTS/Fatigue.AceView.baddy.txt
grep Minus RESULTS/Expression/quasi_unique/AceView/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )next;for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 >> RESULTS/Fatigue.AceView.baddy.txt


grep 'us genes:' RESULTS/Expression/quasi_unique/RefSeq/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 > RESULTS/Fatigue.RefSeq.maybe.txt
grep 'us genes:'  RESULTS/Expression/quasi_unique/AceView/*LibBa*beta.0.txt | gawk '/Group/{if(index($0,"Group1")>0 && index($0,"Group2")>0 )for(i=1;i<=NF;i++)n[$i]++;}END{for(z in n)printf("%s\t%d\n",z,n[z]);}' | sort -k 1,1 > RESULTS/Fatigue.AceView.maybe.txt


cat RESULTS/Fatigue.AceView.maybe.txt ZZZZZ RESULTS/Fatigue.AceView.baddy.txt | gawk '/^ZZZZZ/{zz++;next;}{n[$1]+=zz+1;if(0)print "UU" zz "+" $1 "+" n[$1]}END{for(k in n)if(n[k]==1)print k,n[k];}' | sort > RESULTS/Fatigue.AceView.goody.txt
 cat RESULTS/Fatigue.RefSeq.maybe.txt ZZZZZ RESULTS/Fatigue.RefSeq.baddy.txt | gawk '/^ZZZZZ/{zz++;next;}{n[$1]+=zz+1;if(0)print "UU" zz "+" $1 "+" n[$1]}END{for(k in n)if(n[k]==1)print k,n[k];}' | sort > RESULTS/Fatigue.RefSeq.goody.txt

set toto=RESULTS/$MAGIC.candidate_genes.txt
echo > $toto.1
foreach uu (u nu)
  set qu=""
  if ($uu == nu) set qu="quasi_"
  foreach target (RefSeq AceView)
    foreach pm (Plus Minus)
      foreach typ (GENE MRNAH)
        foreach beta (0 1 2)
          echo "$target\t$uu\t$pm\t$typ\t$beta  RESULTS/Expression/"$qu"unique/$target/$MAGIC.*.$typ.$uu.FC*.beta.$beta.txt"
          cat  RESULTS/Expression/$qu'unique'/$target/$MAGIC.*.$typ.$uu.FC*.beta.$beta.txt | gawk '/# AUC/{ok=0;if($6>62)ok=1;next;}'"/$pm genes/"'{for(i=5 ; i<=NF;i++)n[$i]++;}END{for(k in n)printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\n",k,n[k],target,uu,pm,typ,beta);}' target=$target beta=$beta uu=$uu pm=$pm typ=$typ >> $toto.1
        end
      end
    end
  end
end

echo "# " > $toto
date >> $toto
echo "# Gene\tInstance\tPM\tType\tDubious\tBad\tGene Total\tsubtotal\tRefSeq u no iter\tRefSeq u iter\tRefSeq u reiter\tRefSeq nu no iter\tRefSeq nu iter*tRefSeq nu reiter\tAceView u no iter\tAceView u iter\tAceView u reiter\tAceView nu no iter\tAceView nu iter\tAceView nu reiter" >> $toto

cat RESULTS/Fatigue.*.baddy.txt ZZZZZ $toto.1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;bad[g]="Baddy";gene=g;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}}dubious[gene]="Dubious";next;}}{g=$1;if(substr(g,1,2)=="FC" && length(g)==7)next;if(0+g>0)next;gg[g]+=$2;tt[$6]=1;ggt[g,$5,$6]+=$2;n[g,$3,$4,$5,$6,0+$7]=$2;gene=g;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}}ngene[gene]+=$2;g2gene[g]=gene;}END{for(pmi=0;pmi<2;pmi++){pm="Plus";if(pmi==1)pm="Minus";for(g in gg)for (t in tt)if(ggt[g,pm,t]>0){printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d",g2gene[g],g,pm,t,dubious[g2gene[g]],bad[g],ngene[g2gene[g]],ggt[g,pm,t]);for (ira=0;ira<2;ira++){ra="RefSeq";if(ira>0)ra="AceView";{for(uui=0;uui<2;uui++){uu="u";if(uui>0)uu="nu"; for(iter=0;iter<3;iter++)printf("\t%d",n[g,ra,uu,pm,t,iter]);}}}printf("\n");}}}' > $toto.2

cat $toto.2 |  grep -v Dubious | sort -k 5nr >> $toto
cat $toto.2 |  grep  Dubious | sort -k 5nr  >> $toto


# selected histos
set toto=RESULTS/RefSeq.candidates.histo.txt 
date > $toto
foreach ff ( `ls RESULTS/Expression/*unique/RefSeq/Differential_genes/Gr1_2.RefSeq.GENE.*.1.txt`)
  head -3  $ff | gawk '{print}' >> $toto
  cat RESULTS/RefSeqGenes.txt ZZZZZ $ff | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){ok["X__" $1]=1;next;}if(ok[$1])print}' >>  $toto
end
cat $toto | cut -f 1,2,3,4,5,6,7,8

cat MetaDB/Fatigue/runs.ace  | gawk '/^Runs/{next;}/^Run/{r=$2;next;}/^Sample/{printf("%s\t%s\n",r,$2);}' | sed -e 's/\"//g' | gzip > MetaDB/$MAGIC/run2sample.ace.gz

set toto=RESULTS/RefSeq.candidates.index.txt 
echo -n "# " > $toto
date >> $toto
\rm $toto.1
foreach ff ( `ls RESULTS/Expression/*unique/RefSeq/Gr1_2.RefSeq.GENE.*.ace.gz`)
  gunzip -c  MetaDB/$MAGIC/run2sample.ace.gz ZZZZZ.gz RESULTS/RefSeqGenes.txt.gz ZZZZZ.gz $ff | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){r2s[$1]=$2;next;}}{if(zz<2){ok["X__" $1]=1;next;}}{gsub(/\"/,"",$0);}/^Gene/{g=$2;gok=ok[g];next;}/Run_nU/{if($11=="NAw" || gok<1)next;r=$2;rr[r]=1;gsub("NA/","",$3);ngr[g,r]=$3;}END{printf("# Gene");for (ir = 3991; ir<=4190;ir++) {r="Rhs"ir;printf("\t%s",r);}printf("\n# Sample");for (ir = 3991; ir<=4190;ir++) {r="Rhs"ir;printf("\t%s",r2s[r]);}for(g in ok)if(ok[g]==1){printf("\n%s",g);for (ir = 3991; ir<=4190;ir++){r="Rhs"ir;printf("\t%s",ngr[g,r]);}}printf("\n");}' >>  $toto.1
end
head -1 $toto.1 >> $toto
tail -n +2 $toto.1 | sort >> $toto
cat $toto | cut -f 1,2,3,4,5,6,7,8


 cat MetaDB/ru55555555555e ZZZZZ RESULTS/Expression/quasi_unique/AceView/Fatigue.av.GENE.*.FC*...beta.0.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){r2s[$1]=$2;next;}}/^Run:/{for(i=3;i<=NF;i++)nn[$i]+=i-103;}END{for(r in nn)if(substr(r,1,3)=="Rhs")printf("%d\t%s\t%s\n",nn[r],r,r2s[r]);}' | sort -k 1n > RESULTS/Fatigue.sorted_runs.txt


#########
