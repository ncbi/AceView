#!bin/tcsh -f

if ($1 == FDR1) goto phaseFDR1
if ($1 == FDREXPORT) goto phaseFDREXPORT
if ($1 == GetFDR) goto phaseGetFDR

goto phaseLoop

###################################################################################
### recover the FDR and counts in known cases from MetaDB

phaseGetFDR:

# grab the DEG experiments from the database
bin/tacembly MetaDB <<EOF
  query find Compare COUNT Runs == 2 AND project==$MAGIC
  bql -a -o MetaDB/$MAGIC/DEG_high_low.pretxt -active  select e,r1,r2,u1,u2,type,sense,target,n,threshold from e in @, r1 in e->runs, r2 in e->runs , u1 in r1[1], u2 in r2[1] where u1 < u2, type in e->Counts, sense in type[1], target in sense[1], threshold in target[1],n in threshold[4]
  quit
EOF

cat  MetaDB/$MAGIC/DEG_high_low.pretxt | sort | sed -e 's/\"//g' -e 's/System://' -e 's/Target://'  | gawk -F '\t' '{print $1 "::" $2 "::" $3 "::" $4 "::" $5 "::" $6 "::" $7 "::" $8 "::" $9 "::" $10  ;}'> MetaDB/$MAGIC/DEG_high_low.txt
\rm   MetaDB/$MAGIC/DEG_high_low.pretxt 

goto phaseLoop

##################################################################

phaseFDR1:

scripts/FDR.tcsh GetFDR

echo "FDR.tcsh FDR1 $MAGIC"
echo " Grab the thresholds and number of DEGs from the DEG_FDR_selection  files and parse it into MetaDB"
set oldex=""
echo -n "// " > RESULTS/Expression/$MAGIC.DEG_threshold.ace
date >> RESULTS/Expression/$MAGIC.DEG_threshold.ace
echo -n "// " > RESULTS/Expression/$MAGIC.DEG_AUC.ace
date >> RESULTS/Expression/$MAGIC.DEG_AUC.ace


set nam='2016'
set nam='var'
set nam=''
set nam='_jan29_Seuil0'       # threshold=0 in scoreHistoss, new Agilent norm
set nam='_jan31_SeuilX0'      # 2015_01_31  threshold crossOver + 0 0 in scoreHistos 2015_01_31s, new Agilent norm
set nam='_jan31_Gini'         # 2015_01_31 Gini threshold crossOver + 0 0 in scoreHistos, new Agilent norm
set nam='_feb1_oldNorm_Gini'  # 2015_02_01 Gini, flipped back to old Agilent norm
set nam='_feb1_GiniXm1'       # 2015_02_01 Gini, flipped back to old Agilent norm, seuil crossOver - 1  
set nam='_feb2_X0'            # 2015_02_01 Gini, flipped back to old Agilent norm, seuil crossOver - 1  
set nam='_feb6_X0_weighted'   # 2015_02_06 magic weighted score, 
set nam='_feb10_ali2013'      # 2015_02_10, idem mais ali avec code 2013
set nam='_may19_2015'         # 2015_05_19, we fixed 2 erroneous swapped samples, ali selon code 2013
set nam='_may21_2015_bb'      # 2015_05_19, we fixed 2 erroneous swapped samples, ali selon code 2013, (baddy batchies) subtracted intrinsic noise measured as the average differential score of each gene in 40 random resampling, from the test differential score, beore running the prediction module. This intrinsic noise often corresponds to batch effects or to imbalance of hidden characters, like sex. 

set nam='2016'
set nam='2016_06_05_12000_groups'
set nam='2016_11_29'
if ($?FDR_nam) then
  set nam=$FDR_nam
endif


foreach GM (Gene mRNAH INTRON)
  cat MetaDB/$MAGIC/DEG_high_low.txt  | gawk '{split($1,aa,"::");print aa[1];}' | sort -u | gawk '{printf("Compare %s\n-D Counts Diff_%s \"1>2%s\"\n-D Counts Diff_%s \"2>1%s\"\n\n",$1,GM,nam,GM,nam);}'  GM=$GM nam=$nam  >> RESULTS/Expression/$MAGIC.DEG_threshold.ace
end

foreach target (`ls -d RESULTS/Expression/unique/* | gawk '{print substr($1,27)}' | sort -r `)
  foreach GM (GENE MRNAH MA SNP INTRON SF24A SF24B SF31A SF31B)
   if (0 && $GM != MRNAH && $GM != GENE) continue 
   if (0 && $GM != GENE) continue 
    set target2=$target
    if ($target == snp) set target2=AceView
    if ($target == av) set target2=AceView
     
    foreach mm (`cat MetaDB/$MAGIC/DEG_high_low.txt  | gawk '{split($1,aa,"::");printf("%s::%s::%s::%s::%s\n",aa[1],aa[2],aa[3],aa[4],aa[5]);}'   | sort -u `)
      set ex=`echo "$mm" | gawk '{split($1,aa,"::");print aa[1];}'`
      set r1=`echo "$mm" | gawk '{split($1,aa,"::");print aa[2];}'`
      set r2=`echo "$mm" | gawk '{split($1,aa,"::");print aa[3];}'`
      set u1=`echo "$mm" | gawk '{split($1,aa,"::");print aa[4];}'`
      set u2=`echo "$mm" | gawk '{split($1,aa,"::");print aa[5];}'`

      set ff=tmp/GENEINDEX/Results/$MAGIC.$target2.$GM.u.$ex.$r1.$r2.DEG_FDR_selection.txt
      if (-e $ff) then
        cat $ff | gawk -F '\t' '{if($2 == r1 && $3 == r2 && ($4 == "Threshold" || $5 == "Selected threshold")){deg="DEG";if(GM=="MRNAH")deg="DET";if(GM=="MA")deg="DEP";if(GM==SNP)deg="DES";if(substr(target,1,5)=="RNA_at")deg="DET";printf("Compare %s\nDiff_%s \"%s>%s%s\"  %s %s FDR  %s %s %s\n\n", ex, GM,u1,u2,nam,target,$5, $9,deg, $7);}}' ex=$ex GM=$GM r1=$r1 r2=$r2 u1=$u1 u2=$u2 target=$target2  nam=$nam >> RESULTS/Expression/$MAGIC.DEG_threshold.ace
        #\rm $ff
      endif
      set ff=tmp/GENEINDEX/Results/$MAGIC.$target2.$GM.u.$ex.$r2.$r1.DEG_FDR_selection.txt
      if (-e $ff) then
        cat $ff | gawk -F '\t' '{if($2 == r1 && $3 == r2 && ($4 == "Threshold"|| $5 == "Selected threshold")){deg="DEG";if(GM=="MRNAH")deg="DET";if(GM=="MA")deg="DEP";if(GM==SNP)deg="DES";if(substr(target,1,5)=="RNA_at")deg="DET";printf("Compare %s\nDiff_%s \"%s>%s%s\"  %s %s FDR  %s %s %s\n\n", ex, GM,u1,u2,nam,target,$5, $9,deg, $7);}}' ex=$ex GM=$GM r1=$r2 r2=$r1 u1=$u2 u2=$u1 target=$target2  nam=$nam >> RESULTS/Expression/$MAGIC.DEG_threshold.ace
        #\rm $ff
      endif
      set ff=tmp/GENEINDEX/Results/$MAGIC.$target2.$GM.u.$ex.beta.0.txt
      if (-e $ff) then
        cat $ff | head -20 | gawk -F '\t' '/^##AUC:/{auc=$3;printf("Compare %s\n-D AUC %s%s %s\n\n",ex,GM,nam,target);if(auc>0)printf("Compare %s\nAUC %s.%s %s %.1f\n\n",ex,GM,nam,target,auc);next;}/^##AUC2:/{auc2=$5;if(auc>0 && auc2>0)printf("Compare %s\nAUC %s.%s %s %.1f AUC2 %.1f\n\n",ex,GM,nam,target,auc,auc2);auc=0;auc2=0;}'   ex=$ex GM=$GM target=$target2 nam=$nam >> RESULTS/Expression/$MAGIC.DEG_AUC.ace
        #\rm $ff
      endif
    end
  end
end

ls -ls RESULTS/Expression/$MAGIC.DEG_*.ace
bin/tacembly MetaDB <<EOF
  read-models
  pparse  RESULTS/Expression/$MAGIC.DEG_threshold.ace
  pparse  RESULTS/Expression/$MAGIC.DEG_AUC.ace
  save
  quit
EOF

# rexport the Compare class, now containing the updated thresholds
scripts/FDR.tcsh GetFDR

set titi=RESULTS/Expression/$MAGIC.Compare.score.FDR_selection.$nam.txt
touch $titi.1
\rm $titi.*

echo -n '# ' > $titi
date >> $titi
echo "Histogram of the number of differentially expressed gene per score in different conditions" >> $titi
echo "For each experiment, and each method: RNA_seq AceView, RefSeq, EBI, genes or transcripts ; Agilent microarray, RNA-seq at probe locations ; SNPs." >> $titi
echo "the histogram of the number of genes with a given differential score is presented in 2 columns per comparison: A>B and B>A." >> $titi
echo "As a control, the noise is evaluated in the same conditions, by randomly resampling 40 times the same samples into 2 random classes alpha and beta," >> $titi
echo "of size A and B, each containg the same proportion of A and B, so the alpha-beta split is orthogonal to the A-B split." >> $titi
echo "The control histogram corresponds to the average of the 40 resamplings." >> $titi
echo "The false discovery rate (FDR) in each bin is estimated by dividing the number of controls by the number of observed differential features." >> $titi


echo Score | gawk '{printf("\nScore");for(i=0;i<=200;i++)printf("\t%d",i);printf("\n");}' >> $titi.0
foreach mm (`cat MetaDB/$MAGIC/DEG_high_low.txt  | gawk '{split($1,aa,"::");printf("%s::%s::%s\n",aa[1],aa[2],aa[3]);}' | sort -u`)


  set ex=`echo "$mm" | gawk '{split($1,aa,"::");print aa[1];}'`
  set r1=`echo "$mm" | gawk '{split($1,aa,"::");print aa[2];}'`
  set r2=`echo "$mm" | gawk '{split($1,aa,"::");print aa[3];}'`

  cat $titi.0 >> $titi.1
  cat $titi.0 >> $titi.2
  cat $titi.0 >> $titi.3

  foreach target (`ls -d RESULTS/Expression/unique/* | gawk '{print substr($1,27)}'`)
    foreach GM (GENE MRNAH MA snp)
      
      set target2=$target
      if ($target == snp) set target2=AceView
      if ($target == av) set target2=AceView

      set ff=tmp/GENEINDEX/Results/$MAGIC.$target2.$GM.u.$ex.$r1.$r2.DEG_FDR_selection.txt
      if (! -e $ff) continue
     
      echo -n "$ex $target2 $GM Experiment " >> $titi.1
      cat $ff |  gawk -F '\t' '{if($2 == r1 && $3 == r2 && $4 == "Experiment"){for(i=5;i<=205;i++)printf("\t%s",$i);printf("\n");}}' r1=$r1 r2=$r2 >> $titi.1        
      echo -n "$ex $target2 $GM Resampling " >> $titi.1
      cat $ff |  gawk -F '\t' '{if($2 == r1 && $3 == r2 && match($4,"Resampling")){for(i=5;i<=205;i++)printf("\t%s",$i);printf("\n");}}' r1=$r1 r2=$r2 >> $titi.1 

      echo -n "$ex $target2 $GM TDR " >> $titi.2
      cat $ff |  gawk -F '\t' '{if($2 == r1 && $3 == r2 && $4 == "TDR"){for(i=5;i<=205;i++)printf("\t%s",$i);printf("\n");}}' r1=$r1 r2=$r2 >> $titi.2

      echo -n "$ex $target2 $GM FDR " >> $titi.3
      cat $ff |  gawk -F '\t' '{if($2 == r1 && $3 == r2 && $4 == "FDR"){for(i=5;i<=205;i++)printf("\t%s",$i);printf("\n");}}' r1=$r1 r2=$r2 >> $titi.3

    end
  end
end
cat   $titi.1  $titi.2 $titi.3 | scripts/transpose >> $titi

\rm $titi.?
echo $titi
goto phaseLoop

##################################################################
#### export the DEG_high_low ace files | grep AGLuK | grep  -v RNA_at

phaseFDREXPORT:

scripts/FDR.tcsh GetFDR

set target=`echo $Etargets | gawk '{print $1}'`
      set target2=$target
      if ($target == snp) set target2=AceView
      if ($target == av) set target2=AceView

if (! -d RESULTS/Expression/unique/$target/Diff_genes) goto phaseLoop
set tutu=RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.significant_mrna


if (-e $tutu.txt) \rm $tutu.txt
foreach mm (`cat MetaDB/$MAGIC/DEG_high_low.txt  | grep '::av::' | grep Diff_mRNAH `)
  set ff=`echo "$mm" | gawk '{split($1,aa,"::");GM=substr(aa[6],5);target=aa[8];ex=aa[1];r1=aa[2];r2=aa[3];print "RESULTS/Expression/unique/" target "/Diff_genes/" MAGIC "." target ".MRNAH.u." ex "." r1 "_" r2 ".diffGenes.0.txt";}' MAGIC=$MAGIC`
  ls -ls $ff
  set ex=`echo "$mm" | gawk '{split($1,aa,"::");print aa[1];}'`
  set target=`echo "$mm" | gawk '{split($1,aa,"::");print aa[8];}'`
  set th=`echo "$mm" | gawk '{split($1,aa,"::");th=aa[10];if(0+th<10)th=10;print th;}'`
  set GM=`echo "$mm" | gawk '{split($1,aa,"::");GM=substr(aa[6],6);target=aa[8];if(target=="av")print substr(GM,1,4);else if(GM=="MA")print "Probe"; else print "Sequence";}'`
  set ff=RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.$ex.compare.txt 

    cat $ff | gawk -F '\t' '/^#/{next;}{g=$1;score=$2;if(score<th)next;fc=$12;if(fc>0){oldg=g;olds=score;z1=$13;next}z2=$13;fc=$12;if(g==oldg && score==olds && fc<0){nn++;printf("%s\t%s\t%s\t%1f\t%.1f\t%.1f\t%.1f\t%.1f\n", GM,g,ex,nn,score,-fc,z2,z1);}oldg="";z1=-1;}' ex=$ex th=$th GM=$GM >> $tutu.txt

  endif
end

cat $tutu.txt | gawk -F '\t' '{printf("%s %s\nDEG %s %d score %s fc %s %s %s\n\n",$1,$2,$3,$4,$5,$6,$7,$8);}' > $tutu.ace
set n=`wc -l $tutu.ace`
if ($n < 2) \rm $tutu.ace

goto phaseLoop

RESULTS/Expression/unique/RefSeq/Diff_genes/NB.RefSeq.GENE.u.HR_Died_94_HR_Survived_84.diffGenes.1.txt


goto phaseLoop

##################################################################
####

phaseLoop:
 echo done
exit 0

##################################################################
##################################################################

set RESULTS=RESULTS_2013_09

set toto=$RESULTS/Expression/Rarely_not_expressed_genes.txt
gunzip -c $RESULTS/Expression/unique/av/NB.av.GENE.u.expression_index.txt.gz |  head -100 | gawk '/^#/{printf("\t\t\t");print}' > $toto

gunzip -c $RESULTS/Expression/unique/av/NB.av.GENE.u.expression_index.txt.gz | gawk -F '\t'  '/^#Run/{for(i=1;i<=NF;i++){if($i=="_SumOfAllReadsInProject")n1=i;if(match($i,"Rhs"))n2=i;}next;}{k=0;m=0;for(i=n1+1;i<=n2;i++)if(substr($i,1,1)=="N")k++;else m++;if (k>0 && k<=30){printf("%d\t%d\t%.2f",k,m,$n1);for(i=1;i<=n2;i++)printf("\t%s",$i);printf("\n");}}' |  sort -k 1,1n -k 3,3nr > $toto.1

cat $toto.1 | gawk '{nn[$1]++;}END{for(k in nn){printf("Missing_in\t%s\t%d\n",k,nn[k]);}}' | sort -k 2n | gawk -F '\t' '{printf("Missing_in\t%s\t%d\t%d\n",$2,$3,cumul+$3);if($2>=1)cumul+=$3;}' >> $toto
cat $toto.1 | gawk -F '\t' '{if($3<10)next;for(i=16;i<NF;i++)if(substr($i,1,1)=="N"){k[0,i]++;k[$1,i]++;}else m[i]++;nf=NF;}END{printf("# Exp");for(i=2;i<=nf;i++)printf("\t%d",m[i]);for(j=0;j<=15;j++){if (j==0)printf("\n# Not Expressed in up to 30 samples");else printf("\n# Not Expressed in up to %d samples", j);for(i=2;i<=nf;i++){printf("\t%d",k[j,i]+kk[i]);if(j>0)kk[i]+=k[j,i]}}printf("\n");}' >> $toto
cat $toto.1 | gawk '{if($1 != old) printf("\n\n\n");old=$1;if($3 <10)next;print;}' >> $toto







