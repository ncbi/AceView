#!bin/tcsh -f

echo ZZZZZ > ZZZZZ
set R=RESULTS/HISTOS


# method D3 on AceView

set tt=RefSeq/NB.RefSeq
set tt=AceView/NB.av

set pass=99
set pass=1

set dd=Mix4.nu.min10.max50.NA1
set dd=Mix4.nu.min10.max50.NA1
set dd=Dan73wide
set dd=D9
set dd=D10_3
set dd=D12d
set dd=D12f
set dd=D12m
set dd=D14b

set ss0=RESULTS/Expression.$dd/quasi_unique/$tt.nu.HR_Died_TrainingSet_43_HR_Survived_TrainingSet_42.alpha.pass$pass.txt

set ss0=OTHER_PIPELINES/Results/RefSeq_genes_BGI/R2.0.Shift1.MaxGene300/NB.NB_female_NB_male.alpha.pass1.txt
echo "#Plus" > $R/$dd.list
cat $ss0 | grep Plus | gawk '{for(i=3;i<=NF;i++)print $i}' >> $R/$dd.list
echo "#Minus" >> $R.$dd.list
cat $ss0 | grep Minus | gawk '{for(i=3;i<=NF;i++)print $i}' >> $R/$dd.list

set ss1=RESULTS/Expression.$dd/quasi_unique/$tt.nu.histo.HR_Survived_TrainingSet_42.HR_Died_TrainingSet_43.pass$pass.txt 
set ss2=RESULTS/Expression.$dd/quasi_unique/$tt.nu.histo.HR_Died_TrainingSet_43.HR_Survived_TrainingSet_42.pass$pass.txt

set ss1=OTHER_PIPELINES/Results/RefSeq_genes_BGI/R2.0.Shift1.MaxGene300/NB.histo.NB_female.NB_male.pass1.txt
set ss2=OTHER_PIPELINES/Results/RefSeq_genes_BGI/R2.0.Shift1.MaxGene300/NB.histo.NB_male.NB_female.pass1.txt


cat  $R/$dd.list ZZZZZ $ss1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=1;i<=NF;i++)ok[$i]=1;next;}if(ok[$1]){printf("%s\t%s\t%s\t%s",$1,$13,$17,$18);for (i=25;i<=NF;i++)printf("\t%s",$i);printf("\n");}}'  > $R/$dd.f.hr.histo.txt
cat  $R/$dd.list ZZZZZ $ss2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=1;i<=NF;i++)ok[$i]=1;next;}if(ok[$1]){printf("%s\t%s\t%s\t%s",$1,$13,$17,$18);for (i=25;i<=NF;i++)printf("\t%s",$i);printf("\n");}}'  > $R/$dd.u.hr.histo.txt

cat $R/$dd.f.hr.histo.txt | gawk -f $R/histo_select.awk NA=NA uf=f > $R/_toto_hr
cat $R/$dd.u.hr.histo.txt | gawk -f $R/histo_select.awk NA=NA uf=u >> $R/_toto_hr

cat $R/_toto_hr | sort -u > $R/$dd.HISTO.txt
\rm $R/_*

cat  $R/$dd.HISTO.txt | head -122 | tail -120 | gawk '/NA/{print $2}' | sort -u 
cat $ss1 | gawk '/^#/{next;}{n++;if(n%3==0)print}' | cut -f 1,18 | head -20

#######

cat $R/tutu | gawk '{for (i=1 ; i<=NF ; i++) print $i}' > genePlus.list
cat $R/tutu | gawk '{for (i=1 ; i<=NF ; i++) print $i}' > geneMinus.list

# set ss1=RESULTS/Expression.D12m/quasi_unique/AceView/NB.av.nu.histo.OS_Survived_TrainingSet.OS_Died_TrainingSet.pass1.txt
# set ss2=RESULTS/Expression.D12m/quasi_unique/AceView/NB.av.nu.histo.OS_Died_TrainingSet.OS_Survived_TrainingSet.pass1.txt

cat $ss1 |  gawk -F '\t' '/^#/{next;}{printf("%s\t%s\t%s\t%s",$1,$13,$17,$18);for (i=25;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $ss1.x
cat $ss2 |  gawk -F '\t' '/^#/{next;}{printf("%s\t%s\t%s\t%s",$1,$13,$17,$18);for (i=25;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $ss2.x

cat $ss1.x  | gawk -f $R/histo_select.awk NA=NA uf=f  | sort -k 2nr  > $ss1.y
cat $ss2.x  | gawk -f $R/histo_select.awk NA=NA uf=u  | sort -k 2nr  > $ss2.y

cat $ss1.y | grep MMM | sort -k 3nr | head -30 | cut -f 2 >  geneMinus.list
cat $ss2.y | grep MMM | sort -k 3nr | head -30 | cut -f 2 >  genePlus.list

cat  genePlus.list geneMinus.list > $R/$dd.list

# cat $ss1.y $ss2.y | grep MMM | sort -k 3nr | head -40  | cut -f 2 >  genePlus.list
# cat $ss1.y $ss2.y | grep MMM | sort -k 3nr | head -40  | cut -f 2 >  geneMinus.list


set target=av
set uu=nu

cat genePlus.list geneMinus.list

bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.$target.$uu.with_index.ace.gz  -$uu  -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.info.ace  -o tmp/GENEINDEX/Results/$MAGIC.$target.$uu -gzo  -pA -method Gene_AceView -correl -compare  -genePlus genePlus.list -geneMinus geneMinus.list -keepIndex

cat runs_MYCN_over_16.5.list ZZZZZ tmp/GENEINDEX/Results/$MAGIC.$target.$uu.sampleClassificationBySignature.txt |  gawk -F '\t' -f scripts/NB.roc.awk justLowMYCN=0 parity=1 >  tmp/GENEINDEX/Results/$MAGIC.$target.u.outcome.txt

# construct a ROC curve
scripts/NB.roc.tcsh $target 1 .
#cat _roc.$target | grep MCC
cat _roc.$target | grep uf
head -4 _km

