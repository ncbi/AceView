#!bin/tcsh -f

# 2022_10_1
# Following a discussion with Josh Cherry we look for statistics of substitutions behind any triplet


# cd Fatigue123
set type=all
set type=NotRejected
set type=rejected

set type=$1
set phase=$2


set toto1=RESULTS/SNV/$MAGIC.triplets.$type.txt

if ($type == rejected) then
  set ff=RESULTS/SNV/$MAGIC.'*rejected*'

  if ($MAGIC == Corona) then
    cd ~/CRN
    set ff=RESULTS/SNV/crn.snp_list_with_allele_frequency_per_sample.txt
  endif

endif

if ($type == all) then
  set ff="RESULTS/SNV/$MAGIC.snp_list_with_allele_frequency_per_sample.txt RESULTS/SNV/$MAGIC."'*rejected*'

  if ($MAGIC == Corona) then
    cd ~/CRN
    set ff=RESULTS/SNV/crn.snp_list_with_allele_frequency_per_sample.txt
  endif

endif

if ($type == NotRejected) then

  set ff=RESULTS/SNV/$MAGIC.snp_list_with_allele_frequency_per_sample.txt 

  if ($MAGIC == Corona) then
    cd ~/CRN
    set ff=RESULTS/SNV/crn.snp_list_with_allele_frequency_per_sample.txt
  endif
endif


if ($phase == 1) goto phase1
if ($phase == 2) goto phase2
if ($phase == 3) goto phase3
echo "Unknown pahse"

goto phaseHLoop


phase1:
# col 5: a>g ,    6: target snippet  , 19,20,21,22   2 strand counts of mutant and wild type


# the format changed a little since the Fatigue project, we added the VCF name
if ($MAGIC == Fatigue || $MAGIC == NB) then
  cat tyty2 | cut -f 5,6,7,10,19-22  | gawk -F '\t' '{if(substr($1,2,1)==">"){$1=toupper($1);w=$2;gsub(/[atgc]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($2,n-2,4)"."$1;g=$3;gz[g]=z;mrna=$4;x=substr($1,1,1);w=substr(mrna,length(mrna)-2,1);if(w==x){pp[g,1]+=$5;pp[g,2]+=$6;pp[g,3]+=$7;pp[g,4]+=$8;ng[g]++;if(ng[g]==1)nz[z]++;}}}END{for(g in ng) {printf ("%s\t%s\t%d",g,gz[g],nz[gz[g]]);for(i=1;i<=4;i++)printf("\t%d",pp[g,i]);printf("\n");}}' > tyty3
  cat $ff | cut -f 5,6,7,10,19-22    | gawk -F '\t' '{if(substr($1,2,1)==">"){$1=toupper($1);w=$2;gsub(/[atgc]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($2,n-2,4)"."$1;g=$3;gz[g]=z;mrna=$4;x=substr($1,1,1);w=substr(mrna,length(mrna)-2,1);if(w==x){pp[g,1]+=$5;pp[g,2]+=$6;pp[g,3]+=$7;pp[g,4]+=$8;ng[g]++;if(ng[g]==1)nz[z]++;}}}END{for(g in ng) {printf ("%s\t%d",gz[g],nz[gz[g]]);for(i=1;i<=4;i++)printf("\t%d",pp[g,i]);printf("\n");}}' > totog.$type

  cat totog.$type | gawk -F '\t' '{z=$1;nn[z]=$2;for(i=3;i<=6;i++)nz[z,i-2]+=$i;}END{for(z in nn)printf("%s\t%d\t%d\t%d\t%d\t%d\n",z,nn[z],nz[z,1],nz[z,2],nz[z,3],nz[z,4]);}' > toto4.$type
  cat $ff | cut -f 5,6,7,10,19-22  | gawk -F '\t' '{if(substr($1,2,1)==">"){$1=toupper($1);w=$2;gsub(/[atgc]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($2,n-2,4)"."$1;g=$3;gz[g]=z;mrna=$4;x=substr($1,1,1);w=substr(mrna,length(mrna)-2,1);if(w==x){pp[z,1]+=$5;pp[z,2]+=$6;pp[z,3]+=$7;pp[z,4]+=$8;nnn[z]++;}}}END{for(z in nnn) {printf ("%s\t%d",z,nnn[z]);for(i=1;i<=4;i++)printf("\t%d",pp[z,i]);printf("\n");}}' > toto1.$type
else if ($MAGIC == Corona)
  cat $ff | cut -f 5,8,20-23  | gawk -F '\t' '{if(substr($1,2,1)==">"){w=$2;gsub(/[atgc]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($2,n-2,4)"."$1;for(i=3;i<=6;i++)pp[z,i-2]+=$i;nnn[z]++;}}END{for(z in nnn) {printf ("%s\t%d",z,nnn[z]);for(i=1;i<=4;i++)printf("\t%d",pp[z,i]);printf("\n");}}' > toto1.$type
else
  cat $ff | cut -f 5,8,20-23  | gawk -F '\t' '{if(substr($1,2,1)==">"){w=$2;gsub(/[ATGC]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($2,n-2,4)"."$1;for(i=3;i<=6;i++)pp[z,i-2]+=$i;nnn[z]++;}}END{for(z in nnn) {printf ("%s\t%d",z,nnn[z]);for(i=1;i<=4;i++)printf("\t%d",pp[z,i]);printf("\n");}}' > toto1.$type
  cat $ff | cut -f 5,7,8,20-23  | gawk -F '\t' '{if(substr($1,2,1)==">"){w=$3;gsub(/[ATGC]/,"x",w);k=split(w,aa,"x");if(k<2)next;n=length(aa[1]);if(n<4)next;z=substr ($3,n-2,4)"."$1;g=$3;gz[g]=z;for(i=4;i<=7;i++)pp[g,i-3]+=$i;ng[g]++;if(ng[g]==1)nz[z]++;}}END{for(g in ng) {printf ("%s\t%d",gz[g],nz[gz[g]]);for(i=1;i<=4;i++)printf("\t%d",pp[g,i]);printf("\n");}}' > totog.$type
  cat totog.$type | gawk -F '\t' '{z=$1;nn[z]=$2;for(i=3;i<=6;i++)nz[z,i-2]+=$i;}END{for(z in nn)printf("%s\t%d\t%d\t%d\t%d\t%d\n",z,nn[z],nz[z,1],nz[z,2],nz[z,3],nz[z,4]);}' > toto4.$type

endif



echo -n "### $toto1 :" > $toto1
date >> $toto1
echo "### Are substitutions following a given triplet on the the positive strand of the genes supported by plus or by minus reads, i.e. do we have a sequncing starnd bias indecating a probable systematic error of teh sequencing machine" >> $toto1
echo "# Triplet-substitution\tNumber of SNP in project $MAGIC\tPlus reads support the variant\tMinus reads support the variant \tPlus reads support the reference\tMinusreads support the reference\tStrand ratio in the variant\tStrand ration in the reference\tStrand ratio in the coverage\tAverage allele frequency plus strand\tAverage allele frequency minus strand\tDifference" >> $toto1
cat toto1.$type | gawk -F '\t' '{mp=$3;mm=$4;wp=$5;wm=$6;m=mp+mm;w=wp+wm;cp=mp+wp;cm=mm+wm;n=m+w;printf("%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$1,$2,mp,mm,wp,wm,100*mp/(.001+m),100*wp/(.0001+w),100*cp/(.0001+n),100*mp/(.0001+cp),100*mm/(.0001+cm),100*mp/(.0001+cp)-100*mm/(.0001+cm))}'  | sort >> $toto1

goto phaseLoop

phase2:
pahseHisto:

if (! -d HISTOS.$type) mkdir HISTOS.$type

cat totog.$type | cut -f 1 | sed -e 's/\./_/' -e 's/>/_/'  | sort -u > totog2.$type
foreach gg (`cat totog2.$type`)
  if (-e  HISTOS.$type/titi.$gg.txt) continue
  echo $gg
  cat totog.$type | sed -e 's/\./_/' -e 's/>/_/' | grep $gg | gawk '{printf("%d\t%d\n",$3,$3+$4);}' | bin/histo -smooth -o HISTOS.$type/titi.$gg &
end

goto phaseLoop

phase3:
pahseHisto2:

\rm -rf HISTOS.$type/all
foreach gg (`cat totog2.$type`)
  cat HISTOS.$type/titi.$gg.txt | gawk '/^#/{next;}{if(NF==2) printf("%s\t%s\t%s\n",gg,$1,$2);}' gg=$gg >> HISTOS.$type/all
end


cat HISTOS.$type/all | gawk '{g=$1;i=g2i[g]+0;if(i==0){ii++;i=ii;i2g[i]=g;g2i[g]=i;}k=int($2+0);z[i,k]=$3;}END{for(i=1;i<=ii;i++)printf("\t%s",i2g[i]);for(k=0;k<=100;k++){printf("\n%d",k);for(i=1;i<=ii;i++)printf("\t%.3f",z[i,k]);}printf("\n");for(i=1;i<=ii;i++)printf("\t%s",i2g[i]);printf("\n");}' > toto5.$type


set toto=RESULTS/SNV/$MAGIC.triplet.$type.histo.txt
echo -n "### $toto : " > $toto
date >> $toto
cat toto5.$type | transpose > toto5t.$type
cat toto5t.$type | head -1 | gawk -F '\t' '{for (i=2;i<=102;i++)printf("\t%s",$i);printf("\tVariant\tCumul\t0-10\t20-80\t30-70\t90-100\t1-10%%\t20-80%%\t30-70%%\t90-100%%\n");}' > toto5s.$type
cat toto5t.$type | tail -n +2 | gawk -F '\t' '{w=substr($1,1,4) substr($1,7,2);printf("%s",w);t=0;t1=0;t2=0;t3=0;t4=0;for(i=2;i<=12;i++)t1+=$i;for(i=22;i<=82;i++)t2+=$i;for(i=32;i<=72;i++)t3+=$i;for(i=92;i<=102;i++)t4+=$i;for(i=2;i<=102;i++)t+=$i;t+=.00001;t=t;for(i=2;i<=102;i++)printf("\t%.3f",$i);printf("\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",w,t,t1,t2,t3,t4,100*t1/t,100*t2/t,100*t3/t,100*t4/t);}' | sort -k 112nr >> toto5s.$type
cat toto5s.$type | transpose >> $toto 

goto phaseLoop

phaseLoop:
 echo 'done'
