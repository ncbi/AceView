#!bin/tcsh -f

# 2013_09_21
# given a rarrangement detected in RNA_seq
# grab all hits to the 2 genes and all corresponding fasta files
# realign locally and analyse the proportions

set gene1=$1
set gene2=$2
set run=$3

echo -n "Import the genome zone for $gene1 $gene2 : "
date
set dd="tmp/TranslocGenePair/GenePair.$gene1.$gene2"
set ici=`pwd`
pushd $dd
  mkdir database TABIX
  ln -s $ici/metaData/wspec.aceview_web_site wspec
  echo y | $ici/bin/tacembly .
popd

set map1=`cat tmp/TranslocGenePair/g2g2r.txt | gawk '{if(ok)next;split($1,aa,"::");split($2,bb,"::");if(aa[1]==g){ok++;print $1;}if(bb[1]==g){ok++;print $2;}}' g=$gene1`
set map2=`cat tmp/TranslocGenePair/g2g2r.txt | gawk '{if(ok)next;split($1,aa,"::");split($2,bb,"::");if(aa[1]==g){ok++;print $1;}if(bb[1]==g){ok++;print $2;}}' g=$gene2`

if ($run == ZZZZZ) then

echo $map1 | gawk '{split($1,aa,":");c=aa[3];split(aa[4],bb,"_");a1=bb[1];a2=bb[2];if(a1<a2){a1-=300;a2+=300;}else{a1+=300;a2-=300;cpl=" | bin/dna2dna -complement"; a0=a1;a1=a2;a2=a0;}printf("bin/dna2dna -I fasta -i TARGET/CHROMS/%s.chrom_%s.fasta.gz -leftClipAt %d  -rightClipAt %d -O fasta %s > %s/%s.F.genome.fasta\n", species,c,a1,a2,cpl,dd,gene1)}' species=$species dd=$dd gene1=$gene1 > $dd/getdna1

echo $map2 | gawk '{split($1,aa,":");c=aa[3];split(aa[4],bb,"_");a1=bb[1];a2=bb[2];if(a1<a2){a1-=300;a2+=300;}else{a1+=300;a2-=300;cpl=" | bin/dna2dna -complement"; a0=a1;a1=a2;a2=a0;}printf("bin/dna2dna -I fasta -i TARGET/CHROMS/%s.chrom_%s.fasta.gz -leftClipAt %d  -rightClipAt %d -O fasta %s > %s/%s.F.genome.fasta\n", species,c,a1,a2,cpl,dd,gene1)}' species=$species dd=$dd gene1=$gene2 > $dd/getdna2

if (! -e  $dd/$gene1.F.genome.fasta) then
  source $dd/getdna1
  bin/dna2dna -I fasta -i  $dd/$gene1.F.genome.fasta -O fasta -complement >  $dd/$gene1.R.genome.fasta
endif
  foreach ff (F1F2 F2F1 F1R2 R1F2)
    echo "DNA $gene1.$ff" >>  $dd/genome.ace
    tail -n +2 $dd/$gene1.F.genome.fasta >> $dd/genome.ace
    echo >> $dd/genome.ace
  end

if (! -e  $dd/$gene2.F.genome.fasta) then
  source $dd/getdna2
  bin/dna2dna -I fasta -i  $dd/$gene2.F.genome.fasta -O fasta -complement >  $dd/$gene2.R.genome.fasta
endif
  foreach ff (F1F2 F2F1 F1R2 R1F2)
    echo "DNA $gene2.$ff" >>  $dd/genome.ace
    tail -n +2 $dd/$gene2.F.genome.fasta >> $dd/genome.ace
    echo >> $dd/genome.ace
  end

endif

#### construct the rearranged genomes
# F1F2 Gene1 -->       Gene2--->
# F2F1 Gene2 -->       Gene1--->
# F1R2 Gene1 -->       <----Gene2
# R1F2 Gene1 <--       --->Gene2
echo "export the reaaranged genome : F1.F2 F2.F1 R1.F2 F1.R2"

# export Nn nnn between the 2 genes so that position ND is always between the 2 genes
set N1=`cat $dd/$gene1.F.genome.fasta | gawk '/^>/{next;}{n+=length($1);}END{print n}'`
set N2=`cat $dd/$gene2.F.genome.fasta | gawk '/^>/{next;}{n+=length($1);}END{print n}'`
set ND=`echo $N1 $N2 | gawk '{n=$1;if($2>n)n=$2;print n+50;}'`


cat <<EOF>> $dd/methods.ace

Method "Genefinder"
Colour   MIDBLUE
Show_up_strand
Bumpable
Right_priority   1.8
Show_text
EMBL_feature     "CDS"

Method "HighJump"
Colour   RED
Show_up_strand
Bumpable
Right_priority   2.8
Show_text
EMBL_feature     "CDS"

Method "Jump"
Colour   CYAN
Show_up_strand
Bumpable
Right_priority   3.8
Show_text
EMBL_feature     "CDS"

Method "LowJump"
Colour   YELLOW
Show_up_strand
Bumpable
Right_priority   4.8
Show_text
EMBL_feature     "CDS"

EOF

if ($run == ZZZZZ) then

echo $N1 $N2 $gene1 $gene2 | gawk '{n1=$1;n2=$2; dn1=n2-n1;if(dn1<0)dn1=0;dn2=n1-n2;if(dn2<0)dn2=0;g1=$3;g2=$4;S="Sequence";SS="Subsequence";printf("%s F1.F2\n%s %s.F1F2 %d %d\n%s %s.F1F2 %d %d\n\n",S,SS,g1,dn1+1,n1+dn1,SS,g2,n1+dn1+101,n1+dn1+100+n2);printf("%s F2.F1\n%s %s.F2F1 %d %d\n%s %s.F2F1 %d %d\n\n",S,SS,g2,dn2+1,dn2+n2,SS,g1,n2+dn2+101,n2+dn2+n1+100);printf("%s F1.R2\n%s %s.F1R2 %d %d\n%s %s.F1R2 %d %d\n\n",S,SS,g1,dn1+1,dn1+n1,SS,g2,dn1+n1+n2+100,dn1+n1+101); printf("%s R1.F2\n%s %s.R1F2 %d %d\n%s %s.R1F2 %d %d\n\n",S,SS,g1,dn1+n1,dn1+1,SS,g2,dn1+n1+101,dn1+n1+n2+100);  }' >> $dd/genome.ace

# recover the coordinates of the transcripts
foreach gene ($gene1 $gene2)
  gunzip -c  tmp/METADATA/mrnaRemap.gz | gawk '{if($8==gene)print}' gene=$gene> $dd/_t
  foreach fr (F1F2 F2F1 F1R2 R1F2)
    cat $dd/_t | gawk -F '\t' '{m=$2 "." fr ;mm[m]++;k=mm[m];m1[m,k]=$6;m2[m,k]=$7;if(g1<1){g1=$6;g2=$6;}if($6<g1)g1=$6;if($7<g1)g1=$7;if($6>g2)g2=$6;if($7>g2)g=$7;if($6<$7)isDown=1;}END{printf("Sequence %s.%s\n",gene,fr); for (m in mm) {if(isDown==1){uu1=301+m1[m,1]-g1;uu2=301+m2[m,mm[m]]-g1;}else {uu1=301+g2-m1[m,1];uu2=301+g2-m2[m,mm[m]];}printf("Subsequence %s %d %d\n", m,uu1,uu2);}  for (m in mm) {printf("\nSequence %s\nCDS\nMethod Genefinder\n",m);for(i=1;i<=mm[m];i++){if(isDown==1){uu1=1+m1[m,i]-m1[m,1];uu2=1+m2[m,i]-m1[m,1];}else {uu1=1+m1[m,1]-m1[m,i];uu2=1+m1[m,1]-m2[m,i];}printf("Source_exons %d %d\n",uu1,uu2);}printf("\n");}}'  gene=$gene fr=$fr >> $dd/genome.ace
   end
  \rm $dd/_t
end
 
bin/tacembly $dd <<EOF
  pparse  $dd/genome.ace
  pparse  $dd/methods.ace
  save
  find sequence F1.F2
  dna  $dd/F1.F2.genome.fasta
  find sequence F2.F1
  dna  $dd/F2.F1.genome.fasta
  find sequence F1.R2
  dna  $dd/F1.R2.genome.fasta
  find sequence R1.F2
  dna  $dd/R1.F2.genome.fasta
  quit  
EOF

goto phaseLoop
endif

bin/tacembly $dd <<EOF
  find sequence F1.F2
  dna  $dd/F1.F2.genome.fasta
  find sequence F2.F1
  dna  $dd/F2.F1.genome.fasta
  find sequence F1.R2
  dna  $dd/F1.R2.genome.fasta
  find sequence R1.F2
  dna  $dd/R1.F2.genome.fasta
  quit  
EOF

# identification des clones utiles
echo "Import from tmp/COUNT/$run/_.hits.gz the relevant clones"
if (! -d $dd/$run/fastc.all) then
  if (! -d $dd/$run)  mkdir $dd/$run
  foreach lane (`cat Fastc/$run/LaneList`) 
    if (-e $dd/$lane.fastc) continue
    echo "       import from $lane"
    gunzip -c tmp/COUNT/$lane.hits.gz | gawk -F '\t' '{ if(first<1){first=1;split(map1,aa,":");chrom1=aa[3];split(aa[4],bb,"_");b1=bb[1];b2=bb[2];split(map2,aa,":");chrom2=aa[3];split(aa[4],bb,"_");c1=bb[1];c2=bb[2];if(b1>b2){b0=b1;b1=b2;b2=b0;}if(c1>c2){c0=c1;c1=c2;c2=c0;}}g=$9;if($8==tcl && g==g1 || g==g2)print ;else {a1=$12;a2=$3;if(a1>a2){a0=a1;a1=a2;a2=a0;}if($8 == "Z_genome" && (($11==chrom1 && a1<b2 && a2>b1) || ($11==chrom1 && a1<b2 && a2>b1)))print;}}' g1=$gene1 g2=$gene2 map1=$map1 map2=$map2 >  $dd/$lane.hits
    cat   $dd/$lane.hits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);print c;}' | sort -u > $dd/$lane.list
    bin/dna2dna -i Fastc/$lane.fastc.gz -I fastc -O fastc -select $dd/$lane.list -o $dd/$lane
  end

  # grab the bridging clones 
  cat $dd/$run/*.*.hits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);g=$9;if(g==g1)cc1[c]=1;if(g==g2)cc2[c]=1;}END{for(c in cc1)if(cc2[c]==1)print c}' g1=$gene1 g2=$gene2 | sort > $dd/$run/bridging.list

  # group all fastc reads (they belong to the same run so the ids are distinct)
  cat  $dd/$run/*.fastc > $dd/$run/fastc.all
endif

set step=1
echo "Align the clones on the rearranged genome sections"
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
echo "       align on $fr"
  if (-e  $dd/$run/$fr.genome.hits) continue
  bin/clipalign -i  $dd/$run/fastc.all -t  $dd/$fr.genome.fasta -minAli 16 -splice -seedOffset 1 -seedShift 5 -clipPolyA -o $dd/$run/$fr.genome
end


# select the best composite hits and filter on 90% quality
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  # split the hit table at the boundary and separatelly clean up the 2 sides
  cat  $dd/$run/$fr.genome.hits | gawk -F '\t' '{a1=$12;a2=$13;if (a1<ND && a2<ND) print}' ND=$ND >   $dd/$run/$fr.genome.hits.1
  cat  $dd/$run/$fr.genome.hits | gawk -F '\t' '{a1=$12;a2=$13;if (a1>ND && a2>ND) print}' ND=$ND >   $dd/$run/$fr.genome.hits.2
  cat  $dd/$run/$fr.genome.hits.1 | bin/bestali -exportBest > $dd/$run/$fr.genome.besthits.1
  cat  $dd/$run/$fr.genome.hits.2 | bin/bestali -exportBest > $dd/$run/$fr.genome.besthits.2
  # then merge, sort by score and remove cases where the lower score is not extending the coverage in read coordinates
  cat $dd/$run/$fr.genome.besthits.[12] | sort -k 1,1 -k 2,2nr |  gawk -F '\t' '{r=$1;if(r==old){y1=$6;y2=$7;u1=x1;if(y1>u1)u1=y1;u2=x2;if(y2<u2)u2=y2;du=u2-u1;if(2*du>y2-y1 || 2*du>x2-x1 || du < -25)next;}else{x1=$6;x2=$7;}old=r;if($2>=25)print;}' >  $dd/$run/$fr.genome.betterhits
  # select high single hit or high composite hit quality
  cat  $dd/$run/$fr.genome.betterhits |  gawk -F '\t' '{if(100*$2>90*$4)print;r=$1;if(r==old){y1=$6;y2=$7;b1=$12;u1=x1;if(y1>u1)u1=y1;u2=x2;score+=$2;if(y2<u2)u2=y2;du=u2-u1;if(du>1)score=score-du+1;if ( 100*score > 90*ln){print z; print;}}score=$2;ln=$4;x1=$6;x2=$7;a2=$13;z=$0;old=r;}'   | sort -k 1,1 -k 12,12n >  $dd/$run/$fr.genome.besthits 
end


echo "Synthetize the results, single read jump (4 orientations) or facing pair (2 orientations)"
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  # grab the bridging clones
  cat  $dd/$run/$fr.genome.besthits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);s=substr(r,length(r));cc[c]=1;a1=$12;a2=$13;if(s==">" && a1<a2 && a1<ND)cc1a[c]=1;if(s==">" && a1<a2 && a2>ND)cc2a[c]=1;if(s==">" && a1>a2 && a2<ND)rr1a[c]=1;if(s==">" && a1>a2 && a1>ND)rr2a[c]=1;if(s=="<" && a1<a2 && a1<ND)cc1b[c]=1;if(s=="<" && a1<a2 && a2>ND)cc2b[c]=1;if(s=="<" && a1>a2 && a2<ND)rr1b[c]=1;if(s=="<" && a1>a2 && a1>ND)rr2b[c]=1;}END{for(c in cc)if(cc1a[c]*cc2a[c]+cc1b[c]*cc2b[c]+rr1a[c]*rr2a[c]+rr1b[c]*rr2b[c] + cc1a[c]*rr2b[c]+cc1b[c]*rr2a[c] > 0 )print c}' ND=$ND >  $dd/$run/$fr.bridging_clone.list
  echo ZZZZZ > ZZZZZ
  cat $dd/$run/$fr.bridging_clone.list ZZZZZ  $dd/$run/$fr.genome.besthits |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{r=$1;c=substr(r,1,length(r)-1);if(ok[c]==1)print;}' | sort -k 1,1 -k 12,12n >  $dd/$run/$fr.bridging_clone.hits
 
end

# look in this list for exact introns
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  cat  $dd/$run/$fr.bridging_clone.hits |  gawk -F '\t' '{r=$1;if(r==old){y1=$6;y2=$7;b1=$12;u1=x1;if(y1>u1)u1=y1;u2=x2;score+=$2;if(y2<u2)u2=y2;du=u2-u1;if (du >= -5 && a2<ND && b1 > ND && 100*(score-du+1) > 90*ln){print z; print;}}score=$2;ln=$4;x1=$6;x2=$7;a2=$13;z=$0;old=r;}' ND=$ND >  $dd/$run/$fr.bridging_reads.hits
  cat $dd/$run/$fr.bridging_reads.hits | sort -k 1,1 -k 12,12n |  gawk -F '\t' '{r=$1;b1=$12;b2=$13;if(r==old)nn[a2 "\t" b1]+=$3;old=r;a1=b1;a2=b2;}END{for(k in nn) printf("%s\t%s\t%s\t%s\t%d\n",g1,g2,fr,k,nn[k]);}' g1=$gene1 g2=$gene2 fr=$fr >   $dd/$run/$fr.bestintrons
end

if (-e $dd/$run/all.genome.introns.1) \rm  $dd/$run/all.genome.introns.1
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  cat $dd/$run/$fr.genome.introns | sort -k 11nr |  gawk '{if(($5<ND && $7>ND) || ($7<ND && $5>ND)){if(10*$11>x){printf ("%s\t%s\t",run,fr);print;}if(x==0)x=$11;}}' fr=$fr ND=$ND run=$run >> $dd/$run/all.genome.introns.1
end
cat  $dd/$run/all.genome.introns.1  | sort -k 13nr |  gawk '{if(10*$13>x)print;if(x==0)x=$13;}' > $dd/$run/all.genome.introns
#\rm $dd/$run/all.genome.introns.1

echo ND=$ND
ls -ls $dd/$run/all.genome.introns*

# Find the orientation of the best intron
set bestfr=`cat  $dd/$run/*.bestintrons | gawk '{n[$3]++}END{for(k in n)printf("%s\t%d\n", k,n[k])}' | sort -k 2nr | head -1 | cut -f 1`

# construct the wiggle
if (! -e  $dd/TABIX/$run/wiggle.ace) then
  set step=10
  mkdir $dd/TABIX/$run $dd/TABIX/$run.bridge
  foreach fr ($bestfr)
    if (! -e $dd/$run/$fr.genome.hits.step$step.BV4) then
      cat  $dd/$run/$fr.bridging_clone.hits | bin/wiggle     -strand -I BHIT  -out_step $step -O BV |  gawk '/^variableStep/{zz++;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}' chrom=$fr | bgzip -c > $dd/TABIX/$run.bridge/$fr.u.f.tabix.gz
      cat  $dd/$run/$fr.bridging_clone.hits | bin/wiggle -antistrand -I BHIT  -out_step $step -O BV |  gawk '/^variableStep/{zz++;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}' chrom=$fr | bgzip -c > $dd/TABIX/$run.bridge/$fr.u.r.tabix.gz
      tabix -s 1 -b 2 -e 2  $dd/TABIX/$run.bridge/$fr.u.f.tabix.gz
      tabix -s 1 -b 2 -e 2  $dd/TABIX/$run.bridge/$fr.u.r.tabix.gz
      cat  $dd/$run/$fr.genome.besthits | bin/wiggle     -strand -I BHIT  -out_step $step -O BV |  gawk '/^variableStep/{zz++;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}' chrom=$fr | bgzip -c > $dd/TABIX/$run/$fr.u.f.tabix.gz
      cat  $dd/$run/$fr.genome.besthits | bin/wiggle -antistrand -I BHIT  -out_step $step -O BV |  gawk '/^variableStep/{zz++;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}' chrom=$fr | bgzip -c > $dd/TABIX/$run/$fr.u.r.tabix.gz
      tabix -s 1 -b 2 -e 2  $dd/TABIX/$run/$fr.u.f.tabix.gz
      tabix -s 1 -b 2 -e 2  $dd/TABIX/$run/$fr.u.r.tabix.gz
      echo toto | gawk '{n=n1;if(n2>n1)n=n2;printf("Sequence %s\nIntMap %s 1 %d\n\nMap %s\n-D Wiggle\nWiggle %s %s/%s\nWiggle %s.bridge %s.bridge/%s\n\nRun %s\nW_colour_plus darkgreen\nW_colour_minus lightGreen\n\nRun %s.bridge\nW_colour_plus MAGENTA\nW_colour_minus LIGHTMAGENTA\n\n",fr,fr,n+n2+100,fr,run,run,fr,run,run,fr,run,run);}' fr=$fr run=$run n1=$N1 n2=$N2 dd=$dd step=$step >> $dd/TABIX/$run/wiggle.ace
    endif
  end
endif

cat  $dd/$run/*.bestintrons | gawk -F '\t' '{nn++;mt="Jump";if($6<3)next;if($6<10)mt="LowJump";if($6>=100)mt="HighJump";printf("Sequence %s\nGenomic\nSubsequence %s.%s.%d.jump.%d %d %d\n\n",$3,$3,run,$6,nn,$4-30,$5+30);printf("Sequence  %s.%s.%d.jump.%d\n-D Method\nMethod %s\n-D Source_exons\n",$3,run,$6,nn,mt);u1=1;u2=30;v1=$5-$4+1;v2=$5-$4+30; printf("Source_exons %d %d\nSource_exons %d %d\n\n",u1,u2,v1,v2) ;}' run=$run > $dd/$run/all.genome.introns.ace

laba:


echo "pparse $dd/TABIX/$run/wiggle.ace" >  $dd/$run/_r
echo "pparse $dd/$run/all.genome.introns.ace" >>  $dd/$run/_r
echo "pparse $dd/methods.ace" >>  $dd/$run/_r
echo save >>  $dd/$run/_r
above:
  sleep 1
  if (-e $dd/database/lock.wrm) goto above

bin/tacembly $dd < $dd/$run/_r

goto phaseLoop


phaseLoop:
 echo -n "done "
 date
