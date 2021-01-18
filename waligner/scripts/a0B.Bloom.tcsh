#!bin/tcsh -ef

set task=$1
set run=$2
set type=$3
set type2=$4

echo $task $run

set solid="cat"
foreach run2 (`cat MetaDB/$MAGIC/RunSolidList`)
  if ($run == $run2) set solid=" dna2dna -I csfastc -O fastc "
end

if ($task == telomere) then

  if ($type == any) then
    echo "gunzip -c Fastc/$run/"'*'".fastc.gz | $solid | bin/dnabloom -telomere -run $run -max_threads 8 -o tmp/Bloom/$run/a0B"
    gunzip -c Fastc/$run/*.fastc.gz | $solid | bin/dnabloom -telomere -max_threads 8  -run $run -o tmp/Bloom/$run/a0B
  endif
  if ($type == unaligned) then
    echo "gunzip -c tmp/Unaligned/$run/"'*'".fastc.gz | $solid | bin/dnabloom -telomere -run $run -max_threads 8 -o tmp/Bloom/$run/a0B"
    gunzip -c tmp/Unaligned/$run/*.fastc.gz | $solid | bin/dnabloom -telomere -max_threads  8 -run $run -o tmp/Bloom/$run/a0B.unaligned
  endif
  exit 0
endif

if ($task == telomereReport) then

  # regroup the sublibs
  foreach run (`cat MetaDB/$MAGIC/r2sublib | cut -f 1 | sort -u`)
    if (! -d tmp/Bloom/$run) then
      $mkDir Bloom $run
    endif

    if ($type == any) then
      set ff=tmp/Bloom/$run/a0B.telomere_distribution.txt
    endif
    if ($type == unaligned) then
      set ff=tmp/Bloom/$run/a0B.unaligned.telomere_distribution.txt
    endif
    touch $ff.1
    \rm $ff.1
    foreach sublib (`cat  MetaDB/$MAGIC/r2sublib | gawk '{if($1 == run)print $2;}' run=$run`)
      if ($type == any) then
        set ff1=tmp/Bloom/$sublib/a0B.telomere_distribution.txt
      endif
      if ($type == unaligned) then
        set ff1=tmp/Bloom/$sublib/a0B.unaligned.telomere_distribution.txt
      endif
      cat $ff1 | gawk '/^#/{next;}{print}' > $ff.1 
    end
    cat $ff.1 | gawk -F '\t' '{if(NF>nf)nf=NF;t=$2;n=t2n[t]; if(t2n[t]<1){nt++;t2n[t]=nt;n=nt;n2t[nt]=t;}for(i=3;i<=NF;i++)z[n,i]=$i;}END{for (n=1;n<=nt;n++){printf("%s\t%s",run,n2t[n]);for (i=3 ; i <= nf ; i++) printf("\t%s",z[n,i]);printf("\n");}}' run=$run > $ff
  end

  if ($type == any) then
    set toto=RESULTS/Mapping/$MAGIC.telomere.txt
    set toto=tmp/Bloom/$MAGIC.telomere.txt
  endif
  if ($type == unaligned) then
    set toto=RESULTS/Mapping/$MAGIC.unaligned.telomere.txt
    set toto=tmp/Bloom/$MAGIC.unaligned.telomere.txt
  endif
  touch $toto.1
  \rm $toto.[0123]
  echo -n "## $toto : " > $toto.0
  date >> $toto.0
  
  echo "## Histogram of prevalence of telomeric and poly-A homopolymers of length at least 13 nucleotides in project $MAGIC." >> $toto.0
  if ($type == any) then
    echo "##  This table is counted over all fragments present in the original fastq files." >> $toto.0
  endif
  if ($type == unaligned) then
    echo "##  This table is counted over the unaligned fragments." >> $toto.0
  endif
  echo "## A fragment is either a read-pair (in case of paired end sequencing) or a read (in case of single end sequencing)."  >> $toto.0
  echo "## The characteristic telomeric motif is the homopolymer (TTAGGG/CCCTAA)n, poly-A is (A/T)n" >> $toto.0
  echo "## A G>C variation is allowed in the variant telomeric motifs" >> $toto.0 
  echo "## This analysis is non stranded: the motif and its complement are both counted." >> $toto.0
  echo "## The maximal number of motifs is dictated by the sequencing length. In a pure poly-A 100+100 read pair, there would be 88+88=176 (A)13 motifs" >> $toto.0
  echo "## Therefore fragments with 11 A13 motif could either include 11 isolated polyA[13] motifs, or a single 23-mer polyA motif, or a polyA[19] (counting 7) and a polyA[16] (counting 4) etc." >> $toto.0
  echo "## A fragment with 11 telomeric motifs has necessarilly 11 distinct motifs occupying at least 66 bases, sinc e the telomere motif is non-sliding" >> $toto.0
  echo "## Each very long reads is cut in regular fragments of around 1000 bases, then the counts are normalized per million such fragments. e.g. 7 segments of length 1050 for a read of length 7350" >> $toto.0
  echo "## The number of fragments in a given run with the given number of motifs is reported in the table, as cumuls showing fragments with at least so many motifs." >> $toto.0
  
  echo >> $toto.0
  
  # count the columns
  set nn=0
  foreach run ( `cat MetaDB/$MAGIC/RunList`)
    set ff=tmp/Bloom/$run/a0B.telomere_distribution.txt
    if (-e $ff) then
      set nn=`cat $ff | head -12 | gawk -F '\t' '{if (NF>n) n = NF ; }END{print n-3;}' n=$nn`
    endif
  end
  echo "reporting over $nn columns $type"
  # cumulate the values in a single file
    set ff=tmp/Bloom/$MAGIC.telomere.$type.ace  
    echo " " > $ff
  echo > $toto.1
  foreach run ( `cat MetaDB/$MAGIC/RunList`)
    set ff="x"
    if ($type == any) then
      set ff=tmp/Bloom/$run/a0B.telomere_distribution.txt
    endif
    if ($type == unaligned) then
      set ff=tmp/Bloom/$run/a0B.unaligned.telomere_distribution.txt
    endif
    # echo $ff
    if (-e $ff) then
      cat $ff | gawk -F '\t' '/^$/{next;}{printf("%s",run);for(i=2;i<=NF;i++)printf ("\t%s",$i);printf("\n");}' run=$run | sed -e 's/per_fragment/per_read_or_read_pair/g' | sed -e 's/per_kb/per_1_kb_fragment/g'  >> $toto.1
    endif
  end

  ls -ls $toto.?

  cat $toto.1  | gawk  -F '\t' '/^#/{next}{if($3+0<3)next;n=0;for(i=0;i<=NF;i++)if(0+$i>0)n=i;printf("%s\t%s",$2,$1);x=1000000.0*$6/($3+1);printf("\t%.2f",x);for(i=3;i<=n;i++)printf("\t%s",$i);for(i=n+1;i<=nn;i++)printf("\t");printf("\n");}' nn=$nn > $toto.3

   foreach fp (per_read_or_read_pair per_1_kb_fragment)
    foreach motif (A13-$fp-cumul Telomeric_6-$fp-cumul Tel_C_6-$fp-cumul imagT_6-$fp-cumul imagC_6-$fp-cumul)
      echo $motif.$fp
      echo -n "#Type\tRun\tIndex:FPM $fp for at least 3 motifs\tNumber of fragments\tAt least one motif" >> $toto.2
      echo $nn | gawk '{for(i=2;i<=$1;i++)printf ("\t%d motifs",i);printf("\n");}' >> $toto.2
      cat $toto.3 | gawk -F '\t' '{if($1==motif) print}' motif="$motif"  >> $toto.2
      echo "\n\n\n"  >> $toto.2
    end
  end

  cat $toto.1 | sed -e 's/per_1_kb_fragment/per_million_1_kb_fragments/g' | gawk  -F '\t' '/^#/{next}{if($3+0<3)next;n=0;for(i=0;i<=NF;i++)if(0+$i>0)n=i;printf("%s\t%s",$2,$1);x=1000000.0*$6/($3+1);printf("\t%.2f",x);z=1000000.00/$3;for(i=3;i<=n;i++)printf("\t%.2f",z*$i);for(i=n+1;i<=nn;i++)printf("\t");printf("\n");}' nn=$nn > $toto.3

  foreach fp (per_million_1_kb_fragments)
    foreach motif (A13-$fp-cumul Telomeric_6-$fp-cumul Tel_C_6-$fp-cumul imagT_6-$fp-cumul imagC_6-$fp-cumul)
      echo $motif.$fp
      echo -n "#Type\tRun\tIndex:FPM Normalized for at least 3 motifs\tNumber of fragments\tAt least one motif" >> $toto.2
      echo $nn | gawk '{for(i=2;i<=$1;i++)printf ("\t%d motifs",i);printf("\n");}' >> $toto.2
      cat $toto.3 | gawk -F '\t' '{if($1==motif) { print}}' motif="$motif"  >> $toto.2
      echo "\n\n\n"  >> $toto.2
    end
  end
  cat $toto.0 $toto.2 > $toto

  set fface=tmp/Bloom/$MAGIC.telomere.$type.ace  
  touch $fface
  \rm $fface
  if (! -e $fface) then
  foreach motif (`cat $toto | cut -f 1 | grep cumul | sort -u`)
      echo $motif 
      cat $toto | gawk -F '\t' '{if($1==motif) print}' motif="$motif" | gawk -F '\t' '{run=$2; if($4<1000)next ;printf("Ali %s\nBloom \"%s.%s\" %f %d %d %d %d %d  %d %d %d %d %d\n\n",run,type,$1,$3,$4,$5,$6,$7,$8,$9,$14,$24,$54,$104);}' type=$type  >> $fface
    end
    echo bb
    ls -ls $toto.? $ff
      tbly MetaDB << EOF
        read-models
        parse  $fface
        save
        quit
EOF
  endif

  # \rm $toto.[0123]
  if (! -d RESULTS/Bloom) mkdir RESULTS/Bloom
  foreach run ( `cat MetaDB/$MAGIC/RunListSorted`)
    \cp tmp/Bloom/$run/a0B.telomere_distribution.txt RESULTS/Bloom/a0B.$run.telomere_distribution.txt
  end
  \cp $toto RESULTS/Bloom
  \cp $toto RESULTS/Mapping
  echo $toto
   exit 0
endif
  
###########################################################################################
if ($task == telomereRemap) then
   echo -n  "a0B.Bloom telomereRemap :"
   date
# construct a fasta file
    cat tmp/Bloom/$run/a0B.telomere_hits.txt | gawk -F '\t' '/^#/{next}{k++;z=$2 "__" $3 "__" k "#" $4;printf(">%s\n%s\n",z,$5);}' | gzip > tmp/Bloom/$run/a0m.fastc.gz
    bin/clipalign -i tmp/Bloom/$run/a0m.fastc.gz -t TARGET/Targets/$species.genome.fasta.gz -minAli 30 -maxHit 30 -gzo -o tmp/Bloom/$run/a0m

   exit 0
endif

###########################################################################################
## genome loop

if ($task == METADATA) then

  if (! -d TARGET/Bloom) then
    mkdir TARGET/Bloom
    echo "These files are created automatically by the command : a0B.Bloom.tcsh METADATA " >  TARGET/Bloom/README.Bloom.automatic_creation
  endif
  foreach target ($RNAtargets $DNAtargets)
    if ($target == genome || $target == gdecoy) continue
    if (! -e  TARGET/Bloom/$target.telomere_distribution.txt) then
      echo "Telomere bloom $target"
      gunzip -c TARGET/Targets/$species.$target.fasta.gz  | bin/dnabloom -telomere -max_threads 4  -run $target -I fasta -noPolyA  -o  TARGET/Bloom/$target
    endif
  end
  
  # decompose the genome in sections
  
  if (! -e  TARGET/Bloom/chroms.telomere_distribution.txt) then
    foreach chrom ($chromSetAll) 
      if (-e TARGET/CHROMS/$species.chrom_$chrom.fasta.gz) then
        echo  "Telomere bloom chromosome $chrom"
        gunzip -c TARGET/CHROMS/$species.chrom_$chrom.fasta.gz |  bin/dnabloom -telomere -max_threads 0 -block 300  -run chr$chrom -I global -noPolyA >>  TARGET/Bloom/chroms.telomere_distribution.txt
      endif
    end
  endif

  if (! -e  TARGET/Bloom/chroms.sections.telomere_distribution.txt) then
    foreach chrom ($chromSetAll) 
      if (-e TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  && ! -e TARGET/Bloom/chr$chrom.sections.telomere_distribution.txt) then
        echo  "Telomere bloom chromosome $chrom"
        set n=`gunzip -c TARGET/CHROMS/$species.chrom_$chrom.fasta.gz | wc | gawk '{print $1}'`
        echo $n | gawk '{dixM=200000;r=$1 % dixM;nn=int($1/dixM);for(i=0;i<nn/2;i++)printf("%d %d\n",1+dixM*i,dixM);printf("%d %d\n",1+dixM*i,r);for(;i<nn;i++)printf("%d %d\n",1+dixM*i + r,dixM);}' > TARGET/Bloom/c
        cat TARGET/Bloom/c | gawk '{n++;printf("gunzip -c %s | tail -n +%d | head -%d | bin/dnabloom -telomere -max_threads 0 -run chr%s.%s -noPolyA -block 300 -I global >> %s\n",f,$1,$2,chrom,n,g);}' chrom=$chrom f=TARGET/CHROMS/$species.chrom_$chrom.fasta.gz g=TARGET/Bloom/chroms.sections.telomere_distribution.txt >  TARGET/Bloom/d
        source  TARGET/Bloom/d
      endif
    end
  endif
  
  
    if (! -d RESULTS/Mapping) mkdir RESULTS/Mapping
  
    set toto=tmp/Bloom/$MAGIC.targets.telomere.txt
   echo -n "## $toto : " > $toto.0
  date >> $toto.0
  
  echo "## Histogram of prevalence of telomeric homopolymers of length at least 13 nucleotides in the reference targets." >> $toto.0
  if ($type == any) then
    echo "##  This table is counted over all fragments present in the original fastq files." >> $toto.0
  endif
  if ($type == unaligned) then
    echo "##  This table is counted over the unaligned fragments." >> $toto.0
  endif
  echo "## A fragment is either a read-pair (in case of paired end sequencing) or a read (in case of single end sequencing)."  >> $toto.0
  echo "## The characteristic telomeric motif is the homopolymer (TTAGGG/CCCTAA)n" >> $toto.0
  echo "## This analysis is non stranded: the motif and its complement are both counted." >> $toto.0
  echo "## The maximal number of motifs is equal to the sequence length -12. For example, in a pure poly-A 100 base sequence, there would be 88 (A)13 motifs" >> $toto.0
  echo "## Therefore fragments with 11 motifs could either include 11 isolated 13-mers, or a single 23-mer motif, or a 19-mer (counting 7) and a 16-mer (counting 4) etc." >> $toto.0
  echo "## The number of fragments in a given run with the given number of motifs is reported in the table, as plain counts, then as cumuls showing fragments with at least so many motifs." >> $toto.0
  
  echo >> $toto.0

    foreach target ($RNAtargets $DNAtargets)
      if (-e TARGET/Bloom/$target.telomere_distribution.txt) cat  TARGET/Bloom/$target.telomere_distribution.txt | grep TERRA | grep cumul | sed -e 's/TERRA/telomeric/' -e 's/0-cumul/cumul/' >> $toto.1
    end
    cat TARGET/Bloom/chroms.telomere_distribution.txt | grep TERRA | grep cumul | sed -e 's/TERRA/telomeric/' -e 's/0-cumul/cumul/' >> $toto.1
    cat TARGET/Bloom/chroms.sections.telomere_distribution.txt | grep TERRA | grep cumul | sed -e 's/TERRA/telomeric/' -e 's/0-cumul/cumul/' >> $toto.1
  
  
    cat $toto.1 | gawk -F '\t' '{printf ("%s\t%s\t%s",$2,$1,$3);for(i=4;i<=NF ;i++){if($i>0)printf("\t%s",$i);}printf("\n");}' | scripts/transpose > $toto.2
    cat $toto.2 | gawk 'BEGIN{n=-3;}{n++;if(n<=10 || (n<=100 && (n%10) == 0) ||  (n<=1000 && (n%100) == 0) ||  (n<=10000 && (n%1000) == 0) ) {printf("%d\t",n);print} }' | scripts/transpose > $toto.3
    
    cat $toto.0 $toto.3 >  $toto

    \rm $toto.?
  
  exit 0
endif
  
if ($task == fastq) then

   echo "gunzip -c $type $type2 | dnabloom -telomere -max_threads 4 -run $run -I fastq -o Bloom/$run"
         gunzip -c $type $type2 | dnabloom -telomere -max_threads 4 -run $run -I fastq -o Bloom/$run
 
  exit 0
  foreach ff (`ls *_1.fastq.gz`)
    set run=`echo $ff | gawk '{i = index($1,"_1.fastq.gz");print substr($1,1,i-1);}'`
    set ff2=`echo $ff | sed -e 's/_1\.fastq.gz/_2.fastq.gz/'`

    scripts/submit Bloom/$run "scripts/a0B.Bloom.tcsh  fastq $run $ff $ff2"  32G
    
  end

endif 
  
##############################################
## motif distribution along the chromosomes

cat /home/mieg/SEQC2_2020/TARGET/CHROMS/Centromere_telomeres.ace | gawk '/^Map/{gsub(/\"/,"",$2);chrom=$2;}/^Centromere_telomeres/{t1=$2;t2=$5;c=($3+$4)/2;chLn[chrom]=t2;chCtr[chrom]=c;printf("%s\t%d\t%d\n",chrom,t2,c);}' | gzip > tmp/METADATA/chrom2centre2length.txt.gz
  
foreach run (`cat MetaDB/$MAGIC/RunList`)

  if (-e tmp/Bloom/$run/a0m.hits.gz && ! -e tmp/Bloom/$run/chrom2motif.txt) then
    zcat tmp/METADATA/chrom2centre2length.txt.gz ZZZZZ.gz  tmp/Bloom/$run/a0m.hits.gz | gawk -F '\t' '/^#/{next}/^ZZZZZ/{zz++;next;}{if(zz<1){if ($2+0>0){c=$1;ctr[c]=$2;cLn[c]=$3;}next;}{NN=10;split($1,aa,"__");type=aa[1];chrom=$11;n=$2; a=($12+$13)/2;if(cLn[chrom]+0==0)next;split($1,bb,"#");mult=0+bb[2];s=int(NN*a/cLn[chrom]);cc[chrom]++;nn[chrom,s,type]+=n*mult;types[type]+=n*mult;}END{printf("#Run\tType");for(c in cLn)for(i=0;i<NN;i++)printf("\t%s.%d",c,i+1);for(t in types){printf("\n%s\t%s",run,t);for(c in cLn)for(i=0;i<NN;i++)printf("\t%d",nn[c,i,t]);}printf("\n");}' run=$run > tmp/Bloom/$run/chrom2motif.txt
  endif

end

cat /home/mieg/SEQC2_2020/TARGET/CHROMS/Centromere_telomeres.ace | gawk '/^Map/{gsub(/\"/,"",$2);chrom=$2;}/^Centromere_telomeres/{t1=$2;t2=$5;c=($3+$4)/2;chLn[chrom]=t2;chCtr[chrom]=c;printf("%s\t%d\t%d\n",chrom,t2,c);}' > tmp/METADATA/chrom2centre2length.txt
  
set toto=RESULTS/Bloom/$MAGIC.chrom2motif.txt
echo -n "### $toto :" > $toto
date >> $toto
echo "## The reads containing at least 3 telomeric motifs were selected and remapped"
foreach run (`cat MetaDB/$MAGIC/GroupListSorted MetaDB/$MAGIC/RunsListSorted`)

  echo $run
  echo ' ' > tmp/Bloom/_a
  foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if ($1==run)print $2}END{print run}' run=$run`)
    if (-e tmp/Bloom/$run2/a0m.hits.gz) then
      zcat tmp/Bloom/$run2/a0m.hits.gz | gawk -F '\t' '/^#/{next;}{p=$1;gsub(/>/,"",p);gsub(/</,"",p);if(p==old)next;old=p;print}' >> tmp/Bloom/_a
      echo -n ".. $run2  "
      ls -ls tmp/Bloom/_a
    endif
  end

  cat tmp/METADATA/chrom2centre2length.txt ZZZZZ tmp/Bloom/_a | gawk -F '\t' '/^#/{next}/^ZZZZZ/{zz++;next;}{if(zz<1){if ($2+0>0){c=$1;nC++;CC[nC]=c;ctr[c]=$3;cLn[c]=$2;}next;}{NN=10;split($1,aa,"__");type=aa[1];chrom=$11;n=aa[2]; a=($12+$13)/2;if(cLn[chrom]+0==0)next;split($1,bb,"#");mult=0+bb[2];s=int(NN*a/cLn[chrom]);cc[chrom]++;nn[chrom,s,type]+=n*mult;types[type]+=n*mult;if(0)print chrom,s,type,n,mult;}}END{printf("#\t\t");for(iC=1;iC<=nC;iC++){c=CC[iC];printf("\t");for(i=0;i<NN;i++)printf("\tchrom %s part %d",c,i+1);}printf("#Run\tType\tNumber of motifs in mapped telomeric reads");for(iC=1;iC<=nC;iC++){c=CC[iC];printf("\t");for(i=0;i<NN;i++)printf("\tfrom %d to %d",1+(i)*cLn[c]/10,(i+1)*cLn[c]/10);}for(t in types){printf("\n%s\t%s\t%d",run,t,types[t]);for(iC=1;iC<=nC;iC++){c=CC[iC];printf("\t");for(i=0;i<NN;i++)printf("\t%d",nn[c,i,t]);}}printf("\n");}' run=$run >>  $toto
 ls -ls tmp/Bloom/_a
end
mv $toto _a 
cat _a | gawk '/^#/{print;}' | sort -u > $toto
cat _a | gawk '/^#/{next;}{print;}' |  sort -k 2,2 -k 1,1 >> $toto


##############################################
## preselection des genes A* et B* et des reads qui vont dans A* et B*

## select gene A*

foreach x (A B C )
  gunzip -c TARGET/Targets/hs.av.fasta.gz | gawk '/^>/{ok=0;if(toupper(substr($1,2,1))==x)ok=1;}{if(ok==1)print;}' x=$x > toto.av.$x.fasta
end

## select in run Rhs5337 the reads hitting gene A B C
set run=Rhs5337
foreach x (A)
  echo $x
  echo toto > toto.$run.$x.fastc
  foreach lane (`cat Fastc/$run/LaneList`)
     echo "--- $x $lane"
     gunzip -c tmp/COUNT/$lane.hits.gz | gawk -F '\t' '{if($8=="ET_av" && toupper(substr($11,1,1))==x){z=substr($1,1,length($1)-1);print z;}}' x=$x | sort -u > toto
     dna2dna -i Fastc/$lane.fastc.gz -I fastc -O fastc -select toto >>  toto.$run.$x.fastc
  end
end
