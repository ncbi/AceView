#!/bin/tcsh -f

set phase=$1
set run=$2

echo hello from virus.tcsh $phase
if ($phase == vir1) then

  foreach lane (`cat Fastc/$run/LaneList`)
    if (-e tmp/COUNT/$lane.hits.gz && ! -e  tmp/VIRUS/$lane.v_virus.hits) then
      gunzip -c tmp/COUNT/$lane.hits.gz | gawk -F '\t' '{if($8 == "v_virus" || $8 == "b_bacteria" || $8 == "D_transposon") print > "tmp/VIRUS/" lane "." $8 ".hits" ;}' lane=$lane
    endif
  end
  if (! -e  tmp/VIRUS/$run/vir1.virus.count) then
      ls  tmp/VIRUS/$run/*.v_virus.hits >  tmp/VIRUS/$run/_f
      bin/bestali -autoH -target_class v_virus -inFileList  tmp/VIRUS/$run/_f -run $run > tmp/VIRUS/$run/vir1.virus.count 
      ls  tmp/VIRUS/$run/*.b_bacteria.hits >  tmp/VIRUS/$run/_f
      bin/bestali -autoH -target_class b_bacteria  -inFileList  tmp/VIRUS/$run/_f -run $run > tmp/VIRUS/$run/vir1.bacteria.count
      ls  tmp/VIRUS/$run/*.D_transposon.hits >  tmp/VIRUS/$run/_f
      bin/bestali -autoH -target_class D_transposon  -inFileList  tmp/VIRUS/$run/_f -run $run > tmp/VIRUS/$run/vir1.transposon.count
  endif

  goto phaseLoop
endif


if ($phase == vir2) then

## associate the Accession to its full title
  if (1) then
    gunzip -c TARGET/Targets/$species.virus.fasta.gz     | gawk '/^>/{n=split($1,aa,"|");s=substr($0,2);if(n>=4){printf("%s\t",aa[4]);print s;}else printf("%s\t%s\n",s,s);}' > tmp/VIRUS/virus_names.txt
    gunzip -c TARGET/Targets/$species.bacteria.fasta.gz  | gawk '/^>/{n=split($1,aa,"|");s=substr($0,2);if(n>=4){printf("%s\t",aa[4]);print s;}else printf("%s\t%s\n",s,s);}' > tmp/VIRUS/bacteria_names.txt
    gunzip -c TARGET/Targets/$species.transposon.fasta.gz  | gawk '/^>/{n=split($1,aa,"|");s=substr($0,2);if(n>=4){printf("%s\t",aa[4]);print s;}else {gsub(/\t/," ",s);printf("%s\t%s\n",substr($1,2),s);}}' > tmp/VIRUS/transposon_names.txt
  endif


#####
## create a complete table using geneindex.c

  set toto=tmp/VIRUS/$MAGIC.virus_bacteria
  echo ' ' > $toto.counts.ace
  if (-e $toto.777) \rm $toto.777
  touch  $toto.777
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    foreach target (virus bacteria transposon)
      if (-e  tmp/VIRUS/$run/vir1.$target.count) then
        cat  tmp/VIRUS/$run/vir1.$target.count | gawk  '/^#/{next;}{g=$2;n=split(g,aa,"|");if(n>=4)g=aa[4];printf("Gene \"%s\"\nRun_U \"%s\" 0.00 %d seqs %d tags %d kb %d raw_tags\n\n",g,$1,$3,$3,$3/10,$4);}' >> $toto.counts.ace
        cat  tmp/VIRUS/$run/vir1.$target.count  >> $toto.777
      endif
    end
  end

# We wish to add the names
  echo ZZZZZ > ZZZZZ
  cat tmp/VIRUS/virus_names.txt tmp/VIRUS/bacteria_names.txt tmp/VIRUS/transposon_names.txt ZZZZZ $toto.777 | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){seq=$1;title[seq]=$2;next;}}{if(zz==1){seq=$2;n=split(seq,aa,"|");if(n>=4)seq=aa[4];n1[seq]+=$3;n2[seq]+=$4;}}END{for (seq in title)printf("%d\t%d\t%s\t%s\n",n1[seq],n2[seq],seq,title[seq]);}' | sort -k 1,1nr >  $toto.virus_bacteria_names.sorted

  cat MetaDB/$MAGIC/runs.ace > $toto.info.ace
  cat $toto.virus_bacteria_names.sorted | gawk -F '\t' '{i=index($4," ") ; g=$3;n=split(g,aa,"|");if(n>=4)g=aa[4];if(i<0)i=0;printf("Gene \"%s\"\nTitle \"%s\"\nGeneId \"%s\"\nTargeted\n\n",g,substr($4,i+1),$4);}' >> $toto.info.ace
  gunzip -c TARGET/Targets/$species.bacteria.TM.txt.gz TARGET/Targets/$species.virus.TM.txt.gz |  gawk -F '\t' '{g=$1;n=split(g,aa,"|");if(n>=4)g=aa[4];printf("Gene \"%s\"\nLength %d\n\n", g, $2);}' >>  $toto.info.ace

# we wish to sort the runs, then the group, but we add the sublibs so that their count is parsed, but anyway geneindex.c considers sublibs as private
  cat  MetaDB/$MAGIC/RunsListSorted  MetaDB/$MAGIC/GroupListSorted MetaDB/$MAGIC/RunList  > tmp/VIRUS/$MAGIC.GroupsRunsListSorted
  bin/geneindex -export at -deepGene $toto.counts.ace -u -runList tmp/VIRUS/$MAGIC.GroupsRunsListSorted -runAce $toto.info.ace -o $toto 
  if (-d RESULTS/VIRUS && ! -d RESULTS/Microbiome) mv RESULTS/VIRUS RESULTS/Microbiome
  if (! -d RESULTS/Microbiome) mkdir RESULTS/Microbiome
  set toto2=RESULTS/Microbiome/$MAGIC.virus_bacteria.reads_aligned_per_run.txt
  echo -n "## $MAGIC.virus_bacteria.reads_aligned_per_run.txt " > $toto2
  cat $toto.reads_aligned_per_gene.txt  | gawk '/^#/{print}'  >> $toto2
  echo >> $toto2
  cat $toto.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr > $toto.999
  cat  tmp/VIRUS/virus_names.txt ZZZZZ $toto.999    | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;ok[g]=1;nam[g]=$2;next;}}{if(ok[$1]>=1){ok[$1]=2;print;}}END{for(k in ok){if(ok[k]==1)printf("%s\t%s\n",k,nam[k]) ;}}' | scripts/tab_sort -k 8,8nr >> $toto2
  echo "\n\n" >> $toto2
  cat  tmp/VIRUS/bacteria_names.txt ZZZZZ $toto.999 | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;ok[g]=1;nam[g]=$2;next;}}{if(ok[$1]>=1){ok[$1]=2;print;}}END{for(k in ok){if(ok[k]==1)printf("%s\t%s\n",k,nam[k]) ;}}'  | scripts/tab_sort -k 8,8nr >> $toto2
  echo "\n\n" >> $toto2
  cat  tmp/VIRUS/transposon_names.txt ZZZZZ $toto.999 | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;ok[g]=1;nam[g]=$2;next;}}{if(ok[$1]>=1){ok[$1]=2;print;}}END{for(k in ok){if(ok[k]==1)printf("%s\t%s\n",k,nam[k]) ;}}'  | scripts/tab_sort -k 8,8nr >> $toto2
  \rm $toto.999 $toto.777
endif

######
## create a .ace file : Ali->High_genes Virus/baceria counting the 5 bests
  set toto=tmp/VIRUS/$MAGIC.virus_bacteria.ace
  set toto4=tmp/VIRUS/$MAGIC.virus_bacteria.high_genes.ace
  echo ' ' > $toto4
  cat $toto | gawk '/^Gene/{g=$2;gsub(/\"/,"",g);i=split(g,aa,"|");if (i>=4)g=aa[4];gg[g]=1;next;}{ok=0;}/^Run_U/{ok=1;}/^Group_U/{ok=1;}{if(ok==1){printf("%s\t%s\t%s\t%s\t%s\n",$2,g,$4,$6,$8);}}' | sort -k 1,1 -k 3,3nr > $toto4.1

  foreach target (virus bacteria transposon)
    cat tmp/VIRUS/$target'_'names.txt ZZZZZ $toto4.1 |  gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;nam[$1]=$2;next;}}{if(ok[$2]==1)printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,nam[$2]);}' > $toto4.2
    cat  $toto4.2 | gawk -F '\t' '{if ($1 != old){n=0;old=$1;old4=0;}n++;if(n>5 && $4 < old4)next;old4=$4;}{print}'  > $toto4.3
    cat $toto4.3 ZZZZZ $toto4.3 |  gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){n[$1]++;k[$1]=$3;next;}if(n[$1]<=5||k[$1]<$3)print}'  > $toto4.4

    cat  $toto4.4 | gawk -F '\t' '{if($1 != old){printf("\nAli %s\n-D High_genes %s\n",$1,target);old=$1;}printf("High_genes %s \"%s\" %d \"%s\"\n",target,$2,$4,$6);}END{printf("\n");}' target=$target >> $toto4
  end

 echo " pparse  $toto4" | bin/tacembly MetaDB -noprompt
 \rm $toto4.*

goto phaseLoop
exit 0

cat <<EOF > toto
NC_001806.2 gi|820945227|ref|NC_001806.2| Human herpesvirus 1 strain 17, complete genome
NC_001798.2 gi|820945149|ref|NC_001798.2| Human herpesvirus 2 strain HG52, complete genome
NC_001348.1 gi|9625875|ref|NC_001348.1| Human herpesvirus 3, complete genome
NC_009334.1 gi|139424470|ref|NC_009334.1| Human herpesvirus 4, complete genome
NC_007605.1 gi|82503188|ref|NC_007605.1| Human herpesvirus 4 complete wild type genome
NC_006273.2 gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome
NC_001664.2 gi|224020395|ref|NC_001664.2| Human herpesvirus 6A, complete genome
NC_000898.1 gi|9633069|ref|NC_000898.1| Human herpesvirus 6B, complete genome
NC_001716.2 gi|51874225|ref|NC_001716.2| Human herpesvirus 7, complete genome
NC_009333.1 gi|139472801|ref|NC_009333.1| Human herpesvirus 8, complete genome
EOF
cat toto | gawk '{print $1}' >  tmp/VIRUS/herpes.list

# export the intron support
if (-e tmp/VIRUS/$MAGIC.intron.pre) \rm tmp/VIRUS/$MAGIC.intron.pre
foreach run (`cat MetaDB/$MAGIC/RunList`)
  if (-d tmp/PHITS_virus/$run ) then
    gunzip -c tmp/PHITS_virus/$run/*.introns.gz | gawk -F '\t' '/^#/{next}{n[$3 "\t" $5 "\t" $7 "\t" $9]+=$11;}END{for(k in n) printf ("%s\t%d\n",k,n[k])}' | sort -k 1 -k 2n | gawk -F '\t' '{if($2<$3){dx=$3-$2+1;s="+";}else{dx=$2-$3+1;s="-";}if($5>0)printf("%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,$2+0,$3+0,s,dx,$4,$5);}' >> tmp/VIRUS/$MAGIC.intron.pre
  endif
end
cat  tmp/VIRUS/$MAGIC.intron.pre | gawk -F '\t' '{n[$1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6] += $7;}END{for (k in n)printf("%s\t%d\n",k,n[k]);}' | sort -k 1,1 -k 7nr > tmp/VIRUS/$MAGIC.intron.pre2
 
if (-e tmp/VIRUS/$MAGIC.intron.preB) \rm tmp/VIRUS/$MAGIC.intron.preB
foreach run (`cat MetaDB/$MAGIC/RunList`)
  if (-d tmp/PHITS_bacteria/$run ) then
    gunzip -c tmp/PHITS_bacteria/$run/*.introns.gz | gawk -F '\t' '/^#/{next}{n[$3 "\t" $5 "\t" $7 "\t" $9]+=$11;}END{for(k in n) printf ("%s\t%d\n",k,n[k])}' | sort -k 1 -k 2n | gawk -F '\t' '{if($2<$3){dx=$3-$2+1;s="+";}else{dx=$2-$3+1;s="-";}if($5>0)printf("%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,$2+0,$3+0,s,dx,$4,$5);}' >> tmp/VIRUS/$MAGIC.intron.preB
  endif
end
cat  tmp/VIRUS/$MAGIC.intron.preB | gawk -F '\t' '{n[$1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6] += $7;}END{for (k in n)printf("%s\t%d\n",k,n[k]);}' | sort -k 1,1 -k 7nr > tmp/VIRUS/$MAGIC.intron.pre2B
 
echo -n "## Support for introns or small deletions in viruses in project $MAGIC : " > RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
date >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
echo "## Introns are discovered by mapping discontinuously using Magic on a number of pre-selected viruses" >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
cat  tmp/VIRUS/virus_names.txt ZZZZZ  tmp/VIRUS/$MAGIC.intron.pre2 |  gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){split($2,aa," ");ok[aa[1]]=$2;next}}BEGIN{printf("# Accession\tTitle\tCoordinate of intron first base on accession\tCoordinate of intron last base\tIntron boundary\tIntron length (nt)\tStrand\tNumber of reads supporting the intron\n");}{if($1 != old)printf("\n");old=$1;printf("%s\t%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,ok[$1],$2,$3,$6,$5,$4,$7);}' >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt


# extract the herpesvirus names and export their global coverage
cat  tmp/VIRUS/virus_names.txt | grep herpes | cut -f 1 >  tmp/VIRUS/herpes.list
foreach v (`cat tmp/VIRUS/herpes.list`)
  echo ".... virus.tcsh: cumulating the wiggles : $v"
  echo >  tmp/VIRUS/_h
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    cat tmp/VIRUS/$run/*.vir1.hits | gawk -F '\t' '{split ($11,aa,"|");if(aa[4] == v) print;}' v=$v >>  tmp/VIRUS/_h
  end
  set n=`wc -l  tmp/VIRUS/_h | gawk '{print $1}' `
  if ($n > 1) then
    wiggle -i  tmp/VIRUS/_h -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v -force_unique 
  endif
end

foreach v (`cat tmp/VIRUS/herpes.list`)
  cat tmp/VIRUS/$MAGIC.$v.BV | gawk '/^variabl/{print ; print "1\t1" ; next;}{print}' |  wiggle -I BV -O BF -o  tmp/VIRUS/$MAGIC.$v
end

\rm $MAGIC.toto
foreach v (`cat tmp/VIRUS/herpes.list`)
  cat  tmp/VIRUS/$MAGIC.$v.BF  | scripts/transpose >> $MAGIC.toto
end


 set toto=RESULTS/Microbiome/$MAGIC.herpes.BF
 echo -n '## ' > $toto
 date >> $toto
 echo "## Coverage plots of selected viruses in all samples of project $MAGIC" >> $toto
 cat $MAGIC.toto  | scripts/transpose | gawk -F '\t' '{ln++; if (ln<3 || ln == 4)next;if(ln == 3) {printf ("Coordinate");for(i=1;i<=NF;i++){split($i,aa," ") ;printf("\t%s",substr(aa[2],7));}printf ("\n");next;}printf("%d\t", 10 * (ln-4)) ; print}' >> $toto

 cat $toto | head -2 > $toto.1
 cat $toto | head -3 | tail -1 | transpose > $toto.2
 cat tmp/VIRUS/virus_names.txt ZZZZZ $toto.2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$2]=1;next}}{k0=$1;for(k in ok){if(index(k,$1)>0)k0=k;}print k0}' | transpose  >> $toto.1
 cat $toto | tail -n +4 >> $toto.1
 mv $toto.1 $toto
 \rm $toto.2


\rm  $MAGIC.toto.Rh $MAGIC.toto



# extract the herpesvirus names and export their global coverage

foreach v (NZ_CP007181.1 AL111168 HE978252.1 CP000814.1)
  echo ".... virus.tcsh: cumulating the wiggles : $v"
  echo >  tmp/VIRUS/_h
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    cat tmp/VIRUS/$run/*.vir1.hits | gawk -F '\t' '{split ($11,aa,"|");if($11 == v || aa[4] == v) print;}' v=$v >>  tmp/VIRUS/_h
  end
  set n=`wc -l  tmp/VIRUS/_h | gawk '{print $1}' `
  if ($n > 1) then
    wiggle -i  tmp/VIRUS/_h -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v -force_unique 
  endif
end



phaseLoop:
 echo done
