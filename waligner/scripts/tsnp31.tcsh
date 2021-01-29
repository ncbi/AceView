#!/bin/tcsh

########################################################################################################################################################################

phaseTsnp3:

if (! -e tmp/TSNP/tsnp3.$MAGIC.val.sorted.txt) then
  echo "//" >  tmp/TSNP/tsnp3.$MAGIC.val.txt
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e  tmp/TSNP/$run/tsnp2.val.txt)  then
      cat tmp/TSNP/$run/tsnp2.val.txt | gawk -F '\t' '{split($2,aa,"___");split(aa[1],bb,"__");printf("%s\t%s\t%s\t%d\n", bb[2],bb[1],$1,$3+$4)}' >> tmp/TSNP/tsnp3.$MAGIC.val.txt
    endif
  end
  cat  tmp/TSNP/tsnp3.$MAGIC.val.txt | sort -V >  tmp/TSNP/tsnp3.$MAGIC.val.sorted.txt
  \rm tmp/TSNP/tsnp3.$MAGIC.val.txt 
endif

tace MetaDB <<EOF
  query find run file
  select -o r2st.txt t,r from r in @, t in r->sorting_title
  query find run file
  select -o r2t.txt t,r from r in @, t in r->title
EOF
cat g2st.txt r2st.txt > gr2st.txt
echo ZZZZZ >>  gr2st.txt
cat g2st.txt r2t.txt >> gr2st.txt

date
cat ZZZZZ ZZZZZ tmp/TSNP/tsnp3.$MAGIC.val.sorted.txt | gawk -F '\t' -f scripts/tsnp3.a.awk any=0 | sort -u  > tmp/TSNP/tsnp3.$MAGIC.val.tsf
date
cat ZZZZZ ZZZZZ tmp/TSNP/tsnp3.$MAGIC.val.sorted.txt | gawk -F '\t' -f scripts/tsnp3.a.awk any=1 | sort -u  >> tmp/TSNP/tsnp3.$MAGIC.val.tsf
date
cat gr2st.txt ZZZZZ tmp/TSNP/tsnp3.$MAGIC.val.sorted.txt | gawk -F '\t' -f scripts/tsnp3.a.awk any=2 | sort -u  >> tmp/TSNP/tsnp3.$MAGIC.val.tsf
date

############################
# HACK may 12 for the big deletions
cat tmp/TSNP/tsnp3.$MAGIC.val.tsf | gawk -F '\t' '/:Del/{split($1,aa,":");x=aa[2]+0;if(x < 55 || x > 78)next;print;}' > _a
cat gr2st.txt ZZZZZ _a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){r2st[$2]=$1;next;}}{if(zz==1){r2t[$2]=$1;next;}}{printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\n",$1,r2st[$2],r2t[$2],$2,$4,$8,$9);}' | sort -V > _b
# now we have the del, the sorting title, and the 3 counts
# we want to export for each run the sum of all var, the average of the RefDonor and the precentage 100*v/ sum v  + average ref donor

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=0 > _c1
set toto=$MAGIC.LongDel.frequency.txt
echo -n "### $toto : " > _d1
date >> _d1
cat _c1 | gawk '/^#/{print}' >> _d1
cat _c1 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d1
\cp _d1  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=9 > _c2
set toto=$MAGIC.LongDel.counts.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=1 > _c2
set toto=$MAGIC.LongDel.countsDel.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=2 > _c2
set toto=$MAGIC.LongDel.countsWproximal.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=3 > _c2
set toto=$MAGIC.LongDel.countsWdistal.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk forceX1=76 type=9 > _c2
set toto=$MAGIC.LongDel.tripleCount.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 3n >> _d2
\cp _d2  RESULTS/$toto

###########################
## all other deletions
cat tmp/TSNP/tsnp3.$MAGIC.val.tsf | gawk -F '\t' '/:DelIns/{next;}/:Del/{split($1,aa,":");x=aa[2]+0;if(x < 55 || x > 78)print;next;}' > _a
cat gr2st.txt ZZZZZ _a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){r2st[$2]=$1;next;}}{if(zz==1){r2t[$2]=$1;next;}}{printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\n",$1,r2st[$2],r2t[$2],$2,$4,$8,$9);}' | sort -V > _b
# now we have the del, the sorting title, and the 3 counts
# we want to export for each run the sum of all var, the everage of the RefDonor and the precentage 100*v/(sum(v) + average(ref donor))


cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=0 > _c2
set toto=$MAGIC.Deletions.frequency.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n -k 1,1r | gawk -F '\t' '{if($2 "x" $3 ==old)next;old=$2 "x" $3; print;}' >> _d2
\cp _d2  RESULTS/$toto


cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=1 > _c2
set toto=$MAGIC.Deletions.countDel.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n  -k 1,1r | gawk -F '\t' '{if($2 "x" $3 ==old)next;old=$2 "x" $3; print;}' >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=2 > _c2
set toto=$MAGIC.Deletions.countWproximal.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n -k 1,1r  | gawk -F '\t' '{if($2 "x" $3 ==old)next;old=$2 "x" $3; print;}' >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=3 > _c2
set toto=$MAGIC.Deletions.countWdistal.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n  -k 1,1r | gawk -F '\t' '{if($2 "x" $3 ==old)next;old=$2 "x" $3; print;}' >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=9 > _c2
set toto=$MAGIC.Deletions.tripleCount.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n -k 1,1r | gawk -F '\t' '{if($2 "x" $3 ==old)next;old=$2 "x" $3; print;}' >> _d2
\cp _d2  RESULTS/$toto

set toto=$MAGIC.AllDels.AnyRun.txt
echo -n "### $toto : " > $toto
date >> $toto
cat RESULTS/LongDel.tripleCount.txt | gawk '/^#/{n++;if(n>1)print}' | cut -f 1-10  > $toto
cat RESULTS/LongDel.tripleCount.txt RESULTS/Deletions.tripleCount.txt | gawk '/^#/{next;}{print}' | sort -k 5nr | cut -f 1-10 >> $toto


###########################
## short insertions
cat tmp/TSNP/tsnp3.$MAGIC.val.tsf | gawk -F '\t' '/:Ins/{k=split($1,aa,":");print}' > _a
cat gr2st.txt ZZZZZ _a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){r2st[$2]=$1;next;}}{if(zz==1){r2t[$2]=$1;next;}}{printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\n",$1,r2st[$2],r2t[$2],$2,$4,$5,$6);}' | sort -V > _b
# now we have the del, the sorting title, and the 3 counts

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=0 > _c2
set toto=$MAGIC.Insertions.frequency.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=1 > _c2
set toto=$MAGIC.Insertions.countIns.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=9 > _c2
set toto=$MAGIC.Insertions.tripleCount.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

###########################
## substitutions
cat tmp/TSNP/tsnp3.$MAGIC.val.tsf | gawk -F '\t' '/:Sub/{k=split($1,aa,":");print}' > _a
cat gr2st.txt ZZZZZ _a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){r2st[$2]=$1;next;}}{if(zz==1){r2t[$2]=$1;next;}}{printf("%s\t%s\t%s\t%s\t%d\t%d\t%d\n",$1,r2st[$2],r2t[$2],$2,$4,$5,$6);}' | sort -V > _b
# now we have the del, the sorting title, and the 3 counts
# we want to export for each run the sum of all var, the everage of the RefDonor and the precentage 100*v/(sum(v) + average(ref donor))


cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=100 > _c2
set toto=$MAGIC.Substitutions.frequency.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=101 > _c2
set toto=$MAGIC.Substitutions.count.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=102 > _c2
set toto=$MAGIC.Substitutions.reference.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

cat _b | gawk -F '\t' -f scripts/tsnp3.c.awk type=109 > _c2
set toto=$MAGIC.Substitutions.doubleCount.txt
echo -n "### $toto : " > _d2
date >> _d2
cat _c2 | gawk '/^#/{print}' >> _d2
cat _c2 | gawk '/^#/{next;}{print}' | sort -k 2n -k 3n >> _d2
\cp _d2  RESULTS/$toto

echo -n "### All SNV insertions deletions " >   RESULTS/MMM.any.Frequency.txt 
date >>   RESULTS/MMM.any.Frequency.txt 
cat RESULTS/MMM.*.frequency.txt | gawk '/^#{if(zz!=1)print ; next;}{n=1;exit;}'  >>  RESULTS/MMM.any.Frequency.txt 
cat RESULTS/MMM.*.frequency.txt | gawk '/^#{next;}{n=1;print}' | sort -k 4n -k 1,1  >>  RESULTS/MMM.any.Frequency.txt 
