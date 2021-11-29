#!bin/tcsh -f

set phase=$1
set run=$2
set v1=$3
set v2=$4

# clip the reads, count the tags in range [18,35]
# call high the tags over 2*10^-6 , mid 2*10^-5  low 2*10^-4 (approximately 1000, 100, 10 reads per sequence)

if ($phase == clip) then
echo "$run $v1 $v2"
 if (! -e tmp/ClippedFastc/$run/f.filtered.1.tc) then
   # clip and count the clipped reads, often they merge better than the raw reads
   # filter at minimal coverage 10 and length range [18, 35]
    set toto=tmp/ClippedFastc/$run/f.filtered
   gunzip -c Fastc/$run/*.fastc.gz | gawk '/^>/{n=split($1,aa,"#");mult=aa[2]+0;if(mult==0)mult=1;next}{n=split($1,aa,"><");for(i=1;i<=1;i++)printf("%s\t%d\n",aa[i],mult);}' > $toto.F.tc
   gunzip -c Fastc/$run/*.fastc.gz | gawk '/^>/{n=split($1,aa,"#");mult=aa[2]+0;if(mult==0)mult=1;next}{n=split($1,aa,"><");for(i=2;i<=n;i++)printf("%s\t%d\n",aa[i],mult);}' > $toto.R.tc
   if ($v1 != "X") then
     mv $toto.F.tc $toto.F1.tc
     dna2dna -I tc -O tc -rightClipOn $v1 -i $toto.F1.tc -o $toto.F
   endif
   touch $toto.R.tc
   if ($v2 != "X") then
     mv $toto.R.tc $toto.R1.tc
     dna2dna -I tc -O tc -rightClipOn $v2 -i $toto.R1.tc | dna2dna -I tc -O tc -o $toto.R -complement
   endif
   cat $toto.[FR].tc  | dna2dna -I tc -O tc -minLength 18 -maxLength 35  | sort -k 2nr >   $toto.1.tc
   \rm $toto.[FR]*

   endif
   set n18=`cat tmp/ClippedFastc/$run/f.filtered.1.tc  | gawk -F '\t' '{t += $2}END{print t}'`
   set seuil=`echo $n18 | gawk '{printf("%d",$1/10000);}'`
   cat tmp/ClippedFastc/$run/f.filtered.1.tc  | gawk -F '\t' '{if(1000000*$2 >= t && $2>= 10) print}' t=$n18 > tmp/ClippedFastc/$run/f.filtered.10.tc
   cat tmp/ClippedFastc/$run/f.filtered.10.tc  | gawk -F '\t' '{if(100000*$2 >= t) print}' t=$n18 > tmp/ClippedFastc/$run/f.filtered.100.tc
   cat tmp/ClippedFastc/$run/f.filtered.100.tc  | gawk -F '\t' '{if(10000*$2 >= t) print}' t=$n18 > tmp/ClippedFastc/$run/f.filtered.1000.tc

    cat tmp/ClippedFastc/$run/f.filtered.10.tc | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.10 
    cat tmp/ClippedFastc/$run/f.filtered.100.tc | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.100 
    cat tmp/ClippedFastc/$run/f.filtered.1000.tc  | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.1000 


   echo "Ali $run"  > tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.1.tc | gawk -F '\t' '{if (t>=0) {s++; t+= $2;}} END {printf ("N_18_35 -10 %d %d\n", s, t);}' >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.10.tc | gawk -F '\t' '{s++; t+= $2; } END {printf ("N_18_35 -6 %d %d\n", s, t);}' >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.100.tc | gawk -F '\t' '{s++; t+= $2; } END {printf ("N_18_35 -5 %d %d\n", s, t);}' >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.1000.tc | gawk -F '\t' '{s++; t+= $2; } END {printf ("N_18_35 -4 %d %d\n", s, t);}' >> tmp/ClippedFastc/$run/f.stats.ace

   cat tmp/ClippedFastc/$run/f.filtered.10.tc | gawk -F '\t' '{ln = length($1); n = $2 ; nln[ln] += n ; if(ln > lnMax) lnMax = ln ; jj=0 ; k = 10; while (n >= k) { z[jj]+=n ; jj++; k *= 10 ;}}END{printf("Clipped_multiplicity") ;jj = 0 ; k = 10 ; while (z[jj] > 0) {printf (" %d \"seen %d times\" ", z[jj], k); k *= 10 ;jj++;} printf("\n") ; for(ln=1;ln<=lnMax ; ln++) if (nln[ln] > 0) printf("Preclipped_length %d %d\n", ln, nln[ln]) ;} '  >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.count ZZZZZ tmp/ClippedFastc/$run/f.filtered.100 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^Tags_processed/{nAll = $2; next;}{if(zz < 1)next;}{n=$2;if(50 *n >= nAll)printf("High_short %s %d\n", $1, n);}' | sort -k 3nr >> tmp/ClippedFastc/$run/f.stats.ace
   echo  >> tmp/ClippedFastc/$run/f.stats.ace

 endif
endif

exit 0

if ($phase == leming) then
   foreach run (`cat MetaDB/$MAGIC/RunsList`)
     cat leming.s1.tc  ZZZZZ tmp/ClippedFastc/$run/f.filtered.1.tc | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){nam[$1]=$3;ok[$1]=1;next;}}{if(ok[$1]>0)nn[$1]+=$2;}END{for(k in ok)if(nn[k]>0)printf("%s\t%s\t%d\n",k,run,nn[k]);}' run=$run > tmp/ClippedFastc/$run/leming.s1.counts
   end
   cat tmp/ClippedFastc/*/leming.s1.counts | sort | gawk -F '\t' '{g=$1;if(g!=old)printf("\nGene %s\n",$1);printf("Run_U %s %d\n",$2,$3);}' > tmp/ClippedFastc/leming.s1.ace


endif
exit 0



 
