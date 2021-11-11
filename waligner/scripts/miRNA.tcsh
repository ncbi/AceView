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
   gunzip -c Fastc/$run/*.fastc.gz | gawk '/^>/{n=split($1,aa,"#");mult=aa[2]+0;if(mult==0)mult=1;next}{n=split($1,aa,"><");for(i=1;i<=1;i++)printf("%s\t%d\n",aa[i],mult);}' | dna2dna -I tc -O tc -rightClipOn $v1 | dna2dna -I tc -O tc -minLength 18 -maxLength 35 -o  tmp/ClippedFastc/$run/f.filtered.1
   if ($v2 != "X") then
     gunzip -c Fastc/$run/*.fastc.gz | gawk '/^>/{n=split($1,aa,"#");mult=aa[2]+0;if(mult==0)mult=1;next}{n=split($1,aa,"><");for(i=2;i<=n;i++)printf("%s\t%d\n",aa[i],mult);}' | dna2dna -I tc -O tc -rightClipOn $v2 |  dna2dna -I tc -O tc -minLength 18 -maxLength 35  -o tmp/ClippedFastc/$run/f.filtered.1
   endif
   set t=`cat tmp/ClippedFastc/$run/f.filtered.1.tc  | gawk -F '\t' '{t += $2}END{print t}'`
   
   cat tmp/ClippedFastc/$run/f.filtered.1.tc  | gawk -F '\t' '{if(2*500000*$2 >= t && t>= 100) print}' t=$t > tmp/ClippedFastc/$run/f.filtered.10.tc
   cat tmp/ClippedFastc/$run/f.filtered.10.tc  | gawk -F '\t' '{if(2*50000*$2 >= t) print}' t=$t > tmp/ClippedFastc/$run/f.filtered.100.tc
   cat tmp/ClippedFastc/$run/f.filtered.100.tc  | gawk -F '\t' '{if(2*5000*$2 >= t) print}' t=$t > tmp/ClippedFastc/$run/f.filtered.1000.tc
 endif

 if (! -e  tmp/ClippedFastc/$run/f.filtered.10.fastc) then
    cat tmp/ClippedFastc/$run/f.filtered.10.tc | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.10 
    cat tmp/ClippedFastc/$run/f.filtered.100.tc | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.100 
    cat tmp/ClippedFastc/$run/f.filtered.1000.tc  | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.1000 
 endif

 if (! -e  tmp/ClippedFastc/$run/f.stats.ace) then
   echo "Ali $run"  > tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.10.tc | gawk -F '\t' '{ln = length($1); n = $2 ; nln[ln] += n ; if(ln > lnMax) lnMax = ln ; jj=0 ; k = 10; while (n >= k) { z[jj]+=n ; jj++; k *= 10 ;}}END{printf("Clipped_multiplicity") ;jj = 0 ; k = 10 ; while (z[jj] > 0) {printf (" %d \"seen %d times\" ", z[jj], k); k *= 10 ;jj++;} printf("\n") ; for(ln=1;ln<=lnMax ; ln++) if (nln[ln] > 0) printf("Preclipped_length %d %d\n", ln, nln[ln]) ;} '  >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.count ZZZZZ tmp/ClippedFastc/$run/f.filtered.100 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^Tags_processed/{nAll = $2; next;}{if(zz < 1)next;}{n=$2;if(50 *n >= nAll)printf("High_short %s %d\n", $1, n);}' | sort -k 3nr >> tmp/ClippedFastc/$run/f.stats.ace
   echo  >> tmp/ClippedFastc/$run/f.stats.ace
 endif 

endif

exit 0



 
