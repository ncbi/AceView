#!bin/tcsh -f

set phase=$1
set run=$2
set v1=$3


if ($phase == clip) then

 if (! -e tmp/ClippedFastc/$run/f.tc) then
   # clip
   gunzip -c Fastc/$run/*.fastc.gz | dna2dna -I fastc -O fastc -rightClipOn $v1 -o  tmp/ClippedFastc/$run/f
   # count the clipped reads, often they merge better than the raw reads
   dna2dna -i tmp/ClippedFastc/$run/f.fastc -I fastc -O tc -o tmp/ClippedFastc/$run/f 
 endif

 if (! -e tmp/ClippedFastc/$run/f.count) then
   dna2dna -i tmp/ClippedFastc/$run/f.tc -I tc -count -o tmp/ClippedFastc/$run/f
 endif

 if (! -e  tmp/ClippedFastc/$run/f.filtered.10) then
   # filter at minimal coverage and length range [18, 35]
   cat  tmp/ClippedFastc/$run/f.tc | gawk -F '\t' '{if ($2 < minCount) next; k=length($1);if(k>=18 && k <= 35 ) printf ("%s\t%d\n", $1, $2);}' run=$run minCount=10 | bin/dna2dna -I tc -O tc -minEntropy 12 > tmp/ClippedFastc/$run/f.filtered.10
   cat tmp/ClippedFastc/$run/f.filtered.10  | gawk -F '\t' '{if ($2 >= 100)print}' >  tmp/ClippedFastc/$run/f.filtered.100
   cat tmp/ClippedFastc/$run/f.filtered.100  | gawk -F '\t' '{if ($2 >= 1000)print}' >  tmp/ClippedFastc/$run/f.filtered.1000
  endif
 if (! -e  tmp/ClippedFastc/$run/f.filtered.10.fastc) then
    cat tmp/ClippedFastc/$run/f.filtered.10 | dna2dna -I tc -O fastc -o tmp/ClippedFastc/$run/f.filtered.10 
 endif
 if (! -e  tmp/ClippedFastc/$run/f.stats.ace) then
   echo "Ali $run"  > tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.filtered.10 | gawk -F '\t' '{ln = length($1); n = $2 ; nln[ln] += n ; if(ln > lnMax) lnMax = ln ; jj=0 ; k = 10; while (n >= k) { z[jj]+=n ; jj++; k *= 10 ;}}END{printf("Clipped_multiplicity") ;jj = 0 ; k = 10 ; while (z[jj] > 0) {printf (" %d \"seen %d times\" ", z[jj], k); k *= 10 ;jj++;} printf("\n") ; for(ln=1;ln<=lnMax ; ln++) if (nln[ln] > 0) printf("Preclipped_length %d %d\n", ln, nln[ln]) ;} '  >> tmp/ClippedFastc/$run/f.stats.ace
   cat tmp/ClippedFastc/$run/f.count ZZZZZ tmp/ClippedFastc/$run/f.filtered.100 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^Tags_processed/{nAll = $2; next;}{if(zz < 1)next;}{n=$2;if(50 *n >= nAll)printf("High_short %s %d\n", $1, n);}' | sort -k 3nr >> tmp/ClippedFastc/$run/f.stats.ace
   echo  >> tmp/ClippedFastc/$run/f.stats.ace
 endif 

endif

exit 0

