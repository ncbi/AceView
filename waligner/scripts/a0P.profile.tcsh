#!bin/tcsh -f

set run=$1
set prefix=$2

   set title=`cat MetaDB/$MAGIC/RunTitleList | gawk -F '#@#' '{if($1 == run){ gsub (/:/,"_",$2);gsub (/\(/,"_",$2);gsub (/\)/,"_",$2);gsub (/\;/,"_",$2);gsub (/\,/,"_",$2);print $1; last;}}' run=$run`

   set fastctype=fastc
   foreach run2 (`cat MetaDB/$MAGIC/RunSolidList`)
      if ($run == $run2) set fastctype=csfastc
   end


echo "bin/dna2dna  -I fastc -minEntropy $minEntropy   -gnuplot 10000 -o tmp/Profiles/$run/readsBeforeAlignment  -title $run ATGC profile after removal of low entropy sequences and clipping of terminal N in $MAGIC $title"

# for solid, transform to ccfa and pretent it is fastc

if ($fastctype == csfastc) then
  gunzip -c Fastc/$run/$prefix.*.fastc.gz | head -20 | bin/dna2dna  -I csfastc -O ccfa  > tmp/Profiles/$run/test.csfastc
  gunzip -c Fastc/$run/$prefix.*.fastc.gz | head -20 | bin/dna2dna  -I ccfa -O ccfa > tmp/Profiles/$run/test.ccfa
  set n1=`cat tmp/Profiles/$run/test.csfastc | wc -l`
  set n2=`cat tmp/Profiles/$run/test.ccfa | wc -l`
  if ($n1 > $n2) then
    gunzip -c Fastc/$run/$prefix.*.fastc.gz | bin/dna2dna  -I $fastctype -O ccfa | bin/dna2dna  -I fastc -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run Solid color profile after removal of low entropy sequences in $MAGIC $prefix $title"
  else
    gunzip -c Fastc/$run/$prefix.*.fastc.gz | bin/dna2dna  -I fastc -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run Solid color profile after removal of low entropy sequences in $MAGIC $prefix $title"
  endif
  
else
  gunzip -c Fastc/$run/$prefix.*.fastc.gz | bin/dna2dna  -I $fastctype -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run ATGC profile after removal of low entropy sequences in $MAGIC $prefix $title"
endif

if (-e tmp/Profiles/$run/readsBeforeAlignment.parsed) \rm  tmp/Profiles/$run/readsBeforeAlignment.parsed

