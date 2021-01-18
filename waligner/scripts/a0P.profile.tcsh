#!bin/tcsh -f

set run=$1
set prefix=$2

   set title=`cat MetaDB/$MAGIC/Run2Title.txt | gawk -F '\t' '{gsub(/\"/,"",$0);}{if($1 == run){ gsub (/:/,"_",$3);gsub (/\(/,"_",$3);gsub (/\)/,"_",$3);gsub (/\;/,"_",$3);gsub (/\,/,"_",$3);if(length($3)<2)$3=run;print $3; last;}}' run=$run`

   set fastctype=fastc
   foreach run2 (`cat MetaDB/$MAGIC/RunSolidList`)
      if ($run == $run2) set fastctype=csfastc
   end

set minEntropy=0
if (-e Fastc/$run/Max_probe_length) then
  set m=`cat Fastc/$run/Max_probe_length | gawk '{if($1>0 && mx>$1/3)mx=$1/3;}END{print mx;}' mx=$minEntropy`
  set minEntropy=$m
endif

echo "bin/dna2dna  -I fastc -minEntropy $minEntropy   -gnuplot 1000 -o tmp/Profiles/$run/readsBeforeAlignment  -title $run ATGC profile after removal of low entropy sequences and clipping of terminal N in $MAGIC $title"

# for solid, transform to ccfa and pretent it is fastc

if ($fastctype == csfastc) then
  gunzip -c Fastc/$run/*.*.fastc.gz | head -20 | bin/dna2dna  -I csfastc -O ccfa  > tmp/Profiles/$run/test.csfastc
  gunzip -c Fastc/$run/*.*.fastc.gz | head -20 | bin/dna2dna  -I ccfa -O ccfa > tmp/Profiles/$run/test.ccfa
  set n1=`cat tmp/Profiles/$run/test.csfastc | wc -l`
  set n2=`cat tmp/Profiles/$run/test.ccfa | wc -l`
  if ($n1 > $n2) then
    gunzip -c Fastc/$run/*.*.fastc.gz | bin/dna2dna  -I $fastctype -O ccfa | bin/dna2dna  -I fastc -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run Solid color profile after removal of low entropy sequences in $MAGIC $prefix $title"
  else
    gunzip -c Fastc/$run/*.*.fastc.gz | bin/dna2dna  -I fastc -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run Solid color profile after removal of low entropy sequences in $MAGIC $prefix $title"
  endif
  
else
  gunzip -c Fastc/$run/*.*.fastc.gz | bin/dna2dna  -I $fastctype -minEntropy $minEntropy  -gnuplot 1000 -o tmp/Profiles/$run/$prefix.readsBeforeAlignment.1  -title "$run ATGC profile after removal of very low entropy sequences below 5 in $MAGIC $prefix $title"
endif

if (-e tmp/Profiles/$run/readsBeforeAlignment.parsed) \rm  tmp/Profiles/$run/readsBeforeAlignment.parsed

exit 0

 
