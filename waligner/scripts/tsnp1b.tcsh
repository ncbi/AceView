#!/bin/tcsh -f
set run=$1
set zone=$2 
set minCount=$3
set minCover=$4
set minFrequency=$5

bin/tsf -i  tmp/TSNP/$run/$zone/_r --group $run -o  tmp/TSNP/$run/$zone/tsnp1.2.0.pre_deUno
cat tmp/TSNP/$run/$zone/tsnp1.2.0.pre_deUno.tsf | gawk -F '\t' '/^#/{print}{if ($7 > minCount && $8 > minCount && 200*$7 > $minSnpFrequency * ($9+$10))print}' minCover=$minCover minCount=$minCount minFrequency=$minFrequency > tmp/TSNP/$run/$zone/tsnp1.2.0.deUno.tsf 

\rm   tmp/TSNP/$run/$zone/_r tmp/TSNP/$run/$zone/tsnp1.2.0.pre_deUno.tsf
