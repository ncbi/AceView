#!bin/tcsh -f
set lane=$1

gunzip -c Fastc/$lane.fastc.gz | gawk '{nn++;if(nn%2==0){i=index($1,"><");$1=substr($1,1,i-1);}print}' > tmp/JUNK/$lane.F.fastc
gunzip -c Fastc/$lane.fastc.gz | gawk '{nn++;if(nn%2==0){i=index($1,"><");$1=substr($1,i+2);}print}' > tmp/JUNK/$lane.R.fastc

bin/dna2dna -I fastc -i  tmp/JUNK/$lane.F.fastc -O tc | bin/dna2dna -I tc -O fastc -gzo -o tmp/JUNK/$lane.FF -minEntropy 16
bin/dna2dna -I fastc -i  tmp/JUNK/$lane.R.fastc -O tc | bin/dna2dna -I tc -O fastc -gzo -o tmp/JUNK/$lane.RR -minEntropy 16

