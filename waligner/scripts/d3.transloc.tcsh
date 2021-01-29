#!bin/tcsh -f

set run=$1
gunzip -c tmp/PHITS_genome/$run/*.overhangs.gz | bin/geneelements -newIntrons -minimalSupport 8 -non_gt_ag 8 -overhangLength 8 -transsplicing 8 -run $run > tmp/Transloc/$run/d3.$run.de_duo.txt

 cp  tmp/Transloc/$run/d3.$run.de_duo.txt RESULTS/SNV/TRANSLOC
 zcat tmp/Transloc/$run/d1.transloc.txt.gz > RESULTS/SNV/TRANSLOC/d1.$run.de_uno.txt

