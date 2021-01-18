#!bin/tcsh -f

set run=$1
set lane=$2
set target=$3

echo "bin/clipalign -i tmp/Unaligned/$lane.fastc.gz -t PARASITES/$MAGIC.parasites.fasta  -maxHit 30 -minAli 30 -errRateMax 6 -clipPolyA  -target_class b_bacteria -minEntropy 16 -o tmp/PHITS_$target/$lane -gzo"

  bin/clipalign -i tmp/Unaligned/$lane.fastc.gz -t PARASITES/$MAGIC.parasites.fasta  -maxHit 30 -minAli 30 -errRateMax 6 -clipPolyA  -target_class b_bacteria -minEntropy 16 -o tmp/PHITS_$target/$rlane -gzo
endif

gunzip -c tmp/PHITS_$target/$lane.hits.gz  | gawk -F '\t' '{p=$11;nn[p]+=$3;}END{for (k in nn)if(nn[k]>0)printf("%d\t%s\n",nn[k],k);}'  > tmp/PHITS_$target/$lane.p4


