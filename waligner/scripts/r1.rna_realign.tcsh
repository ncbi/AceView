#!bin/tcsh -f

echo "scripts/r1.rna_realign.tcsh $*"

set run=$1
set lane=$2

set pair=""
foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
  if ($run == $run2) then
    set pair="-pair 500"
  endif
end
  
set geneRemap=" -geneRemap tmp/METADATA/mrnaRemap.gz "

# collect the unaligned reads

gunzip -c tmp/COUNT/$lane.hits.gz | gawk '/ET_av/{if ($22 == -1 || $22 == -10) print substr($1,1,length($1)-1);}' | sort -u >  tmp/RNA_editing/$lane.r1.orphans
dna2dna -I fastc -O fastc -gzo -i Fastc/$lane.fastc.gz -select tmp/RNA_editing/$lane.r1.orphans -o  tmp/RNA_editing/$lane.r1.orphans


 bin/clipalign -t TARGET/Targets/hs.av.fasta.gz -i  tmp/RNA_editing/$lane.r1.orphans.fastc.gz -minAli 35 -A2G -o tmp/RNA_editing/$lane.r1.av.A2G.MRNAH  -gzo -MRNAH  -seedLength 16 -maxHit 10
 bin/clipalign -t TARGET/Targets/hs.av.fasta.gz -i  tmp/Unaligned/$lane.fastc.gz -minAli 35 -A2G -o tmp/RNA_editing/$lane.unali.av.A2G.MRNAH  -gzo -MRNAH  -seedLength 16 -maxHit 10

gunzip -c tmp/RNA_editing/$lane.r1.av.A2G.MRNAH.hits.gz tmp/RNA_editing/$lane.unali.av.A2G.MRNAH.hits.gz  | gawk -F '\t' '{if($2 < 160)next ; n=split($17,aa,"r>g");if (n > 4)print;}' > tmp/RNA_editing/$lane.r1.av.5_A2G.hits
cat  tmp/RNA_editing/$lane.r1.av.5_A2G.hits |  gawk -F '\t' '{n[$9]+=$3;}END{for(g in n)if(n[g]>=1)printf("%s\t%d\n",g,n[g]);}' | sort -k 2nr >  tmp/RNA_editing/$lane.r1.av.5_A2G.genes
head -10 tmp/RNA_editing/$lane.r1.av.5_A2G.genes

touch tmp/RNA_editing/$lane.r1.done
exit 0




# realign in standard mode to clean up the pair assignment errors
 bin/clipalign -t TARGET/Targets/hs.av.fasta.gz -i  tmp/RNA_editing/$lane.r1.orphans.fastc.gz -minAli 35 -o tmp/RNA_editing/$lane.r1.av.normal.MRNAH -gzo -clipPolyA -MRNAH -seedLength 20
 gunzip -c  tmp/RNA_editing/$lane.r1.av.normal.MRNAH.hits.gz | bin/bestali -filter 35 -maxHit 3 -strategy RNA_seq $geneRemap $pair  -exportBest -gzo -o  tmp/RNA_editing/$lane.r1.av.normal.MRNAH.best
 gunzip -c  tmp/RNA_editing/$lane.r1.av.normal.MRNAH.best.hits.gz | gawk '{if($22 == -1 || $22 == -10) print substr($1,1,length($1)-1);}' | sort -u >  tmp/RNA_editing/$lane.r1.orphans2
 dna2dna -I fastc -O fastc -gzo -i  tmp/RNA_editing/$lane.r1.orphans.fastc.gz -select tmp/RNA_editing/$lane.r1.orphans2 -o  tmp/RNA_editing/$lane.r1.orphans2

# realign in RNA-edit mode
 bin/clipalign -t TARGET/Targets/hs.av.fasta.gz -i  tmp/RNA_editing/$lane.r1.orphans.fastc.gz -minAli 35 -A2G -o tmp/RNA_editing/$lane.r1.av.A2G.MRNAH  -gzo -MRNAH  -seedLength 16 -maxHit 10

# grep lines with aligned pairs
   gunzip -c  tmp/RNA_editing/$lane.r1.av.A2G.MRNAH.hits.gz | bin/bestali -filter 35 -maxHit 10 -strategy RNA_seq $geneRemap $pair  -exportBest -gzo -o  tmp/RNA_editing/$lane.r1.av.A2G.MRNAH.best

gunzip -c tmp/RNA_editing/$lane.r1.av.normal.MRNAH.hits.gz | gawk -F '\t' '{if($22 == -10) next; if($2 > 160)print;}' | cut -f 1 | sort -u | gzip >  tmp/RNA_editing/$lane.r1.av.normal.MRNAH.list.gz
gunzip -c  tmp/RNA_editing/$lane.r1.av.normal.MRNAH.list.gz ZZZZZ.gz tmp/RNA_editing/$lane.r1.av.A2G.MRNAH.best.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next}}{if(ok[$1] == 1 || $22 == -10) next; if($2 > 160)print;}' | wc

gunzip -c tmp/RNA_editing/$lane.r1.av.normal.MRNAH.list.gz ZZZZZ.gz tmp/RNA_editing/$lane.r1.av.A2G.MRNAH.best.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next}}{if(ok[$1] == 1 || $22 == -10) next; if($2 > 160)print;}' | gawk -F '\t' '{n=split($17,aa,"r>g");if (n > 4)print;}' | cut -f 9 | tags | sort -k 2n



 bin/clipalign -t TARGET/Targets/hs.genome.fasta.gz -i tmp/Unaligned/Rhs4200/any.fastc -splice -minAli 24 -A2Gf -o tmp/Unaligned/Rhs4200/A2Gf.genome
touch tmp/RNA_editing/$lame.r1.done
exit 0

 bin/clipalign -t TARGET/Targets/hs.genome.fasta.gz -i tmp/Unaligned/Rhs4200/any.fastc -splice -minAli 24 -A2Gf -o tmp/Unaligned/Rhs4200/A2Gf.genome
 bin/clipalign -t TARGET/Targets/hs.genome.fasta.gz -i tmp/Unaligned/Rhs4200/any.fastc -splice -minAli 24 -A2Gf -o tmp/Unaligned/Rhs4200/A2Gf.genome


cat tmp/Unaligned/Rhs4200/[AT]2[GC]?.genome.hits | gawk -F '\t' '{if($1==old)next;old=$1;if($2>160 && $10==1)print}' | sort -k 11,11 -k 12,12n | wc

set toto=tmp/RNA_editing/$MAGIC.run2gene.r2g 

#####################
## histo of rate of a>g editing extracted from BRS, CF :
#      s2r.a2g_collect.tcsh
