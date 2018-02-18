#!bin/tcsh -f


clipalign -i tmp/NEWHITS_snp/all.edited_sequence.m.snp.fasta.gz -t TARGET/Targets/$species.genome.fasta.gz -nonBest -errMax 2 -seedLength 32 -seedOffset 1 -seedShift 20 -o tmp/NEWHITS_snp/m.snp2genome -gzo

set toto=RESULTS/SNV/ambiguous_SNV.txt
date > $toto
echo "Look for sequences centered on a SNV mapping several times to the genome: " >> $toto
echo -n "How many have more than 1 hit " >> $toto
gunzip -c tmp/NEWHITS_snp/m.snp2genome.hits.gz | gawk '{n[$1]++;z[$1]=z[$1] "\n" $0;}END{for (k in n) if (n[k]>1)nn++ ; print nn;}' >> $toto
echo "How many sequences are not realigned" >> $toto
gunzip -c  tmp/NEWHITS_snp/all.edited_sequence.m.snp.fasta.gz ZZZZZ.gz tmp/NEWHITS_snp/m.snp2genome.hits.gz  | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^>/{z[substr($1,2)]=1;next;}{if(zz==1)z[$1]=0}END{for (k in z)if(z[k]==1)print k}' > tutu.$$
cat 

echo "Here is the list " >> $toto
gunzip -c tmp/NEWHITS_snp/m.snp2genome.hits.gz | gawk '{n[$1]++;z[$1]=z[$1] "\n" $0;}END{for (k in n) if (n[k]>1)print z[k];}' >> $toto

