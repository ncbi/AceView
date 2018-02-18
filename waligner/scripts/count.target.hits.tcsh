#!bin/tcsh -ef

set target=$1
set manip=$2
set tissue=$3
set lane=$4

set toto=MRNAS/$species.$target.transcript2gene.txt.gz
if (! -e $toto) set toto=ZZZZZ.gz

#ls -ls tmp/PHITS_$target/$manip/$tissue/$lane.off*.hits.gz
# we lower then rise again the score of the hits on unspliced mRNA, 
# de facto eliminating them if a hit of equal quality on a spliced mRNA exists

# by leaving a penalty, we only keep the tags which jump the introns
if ($target == av) then

  gunzip -c $toto ZZZZZ.gz tmp/PHITS_$target/$manip/$tissue/$lane.off*.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next}{if(zz==0){gsub(/\"/,"",$0);gene[$1]=$2;next;}}{uuBonus=0;}/G_|unspliced/{uuBonus=1;}/HIT/{probe=$1;score=0+$2;mrna=$6;a1=$8;a2=$9;da=a2-a1;snip=$16;ln=$21;ali=$22; if(score<16)next;score -=uuBonus;if(da<0)ali=-ali;gsub(/MRNA:/,"",mrna);g=mrna;if(gene[mrna])g=gene[mrna];printf("%s\t%06d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n",probe,score,g,mrna,ln,ali,a1,a2,snip);}'  | sort -T $TMPDIR -u -k 1,1 -k 2,2r -k 3,3  -k 4,4 -k 5,5 -k 6,6nr | gawk -F '\t' '{uuBonus=0;sBonus=0;}/G_|unspliced/{uuBonus=1;}{p=$1;score=$2;gene=$3;mrna=$4;ln=$5;ali=$6;a1=$7;a2=$8;snip=$9;if(p==oldp && score<oldscore)next;oldp=p;oldscore=score;score+=uuBonus;da=a2-a1;printf("%s\t%06d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n",p,score,gene,mrna,ln,ali,a1,a2,snip);}' | sort -T $TMPDIR -u -k 1,1 -k 2,2r -k 3,3  -k 4,4 -k 5,5 -k 6,6nr | gawk -F '\t' '{p=$1;score=$2;gene=$3;mrna=$4;ln=$5;ali=$6;a1=$7;a2=$8;snip=$9;da=a2-a1;if(p==oldp){if(score<oldscore)next;if(target != "genome" && gene==oldgene && mrna==oldmrna)next;}else{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1; for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));ng=0;ng2F=0;ng2R=0;oldgene=0;}ng++;if(gene!=oldgene){if(da>0)ng2F++;else ng2R++;}z[ng]= p "\t" score "\t" gene "\t" mrna "\t" ln "\t" ali "\t" a1 "\t" a2 "\t" snip ;oldp=p;oldscore=score;oldgene=gene;oldmrna=mrna;next;}END{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1;for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));}'  target=$target | gzip > tmp/COUNT/$manip/$tissue/$lane.$target.list.gz

else

if ($target == OBSOLETEgenome) then
ls -ls tmp/PHITS_chrom*/$manip/$tissue/$lane.off*.hits.gz
  gunzip -c $toto ZZZZZ.gz tmp/PHITS_chrom*/$manip/$tissue/$lane.off*.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next}{if(zz==0){gsub(/\"/,"",$0);gene[$1]=$2;next;}}/HIT/{probe=$1;score=0+$2;mrna=$6;a1=$8;a2=$9;da=$9-$8;snip=$16;ln=$21;ali=$22; if(score<16)next;if(da<0)ali=-ali;gsub(/MRNA:/,"",mrna);g=mrna;if(gene[mrna])g=gene[mrna];printf("%s\t%06d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n",probe,score,g,mrna,ln,ali,a1,a2,snip);}'  | sort -T $TMPDIR -u -k 1,1 -k 2,2r -k 3,3  -k 4,4 -k 5,5 -k 6,6nr | gawk -F '\t' '{p=$1;score=$2;gene=$3;mrna=$4;ln=$5;ali=$6;a1=$7;a2=$8;snip=$9;da=a2-a1;if(p==oldp){if(score<oldscore)next;if(target != "genome" && gene==oldgene && mrna==oldmrna)next;}else{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1;for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));ng=0;ng2F=0;ng2R=0;oldgene=0;}ng++;if(gene!=oldgene){if(da>0)ng2F++;else ng2R++;}z[ng]= p "\t" score "\t" gene "\t" mrna "\t" ln "\t" ali "\t" a1 "\t" a2 "\t" snip ;oldp=p;oldscore=score;oldgene=gene;oldmrna=mrna;next;}END{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1;for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));}' target=$target | gzip > tmp/COUNT/$manip/$tissue/$lane.$target.list.gz

else

  gunzip -c $toto ZZZZZ.gz tmp/PHITS_$target/$manip/$tissue/$lane.off*.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next}{if(zz==0){gsub(/\"/,"",$0);gene[$1]=$2;next;}}/HIT/{probe=$1;score=0+$2;mrna=$6;a1=$8;a2=$9;da=$9-$8;snip=$16;ln=$21;ali=$22; if(score<16)next;if(da<0)ali=-ali;gsub(/MRNA:/,"",mrna);g=mrna;if(gene[mrna])g=gene[mrna];printf("%s\t%06d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n",probe,score,g,mrna,ln,ali,a1,a2,snip);}'  | sort -T $TMPDIR -u -k 1,1 -k 2,2r -k 3,3  -k 4,4 -k 5,5 -k 6,6nr | gawk -F '\t' '{p=$1;score=$2;gene=$3;mrna=$4;ln=$5;ali=$6;a1=$7;a2=$8;snip=$9;da=a2-a1;if(p==oldp){if(score<oldscore)next;if(target != "genome" && gene==oldgene && mrna==oldmrna)next;}else{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1;for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));ng=0;ng2F=0;ng2R=0;oldgene=0;}ng++;if(gene!=oldgene){if(da>0)ng2F++;else ng2R++;}z[ng]= p "\t" score "\t" gene "\t" mrna "\t" ln "\t" ali "\t" a1 "\t" a2 "\t" snip ;oldp=p;oldscore=score;oldgene=gene;oldmrna=mrna;next;}END{sureStrand=1;if(ng2F > 0 && ng2R > 0)sureStrand = -1;for(i=1;i<=ng;i++)printf("%s\t%d\n",z[i],sureStrand*(ng2F+ng2R));}' target=$target | gzip > tmp/COUNT/$manip/$tissue/$lane.$target.list.gz

endif
endif

exit 0

