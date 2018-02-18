#!bin/tcsh -ef

set chrom=$1
set genomeDb=$2
set ici=`pwd`
$tacembly ../ZH$chrom <<EOF
  parse -NOXREF tmp/PHITS_introns/Needed/introns.needed.$chrom.ace
  query  subsequence
  table -active -o $ici/tmp/PHITS_introns/Needed/introns.needed.$chrom.dna.txt -f $ici/bin/chrom2subseq2dna.def
  table -o $ici/tmp/PHITS_introns/Needed/introns.needed.$chrom.mrna2dna.txt -f $ici/bin/mrnaIntron2dna.def
  quit
EOF
date

cat tmp/PHITS_introns/Needed/introns.needed.$chrom.mrna2dna.txt  | gawk -F '\t' '/^\"/{gsub(/\"/,"",$0);gsub(/CHROMOSOME_/,"",$0);d1=$13;d2=$14;a1=$5;a2=$6;x1=$7;x2=$8;if(a1<a2){u1=a1+x1-1;u2=a1+x2-1;}else{u1=a1-x1+1;u2=a1-x2+1;}printf(">%s__%d_%d\t",$4,u1,u2);for(i=1;i<100 - length(d1);i++)printf("n");printf("%s%s\n",d1,d2);}' | sort -u >  tmp/PHITS_introns/Needed/introns.needed.mrna2dna.$chrom.prefasta
cat tmp/PHITS_introns/Needed/introns.needed.$chrom.mrna2dna.txt  | gawk -F '\t' '/\"/{gsub(/\"/,"",$0);gsub(/CHROMOSOME_/,"",$0);d1=$13;d2=$14;a1=$5;a2=$6;x1=$7;x2=$8;if(a1<a2){u1=a1+x1-1;u2=a1+x2-1;}else{u1=a1-x1+1;u2=a1-x2+1;}foot=$12;if (foot=="gt_ag" || foot=="gc_ag"|| foot=="at_ac")other="";else other="Other "; printf("Element %s__%d_%d\n%s %s\nIn_mRNA \"%s\"\nFrom_gene \"%s\"\n\n",$4,u1,u2,other,foot, $1, $3);}' >  tmp/PHITS_introns/Needed/introns.needed.mrna2dna.$chrom.ace


cat  tmp/PHITS_introns/Needed/introns.needed.$chrom.dna.txt | gawk -F '\t' '{gsub(/\"/,"",$0);gsub(/CHROMOSOME_/,"",$0);if($2~/XI_/)print;}'  | gawk -F '\t' '{s=$2;gsub(/XI_/,"",s);i=length(s);d=$5;j=length(d);z=substr(s,i-1);if(z=="_A"){sA=substr(s,1,i-2);dA=substr(d,1,j-2);iA=substr(d,j-1,2);next;}if(z=="_B"){sB=substr(s,1,i-2);dB=substr(d,3);iB=substr(d,1,2);iAB = iA "_" iB ; if(iAB=="gt_ag" || iAB == "gc_ag" || iAB == "ct_ac"|| iAB == "at_ac") other = "" ; else other="Other " ; if(sA==sB){printf("Sequence %s\ncDNA_clone %s\nForward\nIs_read\nComposite\nColour GRAY\n%s %s\n\n",sB,sB,other,iAB); printf(">%s\n%s%s\n\n",sB,dA,dB);}}}' >  tmp/PHITS_introns/Needed/introns.needed.$chrom.out.ace

cat  tmp/PHITS_introns/Needed/introns.needed.$chrom.dna.txt | gawk -F '\t' '{gsub(/\"/,"",$0);gsub(/CHROMOSOME_/,"",$0);if($2~/XI_/)print;}'  | gawk -F '\t' '{s=$2;gsub(/XI_/,"",s);i=length(s);d=$5;j=length(d);z=substr(s,i-1);if(z=="_A"){sA=substr(s,1,i-2);dA=substr(d,1,j-2);iA=substr(d,j-1,2);next;}if(z=="_B"){sB=substr(s,1,i-2);dB=substr(d,3);iB=substr(d,1,2);iAB = iA "_" iB ; if(iAB=="gt_ag" || iAB == "gc_ag" || iAB == "ct_ac"|| iAB == "at_ac") other = "" ; else other="Other " ; if(sA==sB){ printf(">%s\t%s%s\n",sB,dA,dB);}}}' >  tmp/PHITS_introns/Needed/introns.needed.$chrom.prefasta

cat  tmp/PHITS_introns/Needed/introns.needed.mrna2dna.$chrom.prefasta  tmp/PHITS_introns/Needed/introns.needed.$chrom.prefasta | sort -u | gawk -F '\t' '{printf("%s\n%s\n",$1,$2);}' >  tmp/PHITS_introns/Needed/introns.needed.$chrom.fasta

exit 0

# split code works better for chrom 1

foreach ii (00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15)
$tacembly  ~/aa/zoo/human/NCBI_36a_R454/ZH1 <<EOF
  parse -NOXREF tmp/PHITS_introns/Needed/introns.needed.1.ace_split$ii
  query  subsequence
  table -active -o $ici/tmp/PHITS_introns/Needed/introns.needed.$chrom.dna.$ii.txt -f $ici/bin/chrom2subseq2dna.def
  quit
EOF
end
