#!bin/tcsh -f

mkdir ~/aaa4/zoo/human/NB_local
cd NB_local

setenv MAGIC NB_Local
mkdir TARGET
mkdir TARGET/Targets

cp ~/NB/TARGET/LIMITS TARGET

~/ace/waligner/scripts/MAGIC int RNA

cat CancerGenesAbove80DE.txt ZZZZZ ~/NB/TARGET/GENES/av.metadata.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{g=$1;if(ok[g]==1){split($2,aa,":");split(aa[2],bb,"-");a1=bb[1];;a2=bb[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}d=500000;printf("%s\t%d\t%d\t%s\n",aa[1],a1 -d,a2+d,g);}}' | sort -k 1,1 -k 2,2n -k 3,3n > toto1
cat toto1 ZZZZZ | gawk -F '\t' '{c=$1;a1=$2;a2=$3;g=$4;if(c != old || a1 > b2 + 1000000){if(old)printf("%s\t%d\t%d\t%s\t%d\n",old,b1,b2,substr(oldg,3),b2-b1+1);b1=a1;oldg="";old=c;} oldg=oldg "_" g;b2=a2;}' > toto2

cat toto2 | gawk -F '\t' '{printf("chr_%s\t1\t%d\t%s\t%d\t%d\n",$4,$5-1,substr($1,4),$2,$3);}' > pseudo_chroms.shadows

bin/dna2dna -i ~/NB/TARGET/Targets/hs.genome.fasta.gz -shadow pseudo_chroms.shadows -O fasta -maxLineLn 60 -gzo > pseudo_chroms.fasta.gz

mv pseudo_chroms.fasta.gz TARGET/Targets/hs.genome.fasta.gz

cat pseudo_chroms.shadows | gawk 'BEGIN{printf("setenv chromSetAll \"");}{printf(" %s",$1);}END{printf("\"\n");}' >> LIMITS
cat << EOF >> LIMITS
setenv DNAtargets "genome"
setenv RNAtargets ""

EOF

source LIMITS
foreach chrom ($chromSetAll)
  bin/dna2dna -i TARGET/Targets/hs.genome.fasta.gz -I fasta -O fasta -maxLineLn 60 -get $chrom -gzo > TARGET/CHROMS/hs.chrom_$chrom.fasta.gz
end


ln -s ../NB/Fastc
tbly ~/NB/MetaDB << EOF
  find project NB
  follow run
  spush
  follow sample
  sor
  undo
  follow DEG
  sor
  spop
  show -a -f NB.runs.ace
  quit
EOF

echo "-R Project NB Local\n" >>  NB.runs.ace
tbly MetaDB << EOF
  pparse  NB.runs.ace
  save
  quit
EOF


MAGIC a0C wait a0D wait a123 wait c1 c2 wait c3 wait c4 wait c5 wait c6
qusage 5
MAGIC WIGGLE INTRON_DISCOVERY
qusage  5
MAGIC  f1 wait f2 wait f3 wait f4
qusage 5




