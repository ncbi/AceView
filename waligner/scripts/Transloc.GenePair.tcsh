#!bin/tcsh -f

# 2013_09_21
# given a rarrangement detetcted in RNA_seq
# grab all hits to the 2 genes and all corresponding fasta files
# realign locally and analyse the proportions

set run=$1
set gene1=$2
set gene2=$3
set inversion=$3

set dd="tmp/TranslocGenePair/GenePair.$gene1.$gene2"

bin/dna2dna -I fasta -i TARGET/Targets/hs.av.fasta.gz -get $gene1.aAug10 -O fasta > $dd/$gene1.F.fasta
bin/dna2dna -I fasta -i TARGET/Targets/hs.av.fasta.gz -get $gene2.aAug10  -O fasta > $dd/$gene2.F.fasta

bin/dna2dna -I fasta -i  $dd/$gene1.F.fasta -O fasta -complement > $dd/$gene1.R.fasta
bin/dna2dna -I fasta -i  $dd/$gene2.F.fasta -O fasta -complement > $dd/$gene2.R.fasta

# export Nn nnn between the 2 genes so that position ND is always between the 2 genes
set N1=`cat $dd/$gene1.F.fasta | gawk '/^>/{next;}{n++;if(n==1)print length($1);}'`
set N2=`cat $dd/$gene2.F.fasta | gawk '/^>/{next;}{n++;if(n==1)print length($1);}'`
set Nn=`echo $N1 $N2 | gawk '{n=$1-$2;if(n<0)n=-n;print n+100;}'`
set ND=`echo $N1 $N2 $Nn | gawk '{n=$1+$2+$3;print int(n/2);}'`

# Gene1 -->       Gene2--->
echo ">GenePair.$gene1.$gene2" >  $dd/F1.F2.fasta
cat $dd/$gene1.F.fasta | tail -n +2 >>  $dd/F1.F2.fasta
echo toto | gawk '{for(i=0;i<Nn;i++)print "n" ;}' Nn=$Nn >>  $dd/F1.F2.fasta
cat $dd/$gene2.F.fasta | tail -n +2 >>  $dd/F1.F2.fasta

# Gene2 -->       Gene1--->
echo ">GenePair.$gene1.$gene2" >  $dd/F2.F1.fasta
cat $dd/$gene2.F.fasta | tail -n +2 >>  $dd/F2.F1.fasta
echo toto | gawk '{for(i=0;i<Nn;i++)print "n" ;}' Nn=$Nn>>  $dd/F2.F1.fasta
cat $dd/$gene1.F.fasta | tail -n +2 >>  $dd/F2.F1.fasta

# Gene1 -->       <----Gene2
echo ">GenePair.$gene1.$gene2" >  $dd/F1.R2.fasta
cat $dd/$gene1.F.fasta | tail -n +2 >>  $dd/F1.R2.fasta
echo toto | gawk '{for(i=0;i<Nn;i++)print "n" ;}' Nn=$Nn>>  $dd/F1.R2.fasta
cat $dd/$gene2.R.fasta | tail -n +2 >>  $dd/F1.R2.fasta

# Gene1 <--       --->Gene2
echo ">GenePair.$gene1.$gene2" >  $dd/R1.F2.fasta
cat $dd/$gene1.R.fasta | tail -n +2 >>    $dd/R1.F2.fasta 
echo toto | gawk '{for(i=0;i<Nn;i++)print "n" ;}' Nn=$Nn>>  $dd/R1.F2.fasta
cat $dd/$gene2.F.fasta | tail -n +2 >> $dd/R1.F2.fasta


# identification des clones utiles

if (! -d $dd/$run) mkdir $dd/$run
foreach lane (`cat Fastc/$run/LaneList`)
  gunzip -c tmp/COUNT/$lane.hits.gz | gawk -F '\t' '{ g=$9;if($8==tcl && g==g1 || g==g2)print ;}' tcl=ET_av g1=$gene1 g2=$gene2 >  $dd/$lane.hits
  cat   $dd/$lane.hits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);print c;}' | sort -u > $dd/$lane.list
  bin/dna2dna -i Fastc/$lane.fastc.gz -I fastc -O fastc -select $dd/$lane.list -o $dd/$lane
end

# grab the bridging clones 
cat $dd/$run/*.*.hits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);g=$9;if(g==g1)cc1[c]=1;if(g==g2)cc2[c]=1;}END{for(c in cc1)if(cc2[c]==1)print c}' g1=$gene1 g2=$gene2 | sort > $dd/$run/bridging.list

# group all fastc reads (they belong to the same run so the ids are distinct)
cat  $dd/$run/*.fastc > $dd/$run/fastc.all
 

foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  bin/bin/clipalign -i  $dd/$run/fastc.all -t  $dd/$fr.fasta -minAli 16 -splice -seedOffset 1 -seedShift 5 -clipPolyA -o $dd/$run/$fr
  cat $dd/$run/$fr.hits | bin/bin/wiggle -I BHIT -O BV > $dd/$run/$fr.hits.BV
end
# grab the bridging reads 
cat  $dd/$run/$fr.hits | gawk -F '\t' '{r=$1;c=substr(r,1,length(r)-1);c=r;a1=$12;if(a1<ND)cc1[c]=1;if(a1>ND)cc2[c]=1;}END{for(c in cc1)if(cc2[c]==1)print c}' ND=$ND >  $dd/$run/$fr.bridge.read.list
echo ZZZZZ > ZZZZZ
cat $dd/$run/$fr.bridge.list ZZZZZ  $dd/$run/$fr.hits |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{if(ok[$1]==1)print;}' >  $dd/$run/$fr.bridge.read.hits


\rm  $dd/$run/all.introns
foreach fr (F1.F2 F2.F1 F1.R2 R1.F2)
  cat $dd/$run/$fr.introns | sort -k 11n | tail -2 | gawk '{if(($5<ND && $7>ND)|| ($7<ND && $5>ND)){printf ("%s\t%s\t",run,fr);print;}}' fr=$fr ND=$ND run=$run >> $dd/$run/all.introns
end


goto phaseLoop

mkdir tmp/GeneLinks/$run/HITS2
bin/dna2dna  -i tmp/GeneLinks/$run/fastc.all -I fastc  -O fastc -select tmp/GeneLinks/$run/list.all > tmp/GeneLinks/$run/fastc.all2
bin/clipalign -i  tmp/GeneLinks/$run/fastc.all2 -t tmp/GeneLinks/$run/toto654.pair.fasta -minAli 16 -splice -seedOffset 1 -seedShift 5 -clipPolyA -o  tmp/GeneLinks/$run/HITS2 
 
# get the best donor acceptor
cat tmp/GeneLinks/Rhs1105/HITS2.overhangs | gawk '{n[$8 "\t" $9]++;}END{for (k in n)print k,n[k];}' | sort -k 3nr | head
set x1=`cat tmp/GeneLinks/Rhs1105/HITS2.overhangs | gawk '{n[$8 "\t" $9]++;}END{for (k in n)print k,n[k];}' | sort -k 3nr | head -1 | gawk '{print $2}'`
set x2=`cat tmp/GeneLinks/Rhs1105/HITS2.overhangs | gawk '{n[$8 "\t" $9]++;}END{for (k in n)print k,n[k];}' | sort -k 3nr | head -2 | tail -1 | gawk '{print $2}'`
# export the corresponding hits
cat tmp/GeneLinks/Rhs1105/HITS2.hits | gawk -F '\t' '{if($13==x1 || $13==x2)print}' x1=$x1 x2=$x2

cat << EOF > tmp/GeneLinks/Rhs1105/HITS2.list
  n.33468364#1
  n.36423699#1
EOF
bin/dna2dna -i tmp/GeneLinks/$run/fastc.all -I fastc  -O fastc -select tmp/GeneLinks/Rhs1105/HITS2.list  > tmp/GeneLinks/Rhs1105/HITS2.fastc


bin/clipalign -t  tmp/GeneLinks/Rhs1105/toto654.pair.fasta -i  tmp/GeneLinks/Rhs1105/HITS2.fastc -minAli 20

mkdir tmp/GeneLinks/database
mkdir tmp/GeneLinks/wspec
cp tmp/XHY/wspec/* tmp/GeneLinks/wspec
echo y | tbly  tmp/GeneLinks

tbly  tmp/GeneLinks << EOF
  parse tmp/GeneLinks/Rhs1105/HITS2.ace
  parse tmp/GeneLinks/Rhs1105/toto654.pair.fasta
  find sequence G*
  edit genomic
  acem
    cdna_1
    quit
  save
EOF


phaseLoop:
 echo -n "done "
 date
