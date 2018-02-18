#!bin/tcsh -f

set ff=$1
if (-e $ff.1.fasta.gz) exit 0

setenv LANG C

echo -n "... samtools: bam -> raw "
date
samtools view $ff -X -f 16 -F 256  | bin/dna2dna -I raw10 -complement -O raw -o $ff.r -gzo -keepNameBam 
samtools view $ff -X -F 272  | bin/dna2dna -I raw10 -O raw -o $ff.f -gzo  -keepNameBam 


echo -n "... alphabetical sort by sequence identifer: raw -> raw.sorted "
date
gunzip -c $ff.f.raw.gz | grep '/1' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.f.1.raw.sorted 
gunzip -c $ff.f.raw.gz | grep '/2' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.f.2.raw.sorted 
gunzip -c $ff.r.raw.gz | grep '/1' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.r.1.raw.sorted 
gunzip -c $ff.r.raw.gz | grep '/2' | gawk '{printf("%s\t%s\n",$2,$1);}' | sort > $ff.r.2.raw.sorted 

echo -n "... merge the .f and the .r reads [f  r].raw.sorted -> fr.raw.sorwted "
date
cat  $ff.f.1.raw.sorted $ff.r.1.raw.sorted | sort >  $ff.fr.1.raw.sorted
cat  $ff.f.2.raw.sorted $ff.r.2.raw.sorted | sort >  $ff.fr.2.raw.sorted

echo -n "... verif against the eventual orphans: raw.sorted -> verified.raw "
date
# it may happen that some reads are absent on one side
echo ZZZZZ > ZZZZZ
cat  $ff.fr.1.raw.sorted  $ff.fr.2.raw.sorted ZZZZZ $ff.fr.1.raw.sorted | gawk '/^ZZZZZ/{zz++;next;}{z=substr($1,1,length($1)-1);if(zz<1){nn[z]++;next;}if(nn[z]==2)print;}' >  $ff.fr.verified.1.raw
cat  $ff.fr.1.raw.sorted  $ff.fr.2.raw.sorted ZZZZZ $ff.fr.2.raw.sorted | gawk '/^ZZZZZ/{zz++;next;}{z=substr($1,1,length($1)-1);if(zz<1){nn[z]++;next;}if(nn[z]==2)print;}' >  $ff.fr.verified.2.raw

echo -n "... raw -> fasta.gz " 
date
cat  $ff.fr.verified.1.raw | gawk '{printf(">%s\n%s\n",$1,$2);}' | gzip > $ff.1.fasta.gz 
cat  $ff.fr.verified.2.raw | gawk '{printf(">%s\n%s\n",$1,$2);}' | gzip > $ff.2.fasta.gz 

\rm $ff.*raw*
ls -ls  $ff.*.fasta.gz 
echo -n "... done "
date


