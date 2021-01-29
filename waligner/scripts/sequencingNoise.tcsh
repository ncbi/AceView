#!/bin/tcsh

## 2020_06_06

## The code variant_caller allows to search N words inside all reads of an Illuminarun or in any other large sequence set
################################################################################################################################################
## 1: Apply this method to measure the substitution sequencing noise
#    1a: create a set of words extracted from the reference geome (or highly expresed transcripts) and add a central substitution xxAxx -> xxGxx
#    1b: count these exact words in the run. Naive expectation: xxAxx = local coverage, xxGxx = zero, but really xxGxx = sequencing noise
#    1c: collect the counts and sort it by kmers, this gives the noise for every xxNxx, independently of the mapping

######
# 1a:

date
set target=TARGET/Targets/$species.genome.fasta.gz
set N=1000  # number of words
if (! -d tmp/SEQNOISE) mkdir  tmp/SEQNOISE
foreach s (`seq 0 30`)
  dna2dna -i $target -createTestSet $N -ctsL 31 -ctsShift $s -ctsStep 31 -O raw -o tmp/SEQNOISE/shift.$s.motifs
end
cat  tmp/SEQNOISE/shift.*.raw >  tmp/SEQNOISE/all.motifs
wc  tmp/SEQNOISE/all.motifs
foreach run (`cat MetaDB/$MAGIC/RunsList`)
  if (! -d tmp/SEQNOISE/$run) mkdir  tmp/SEQNOISE/$run
  if (! -e  tmp/SEQNOISE/$run/sqn.val.txt) then
    scripts/submit  tmp/SEQNOISE/$run/sqn  "bin/variant_caller -count -wLn 31 -wordFile  tmp/SEQNOISE/all.motifs -run $run -o   tmp/SEQNOISE/$run/sqn -antistrand"
  endif
end
qusage 1

cat  tmp/SEQNOISE/*/sqn.val.txt | gawk -F '\t' '{split($2,aa,".");v=aa[1];mut=0;if(aa[2]=="sub")mut=1;vv[v]=1;nn[v,mut]+=$3;seq[v,mut]=substr($7,15,3);}END{for(v in vv){m=nn[v,1]+0;w=nn[v,0]+0;if(m==0)w++;k=int(1000.0*m/(m+w));z=seq[v,0]">"seq[v,1];zw[z]+=w;zm[z]+=m;if(k>100)k=100;kk[z,k]+=m;;}for(z in zw){printf("%s\t%d\t%d\t%f",z,zw[z],zm[z],100*(zm[z]/(.0001+zm[z]+zw[z])));for(k=0;k<=100;k++)printf("\t%d",kk[z,k]);printf("\n");}}' | sort > titi1

set toto=RESULTS/Mismatches/$MAGIC.SequencingNoise.txt
echo -n "### $toto : " > $toto
date >> $toto
echo "## First we construct large set of 31-mers occuring in the sequenced target, and add a random susbtitution of the central base." >> $toto
echo "## Then we count the number of occurences of these reference and mutated 31-mers in all sequence files, before alignemnt." >> $toto
echo "## Some of the variants really occur as SNPs is some runs where they are seen at a high allele frequency, possibly 100%." >> $toto
echo "## Finally, we merge all occurences of the central triplet variations and cumul the counts, isolating the real SNPs which also contibute in a high frequency column." >> $toto 
echo "## For each triplet varication (192 = 3 * 64 cases) the table reports the total counts, and the counts in the first 100 bins per-thousand. The last bin cumulates all variant counts with allele frequeency over 10%." >> $toto 
echo toto | gawk '{printf("# Mismatch\tReference\tVariant\t%% Variant : 100*v/(v+r)\tVariant in per 1000 range ");for(i=1;i<100;i++)printf("\t%d-%d",i-1,i);printf("\t>=10 %%\n");}' >> $toto
cat titi1 >> $toto

date
echo done

#### 
