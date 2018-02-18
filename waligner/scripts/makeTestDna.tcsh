#!bin/tcsh -f

#####################################################################################
# Author mieg@ncbi.nlm.nih.gov
# 2010_11_08
# Aim: construct a set of 50mers differing in a known way from a chosen 99bp reference genome
#      to verify that the aligners correctly report the mismatches
#####################################################################################

#####################################################################################
# Select as you like a "genomic sequence" of length 99
# so each 50-mer extracted from the "genome" covers the central base
set chromW=gtgctgttatgcgatattctcttatagcgcgtctgatctagttataatAActccctcgtatatttaaagagagggccctctctgatcgtatgatatagc


#####################################################################################
# You do not need to modify anything below this line
# when you run this script it generates 3 files
# chrom.wild.fasta : a fasta file with the 99bp reference genome you chose
# test.fasta : a fasta file with 2884 50bp sequences differing from the genome around position 50 
# test.csfasta : the same 2884 sequences, but now translated in SOLiD csfasta format

# You may then align the test file against the genome with you favorite aligner
# and verify how the mismatches are reported
#####################################################################################

# construct a fasta file containing this genome
echo $chromW | gawk '{s=$1;printf(">CHROM\n%s\n",s);}' > chrom.wild.fasta

# construct a set of test sequence altering the central base in all possible ways
echo $chromW | gawk '{s=$1;printf(">testF.wild\n%s\n",s);}' > chrom.modified.f.fasta
foreach aa (A T G C)
  echo $chromW | gawk '{if (aa == substr($1,50,1))next;s1=substr($1,1,49);s2=substr($1,51,49);printf(">testF.subst%s\n%s%s%s\n",aa,s1,aa,s2);}' aa=$aa >> chrom.modified.f.fasta
end
# add test sequences inserting 1 to 3 letters
foreach aa (A T G C AA AT AG AC TA TT TG TC GA GT GG GC CA CT CG CC AAA TTT GGG CCC) 
  echo $chromW | gawk '{s1=substr($1,1,49);s2=substr($1,50,50);printf(">testF.insert%s\n%s%s%s\n",aa,s1,aa,s2);}' aa=$aa >> chrom.modified.f.fasta
end
# add test sequences deleting 1 to 3 letters
  echo $chromW | gawk '{s1=substr($1,1,49);s2=substr($1,51,49);for (i =1;i<=3;i++)printf(">testF.del%d\n%s%s\n",i,s1,substr(s2,i));}' aa=$aa >> chrom.modified.f.fasta

# complement
cat  chrom.modified.f.fasta | bin/dna2dna -complement | sed -e 's/testF/testR/' >  chrom.modified.r.fasta

# transform into solid sequences
cat  chrom.modified.f.fasta | bin/dna2dna -O csfasta > chrom.modified.f.csfasta
cat  chrom.modified.r.fasta | bin/dna2dna -O csfasta > chrom.modified.r.csfasta

# extract all possible 50 mers
cat chrom.modified.f.fasta chrom.modified.r.fasta | gawk '/^>/{s=$1;next;}{n=length($1)-49;for(i=1;i<=n;i++)printf("%s.%d\n%s\n",s,i,substr($1,i,50));}' > test.fasta
cat chrom.modified.f.csfasta chrom.modified.r.csfasta | gawk '/^>/{s=$1;next;}{n=length($1)-49;for(i=1;i<=n;i++)printf("%s.%d\n%s\n",s,i,substr($1,i,50));}' > test.csfasta

#####################################################################################
###########################  end of file  ###########################################
#####################################################################################

## example
#clipalign -t chrom.wild.fasta -i test.csfasta -solidC -seedLength 18 -seedOffset 1 -seedShift 5 -minAli 12 -o toto -silent -splice -strategy RNA_seq
 $bin/clipalign -t chrom.wild.fasta -i test.fasta           -seedLength 18 -seedOffset 1 -seedShift 5 -minAli 12 -o toto -silent -splice
snp -count -run test -fasta chrom.wild.fasta -i toto.hits -o toto.BRS -minFrequency 1
snp -count2mrnas -i toto.BRS -run test -o toto.count
cat toto.count | gawk -F '\t' '/^#/{next}{printf("%s:%s%s\n",$1,$2,$3);}' > toto.snp_list
snp  -fasta chrom.wild.fasta -i toto.hits -aliExtend -snp_list toto.snp_list 
cat tutu.extend_snp 

 




