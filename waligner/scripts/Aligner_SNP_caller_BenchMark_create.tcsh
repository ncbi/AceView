#!bin/tcsh -f
# Authors Danielle and Jean Thierry-Mieg, NCBI
#     mieg@ncbi.nlm.nih.gov

set phase=$1
set TestLength=50
set TestStep=1

if ($3 != "") set TestLength=$3
if ($4 != "") set TestStep=$4


if ($phase == createRNATestSet) goto phasecreateRNATestSet
if ($phase == createMitoTestSet) goto phasecreateMitoTestSet
if ($phase == createGenomicTestSet) goto phasecreateGenomicTestSet

echo "Usage $0 [createRNATestSet | createMitoTestSet |  createGenomicTestSet] Nseq len step
echo "example  $0 createRNATestSet 7 50 3"
echo "      create 7M read of length 50 scanning every 3 bases a set of RefSeq genes" 
exit 1

###################################################################################
###################################################################################
## phase echo '    createTestSet: Make a test set of 10M 50mers, matching the RefSeq, with artificial SNPs at position 50,100,150,200,..."

phasecreateRNATestSet:
phasecreateMitoTestSet:
phasecreateGenomicTestSet:

set xMillion=240000
if ($2 != "") then
  set n1=`echo $2 | gawk '{printf("%d",0+$1);}'`
  if ($n1 < 1 || $n1 > 10) then
    echo "MAGIC createTestSet <n> expects n to be a number of millions between 1 and 10 and will create 6n million reads"
    echo "Invalid <n> parameter: MAGIC createTestSet $2"
    exit 1
  endif
  @ xMillion = 1000000 * $n1
endif

echo "createTestSet:  6 times $xMillion reads (forward, reverse, SOLiD) times (exact, SNP)"

if (! -d  TestData) mkdir  TestData
set myselect=""

if ($phase == createRNATestSet) then
  if (! -e  TARGET/Targets/$species.RefSeq.fasta.gz) then
    echo "Missing file TARGET/Targets/$species.RefSeq.fasta.gz"
    echo "I quit"
    goto phaseLoop
  endif
  set mytarget=TARGET/Targets/$species.RefSeq.fasta.gz
  gunzip -c $mytarget | gawk '/^>/{if(g){n1=nn[g]+0;if(n+0>n1+0){nn[g]=n;g2m[g]=m;}}m=substr($1,2);split(m,aa,"|");g=aa[3];n=0;next;}{n+=length($1);}END{for (g in g2m)print g2m[g];}' >  TestData/RefSeq.longest.list  
  set myselect=" -select TestData/RefSeq.longest.list "
  echo '    createRNATestSet: Make a 2 test sets of x Million 50mers, covering 20 times the RefSeq, one exact one with with artificial SNPs at position 100,200,...'
endif
if ($phase == createGenomicTestSet) then
  if (! -e  TARGET/Targets/$species.genome.fasta.gz) then
    echo "Missing file TARGET/Targets/$species.genome.fasta.gz"
    echo "I quit"
    goto phaseLoop
  endif
  if (1)
    # use CYP2C18andCYP2C19: 212kb
    set chrom=10
    set x1=96400510
    set x2=96650000
  endif
  if (! -e TestData/genome_section.fasta.gz) then
    bin/dna2dna -i TARGET/Targets/$species.genome.fasta.gz -get $chrom -leftClipAt $x1 -rightClipAt $x2 -I fasta -O fasta -gzo -o TestData/genome_section
  endif
  set mytarget=TestData/genome_section.fasta.gz
  echo '    createGenomicTestSet: Make a 2 test sets of x Million 50mers, covering 20 times a section of the genome, one exact one with with artificial SNPs at position 100,200,...'
endif
if ($phase == createMitoTestSet) then
  if (! -e  TARGET/Targets/$species.mito.fasta.gz) then
    echo "Missing file TARGET/Targets/$species.mito.fasta.gz"
    echo "I quit"
    goto phaseLoop
  endif
  if (! -e TestData/genome_section.fasta.gz) then
    gunzip -c  TARGET/Targets/$species.mito.fasta.gz  TARGET/Targets/$species.rrna.fasta.gz | bin/dna2dna -I fasta -O fasta -gzo -o TestData/genome_section
  endif
  set mytarget=TestData/genome_section.fasta.gz
  echo '    createMitoTestSet: Make a 2 test sets of x Million 50mers, covering 20 times the mito, one exact one with with artificial SNPs at position 100,200,...'
endif

# construct a fastc file with a set of stranded 50 mers
# construct the fasta forward reads
if (! -e TestData/test.exact.fasta.gz) then
# construct a fasta file ordered by length

   bin/dna2dna $myselect -i $mytarget -I fasta -makeTest -$xMillion -makeTestLength $TestLength -makeTestStep $TestStep -O fasta -gzo -o TestData/exact_ForwardStrand
   bin/dna2dna $myselect -i $mytarget -I fasta -makeTest  $xMillion -makeTestLength $TestLength -makeTestStep $TestStep  -O fasta -gzo -o TestData/SNP_ForwardStrand
endif
# construct the fasta reverse reads
   bin/dna2dna -i TestData/exact_ForwardStrand.fasta.gz -O fasta -gzo -complement -o TestData/exact_ReverseStrand
   bin/dna2dna -i TestData/SNP_ForwardStrand.fasta.gz  -O fasta -gzo -complement -o TestData/SNP_ReverseStrand

# construct the SOLiD reads
   bin/dna2dna -i TestData/exact_ForwardStrand.fasta.gz -O ccfa -gzo -o TestData/exact_SOLiD
   bin/dna2dna -i TestData/SNP_ForwardStrand.fasta.gz  -O ccfa -gzo  -o TestData/SNP_SOLiD

