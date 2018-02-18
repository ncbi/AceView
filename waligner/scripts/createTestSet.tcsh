#!bin/tcsh -f

set test=$1
set N=$2


if ($N == "") then
  createTestSet.tcsh requires 2 paramaters, a fasta file name and an out put format fasta | fastq
endif

if (! -e $test.fasta) then
  echo "missing file $test.fasta"
  exit 1
endif

# create the variant forward fasta file
bin/dna2dna -createTestSet $N -i $test.fasta -ctsL 50 -ctsStep 1  -gzo -o $test

# reexport as csfasta
bin/dna2dna -i $test.variant.forward.fasta.gz  -I fasta -O csfasta -gzo -o $test.variant.forward
bin/dna2dna -i $test.variant.reverse.fasta.gz  -I fasta -O csfasta -gzo -o $test.variant.reverse
bin/dna2dna -i $test.exact.forward.fasta.gz    -I fasta -O csfasta -gzo -o $test.exact.forward
bin/dna2dna -i $test.exact.reverse.fasta.gz    -I fasta -O csfasta -gzo -o $test.exact.reverse

