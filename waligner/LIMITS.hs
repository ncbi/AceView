#!/bin/tcsh

# qusage 10 ; README.ali a1 ;  qusage 10 ; README.ali a2 ; qusage 10 ; README.ali a4 ;
setenv species hs
setenv seedLength 18
setenv minAli 24
setenv minEntropy 20
setenv minFastqQuality 0
setenv overhangLength 8
setenv bonusStrand 10

setenv project unified
setenv TMPDIR /export/home/TMP

# LIF: CGCCTTGGCCGTACAGCAGCCTCTTACAC se traduit par ttagacatatctccgtcgtagggatccc
setenv exitVector ATmTCGTATGCCGTCTTCTGCTTGAAAAAA,CGCCTTGGCCGTACAGCAGCCTCTTACAC
setenv exitVectorRaw ATmTCGTATGCCGTCTTCTGCTTGAAAAAA,ttagacatatctccgtcgtagggatccc

setenv chromSetAll "X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 Y"


# ATTENTION the pivot offset is the most likely to give a full alignment 
setenv pivot 11
setenv offSets "1 19 27 34 50 80 120 200 300 400 500 600 700 800 900 1000"

# Nov 2009, discard the cloud
setenv targets "rnaGene mito SpikeIn av RefSeq EBI HINV introns genome Line gtweek"
setenv targets "rnaGene mito SpikeIn av RefSeq EBI HINV genome.A genome.B genome.C genome.D genome.E genome.F"
# setenv targets "rnaGene mito SpikeIn av RefSeq EBI HINV genome.A genome.B genome.C"
# setenv targets "rnaGene mito SpikeIn "

# setenv targets " av RefSeq"
 setenv targets "chrom_20"
# setenv targets " EBI HINV "
# setenv targets "genome.A genome.B genome.C genome.D genome.E genome.F"
#setenv targets "genome.D genome.E genome.F "

setenv Etargets "av RefSeq EBI HINV introns genome"
setenv tissues "Brain UHR kidney"
setenv tissues "Brain UHR"
# setenv tissues "kidney"
setenv ntruetissues 0

setenv manips_454  "R454_O R454_A R454_GE R454_SW R454_Ti"
setenv manips_fastc_stranded "ILM_S ILM_RS S_brain HELdge"
setenv manips_fastc_non_stranded "$manips_454 HEL ILM_RnS ILM_nS KIDNEY"
# setenv manips_fastc_non_stranded "HEL ILM_RnS ILM_nS "

setenv manips_ccfa "LIF_S LIF_R"

# fastq with quality, for evaluation of the SNPs
setenv manips_fastq_stranded  "ILM_100"
setenv manips_fastq_non_stranded "ILM_35 ILM_nS100"

setenv genomeDb "/panfs/pan1.be-md.ncbi.nlm.nih.gov/aceview/zoo/human/NCBI_36a_R454"
setenv manips "$manips_fastc_non_stranded  $manips_fastc_stranded $manips_ccfa $manips_fastq_stranded $manips_fastq_non_stranded"
setenv manips "$manips_fastc_non_stranded  $manips_fastc_stranded $manips_ccfa"
#setenv manips "$manips_fastc_stranded  $manips_fastc_non_stranded"
setenv manips_not_amplified "$manips_fastc_stranded  $manips_fastc_non_stranded"
#setenv manips "$manips_fastq_stranded $manips_fastq_non_stranded"
#setenv manips "HELdge $manips_454"
#setenv manips "KIDNEY R454_Ti ILM_S LIF_R"
#setenv manips "R454_Ti ILM_S"

setenv manipsPolyA "$manips_fastc_non_stranded  $manips_fastc_stranded $manips_ccfa $manips_fastq_stranded $manips_fastq_non_stranded"
setenv manips_truly_stranded "ILM_S ILM_RS S_brain HELdge $manips_ccfa $manips_fastq_stranded"

setenv manips_danielle "$manips_454 ILM_nS ILM_S ILM_RS ILM_RnS S_brain HEL HELdge $manips_ccfa $manips_fastq_non_stranded $manips_fastq_stranded"
setenv manips_danielle "$manips_454 ILM_nS ILM_S ILM_RS ILM_RnS S_brain HEL HELdge $manips_ccfa"
setenv manips_danielle "KIDNEY"

echo toto | gawk '{for (i=1;i<=200;i++)printf("f.%d\n",i);}' > Fastc/lanes
#echo toto | gawk '{for (i=1;i<=3;i++)printf("f.%d\n",i);}' > Fastc/lanes

echo "manips   = $manips"
echo "tissues = $tissues"
echo "targets= $targets"
echo "chroms= $chromSetAll"


