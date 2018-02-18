#!bin/tcsh -f

set run=$1
set zone=$2
set ff=$3

set fasta=tmp/SNP_ZONE/$zone.fasta.gz
set qual=""

if (-e tmp/SNP_ZONERUN/$run/$zone.fasta.gz)  set fasta=tmp/SNP_ZONERUN/$run/$zone.fasta.gz
if (0  && -e tmp/SNP_BRS/$run/Quality_profile.txt)   set qual=" --runQuality tmp/SNP_BRS/$run/Quality_profile.txt"


echo "gunzip -c  cat tmp/SNP/$run/s2c.hits_file.list   | bin/snp -aliExtend -unique -snp_list $ff  -fasta $fasta $qual  -gzo -o tmp/SNP/$run/$zone -minCover $minSnpCover -run $run -minAliPerCent 90"

      gunzip -c `cat tmp/SNP/$run/s2c.hits_file.list`  | bin/snp -aliExtend -unique -snp_list $ff  -fasta $fasta $qual  -gzo -o tmp/SNP/$run/$zone -minCover $minSnpCover -run $run -minAliPerCent 90

