#!bin/tcsh -f

set run=$1
set zone=$2

echo "  " > tmp/SNP2/$run/_s13
foreach lane (`cat Fastc/$run/LaneList`)
  echo "  tmp/COUNT/$lane.hits.gz " >> tmp/SNP2/$run/_s13
end

if (-e  tmp/SNP_ZONE/$zone.fasta.gz) then
  set ff=tmp/SNP_ZONE/$zone.fasta.gz
else
  set ff=TARGET/CHROMS/$species.chrom_$zone.fasta.gz
endif

echo "gunzip -c  (cat tmp/SNP2/$run/_s13) | bin/snp -aliExtend -snp_list tmp/SNP/$run/$zone.count.u.snp.gz  -fasta $ff -gzo -o tmp/SNP/$run/$zone"
gunzip -c `cat tmp/SNP2/$run/_s13`  | bin/snp -aliExtend -snp_list  tmp/SNP/$run/$zone.count.u.txt.gz  -fasta $ff -gzo -o tmp/SNP/$run/$zone

echo "s13.done"


exit 0

set zone=mito
set run=NA12878MOD


cat RESULTS/SNV/FDA.snp_list_with_allele_frequency_per_sample.txt | grep 599969 | gawk -F '\t' '{printf("%s\t%s\t%s\n",$2,$3,$5);}' > tyty.list

zcat tmp/COUNT/NA12878MOD/f2.3*.hits.gz | gawk -F '\t' '{if($11== 3 && $13 > 59996600 && $12 < 59997200)print}' > tyty.hits
set ff=TARGET/CHROMS/$species.chrom_3.fasta.gz

cat tyty.hits |  ~/ace/bin.LINUX_4/snp -phasing -snp_list tyty.list -fasta $ff -run $run  -o tyty

