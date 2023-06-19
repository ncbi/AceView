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

