#!/bin/tcsh -f
# export the snps/tsnps using snpsummary

foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/DanLi.$zone -snpType 3 
end

foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/BRS.$zone -snpType 0 -minSnpFrequency 20 -minSnpCover 20
end

