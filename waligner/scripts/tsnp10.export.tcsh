#!/bin/tcsh -f
# export the snps/tsnps using snpsummary

foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  if (-e tmp/TSNP_DB/$zone/mainClone.ace) continue
  cat <<EOF > tmp/TSNP_DB/$zone/mainClone.ace
Clone R
Main_clone
MainTitle  "SNPs $zone"
Species $species

EOF
  bin/tace tmp/TSNP_DB/$zone <<EOF
    read-models
y
    pparse tmp/TSNP_DB/$zone/mainClone.ace
    save
    quit
EOF
end



foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/DanLi.$zone -snpType 3 
end

foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/BRS.$zone -snpType 0 -minSnpFrequency 20 -minSnpCover 20
end

