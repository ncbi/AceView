#!/bin/tcsh
# may 2023   analyse the tsmps and the BRS snps using snpsummary for the SEQC2_capture project

set phase=$1
set  zone=$2



if ($phase == 1) then 
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/$MAGIC.snp3a.$zone --snpType 0  -minSnpFrequency 5 -minSnpCover 20
  goto done
endif

if ($phase == 999) then 
  bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/BRS.$zone --snpType 0 -minSnpFrequency 5 -minSnpCover 20
  goto done
endif

if ($phase == dc) then 
 foreach stype (3)
  foreach pp (SnpA2_A  SnpI2_goodSeq_A  SnpI2_lowQ_A  SnpR2_A SnpA2_B  SnpI2_goodSeq_B  SnpI2_lowQ_B  SnpR2_B SnpA2_C  SnpI2_goodSeq_C  SnpI2_lowQ_C  SnpR2_C)
    bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.any2.$stype.$zone --snpType $stype -e VgGdD --doubleDetect -p $pp
  end
  set pp=SnpA2R2
  # bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$stype.$zone  --snpType $stype -e VISgrdGRD  -p $pp 
 end
 goto done
endif



if ($phase == p) then 
  set toto=toto
  echo "pparse Snp.projects_A.ace" > $toto
  echo "pparse MetaDB/$MAGIC/runs.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/groups.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/samples.ace" >> $toto
  echo "save" >> $toto
  echo "quit" >> $toto

  bin/tace tmp/TSNP_DB/$zone -no_prompt < $toto
  goto done
endif

if ($phase == bed) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF
    query find variant VCF
    select -o tmp/TSNP_DB/$zone/vcf.hg37.txt v,m,a1,b,bb from v in @, m in v->VCF, a1 in m[1], b in m[2], bb in m[3] 
EOF
    cat tmp/TSNP_DB/$zone/vcf.hg37.txt  | gawk -F '\t' '{printf("%s\t%d\t%d\n", $2, $3, $3+1);}' | sort -u >  tmp/TSNP_DB/$zone/vcf.hg37.bed
  goto done
endif



done:
  echo done
  exit 0



  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    scripts/submit tmp/TSNP_DB/$zone "scripts/snp3.tcsh bed $zone"
  end


qusage 1

cat tmp/TSNP_DB/zoner.*/Snp*_A.3.double_counts.txt | gawk '/libs/{next;}/^RNA/{nn[$1 "\t" $2]+=$3;}END{for(k in nn)printf("%s\t%d\n",k,nn[k]);}' | sort

set toto=RESULTS/SNV/SnpA2R2.DanLi.SNP_summary.txt 
 cat tmp/TSNP_DB/zoner.*/SnpA2R2.3.zoner.*.SNP_summary.txt  | head -12 | gawk '/^#/{print}' > $toto
 cat tmp/TSNP_DB/zoner.*/SnpA2R2.3.zoner.*.SNP_summary.txt  | gawk '/^#/{next;}{print}' | sort -k 2,3 -k 3,3n  >> $toto
# histo of our frequencies
 cat $toto | gawk -F '\t' '/^#/{next;}{n=int($21/10);print n;}' | tags | sort -k 1n

# histo of dan li frequencies
cat DanLi/DanLi.zoner.*.ace | gawk '/^DanLi_counts AGLR2/{x=int($7/10);print x}' | tags | sort -k 1n


# compare the coverage of magic/danli
cat tmp/TSNP_DB/zoner.*/SnpA2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($23,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n
cat tmp/TSNP_DB/zoner.*/SnpR2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($25,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n
