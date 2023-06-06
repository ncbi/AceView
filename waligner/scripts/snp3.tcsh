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

if ($phase == bed) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF
    query find variant VCF
    select -o tmp/TSNP_DB/$zone/vcf.hg37.txt v,m,a1,b,bb from v in @, m in v->VCF, a1 in m[1], b in m[2], bb in m[3] 
EOF
    cat tmp/TSNP_DB/$zone/vcf.hg37.txt  | gawk -F '\t' '{printf("%s\t%d\t%d\n", $2, $3, $3+1);}' | sort -u >  tmp/TSNP_DB/$zone/vcf.hg37.bed
  goto done
endif



if ($phase == p) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    read-models
y
    pparse Snp.projects_A.ace
    pparse MetaDB/$MAGIC/runs.ace
    pparse MetaDB/$MAGIC/groups.ace
    pparse MetaDB/$MAGIC/samples.ace
    save
    quit
EOF
  goto done
endif


if ($phase == R) then 
    set remap2g=remap2genome
# if -o filename is not provided, tsnp directly edits the database
    bin/tsnp --db_$remap2g  tmp/METADATA/mrnaRemap.gz  --db tmp/TSNP_DB/$zone --force 
    bin/tsnp --db_translate --db tmp/TSNP_DB/$zone  -p $MAGIC

  goto done
endif

if ($phase == E) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    s -o Global_SNP_LIST/brs.$zone.list v from v in ?Variant where v#VCF 
    quit
EOF
  goto done
endif

if ($phase == W) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    query find run PacBio.ROC
R3.*
    kill
    query find run PacBio.*flnc*
    kill
    query find variant BRS_counts
    // select -o tmp/TSNP_DB/$zone/brs.fix v,r from v in @, r in v->brs_counts where ! r#project
    pparse tmp/TSNP_DB/$zone/brs.fix.ace
    save
    quit
EOF
  goto done
endif

if ($phase == w) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    read-models
y
    pparse DanLi/DanLi.hg37.$zone.ace

    query find variant vcf_hg38 
    spush
    pparse DanLi/brs.vcf38.$zone.ace
    sminus
    spop
    edit -D vcf_hg38

    query find variant  WTrue
    spush
    pparse DanLi/wendell.true_positives.$zone.ace
    sminus
    spop
    edit -D Wtrue

    query find variant  WTrue2
    spush
    pparse DanLi/wendell.true_positives.$zone.ace2
    sminus
    spop
    edit -D Wtrue2

    query find variant  Wfalse
    spush
    key  DanLi/brs.false_positive.$zone.list
    query ! Wfalse
    edit Wfalse
    key  DanLi/brs.false_positive.$zone.list
    sminus
    spop
    edit -D Wfalse

    query find variant  Wfalse2
    spush
    key  DanLi/brs.false_positive.$zone.list2
    query ! Wfalse2
    edit Wfalse2
    key  DanLi/brs.false_positive.$zone.list2
    sminus
    spop
    edit -D Wfalse2

    key DanLi/Fatigue.monomodal.$zone.list
    edit monomodal Fatigue
    save
    quit
EOF
  goto done
endif

# SnpA2_A  SnpI2_goodSeq_A  SnpI2_lowQ_A  SnpR2_A SnpA2_B  SnpI2_goodSeq_B  SnpI2_lowQ_B  SnpR2_B SnpA2_C  SnpI2_goodSeq_C  SnpI2_lowQ_C  SnpR2_C
if ($phase == dc) then 
 foreach stype (1 5 6 7 8)
  foreach pp (SnpA2_A  SnpR2_A SnpA2_B  SnpR2_B)
    bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.any22.$stype.$zone --snpType $stype -e VgGdD --doubleDetect -p $pp
  end
  set pp=SnpA2R2
  set stype=3
  # bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$stype.$zone  --snpType $stype -e VISgrdGRD  -p $pp 
 end
 goto done
endif



done:
  echo done
  exit 0



  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    scripts/submit tmp/TSNP_DB/$zone "scripts/snp3.tcsh dc  $zone"
  end


qusage 1

cat tmp/TSNP_DB/zoner.*/SnpA2_A.any22.5.zone*double_counts.txt | gawk '/libs/{next;}/^RNA/{nn[$1 "\t" $2]+=$3;}END{for(k in nn)printf("%s\t%d\n",k,nn[k]);}' | sort

set toto=RESULTS/SNV/SnpA2R2.DanLi.SNP_summary.txt 
 echo -n "## Whole geneome : " > $toto
 date >> $toto 
 cat tmp/TSNP_DB/zoner.*/SnpA2R2.3.zoner.*.SNP_summary.txt  | head -12 | tail -11 | gawk '/^#/{print}'  >> $toto
 cat tmp/TSNP_DB/zoner.*/SnpA2R2.3.zoner.*.SNP_summary.txt  | gawk '/^#/{next;}{print}' | tab_sort -k 3,3n -k 4,4n  >> $toto
# histo of our frequencies
 cat $toto | gawk -F '\t' '/^#/{next;}{n=int($21/10);print n;}' | tags | sort -k 1n

# histo of dan li frequencies
cat DanLi/DanLi.zoner.*.ace | gawk '/^DanLi_counts AGLR2/{x=int($7/10);print x}' | tags | sort -k 1n


# compare the coverage of magic/danli
cat tmp/TSNP_DB/zoner.*/SnpA2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($23,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n
cat tmp/TSNP_DB/zoner.*/SnpR2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($25,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n


  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    cat tmp/TSNP_DB/$zone/brs.fix | gawk '{printf ("Variant %s\n-D BRS_counts %s\n\n", $1,$2);}' > tmp/TSNP_DB/$zone/brs.fix.ace
  end
