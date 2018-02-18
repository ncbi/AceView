#!bin/tcsh -f

set Strategy=$1
set zone=$2
set snp_type=".homozygous"
set snp_type=".filtered"
set snp_type=".differential"
set snp_type=""


if (! -e tmp/SNP_DB/$zone/snp.parseChrom.done) then

  if (! -d tmp/SNP_DB/$zone) mkdir tmp/SNP_DB/$zone
  if (! -d tmp/SNP_DB/$zone/database) then
    pushd tmp/SNP_DB/$zone
        mkdir database
        ln -s ../../../MetaDB/wspec wspec
        ln -s ../../../metaData/tables.VariantDB tables
        whoami >> wspec/passwd.wrm
        echo y | ../../../bin/tacembly .
    popd
    if (-e TARGET/CHROMS/Centromere_telomeres.ace) then
      echo "pparse TARGET/CHROMS/Centromere_telomeres.ace" | bin/tacembly tmp/SNP_DB/$zone -no_prompt
    endif
  endif
  echo "read-models" > tmp/SNP_DB/$zone/_rs3.$zone

  foreach target (mito SpikeIn $Etargets)
    if (-e  tmp/METADATA/$target.MRNA.info.ace) then
      # gene gene_type geneId mrna mrna-length GC
      echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/SNP_DB/$zone/_rs3.$zone
      echo "pparse  tmp/METADATA/$target.MRNA.info.ace" >> tmp/SNP_DB/$zone/_rs3.$zone
    endif 
    if (-e tmp/METADATA/gtf.$target.goodProduct.ace) then
      # mrna -> product -> good/best product 
      echo "pparse   tmp/METADATA/gtf.$target.goodProduct.ace"  >> tmp/SNP_DB/$zone/_rs3.$zone
    endif
    if (-e  tmp/SNP_DB/$species.$target.mrna.fasta.gz) then
      echo "pparse tmp/SNP_DB/$species.$target.mrna.fasta.gz"  >> tmp/SNP_DB/$zone/_rs3.$zone
    endif 
  end

  # gene mrna product mrna-splicing intmap
  echo "parse TARGET/MRNAS/mrnaStructure.ace" >> tmp/SNP_DB/$zone/_rs3.$zone
  echo "pparse tmp/SNP_DB/RefSeq.cds.ace.gz" >> tmp/SNP_DB/$zone/_rs3.$zone
  # product quality should be in the gff file
  echo "pparse  TARGET/MRNAS/very_good_product.ace" >> tmp/SNP_DB/$zone/_rs3.$zone
  echo "pparse  TARGET/MRNAS/good_product.ace" >> tmp/SNP_DB/$zone/_rs3.$zone
  # mRNA-dna extracted from TARGET/Targets
  echo "pparse tmp/SNP_DB/$species.av.mrna.fasta.gz" >> tmp/SNP_DB/$zone/_rs3.$zone
  if (-e  tmp/SNP_ZONE/$zone.fasta.gz && $Strategy != RNA_seq)   echo "pparse tmp/SNP_ZONE/$zone.fasta.gz"  >> tmp/SNP_DB/$zone/_rs3.$zone
  # Gene titles and intmap
  echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/SNP_DB/$zone/_rs3.$zone

  echo save >> tmp/SNP_DB/$zone/_rs3.$zone
  echo quit >> tmp/SNP_DB/$zone/_rs3.$zone

  echo -n "snp.parseChrom "
  date
  bin/tacembly tmp/SNP_DB/$zone < tmp/SNP_DB/$zone/_rs3.$zone
  touch  tmp/SNP_DB/$zone/snp.parseChrom.done

endif

if (! -e   tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace) then
  echo -n "snp.prepare file  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace "
  date

   # split the file because gawk chokes on very large files and acedb can reconstitute the pieces
    cat  tmp/SNPH/$zone/$MAGIC.snp$snp_type.sorted   | sort -u > tmp/SNPH/toto.$$.s
    split -l 20000  -a 6 tmp/SNPH/toto.$$.s tmp/SNPH/toto.$$.s.p.
    echo ' ' >  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace

    set targetType=IntMap
    if ($Strategy == RNA_seq) set targetType=mRNA
    foreach ff (`ls  tmp/SNPH/toto.$$.s.p.*`)
      echo $ff
      cat $ff | gawk -F '\t' -f scripts/snp2ace.awk format=txt targetType=$targetType >>   tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace
    end

    \rm  tmp/SNPH/toto.$$  tmp/SNPH/toto.$$.*  
endif

# update the Runs
bin/tacembly tmp/SNP_DB/$zone <<EOF
       read-models
       query find run project
       edit -D project
       query find run group
       edit -D group 
       pparse  MetaDB/$MAGIC/runs.ace
       pparse  MetaDB/$MAGIC/samples.ace
       save
       quit
EOF

if (! -e tmp/SNP_DB/$zone/snp.parseSnp.done) then
  if (-e   tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace) then
   echo -n "snp.parse file  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace "
   date

    bin/tacembly tmp/SNP_DB/$zone <<EOF
       pparse  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace
       save
       quit
EOF

    touch  tmp/SNP_DB/$zone/snp.parseSnp.done
  endif
endif

if (0) then
    if (-e  tmp/SNP_DB/$zone/snp.remap.done) \rm tmp/SNP_DB/$zone/snp.remap.done
    if (-e  tmp/SNP_DB/$zone/snp.translate.done) \rm tmp/SNP_DB/$zone/snp.translate.done
    bin/tacembly tmp/SNP_DB/$zone <<EOF
       read-models
       query find  variant observed_sequence
       edit -D  observed_sequence
       query find  variant mrna
       edit -D  mrna
       query find  variant genebox
       edit -D genebox
       query find  variant coding
       edit -D coding
       query find  variant
       show -a -f tmp/SNP_DB/$zone/dumpvar.txt
       parse  tmp/SNP_DB/$zone/dumpvar.txt
       save
       quit
EOF
endif

if (! -e  tmp/SNP_DB/$zone/snp.remap.done) then
  if ($zone == mito) then
    # bin/snp -db_remap2genes tmp/METADATA/mrnaRemap.mito.txt.gz  -db tmp/SNP_DB/$zone
  endif
  set remap2g=remap2genes
  if ($Strategy == RNA_seq) set remap2g=remap2genome

  echo -n "snp.Remap $zone $remap2g against  tmp/METADATA/mrnaRemap.gz  "
  date
  bin/snp -db_$remap2g  tmp/METADATA/mrnaRemap.gz  -db tmp/SNP_DB/$zone

  touch tmp/SNP_DB/$zone/snp.remap.done 
endif

if (! -e  tmp/SNP_DB/$zone/snp.translate.done) then
  echo -n "snp.translate  $zone "
  date
  bin/tacembly tmp/SNP_DB/$zone <<EOF
    query find  variant Translation
    edit -D Translation
    query find  variant IntMap
    show -a -f  tmp/SNP_DB/$zone/toto.intmap.ace IntMap
    edit -D IntMap
    read-models
    parse  tmp/SNP_DB/$zone/toto.intmap.ace
    query find  variant observed_sequence
    edit -D  observed_sequence
    parse  tmp/SNP_DB/RefSeq.fasta  // used in E.coli
    parse tmp/SNP_DB/RefSeq.ace // used in E.coli
    save
    quit
EOF

  bin/snp -db_translate -db tmp/SNP_DB/$zone
  touch tmp/SNP_DB/$zone/snp.translate.done
endif

if (! -e  tmp/SNP_DB/$zone/snp.differential.done2 && -e tmp/SNPH/$zone/$MAGIC.differential_snps ) then
  cat  tmp/SNPH/$zone/$MAGIC.differential_snps | gawk -F '\t' '/^#/{next}{split($1,aa,":"); f1=$5;f2=$12;if (f1<f2){if(f1<1)typ="Gained_in"; else typ="Up_in";}else{if(f2<1)typ="Lost_in";else typ="Down_in";}printf ("Variant \"%s:%s:%s\"\nDifferential %s \"%s\" \"%s\" %s %s\n\n",aa[1],aa[2],aa[3], typ, $9, $2, f1, f2);}' | grep -v "D_HR_surv_noMYCNA_NB399_exome_Normal" | grep -v Ghs3987 | grep -v Ghs4007 >  tmp/SNP_DB/$zone/snp.differential.ace
  bin/tacembly tmp/SNP_DB/$zone <<EOF
   query find variant Differential
   edit -D Differential
   read-models
   pparse  tmp/SNP_DB/$zone/snp.differential.ace
   query find variant differential ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.list
   query find variant Gained_in ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Gained.list
   query find variant Up_in  ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Up.list
   query find variant Down_in ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Down.list
   query find variant Lost_in  ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Lost.list
  save
   quit
EOF
  touch  tmp/SNP_DB/$zone/snp.differential.done
endif

laba:
if (-e RESULTS/StLouisSNP/zone.stl.ace) then
  echo "pparse  RESULTS/StLouisSNP/zone.stl.ace"  | bin/tacembly tmp/SNP_DB/$zone -no_prompt 
endif
echo "bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome  $Reference_genome"
      bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome  "$Reference_genome" 
      bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -differential -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC.differential -Reference_genome  "$Reference_genome "


echo 'bin/snp -i tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt              -db tmp/SNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome   "$Reference_genome "' 
      bin/snp -i tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt              -db tmp/SNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome   "$Reference_genome " 
      bin/snp -i tmp/SNP_DB/$zone/$MAGIC.differential.snp_list_with_allele_frequency_per_sample.txt -db tmp/SNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC.differential -Reference_genome   "$Reference_genome" 

mv  tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt     tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1
cat tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1 | gawk '/^#/{print;next;}' > tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt
cat tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1 | gawk '/^#/{next;}{print;next;}' >> tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt
\rm tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1

# extract the list of genomic SNPs so i can reinject them in the exome study
cat   tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt | gawk -F '\t' '/^#/{next}{chrom=$1;pos=$2;typ=$4; gsub(">","2",typ);printf("%s:%d_%s\n",chrom,pos,typ);}' | sort -u >  tmp/SNP_DB/$zone/$MAGIC.genomic.snp.list

#### snp per gene should we do this here or after filtering
## by the way do we filter here or after gathering all zones 
#echo "## For each gene, we count all SNPs affecting this gene and report the number of observations as wild-type, low, mid, high, pure and the number of samples with at least 0 to 10 doses of any of these SNPs" >> $toto
# cat tmp/SNP_DB/zone*/$MAGIC.snp_prevalence_per_gene.txt 

exit 0

cat tmp/METADATA/RefSeq.MRNA.ln.ace | gawk '/^Transcript/{gsub(/\"/,"",$2);g=$2;next;}/^IntMap/{m=$2;a1=$3;a2=$4;if(a1<a2){ln=a2-a1+1;}else{ln=a1-a2+1;};printf("Product %s\nBest_product\n\nmRNA %s\nGene %s\nProduct %s 1 %d\nSplicing 1 %d 1 %d Exon 1\nDNA mRNA:%s %d\n",g,g,g,g,ln,ln,ln,g,ln);print;printf("\n");}' >  tmp/SNP_DB/RefSeq.ace
gunzip -c  TARGET/Targets/ec.RefSeq.fasta.gz  | gawk '/^>/{split($1,aa,"|");m=substr(aa[1],2);printf(">mRNA:%s\n",m);next;}{print}' >  tmp/SNP_DB/RefSeq.fasta

# analyse the concensus around a particular error : T>G is abundant behing GGt, the worst beeing GGTGG then GGTGn in E.coli or human exome or human RNA-seq, but disappears in LifeTech SOLiD
cat RESULTS/SNV/*.snp_list_with_allele_frequency_per_sample.txt | gawk -F '\t' '{if($28+$29+$30>20 && $11 == motif) {i=index($6,"t");print substr($6,i-2,5);}}' motif="T>G" | tags | sort -k 2n | tail -10


#########################################
#########################################
## identify the rsnumbers
## this script recovers a 3 column table : chrom  coord rsNumber for each chromosome in  tmp/SNP_RS
# info is in /am/ftp-snp/00README
if (! -d  tmp/SNP_RS) mkdir tmp/SNP_RS
set rsdir=/am/ftp-snp/organisms/human_9606_b141_GRCh37p13/chr_rpts/
set ref=GRCh37.p13
if (! -d  tmp/SNP_RS/$ref) mkdir tmp/SNP_RS/$ref
foreach chrom ($chromSetAll)
  if (-e tmp/SNP_RS/pos2rs.$chrom.txt2) continue
  set ff=$rsdir/chr_$chrom.txt.gz
  gunzip -c $ff | gawk -F '\t' '{rs=$1;c=$7;x=$12;withdrawn=$3;suspect=$23;clinvar=$24; if($22 == ref && c==chrom && x != " " && withdrawn == 0) printf ("%s\t%s\t%s\t%s\t%s\n",c,x,rs,suspect,clinvar);}' ref=$ref chrom=$chrom | sort -u | sort -k 1,1 -k2,2n > tmp/SNP_RS/$ref/pos2rs.$chrom.txt &
end

#########################################
#########################################

if (0) bin/snp -pheno -db tmp/SNP_DB/$zone

laba2:
bin/tacembly tmp/SNP_DB/$zone <<EOF
  query find variant
  show -a -f  tmp/SNP_DB/$zone/snp.translated.ace
  list -a -f tmp/SNP_DB/$zone/v.list
  query find variant mm
  table -active -o tmp/SNP_DB/$zone/s20.mm_snps.txt -f tables/s20.mm_snps.def
  query find variant FQ
  list -a -f tmp/SNP_DB/$zone/vfq.list
  query find variant mrna 
  list -a -f tmp/SNP_DB/$zone/vm.list
  query find variant translate
  list -a -f tmp/SNP_DB/$zone/vt.list
  query find Variant coding
  table -active -o tmp/SNP_DB/$zone/snp2coding.txt -f tables/snp2coding.def
  // table -o tmp/SNP_DB/$zone/snp.genotype.$zone -f tables/snp.genotype.run.def
  // table -o tmp/SNP_DB/$zone/snp.pos2pheno.$zone -f tables/snp.pos2pheno.def
  query find variant phenotype
  // table -active -o PHENO/chgMaxAggrColl2.chrom.$zone.txt  -f ../PHENO/variant2map2pheno.def CATEG chgMaxAggrColl2
  quit
EOF
cat  tmp/SNP_DB/$zone/snp2coding.txt  | gawk -F '\t' '{gsub(/\"/,"",$0);snp=$1;if(snp==old)next;old=snp;typ=$2;mrna=$7;if(substr(mrna,1,1)=="_");mrna=substr(mrna,2);x=1+(($8-$9)/3);a=$3;ap=$4;b=$5;bp=$6;printf("%s\t%s:%s:%s:%s:%d:%s:%s\n",snp,mrna,typ,a,ap,x,b,bp);}' > tmp/SNP_DB/$zone/snp2coding.nickname.txt

touch tmp/SNP_DB/$zone/snp.parse.done


echo -n "snp.done"
date

exit 0

set toto=RESULTS/SNV/$MAGIC.SNP_list_with_translation.txt
echo -n "# File $toto : " > $toto
date >> $toto
echo "# Only read pairs aligning uniquely, over at least 90% of their length and as compatible pairs (in close proximity and facing each other) are considered." >> $toto
echo "# To call a SNP or a 1 to 3 bases insertion or deletion, the variant has to be seen in at least one sample, where the variant allele frequency is above $minSnpFrequency %, the coverage above $minSnpCover fold and the variant has at least $minSnpCount supporting reads." >> $toto


cat MetaDB/$MAGIC/RunListSorted ZZZZZ MetaDB/$MAGIC/gtitle.txt ZZZZZ  tmp/SNP_DB/*/snp.translated.ace | gawk -f scripts/snp.translated.awk minSnpCover=$minSnpCover  >  $toto.1
cat $toto.1 | head -20 | gawk '/^#/{print}' >> $toto
cat  $toto.1 | gawk '/^#/{next;}{print}' | sort -k 14,14nr -k 15,15nr -k 16,16nr | sed -e 's/UUUU//g'  >> $toto
\rm $toto.1

###################
## SNP_WIGGLE tracks for UCSC

# 7 Synonymous
# 8 AA_substitution
# 9 Length_variation
# 10 Frameshift
# 11 UTR_5prime
# 12 UTR_3prime
# 13 Non_coding_transcript
#   ......   Kent has splicing i.e. near the intron boundary

if (! -d tmp/SNPW) mkdir tmp/SNPW
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (! -e tmp/SNPW/$group.bed ) then
    cat  tmp/SNP_DB/*/s20.mm_snps.txt | gawk -F '\t' '{gsub(/\"/,"",$0);snp=$1;g=$2;z=$3;c=$4;a1=$5;a2=$6;if(a1>a2){a0=a1;a1=a2;a2=a0};a1--;a2--;if (g==group)printf("chr%s\t%d\t%d\t%s\n",c,a1,a2,z);}' group=$group | sort -k 1,1 -k 2,2n >  tmp/SNPW/$group.bed &
  endif
end
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (-e  tmp/SNPW/$group.bed     && ! -e  tmp/SNPW/$group.bb    ) then
    bedToBigBed  -type=bed4 -tab tmp/SNPW/$group.bed  tmp/WIGGLE/$species.chrom.sizes  tmp/SNPW/$group.bb &
  endif
end

\rm toto
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  cat <<EOF >> toto
track $group'_SNP'
type bigBed
visibility dense
parent Primates_SNPS
maxHeightPixels 500:150:1
bigDataUrl ftp://ftp.ncbi.nlm.nih.gov/repository/acedb/human/bigwig.hg19/Tissue/$group.snp.bb
priority 4

EOF

end


###################
#########################################
#########################################

