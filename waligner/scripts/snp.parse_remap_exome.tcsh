#!bin/tcsh -f

set Strategy=$1
set zone=$2
set snp_type=".homozygous"
set snp_type=".filtered"
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

  foreach target (av mito SpikeIn)
    if (-e  tmp/METADATA/$target.MRNA.info.ace2) then
      # gene gene_type geneId mrna mrna-length GC
      echo "pparse  tmp/METADATA/$target.MRNA.info.ace2" >> tmp/SNP_DB/$zone/_rs3.$zone
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
  touch touch  tmp/SNP_DB/$zone/snp.parseChrom.done

endif

if (! -e   tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace) then
  echo -n "snp.prepare file  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace "
  date

   # split the file because gawk chokes on very large files and acedb can reconstitute the pieces
    cat  tmp/SNPH/$zone/$MAGIC.snp$snp_type.sorted   | sort -u > tmp/SNPH/toto.$$.s
    split -l 20000  -a 6 tmp/SNPH/toto.$$.s tmp/SNPH/toto.$$.s.p.
    echo ' ' >  tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace

    foreach ff (`ls  tmp/SNPH/toto.$$.s.p.*`)
      echo $ff
      cat $ff | gawk -F '\t' -f scripts/snp2ace.awk format=txt targetType=$targetType >>   tmp/SNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace
    end

    \rm  tmp/SNPH/toto.$$  tmp/SNPH/toto.$$.*  
endif

# update the Runs
bin/tacembly tmp/SNP_DB/$zone <<EOF
       read-models
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
       save
       quit
EOF
endif

if (! -e  tmp/SNP_DB/$zone/snp.remap.done) then
  if ($zone == mito) then
    # bin/snp -db_remap2genes tmp/METADATA/mrnaRemap.mito.txt.gz  -db tmp/SNP_DB/$zone
  endif
  if (1) then 
    echo -n "snp.Remap $zone remapping against  tmp/METADATA/mrnaRemap.gz  "
    date
    bin/snp -db_remap2genes tmp/METADATA/mrnaRemap.gz  -db tmp/SNP_DB/$zone
  endif
  touch tmp/SNP_DB/$zone/snp.remap.done 
endif

if (! -e  tmp/SNP_DB/$zone/snp.translate.done) then
  echo -n "snp.translate  $zone "
  date
  bin/tacembly tmp/SNP_DB/$zone <<EOF
    read-models
    save
    quit
EOF

  bin/snp -db_translate -db tmp/SNP_DB/$zone
  touch tmp/SNP_DB/$zone/snp.translate.done
endif

laba:
if (-e RESULTS/StLouisSNP/zone.stl.ace) then
  echo "pparse  RESULTS/StLouisSNP/zone.stl.ace"  | bin/tacembly tmp/SNP_DB/$zone -no_prompt 
endif
echo "bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -o  tmp/SNP_DB/$zone/$MAGIC       -Reference_genome  $Reference_genome"
      bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome  $Reference_genome 

exit 0


 gunzip -c RESULTS/StLouisSNP/*.txt.gz | gawk '/^##/{next}/^#/{for(i=1;i<=NF;i++)tt[i]=$i;next;}{n=length($4);if(n==3){typ= substr($4,1,1) "2" substr($4,3,1); nam=$1 ":" typ ; printf("Variant %s\nTyp %s\nSaint_Louis\nIntMap %s %d 1\n", nam, typ,substr(nam,1,1), substr(nam,3));for(i=2;i<=NF;i++)printf("Saint_Louis \"%s\" \"%s\"\n",tt[i],$i);printf("\n");}}' >  RESULTS/StLouisSNP/zone.stl.ace 
 gunzip -c RESULTS/StLouisSNP/*.txt.gz | gawk '/^##/{next}/^#/{for(i=1;i<=NF;i++)tt[i]=$i;next;}{n=length($4);if(n>3){split($4,aa,";");split(aa[1],bb,">");n1=length(bb[1]);n2=length(bb[2]);if(n1>n2)typ="Del" substr(bb[1],n2+1);else typ="Ins" substr(bb[2],n1+1); nam=$1 ":" typ ; printf("Variant %s\nTyp %s\nSaint_Louis\nIntMap %s %d 1\n", nam, typ,substr(nam,1,1), substr(nam,3));for(i=2;i<=NF;i++)printf("Saint_Louis \"%s\" \"%s\"\n",tt[i],$i);printf("\n");}}' >>  RESULTS/StLouisSNP/zone.stl.ace 

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

#### Retina, export a particular zone

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
## histo des snps
cat tmp/NEWHITS_snp/*/s590.multi.countEdited |  gawk -F '\t' '/^#/{next}{if ($8 > 10 && $9 >= minCover && $10 >= 0) print $9,$10;}' minMut=$minSnpCount minCover=$minSnpCover minF=$minSnpFrequency | bin/histo -smooth -plot -o RESULTS/SNV/$MAGIC.snp_s590_plot

# construct a 4 column file: snp,run,m,w
gunzip -c tmp/SNPH/zoner.[56]/$MAGIC.snp.sorted | gawk -F '\t' '/^#/{next;}{snp=$1 ":" $2 ":" $3;run=$6;c=$9;m=$10;w=$11;if(c>=minCover)printf("%s\t%s\t%d\t%d\n",snp,run);}' minCover=10 > tmp/SNPH/$MAGIC.sorted.distrib

mkdir  tmp/SNPH/Histo  tmp/SNPH/Histo/$MAGIC
cat tmp/SNPH/$MAGIC.sorted.distrib | bin/histo -KL $minSnpFrequency -gzo -o   tmp/SNPH/Histo/$MAGIC"

