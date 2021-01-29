#!bin/tcsh -f

set Strategy=$1
set zone=$2

  if (! -d tmp/TSNP_DB/$zone/database) then
    pushd tmp/TSNP_DB/$zone
        mkdir database
        ln -s ../../../MetaDB/wspec wspec
        ln -s ../../../metaData/tables.VariantDB tables
        whoami >> wspec/passwd.wrm
        echo y | ../../../bin/tacembly .
    popd

    if (-e TARGET/CHROMS/Centromere_telomeres.ace) then
      echo "pparse TARGET/CHROMS/Centromere_telomeres.ace" | bin/tacembly tmp/TSNP_DB/$zone -no_prompt
    endif

    if (! -e TARGET/CHROMS/Centromere_telomeres.ace && -e TARGET/CHROMS/GCF_report.txt) then
       cat GCF_000001405.22_GRCh37.p10_assembly_report.txt | grep assembled-molecule | cut -f 1,7 | sed -e 's/MT/mito/' | gawk '{printf("Map %s\nTitle %s\n\n",$1,$2);}' >  Chromosome_official_name.ace
       # find Map ; list -a -f titi.list
       cat titi.list | grep NT_ | gawk '/NT_/{print ;split($2,aa,"|");printf("Title \"%s\n\n",aa[2]);}' > titi2
       cat titi2 >> Chromosome_official_name.ace

    endif

    if (-e TARGET/CHROMS/Chromosome_official_name.ace) then
      # expecting  Map chr11 \n  Title NC_0012345.11 \n\n
      # this title will be used in the VCF and the HGVS Human Genome Variation Society full name
      echo "pparse TARGET/CHROMS/Chromosome_official_name.ace" | bin/tacembly tmp/TSNP_DB/$zone -no_prompt
    endif
  endif
  echo "read-models" > tmp/TSNP_DB/$zone/_rs3.$zone

 # retrofit the mRNA class name into the mRNA fasta file
  foreach target ( mito SpikeIn $Etargets $Ttargets)
    if (-e TARGET/Targets/$species.$target.fasta.gz && ! -e tmp/TSNP_DB/$species.$target.mrna.fasta.gz) then
      gunzip -c TARGET/Targets/$species.$target.fasta.gz | gawk '/^>/{i=index($1,"|");if(i>0)s=substr($1,2,i-2);else s=substr($1,2);printf(">mRNA:%s\n",s);next;}{print}' | gzip >  tmp/TSNP_DB/$species.$target.mrna.fasta.gz
    endif
  end
   

  foreach target (mito SpikeIn $Etargets $Ttargets)
    if (-e  tmp/METADATA/$target.MRNA.info.ace) then
      # gene gene_type geneId mrna mrna-length GC
      echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
      echo "pparse  tmp/METADATA/$target.MRNA.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
      echo "pparse  tmp/METADATA/$target.MRNA.splicing.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif 
    if (-e tmp/METADATA/gtf.$target.goodProduct.ace) then
      # mrna -> product -> good/best product 
      echo "pparse   tmp/METADATA/gtf.$target.goodProduct.ace"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif
    if (-e TARGET/GTF/$species.magic.good_product.ace.gz) then
      # mrna -> product -> good/best product 
      echo "pparse TARGET/GTF/$species.magic.good_product.ace.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif
    if (-e  tmp/TSNP_DB/$species.$target.mrna.fasta.gz) then
      echo "pparse tmp/TSNP_DB/$species.$target.mrna.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif 
  end

  # gene mrna product mrna-splicing intmap
  echo "parse TARGET/MRNAS/mrnaStructure.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse tmp/TSNP_DB/RefSeq.cds.ace.gz" >> tmp/TSNP_DB/$zone/_rs3.$zone
  # product quality should be in the gff file
  echo "pparse  TARGET/MRNAS/very_good_product.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse  TARGET/MRNAS/good_product.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  # mRNA-dna extracted from TARGET/Targets
  echo "pparse tmp/TSNP_DB/$species.av.mrna.fasta.gz" >> tmp/TSNP_DB/$zone/_rs3.$zone
  if (-e  tmp/TSNP_ZONE/$zone.fasta.gz && $Strategy != RNA_seq)   echo "pparse tmp/TSNP_ZONE/$zone.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
  # Gene titles and intmap
  echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse tmp/SNP_ZONE/$zone.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo save >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo quit >> tmp/TSNP_DB/$zone/_rs3.$zone

  echo -n "snp.parseChrom "
  date
  bin/tacembly tmp/TSNP_DB/$zone < tmp/TSNP_DB/$zone/_rs3.$zone

  touch  tmp/TSNP_DB/$zone/snp.parseChrom.done
  touch  tmp/TSNP_DB/$zone/tsnp0.done

endif
