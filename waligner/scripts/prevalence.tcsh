#!bin/tcsh -f

set target=$1
set GM=$2
set ccc0=$3

      set target2=$target
      if ($target == av) set target2=AceView

set toto=tmp/Prevalence/$MAGIC.$GM.toto

echo "Prevalence statistics"

echo '#' > $toto.anyComparison

foreach ccc (anyGene `cat  MetaDB/$MAGIC/ccc.list` anyComparison)
  if ($ccc0 != anyComparison && $ccc0 != $ccc) continue

  echo "------ $ccc"
  set compare=`echo $ccc | gawk '{split($1,aa,",");print aa[1];}'`
  if ($compare == anyComparison) then
    touch toto.$ccc
   else if ($compare == anyGene) then
    cat tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt | cut -f 5 | sort -u > $toto.anyGene
  else

    echo '#' > $toto.$ccc

    set compare=`echo $ccc | gawk '{split($1,aa,",");print aa[1];}'`
    ######  export the list of up and down differential genes or transcripts in  $toto.$ccc
    foreach sign (1 -1)
      if ($sign == 1) then
        set ff1=`echo $ccc | gawk '{split($1,aa,",");print aa[2];}'`
        set ff2=`echo $ccc | gawk '{split($1,aa,",");print aa[3];}'`
      else
        set ff1=`echo $ccc | gawk '{split($1,aa,",");print aa[3];}'`
        set ff2=`echo $ccc | gawk '{split($1,aa,",");print aa[2];}'`
      endif

      if (-e RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$compare.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.1000.txt && ! -e RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$compare.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.txt) then 
        mv RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$compare.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.1000.txt  RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$compare.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.txt
      endif
      set ff=RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$compare.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.txt 
 
      cat $ff |  gawk -F '\t' '/^#/{next}{print $1}'   >>  $toto.$ccc0
    end

  endif

  if ($ccc0 == anyComparison && $ccc0 != $ccc) continue
 
  ##### In how many individual are these genes expressed for the union
  ## the expected format is
  ## selected genes/mrnas ZZZZZ info about mrna<->gene ZZZZZ info about gene-types ZZZZZ index->expressed ZZZZZ read-count->touched
  
  ## we already have the selected list: append the other files and zip these metadata

  echo ZZZZZ                                         >>  $toto.$ccc
  cat tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt >>  $toto.$ccc
  echo ZZZZZ                                         >>  $toto.$ccc
  cat tmp/METADATA/$target.metadata.txt              >>  $toto.$ccc
 
  if (-e  $toto.$ccc.gz ) \rm   $toto.$ccc.gz
  gzip  $toto.$ccc

  ## Search the expression tables and the read counts (for the touch counts)

    foreach hasGid (any with no)
      foreach GM2 (self m2g)
        if ($GM == GENE && $GM2 == m2g) continue
        set GM3=$GM
        if ($GM2 == m2g) set GM3=GENEviaMRNAH
        set tutu=RESULTS/PrevalenceProfiles/$MAGIC.$target2.$compare.$GM3.$hasGid'_'geneid.txt
        echo -n "# " > $tutu
        date >> $tutu
        echo "# $target $GM : Genes expressed per gene type in at least N individuals, % of annotated genes" >> $tutu
        if ($compare == anyComparison) then
          echo "# resticted to the union of all differential genes annotated in AceView for the $MAGIC project" >> $tutu
        else if ($compare == anyGene) then
          echo "# considering all genes annotaged in AceView" >> $tutu
        else
          echo "# restricted to the $compare differential genes between $ff2 and $ff1 " >> $tutu
        endif
        echo "# statistics derived from the file Expression/unique/$MAGIC.$target2.$GM.u.expression_index" >> $tutu

        gunzip -c $toto.$ccc.gz ZZZZZ.gz RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.expression_index.txt.gz ZZZZZ.gz RESULTS/Expression/unique/$target/$MAGIC.$target2.$GM.u.reads_aligned_per_gene.txt.gz | gawk -F '\t' -f scripts/prevalence.awk  hasGid=$hasGid GM=$GM GM2=$GM2 title=$compare | sort >> $tutu
      end
    end

  \rm $toto.$ccc.gz

end

touch tmp/Prevalence/$target.$GM.$ccc0.done

exit 0




