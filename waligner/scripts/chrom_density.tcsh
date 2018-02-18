#!bin/tcsh -f

set target=$1
set GM0=$2
set ccc0=$3

set GM=$GM0
echo "Genes and DEG chromososomal density"

set toto=tmp/ChromDensity/$MAGIC.$GM.toto

    set target2=$target
    if ($target == av) set target2=AceView

###############################################
## for each ccc

echo '#' > $toto.anyComparison.1
echo '#' > $toto.anyComparison.-1

foreach ccc (`cat  MetaDB/$MAGIC/ccc.list` anyComparison)
  if ($ccc0 != anyComparison && $ccc0 != $ccc) continue

  set compare=`echo $ccc | gawk '{split($1,aa,",");print aa[1];}'`

  set ff1=`echo $ccc | gawk '{split($1,aa,",");print aa[2];}'`
  set ff2=`echo $ccc | gawk '{split($1,aa,",");print aa[3];}'`

  foreach GM ($GM0)
    foreach GM2 (self m2g)
      if ($GM == GENE && $GM2 == m2g) continue
      set GM3=$GM
      if ($GM2 == m2g) set GM3=GENEviaMRNAH
 
      set tutu=RESULTS/ChromDensity/$MAGIC.$target2.$compare.$GM3.txt

      if (-e $tutu.999) \rm $tutu.999
      echo "# " > $tutu
      date >> $tutu
      echo "# $target $GM3 : Density distribution along the genome" >> $tutu

      if ($compare == anyComparison) then
        echo "# restricted to the union of all differential $GM3 annotated in AceView for the $MAGIC project" >> $tutu
      else
        echo "# restricted to the $compare differential genes between $ff2 and $ff1 " >> $tutu
      endif
      echo "# Class\tChromosome\tSection start\tSection end\tExperiment\tUp\tDown\tFeatures in section\t% up\t% down\t% up expected\t%d expected\tchi2 up\tchi2 down\tn-score up\tn-score down\tHighest cases\tLowest cases\tCancer Census genes " >> $tutu

      foreach sign (1 -1)

        echo "------ $ccc up/down = $sign $GM $GM2"

        if ($compare == anyComparison) then
          # mv $toto.$ccc.$sign.$GM2 $toto.$ccc.$sign
	  touch toto.$ccc.$sign.$GM2 
        else
          echo '#' > $toto.$ccc.$sign

    ######  export the list of up and down differential genes or transcripts in  $toto.$ccc.$sign

          if ($sign == 1) then
            set ff1=`echo $ccc | gawk '{split($1,aa,",");print aa[2];}'`
            set ff2=`echo $ccc | gawk '{split($1,aa,",");print aa[3];}'`
          else
            set ff1=`echo $ccc | gawk '{split($1,aa,",");print aa[3];}'`
            set ff2=`echo $ccc | gawk '{split($1,aa,",");print aa[2];}'`
          endif

          set ff=RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target2.$GM.u.$ff1'_'$ff2.diffGenes.0.txt.with_metadata.txt 
 
          cat $ff |  gawk -F '\t' '/^#/{next}{if (0+$2 > 5)printf ("%s\t%s\n",$1,$2);}'   >>  $toto.$ccc0.$sign.$GM2
          # cat $ff |  gawk -F '\t' '/^#/{next}{if (0+$2 > 5)printf ("%s\t%s\n",$1,$2);}'  >> $toto.anyComparison.$sign.$GM2
        endif

        if ($ccc0 == anyComparison && $ccc0 != $ccc) continue
        mv $toto.$ccc.$sign.$GM2 $toto.$ccc.$sign

          ##### Compute the density of the selected objects along the genome
          ## we attribute each gene to its most 3' potition in the orientation of the genome
          ## this is an approximation, but sufficient for computing a density
          ## selected genes/mrnas ZZZZZ info about mrna/gene mapping<->gene
  
          ## we already have the selected list: append the other files and zip these metadata

        echo ZZZZZ                                         >>  $toto.$ccc.$sign
        cat tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt >>  $toto.$ccc.$sign
        echo ZZZZZ                                         >>  $toto.$ccc.$sign
        cat TARGET/GENES/CancerCensus.txt                  >>  $toto.$ccc.$sign
        echo ZZZZZ                                         >>  $toto.$ccc.$sign

        cat $toto.$ccc.$sign | gawk -F '\t' -f scripts/chrom_density.awk  target=$target GM=$GM GM2=$GM2 title=$compare sign=$sign | sort >> $tutu.999
        \rm $toto.$ccc.$sign
      end
      cat  $tutu.999 |  sort -k 1,1 -k 2,2 -k 3,3n -k 4,4n -k 5,5 -k 6,6nr | gawk -F '\t' '{nn++;nn=nn%2;if(nn==1){for(i=1;i<=NF;i++)a[i]=$i;next;}printf("%s\t%s\t%d\t%d\t%s",$1,$2,$3,$4,$5);for(i=7;i<=NF;i++){if(i==7)printf("\t%d\t%d\t%d",a[i],$i,$11);if(i==12 || i==13 || i==14)printf("\t%.2f\t%.2f",a[i],$i);if(i==15)printf("\t%d\t%d",a[i],$i);if(i==17)printf("\t%s\t%s\t%s",a[i],$i,$16);}printf("\n");}' >  $tutu.888
      cat   $tutu.888 | sort -k 1,1 -k 5,5  -k 2,2n -k 3,3n | gawk '{if ($2 != $2 + 0)next;if($2!=old)printf("\n\n\n\n\n\n\n\n\n\n");old=$2;print}' >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      cat   $tutu.888 | sort -k 1,1 -k 5,5  -k 2,2n -k 3,3n | gawk '{if ($2=="X")print}' >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      cat   $tutu.888 | sort -k 1,1 -k 5,5  -k 2,2n -k 3,3n | gawk '{if ($2=="Y")print}' >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      echo >> $tutu
      cat   $tutu.888 | gawk '{if ($2 == $2 + 0 || $2=="X" || $2=="Y")next;print}'  | sort -k 1,1 -k 5,5  -k 2,2  -k 3,3n >> $tutu
      if (-e $tutu.888) rm $tutu.888  
      if (-e $tutu.999) rm $tutu.999
    end
  end
end

##################################################################
touch tmp/ChromDensity/$target.$GM.$ccc0.done

