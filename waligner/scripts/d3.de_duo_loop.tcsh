#!bin/tcsh -f

set chrom=$1
set run=$2
set section=$3
set minX=$4
set maxX=$5
set Strategy=$6


      foreach ii (50 20 10 8 5 2 1)
        set oL=$overhangLength
        if ($species == hs && $ii >= 20) then
          set maxI=100000
          @ oL=$oL + 2
        endif
        set isany=0
# $tissues any
        foreach dd (OR Transloc) 

          if ($Strategy != RNA_seq && $dd != OR) continue 
          if ($Strategy != RNA_seq && $ii != 8) continue 
          if ($species != worm && $dd =~ SL*) continue
          if ($ii > 1 && $dd != OR) continue
          if ($run != RefSeq && $run != Mrna && $ii < 5 && $dd == OR) continue
          if (! -d tmp/$dd) mkdir tmp/$dd

          if ($dd == Transloc) then
            if ($ii != 8) continue
            if ($chrom != $chrom2) continue 
            if (-e tmp/$dd/$run/$run.txt) continue
          else
            if (! -d tmp/$dd/$run/$chrom) mkdir tmp/$dd/$run/$chrom
            if (-e tmp/$dd/$run/$chrom/$run.$ii.$section.txt) continue
          endif

          echo "scripts/d3.de_duo.tcsh $dd $chrom $run $ii $intronMaxLength  $minX $maxX  $oL $section $Strategy"
          scripts/d3.de_duo.tcsh $dd $chrom $run $ii $intronMaxLength  $minX $maxX  $oL $section $Strategy
        end
      end

  if (-e tmp/introns/$run/d4.discovery.counts) \rm  tmp/introns/$run/d4.discovery.counts*
  touch  tmp/OR/$run/$chrom/d3.de_duo_loop.$run.$section.all.done

