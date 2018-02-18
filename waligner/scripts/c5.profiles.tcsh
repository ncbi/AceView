#!bin/tcsh -f

set run=$1

set title=`cat MetaDB/$MAGIC/_q2 | gawk '{gsub(/NULL/," ",$7);if($1 == g)print $7;}' g=$run | sort -u` 

set mxLn=`cat Fastc/$run/Max_probe_length | gawk '{print $1;}'`
if ($mxLn == 0) set mxLn=100
foreach prefix (f f1 f2)
  foreach readNumber (0 1 2)
echo "$prefix $readNumber"
    set n=`cat Fastc/$run/LaneList | sed -e "s/$run\///" | gawk '{if (NF >= 1) {split($1,aa,".");if(aa[1]==p)n++}}END{print n+0}' p=$prefix`
    if ($n == 0)  continue
echo "$prefix $readNumber .. $n"

    foreach type (prefix suffix target_prefix target_suffix clipped_adaptor unaligned)
echo "$prefix $readNumber .. $n .. $type"

      if ($type == unaligned && ! -e tmp/Unaligned/$run/$prefix.1.fastc.gz) continue
      if ($type != unaligned && ! -e tmp/COUNT/$run/$prefix.1.$type.$readNumber.gz) continue

      if ($type == prefix) set tt="complement of the unaligned 5-prime part of the sequences"
      if ($type == suffix) set tt="unaligned 3-prime part of the sequences"
      if ($type == target_prefix) set tt="Target sequence immediatly upstream of the alignments (which starts at 31)"
      if ($type == target_suffix) set tt="Target sequence immediatly downstream of the alignments (which ends at 0)"
      if ($type == clipped_adaptor) set tt="3-prime part of the reads recognized as belonging to the adaptor and clipped"
      if ($type == unaligned) set tt="unaligned sequences"
      ls -ls  tmp/Profiles/$run/$prefix.$type.$readNumber.profile.txt
      if (! -e tmp/Profiles/$run/$prefix.$type.$readNumber.profile.txt) then
        if ($type == unaligned) then
          if (-e tmp/Unaligned/$run/$prefix.1.fastc.gz || -e tmp/Unaligned/$run/$prefix.2.fastc.gz) then
            if ($readNumber == 1) then
               touch  tmp/Profiles/$run/$prefix.$type.$readNumber tmp/Profiles/$run/$prefix.$type.$readNumber.tcsh tmp/Profiles/$run/$prefix.$type.$readNumber.ps
               gunzip -c  tmp/Unaligned/$run/$prefix.*.fastc.gz  | bin/dna2dna -minEntropy 2 -gnuplot $mxLn  -I fastc -title "Unaligned sequences in $run $prefix" -o tmp/Profiles/$run/$prefix.$type.$readNumber
            endif
          endif
          echo "$prefix $readNumber .. $n .. $type ----" 

        else
          if (-e tmp/COUNT/$run/$prefix.1.$type.$readNumber.gz || -e tmp/COUNT/$run/$prefix.2.$type.$readNumber.gz) then
            touch  tmp/Profiles/$run/$prefix.$type.$readNumber tmp/Profiles/$run/$prefix.$type.$readNumber.tcsh tmp/Profiles/$run/$prefix.$type.$readNumber.ps
            gunzip -c  tmp/COUNT/$run/$prefix.*.$type.$readNumber.gz | bin/dna2dna -minEntropy 2 -gnuplot 40  -I raw -title "$type in $run $prefix : $title $tt " -o tmp/Profiles/$run/$prefix.$type.$readNumber
            if ($type ==  target_prefix || $type == target_suffix)  \rm tmp/COUNT/$run/$prefix.*.$type.$readNumber.gz 
          endif
        endif
      endif
    end
  end 
end

# echo "hello tmp/Profiles/$run/c5.profile.done"
# echo "touch       tmp/Profiles/$run/c5.profile.done"
touch       tmp/Profiles/$run/c5.profile.done
