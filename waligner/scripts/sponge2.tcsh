#!bin/tcsh -f
set run=$1
set limit=$2
set mask="$3"
set mainTarget=$4

set uu=u

      foreach chrom (mito $chromSetAll)
        if ($chrom == Un) continue
        foreach fr1 (frns)
          if (! -e tmp/SPONGE/$run/Total.$chrom.$uu.$fr1.$limit.txt ||  ! -e tmp/SPONGE/$run/$mainTarget.gene.$chrom.$uu.ns.1) then
            echo "scripts/sponge_chrom.tcsh $run $chrom $uu $fr1 $limit $mask $mainTarget"
                  scripts/sponge_chrom.tcsh $run $chrom $uu $fr1 $limit $mask $mainTarget
          endif
        end
      end
  
exit 0
   
