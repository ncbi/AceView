#!bin/tcsh -f
set run=$1
set limit=$2
set mask="$3"
set mainTarget=$4

set uu=u

      foreach chrom (mito $chromSetAll)
        if ($chrom == Un) continue
        if (! -e tmp/SPONGE/$run/Total.$chrom.$uu.frns.$limit.txt) then
          echo "scripts/sponge_chrom.tcsh stats $run $chrom $uu frns $limit $mask $mainTarget"
                scripts/sponge_chrom.tcsh stats $run $chrom $uu frns $limit $mask $mainTarget
        endif
        if ($limit == 1) then
          foreach fr (f r)
            if (! -e tmp/SPONGE/$run/$mainTarget.gene.v4.$chrom.$uu.$fr.1) then
              echo "scripts/sponge_chrom.tcsh pressGene $run $chrom $uu $fr $limit $mask $mainTarget"
                    scripts/sponge_chrom.tcsh pressGene $run $chrom $uu $fr $limit $mask $mainTarget
            endif
          end
        endif
      end
touch  tmp/SPONGE/$run/sponge2.gene.$limit.done5
exit 0
   
