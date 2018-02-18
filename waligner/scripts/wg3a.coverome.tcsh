#!bin/tcsh -ef

set group=$1
set chrom=$2
set minCover=$3

set uu=$4

set out=tmp/WIGGLEGROUP/$group/$chrom/wg3a.coverome.$uu
if (-e $out) \rm $out

echo verify there is no error reported 
# grep Status tmp/WIGGLEGROUP/*/*.err | grep -v  '= 0' 


  if($chrom == Un) continue
#group case

  set ok=0
  if (-e $out.txt) \rm $out.txt
  foreach fr (f r ) 
    if (! -e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz) continue
    @ ok = $ok + 1
    echo "ok=$ok"
    ls -ls tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.cumul.gz && $chrom != SpikeIn) then
      gunzip -c tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.cumul.gz | gawk -F '\t' '/^Cumul_at_level/{split ($1,aa,"_");c=chrom;if (c == "SpikeIn") c = "spikeIn"; printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n", group, c, fr,aa[4],$2,$4,$6);if(c != "mito" && c != "spikeIn"){c="genome"; printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n", group, c, fr,aa[4],$2,$4,$6);}}' fr=$fr chrom=$chrom group=$group >> $out.txt
    endif
  end 

  if ($ok > 0) then
    if (! -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu.cumul.gz) then
      gunzip -c tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.[fr].BF.gz | bin/wiggle -I BF -O BF -gzo -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu -cumul
    endif
    gunzip -c tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu.cumul.gz | gawk -F '\t' '/^Cumul_at_level/{split ($1,aa,"_");c=chrom;if (c == "SpikeIn") c = "spikeIn"; printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n", group, c, fr,aa[4],$2,$4,$6);if(c != "mito" && c != "spikeIn"){c="genome"; printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n", group, c, fr,aa[4],$2,$4,$6);}}' fr=ns chrom=$chrom group=$group >> $out.txt
  endif
# evaluate the coverage per threshold

  if (-e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu.BF.gz) then

    foreach cover ($minCover 1 2 5 10 20 50 100 200 500 1000)
      if (! -e  tmp/WIGGLEGROUP/$group/$chrom/coverome.$cover.$uu.txt) then
        gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu.BF.gz | bin/wiggle -I BF -O BV -gauss 20 -minCover $cover -peaks -o tmp/WIGGLEGROUP/$group/$chrom/coverome.$cover.$uu
    end

  endif

  foreach fr (f r )
    if (-e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz) then
      if (! -e  tmp/WIGGLEGROUP/$group/$chrom/coverome_$fr.1.$uu.peaks) then
        foreach cover ($minCover 1 2 5 10 20 50 100 200 500 1000)
          gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz | bin/wiggle -I BF -O BV -gauss 20 -minCover $cover -peaks -o tmp/WIGGLEGROUP/$group/$chrom/coverome_$fr.$cover.$uu
          \rm tmp/WIGGLEGROUP/$group/$chrom/coverome_$fr.$cover.$uu.BV
        end
      endif
    endif
  end

# compute the strand cross correlation, using 

  if (! -e tmp/WIGGLEGROUP/$group/$chrom/strand_shift.txt) then
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.u.EL.BF.gz && -e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.u.ER.BF.gz) then
      bin/wiggle -strand_shift 500 -ssf tmp/WIGGLEGROUP/$group/$chrom/R.chrom.u.EL.BF.gz -ssr tmp/WIGGLEGROUP/$group/$chrom/R.chrom.u.ER.BF.gz -I BF -o tmp/WIGGLEGROUP/$group/$chrom/$uu
    endif
  endif

  touch tmp/WIGGLEGROUP/$group/$chrom/wg3a.coverome.$uu.done


touch $out.done

exit 0
