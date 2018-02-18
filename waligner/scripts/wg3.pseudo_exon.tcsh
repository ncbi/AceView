#!bin/tcsh -f
echo "wg3.pseudo_exon.tcsh is obolete"
exit 0

set chrom=$1
set run=$2
set minCover=$3
 
# set WnewExon=`cat MetaDB/$MAGIC/GroupW_new_exonList | head -1`

# contruct the double stranded frns wiggle for pseudo exon discovery

foreach uu ( u nu)

  if ($uu == nu) continue
#group case
  if (-e tmp/WIGGLEGROUP/$run/$chrom.$uu.ns.BF.gz  || -e tmp/WIGGLEGROUP/$run/$chrom.$uu.f.BF.gz) then
    if (! -e tmp/WIGGLEGROUP/$run/$chrom.frns.$uu.BF.gz) then
      gunzip -c tmp/WIGGLEGROUP/$run/$chrom.$uu.*.BF.gz | bin/wiggle -I BF -O BF -gzo -o tmp/WIGGLEGROUP/$run/$chrom.frns.$uu
    endif
  endif

#run case
  if (-e tmp/WIGGLERUN/$run/$chrom.$uu.ns.BF.gz  || -e tmp/WIGGLERUN/$run/$chrom.$uu.f.BF.gz) then
    if (! -e tmp/WIGGLERUN/$run/$chrom.frns.$uu.BF.gz) then
      gunzip -c tmp/WIGGLERUN/$run/$chrom.$uu.*.BF.gz | bin/wiggle -I BF -O BF -gzo -o tmp/WIGGLERUN/$run/$chrom.frns.$uu
    endif
  endif
# evaluate the coverage per threshold

  if (-e  tmp/WIGGLEGROUP/$run/$chrom.frns.$uu.BF.gz) then

    foreach cover ($minCover 1 2 5 10 20 50 100 200 500 1000)

      if (! -e  tmp/WIGGLEGROUP/$run/$chrom.pseudoexon.$cover.$uu.txt) then
        bin/wiggle -i  tmp/WIGGLEGROUP/$run/$chrom.frns.$uu.BF.gz  -I BF -O BV -minCover $cover | gawk -F '\t' '/^#/{next;}/^track/{next;}/^variableStep/{if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;oldx1=0;oldx2=0;oldy=0;if(x2==0)x2=1;printf("#Chromosome %s Length %.3f Mb, Usable %.3f Mb Covered over %d : %.3f%%\n",chrom, x2/1000000,nn/1000000,minCover,100*nn/x2);nn=0;x2=0;chrom=$1;oldy=0;gsub(/variableStep chrom=/,"",chrom);next;} {n++;if(n<4)next;x1=$1-10;x2=$1+10;if($2>oldy)oldy=$2;if(x1>oldx2+50){if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;oldx1=x1;oldy=0;}oldx2=x2;if($2>oldy)oldy=$2;}END{if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;if(x2==0)x2=1;printf("#Chromosome %s Length %.3f Mb, Usable %.3f Mb Covered over %d : %.3f%%\n",chrom, x2/1000000,nn/1000000,minCover,100*nn/x2);}'  minCover=$cover >   tmp/WIGGLEGROUP/$run/$chrom.pseudoexon.$cover.$uu.txt
      endif
    end

  endif

end

# compute the strand cross correlation

if (! -e tmp/WIGGLEGROUP/$run/$chrom.strand_shift.txt) then
  if (-e tmp/WIGGLEGROUP/$run/$chrom.u.f.BF.gz) then
    bin/wiggle -strand_shift 500 -ssf tmp/WIGGLEGROUP/$run/$chrom.u.f.BF.gz -ssr tmp/WIGGLEGROUP/$run/$chrom.u.r.BF.gz -I BF -o tmp/WIGGLEGROUP/$run/$chrom 
  else 
    if (-e tmp/WIGGLEGROUP/$run/K.$chrom.u.f.BF.gz) then
      bin/wiggle -strand_shift 500 -ssf tmp/WIGGLEGROUP/$run/K.$chrom.u.f.BF.gz -ssr tmp/WIGGLEGROUP/$run/K.$chrom.u.r.BF.gz -I BF -o tmp/WIGGLEGROUP/$run/$chrom 
    endif
  endif
endif



touch tmp/WIGGLEGROUP/$chrom/pseudoexon.done
exit 0


