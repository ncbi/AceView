#!bin/tcsh -f

set RG=$1
set rg=$2
set fr=$3

set uu=u 
set FR=$fr

if ($fr == pp) then
  set FR=frns 
  set uu=pp
endif
if ($fr == pp.f) then
  set FR=pp
  set uu=f
endif
if ($fr == pp.r) then
  set FR=pp
  set uu=r
endif
if ($fr == nu) then
  set FR=frns 
  set uu=nu
endif

echo -n "wg7.make_bigwig $rg any $fr start : "
date

if ($species == corona) then
  set chrom1="chrom=NC_045512"
  set chrom2="chrom=NC_045512v2"
else
  set chrom1="chrom="
  set chrom2="chrom=chr"
endif

echo "construct the wig file tmp/WIGGLE$RG/$rg/any.$fr.wig"

       if (-d tmp/WIGGLE$RG/$rg && ! -e  tmp/WIGGLE$RG/$rg/any.$fr.wig) then
         # if (-e tmp/WIGGLE$RG/$rg/any.$fr.wig) \rm tmp/WIGGLE$RG/$rg/any.$fr.wig
         if ($fr == frns || $fr == nu || $fr == pp || $fr == pp.f || $fr == pp.r ) then
           foreach chrom ($chromSetAll)
             if (! -e  tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz) continue
             set n=`gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz | wc | gawk '{printf("%d", $1 -1); }'`
             if ($n > 0) then
               gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz | head -$n  | gawk '/^#/{next}/^track/{next;}/^fixedStep/{if(index($0,"start=10 step=10")>0){printf("%s %s start=5 %s\n",$1,$2,$4);jump=1;next;}}{if(jump==1){jump=0;next;}}{print}' | sed -e "s/$chrom1/$chrom2/"  | grep -v chrET_av >> tmp/WIGGLE$RG/$rg/any.$fr.wig
>             endif
           end
         else
           foreach chrom ($chromSetAll) 
             if (! -e  tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz) continue
             set n=`gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz | wc | gawk '{printf("%d", $1 -1); }'`
             if ($n >0) then
                gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz | head -$n  | gawk '/^#/{next}/^track/{next;}/^fixedStep/{if(index($0,"start=10 step=10")>0){printf("%s %s start=5 %s\n",$1,$2,$4);jump=1;next;}}{if(jump==1){jump=0;next;}}{print}' | sed -e "s/$chrom1/$chrom2/"  | grep -v chrET_av >> tmp/WIGGLE$RG/$rg/any.$fr.wig
             endif
           end
         endif
       endif

       if (! -d  ~/ftp/$species/bigwig.$UCSCgenomeRelease) mkdir  ~/ftp/$species/bigwig.$UCSCgenomeRelease
       if (! -d  ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC) mkdir  ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC
       if (! -d  ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg) mkdir  ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg

echo "construct the bw file tmp/WIGGLE$RG/$rg/any.$fr.bw"       

       if (-e tmp/WIGGLE$RG/$rg/any.$fr.wig && -e  tmp/WIGGLE/$species.chrom.sizes  && ! -e tmp/WIGGLE$RG/$rg/any.$fr.bw) then
         echo "wigToBigWig tmp/WIGGLE$RG/$rg/any.$fr.wig  tmp/WIGGLE/$species.chrom.sizes  tmp/WIGGLE$RG/$rg/any.$fr.bw"
               wigToBigWig tmp/WIGGLE$RG/$rg/any.$fr.wig  tmp/WIGGLE/$species.chrom.sizes  tmp/WIGGLE$RG/$rg/any.$fr.bw
      endif
# allow propagation time for the file name
      sleep 5  
      if (-e tmp/WIGGLE$RG/$rg/any.$fr.bw && ! -e ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg/any.$fr.bw )  then
    echo "mv tmp/WIGGLE$RG/$rg/any.$fr.bw         ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg/any.$fr.bw"
          mv tmp/WIGGLE$RG/$rg/any.$fr.bw         ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg/any.$fr.bw
      endif

ls -ls  ~/ftp/$species/bigwig.$UCSCgenomeRelease/$MAGIC/$rg
echo -n "done : "
date
attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaa

 attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttc          65 tctaaacgaac tttaa
                                            28240 ttgttttagatttca       28255 tctaaacgaac aaactaaaatgtctg
 21532                                        agtgatgttcttgttaacaa      21552  ctaaacgaac aatgtttgttt
                                                   tcaaattacattaca      25380 cataaacgaac ttatggatt
                                                                gatgagt 26237      acgaac ttatgtac
                                                gagttcctgatcttctgg      26468 tctaaacgaac taaatattatat
                                                             ttgctacatc 27041      acgaac gctttcttattacaaattgggagcttcgcagcgtg
                                                  agcaaccaatggagattgat  27385   taaacgaac atgaaaattattcttttcttgg
                                              cataatgaaacttgtcacgc      27884  ctaaacgaac atgaaa


