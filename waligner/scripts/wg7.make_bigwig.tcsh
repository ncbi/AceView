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
if ($fr == nu) then
  set FR=frns 
  set uu=nu
endif

echo -n "wg7.make_bigwig $rg any $fr start : "
date

echo "construct the wig file tmp/WIGGLE$RG/$rg/any.$fr.wig"

       if (-d tmp/WIGGLE$RG/$rg && ! -e  tmp/WIGGLE$RG/$rg/any.$fr.wig) then
         # if (-e tmp/WIGGLE$RG/$rg/any.$fr.wig) \rm tmp/WIGGLE$RG/$rg/any.$fr.wig
         if ($fr == frns || $fr == nu || $fr == pp ) then
           foreach chrom ($chromSetAll)
             if (! -e  tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz) continue
             set n=`gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz | wc | gawk '{printf("%d", $1 -1); }'`
             if ($n > 0) then
               gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.$FR.$uu.BF.gz | head -$n  | gawk '/^#/{next}/^track/{next;}/^fixedStep/{if(index($0,"start=0")>0){printf("%s %s start=0 %s\n",$1,$2,$4);jump=1;next;}}{if(jump==1){jump=0;next;}}{print}' | sed -e 's/chrom=/chrom=chr/'  | grep -v chrET_av >> tmp/WIGGLE$RG/$rg/any.$fr.wig
             endif
           end
         else
           foreach chrom ($chromSetAll) 
             if (! -e  tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz) continue
             set n=`gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz | wc | gawk '{printf("%d", $1 -1); }'`
             if ($n >0) then
                gunzip -c tmp/WIGGLE$RG/$rg/$chrom/R.chrom.u.$fr.BF.gz | head -$n  | gawk '/^#/{next}/^track/{next;}/^fixedStep/{if(index($0,"start=0")>0){printf("%s %s start=0 %s\n",$1,$2,$4);jump=1;next;}}{if(jump==1){jump=0;next;}}{print}' | sed -e 's/chrom=/chrom=chr/'  | grep -v chrET_av >> tmp/WIGGLE$RG/$rg/any.$fr.wig
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
