#!bin/tcsh -f

set compare=$1
set bigGroup=$2
set bigCover=$3
set controlGroup=$4
set phase=$5

echo "scripts/cnv.compare_coverons.tsch $1 $2 $3 $4 $5"

if ($phase == getCoverons) then
  # given a big group
  # extract the geometry of the coverons at large threshold

  # extract a sponge file

  set ok=0
  echo " tmp/WIGGLEGROUP/$bigGroup/coverome.$bigCover.peaks"

  foreach chrom ($chromSetAll)
    # compute the cumulated u + nu peaks rather than using tmp/WIGGLEGROUP/$bigGroup/$chrom/coverome.$bigCover.peaks
    if (! -e  tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons) then
      gunzip -c tmp/WIGGLEGROUP/$bigGroup/$chrom/R.chrom.frns.u.BF.gz tmp/WIGGLEGROUP/$bigGroup/$chrom/R.chrom.frns.nu.BF.gz | bin/wiggle -peaks  -O BF -I BF -minCover $bigCover -o tmp/WIGGLE_CNV/$chrom/toto.$$
      cat  tmp/WIGGLE_CNV/$chrom/toto.$$.peaks  | gawk -F '\t' '/^#/{next;}{chrom=$1; a1=$2;a2=$3; ln = a2 - a1+1;n1=int(1+ln/2000); da=ln/n1;for(i=0;i<n1;i++){b1=a1+i*da-1;b2=b1+da;nn++;nam= "coveron." chrom "." nn ; printf ("%s\t1\t%s\t%d\t%d\t%s\n",nam,chrom,b1,b2,nam);}}' > tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons
      \rm tmp/WIGGLE_CNV/$chrom/toto.$$.*
    endif
    cat tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons >>  tmp/WIGGLE_CNV/$bigGroup.coverons
  end
  goto done
endif

foreach uu (u nu)
 foreach chrom ($chromSetAll)
  foreach run (`cat tmp/WIGGLE_CNV/cnv.runs1  tmp/WIGGLE_CNV/cnv.runs2 | sort -u`) 
     if  (-e  tmp/WIGGLE_CNV/$chrom/$run.cnv.$uu.sponge) continue

     set f=xx
     set fu=xx
     set fnu=xx
     set funu=xx
     if (-e tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.u.BF.gz) set fu=tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.u.BF.gz
     if (-e tmp/WIGGLEGROUP/$run/$chrom/R.chrom.frns.u.BF.gz) set fu=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.frns.u.BF.gz
     if (-e tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.nu.BF.gz) set fnu=tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.nu.BF.gz
     if (-e tmp/WIGGLEGROUP/$run/$chrom/R.chrom.frns.nu.BF.gz) set fnu=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.frns.nu.BF.gz

     if ($uu == u) set f = $fu
     if ($uu == nu) set f = $fnu

     # if ($fu == xx && $fnu == xx ) continue  # these 3 lines plus the if (0) case were written for the unu mixed case 
     # if ($fu == xx) set f=$fnu
     # if ($fnu == xx) set f=$fu
     # if (1) then
     #    set funu=tmp/WIGGLE_CNV/$chrom/$run.unu
     #    gunzip -c $fu $fnu  | bin/wiggle -O BF -I BF -o $funu -gzo
     #    set f=$funu.BF.gz
     #  endif
     if ($f == xx) continue
     echo "scripts/cnv.compare_coverons_do.tcsh tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons  $f  tmp/WIGGLE_CNV/$chrom/$run.cnv.$uu.sponge"
     scripts/submit  tmp/WIGGLE_CNV/$chrom/$run "scripts/cnv.compare_coverons_do.tcsh tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons  $f  tmp/WIGGLE_CNV/$chrom/$run.cnv.$uu.sponge"
     @ ok = $ok + 1
  end
 end
end

if ($ok > 0) then
  echo "# submitted $ok sponge files, please wait and rerun MAGIC cnv "
  goto done
endif

# measure for the .u. coverons, their coverage in the BigGroup u and nu wiggles, to obtain the unicity index of each coveron index=u+100/u+nu+100 using column 11 = base coverage
if (-e  tmp/WIGGLE_CNV/$bigGroup.u2u.sponge) \rm tmp/WIGGLE_CNV/$bigGroup.u2u.sponge tmp/WIGGLE_CNV/$bigGroup.u2nu.sponge
foreach chrom ($chromSetAll)
  if (-e tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons && -e tmp/WIGGLEGROUP/$bigGroup/$chrom/R.chrom.frns.nu.BF.gz && ! -e tmp/WIGGLE_CNV/$chrom/$bigGroup.u2nu.sponge) then
    bin/geneelements -sponge 10 -spongeFile  tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons -wiggle tmp/WIGGLEGROUP/$bigGroup/$chrom/R.chrom.frns.nu.BF.gz  >>  tmp/WIGGLE_CNV/$chrom/$bigGroup.u2nu.sponge
    bin/geneelements -sponge 10 -spongeFile  tmp/WIGGLE_CNV/$chrom/$bigGroup.coverons -wiggle tmp/WIGGLEGROUP/$bigGroup/$chrom/R.chrom.frns.u.BF.gz  >>  tmp/WIGGLE_CNV/$chrom/$bigGroup.u2u.sponge
  endif
  if (-e  tmp/WIGGLE_CNV/$chrom/$bigGroup.u2nu.sponge && -e tmp/WIGGLE_CNV/$chrom/$bigGroup.u2u.sponge) then
    cat tmp/WIGGLE_CNV/$chrom/$bigGroup.u2nu.sponge >>  tmp/WIGGLE_CNV/$bigGroup.u2nu.sponge
    cat tmp/WIGGLE_CNV/$chrom/$bigGroup.u2u.sponge >>  tmp/WIGGLE_CNV/$bigGroup.u2u.sponge
  endif
end

set u2u=""
if (-e tmp/WIGGLE_CNV/$bigGroup.u2nu.sponge  && -e  tmp/WIGGLE_CNV/$bigGroup.u2u.sponge) then
  set u2u="--u2u tmp/WIGGLE_CNV/$bigGroup.u2u.sponge --u2nu tmp/WIGGLE_CNV/$bigGroup.u2nu.sponge"
endif

set gg=""

if (-e  tmp/METADATA/av.ns.gene.sponge) then
  set gg="--geneSponge tmp/METADATA/av.ns.gene.sponge"
else 
  if (-e  tmp/METADATA/RefSeq.ns.gene.sponge) then
    set gg="--geneSponge tmp/METADATA/RefSeq.ns.gene.sponge"
  endif
endif

foreach uu (u nu)
  echo "bin/cnv --runs MetaDB/$MAGIC/RunListSorted --runsAll tmp/WIGGLE_CNV/cnv.runs1 --runsControl tmp/WIGGLE_CNV/cnv.runs2 --gtitle MetaDB/$MAGIC/gtitle.txt  $gg --coveronSponge tmp/WIGGLE_CNV/$bigGroup.coverons --coverage tmp/WIGGLE_CNV $u2u --$uu -o  tmp/WIGGLE_CNV/$MAGIC.$uu   --db MetaDB --project $MAGIC"
        bin/cnv --runs MetaDB/$MAGIC/RunListSorted --runsAll tmp/WIGGLE_CNV/cnv.runs1 --runsControl tmp/WIGGLE_CNV/cnv.runs2 --gtitle MetaDB/$MAGIC/gtitle.txt  $gg --coveronSponge tmp/WIGGLE_CNV/$bigGroup.coverons --coverage tmp/WIGGLE_CNV $u2u --$uu -o  tmp/WIGGLE_CNV/$MAGIC.$uu   --db `pwd`/MetaDB --project $MAGIC
end

if (! -d  RESULTS/Coverage_and_exons) mkdir  RESULTS/Coverage_and_exons
foreach chrom ( $chromSetAll )
   cat tmp/WIGGLE_CNV/$MAGIC.u.cnv_compare_coverons.scaled.txt | gawk -F '\t' '/^#/{print;next;}{if($2==chrom)print}' chrom=$chrom > RESULTS/Coverage_and_exons/$MAGIC.u.cnv_compare_coverons.scaled.$chrom.txt
   cat tmp/WIGGLE_CNV/$MAGIC.u.cnv_compare_coverons.raw.txt | gawk -F '\t' '/^#/{print;next;}{if($2==chrom)print}' chrom=$chrom > RESULTS/Coverage_and_exons/$MAGIC.u.cnv_compare_coverons.scaled.$chrom.raw.txt
end
\cp tmp/WIGGLE_CNV/$MAGIC.nu.cnv_compare_coverons.*.txt    RESULTS/Coverage_and_exons
\cp tmp/WIGGLE_CNV/$MAGIC.u.cnv_compare_coverons.*lost.txt  RESULTS/Coverage_and_exons


foreach uu (u nu)
 foreach sr (scaled raw)
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.$sr.txt | head -100 | gawk -F '\t' '/^#/{print}' > RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.$sr.txt | gawk -F '\t' '/^#/{next;}{if (0+$2 >0)print}' | sort -k 2,2n  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.$sr.txt | gawk -F '\t' '/^#/{next;}{if ($2 == "X")print}' | sort  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.$sr.txt | gawk -F '\t' '/^#/{next;}{if ($2 == "Y")print}' | sort  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.$sr.txt

  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt | head -100 | gawk -F '\t' '/^#/{print}' > RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt | gawk -F '\t' '/^#/{next;}{if (0+$2 >0)print}' | sort -k 2,2n  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt | gawk -F '\t' '/^#/{next;}{if ($2 == "X")print}' | sort  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt
  cat tmp/WIGGLE_CNV/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt | gawk -F '\t' '/^#/{next;}{if ($2 == "Y")print}' | sort  -k 3,3n | gawk -F '\t' '{c=$2;if(c != old)for(i=0;i<300;i++)printf("\n");old=c;print;}' >> RESULTS/Coverage_and_exons/$MAGIC.$uu.cnv_compare_genes.lost.$sr.txt
end

\cp tmp/WIGGLE_CNV/$MAGIC.*.cnv_compare_*_histo.*txt RESULTS/Coverage_and_exons/


done:
  exit 0

