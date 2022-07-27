#!bin/tcsh -ef

set run=$1
set chrom=$2
set uu=$3
set fr=$4
set limit=$5
set mask="$6"
set mainTarget=$7

# genome gene_and_new_exon new_and_known_exon known_exon RefSeq
echo "sponge_chom.tcsh $run $chrom $uu $fr $limit $mask"


# prepare the frns wiggles if not available, now always done in wg2a for runs and wg2b for groups 
set ww1=toto
set ww2=toto
if (-e  tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz || -e  tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.pp.BF.gz) then
  set ww1=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz 
  set ww2=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.pp.BF.gz
else if (-e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.pp.BF.gz) then
  set ww1=tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz
  set ww2=tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.pp.BF.gz
endif

if (! -e $ww1  && ! -e $ww2) then
  echo "tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.[u|pp].BF.gz not found"
  exit 0
endif
if (! -e $ww1) then
  set ww=$ww2
else if (! -e $ww2) then
  set ww=$ww1
else
  set ww="$ww1,$ww2"
endif

if ($mainTarget != introns) then        
  set toto=tmp/SPONGE/$run/Total.$chrom.$uu.$fr.$limit
  if (! -e $toto) then
    echo Total > $toto
    echo "... bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto"
              bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto
  endif

  grep Total $toto |   sed -e "s/^Total/$chrom\t$fr/" | gawk -F '\t' '{r=run;}/^#/{r="#Run\tChrom"}{printf("%s\t",r);print;}' run=$run > $toto.txt
  \rm $toto
endif

set fr=ns
foreach mainTarget ($Etargets)
  foreach GM (gene mrna)
    set geneMask=tmp/METADATA/$mainTarget.$fr.$GM.sponge
    ls -ls $geneMask
    if (-e $geneMask && $limit == 1) then
      set toto=tmp/SPONGE/$run/$mainTarget.$chrom.$uu.$fr.$limit
      if (-e $toto)  \rm $toto
      set toto=tmp/SPONGE/$run/$mainTarget.$GM.$chrom.$uu.$fr.$limit
      if (-e $toto)  \rm $toto
      set toto=tmp/SPONGE/$run/$mainTarget.$GM.v2.$chrom.$uu.$fr.$limit
      if (! -e $toto) then
        echo "bin/geneelements -sponge $limit -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto"
              bin/geneelements -sponge $limit -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto
      endif
    endif
  end
end



