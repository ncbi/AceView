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
set ww=toto
if (-e  tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz) then
  set ww=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz
else if (-e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz) then
  set ww=tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz
endif

if ($ww == toto) then
  echo "tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz not found"
  exit 0
endif
            
  set toto=tmp/SPONGE/$run/Total.$chrom.$uu.$fr.$limit
  if (! -e $toto) then
    echo Total > $toto
    echo "... bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto"
              bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto
  endif

  grep Total $toto |   sed -e "s/^Total/$chrom\t$fr/" | gawk -F '\t' '{r=run;}/^#/{r="#Run\tChrom"}{printf("%s\t",r);print;}' run=$run > $toto.txt
\rm $toto

if ($fr == frns) set fr=ns
foreach mainTarget ($Etargets)
  set geneMask=tmp/METADATA/$mainTarget.$fr.gene.sponge
  ls -ls $geneMask
  if (-e $geneMask && $limit == 1) then
    set toto=tmp/SPONGE/$run/$mainTarget.$chrom.$uu.$fr.$limit
    if (! -e $toto) then
      echo "bin/geneelements -sponge $limit -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto"
            bin/geneelements -sponge $limit -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto
    endif
  endif
end



