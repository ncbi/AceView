#!bin/tcsh -f

set run=$1
set chrom=$2

echo -n "f2.flagTranscritEnds $run $chrom start : "
date

set wig=tmp/WIGGLERUN/$run/$chrom/R.chrom.u

set ok=1
foreach fr (ELF ELR ERF ERR)
  if (! -e $wig/R.chrom.u.$fr.BF.gz) then
    echo "Missing file   $wig/R.chrom.u.ELF.BF.gz"
    set ok=0
  endif
end

if ($ok == 1) then
  bin/geneelements -run $run  -flagTranscriptEnds -wiggle wig  -minimalSupport $min -o tmp/EHITS.$MAGIC/$chrom/f2e.$run.txt
endif

echo -n "done : "
date




