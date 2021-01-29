#!/bin/tcsh

set phase=$1
set run=$2
set lane=$3

echo "a123MH.matchhits.tcsh $*"
if (! -d TARGET/Targets) exit 1

if (1 && $phase == makeDB) then
  pushd TARGET
  if (! -d MATCHHITSDB) mkdir MATCHHITSDB

  foreach target ($RNAtargets $DNAtargets )
    set ff=$species.$target
    if (-e Targets/$ff.fasta.gz && ! -e MATCHHITSDB/$ff.idx) then
      gunzip -c Targets/$ff.fasta.gz | gawk '/^>/{print;next;}{gsub("-","n",$0);print;}' > MATCHHITSDB/$ff.fasta
      pushd MATCHHITSDB
        matchhits -A mkdb -i $ff.fasta -o $ff --mkidx
        if ($target == genome) then
          dna2dna -decoy -i $ff.fasta -O fasta | gawk '/^>/{print;next;}{gsub("-","n",$0);print;}' > $species.gdecoy.fasta
          matchhits -A mkdb -i $species.gdecoy.fasta -o $species.gdecoy --mkidx
          \rm $species.gdecoy
        endif
      popd
    endif
  end 

  popd
  exit 0
endif

set inok=0
# make fasta files from the fastc files
if (1) then
  foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
    set inok=0
    if ($run == $run2) then
      if (-e Fastc/$lane.fastc.gz && ! -e Fastc/$lane'_1'.fasta) then
        echo "dna2dna -I fastc -O fasta  -splitPairs -i Fastc/$lane.fastc.gz -o Fastc/$lane"
              dna2dna -I fastc -O fasta  -splitPairs -i Fastc/$lane.fastc.gz -o Fastc/$lane
      endif
      set infile="-1 Fastc/$lane"_1".fasta -2 Fastc/$lane"_2".fasta "
      set inok=1
      break
    endif
  end
    if ($inok == 0) then
      if (-e Fastc/$lane.fastc.gz && ! -e Fastc/$lane.fasta.gz) then
        echo "dna2dna -I fastc -O fasta  -i Fastc/$lane.fastc.gz -o Fastc/$lane"
              dna2dna -I fastc -O fasta  -i Fastc/$lane.fastc.gz -o Fastc/$lane
      endif
      ls -ls  Fastc/$lane.fasta Fastc/$lane'_'?.fasta
      set infile="-i Fastc/$lane.fasta "
      set inok=1

    endif

endif

echo "Fasta contruction done, starting alignments"

foreach target ($RNAtargets)
  set subject="--db TARGET/MATCHHITSDB/$species.$target"

  if (-e  tmp/MATCHHITS/$lane.$target.sam) continue
  if (-e  tmp/MATCHHITS/$lane.$target.sam.gz) continue
  if (-e  tmp/MATCHHITS/$lane.$target.sorted_hits.gz) continue
  echo "time matchhits $subject $infile --no-splice -t 1 -o tmp/MATCHHITS/$lane.$target.sam"
        time matchhits $subject $infile --no-splice -t 1 -o tmp/MATCHHITS/$lane.$target.sam
  gzip tmp/MATCHHITS/$lane.$target.sam
end

foreach target ($DNAtargets)
  set subject="--db TARGET/MATCHHITSDB/$species.$target"

  if (-e  tmp/MATCHHITS/$lane.$target.sam) continue
  if (-e  tmp/MATCHHITS/$lane.$target.sam.gz) continue
  if (-e  tmp/MATCHHITS/$lane.$target.sorted_hits.gz) continue
  echo "time matchhits $subject $infile -t 1 -o tmp/MATCHHITS/$lane.$target.sam"
        time matchhits $subject $infile -t 1 -o tmp/MATCHHITS/$lane.$target.sam

# ATTENTION, we assume that the run is stranded plus
# because otherwise we of not know the strand of the intron
if (-e tmp/MATCHHITS/$lane.$target.sam) then
  if (! -e  tmp/MATCHHITS/$lane.$target.introns.gz) then 
    cat tmp/MATCHHITS/$lane.$target.sam | gawk -F '\t' -f scripts/sam2introns.awk | gzip > tmp/MATCHHITS/$lane.$target.introns.gz
  endif
  if (! -e  tmp/MATCHHITS/$lane.$target.aliqc.tsv && -e ~/SRC/CODE/PYTHON_2.7) then
    setenv vdir ~/SRC/CODE/PYTHON_2.7
    bash
    .  $vdir/bin/activate
    python scripts/aliqc.py --SAM -r $lane -i tmp/MATCHHITS/$lane.$target.sam  -f TARGET/Targets/$species.$target.fasta.gz -o  tmp/MATCHHITS/$lane.$target
    exit
  endif
  gzip tmp/MATCHHITS/$lane.$target.sam
endif

touch  tmp/MATCHHITS/$lane.done
exit 0

