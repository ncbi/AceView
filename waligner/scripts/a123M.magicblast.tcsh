#!/bin/tcsh

set phase=$1
set run=$2
set lane=$3

if (! -d TARGET/Targets) exit 1

if ($phase == makeDB) then
  cd TARGET/Targets
  if (! -d BLASTDB) mkdir BLASTDB
  cd BLASTDB

  foreach target ($RNAtargets  $DNAtargets)
    set ff=$species.$target
    if (-e ../$ff.fasta.gz && ! -e $ff.nsq) then
      gunzip -c ../$ff.fasta.gz > $ff
      makeblastdb -in $ff -dbtype nucl
      if ($target == genome) then
        dna2dna -decoy -i ../$ff.fasta.gz -O fasta > $species.gdecoy
        makeblastdb -in  $species.gdecoy -dbtype nucl
        \rm $species.gdecoy
      endif
      \rm $ff
    endif
  end 

  exit 0
endif

set informat=fasta
foreach run2 (`cat MetaDB/RunPairedList`)
  if ($run == $run2) then
    set informat=fastc
    break
  endif
end

foreach target ($RNAtargets)
  if (-e  tmp/MAGICBLAST/$lane.$target.hits.gz) continue
  if (-e  tmp/MAGICBLAST/$lane.$target.sorted_hits.gz) continue
  echo "time magicblast -db TARGET/Targets/BLASTDB/$species.$target -infmt $informat -outfmt tabular -gzo -out tmp/MAGICBLAST/$lane.$target.hits.gz -reftype transcriptome -query Fastc/$lane.fastc.gz -splice F and -limit_lookup F"
        time magicblast -db TARGET/Targets/BLASTDB/$species.$target -infmt $informat -outfmt tabular -gzo -out tmp/MAGICBLAST/$lane.$target.hits.gz -reftype transcriptome -query Fastc/$lane.fastc.gz -splice F and -limit_lookup F"
end

foreach target ($DNAtargets)
  if (-e  tmp/MAGICBLAST/$lane.$target.hits.gz) continue
  if (-e  tmp/MAGICBLAST/$lane.$target.sorted_hits.gz) continue
  echo "time magicblast -db TARGET/Targets/BLASTDB/$species.$target -infmt $informat -outfmt tabular -gzo -out tmp/MAGICBLAST/$lane.$target.hits.gz  -reftype genome -query Fastc/$lane.fastc.gz"
        time magicblast -db TARGET/Targets/BLASTDB/$species.$target -infmt $informat -outfmt tabular -gzo -out tmp/MAGICBLAST/$lane.$target.hits.gz  -reftype genome -query Fastc/$lane.fastc.gz
end

touch  tmp/MAGICBLAST/$lane.done

## a kind of bestali
foreach target ($DNAtargets $RNAtargets) 
  if (-e  tmp/MAGICBLAST/$lane.$target.hits.gz && ! -e tmp/MAGICBLAST/$lane.$target.sorted_hits.gz) then 
    if (! -e  tmp/MAGICBLAST/caption) then
      gunzip -c  tmp/MAGICBLAST/$lane.$target.hits.gz  | head -3 | tail -1 >  tmp/MAGICBLAST/caption
    endif
    gunzip -c  tmp/MAGICBLAST/$lane.$target.hits.gz  | gawk -F '\t' '/^#/{next;}{printf("%s\t%09d\t%s\t", $1,$25,target); print;}' | sort | gzip >  tmp/MAGICBLAST/$lane.$target.sorted_hits.gz
  endif
end

zcat  tmp/MAGICBLAST/$lane.*.sorted_hits.gz | sort | gawk -F '\t' '{if($1 == p && $2 < s)next ; p = $1;  s = $2; print}' | gzip >  tmp/MAGICBLAST/$lane.bestali.gz

exit 0

