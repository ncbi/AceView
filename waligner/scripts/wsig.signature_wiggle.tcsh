#!bin/tcsh -f

set RG=$1
set target=$2
set run=$3

source scripts/target2target_class.txt
set err=0

if ($1 == RUN) then

  foreach lane (`cat Fastc/$run/LaneList`)
    if (! -e tmp/COUNT/$lane.hits.gz) then
      exit 1 
    endif
  end

  foreach lane (`cat Fastc/$run/LaneList`)
    if (-e tmp/COUNT/$lane.sig_hits.gz) then
      gunzip -c  tmp/COUNT/$lane.sig_hits.gz > tmp/WIGGLESIG/$lane.sig_hits
    else
      bin/bestali -i tmp/COUNT/$lane.hits.gz -sigTargets  TARGET/Targets/$species.$target.signature_transcripts_list.txt -o tmp/WIGGLESIG/$lane -target_class ET_av -pair 500
    endif
  end

  cat tmp/WIGGLESIG/$run/*.sig_hits | bin/wiggle -I BHIT -O BV -out_step 10 -gzo -o tmp/WIGGLESIG/$run/sig

  \rm   tmp/WIGGLESIG/$run/*.too_many_hitsHUM
  \rm   tmp/WIGGLESIG/$run/*.geneLinks
  \rm   tmp/WIGGLESIG/$run/*.pairStats
  # \rm  tmp/WIGGLESIG/$run/*.sig_hits
  touch tmp/WIGGLESIG/$run/wsig1.done 

endif

if ($1 == GROUP) then 

  set group=$run
  set err=0
  echo ' ' > tmp/WIGGLESIG/$group/_BV
  foreach run (`cat  MetaDB/$MAGIC/g2r | gawk '{if($1==g)print $2}' g=$group`)
    if (-e  tmp/WIGGLESIG/$run/sig.BV.gz) then
      gunzip -c  tmp/WIGGLESIG/$run/sig.BV.gz >> tmp/WIGGLESIG/$group/_BV
    else
      set err=1
    endif
  end
  if ($err == 0) then
    bin/wiggle -i tmp/WIGGLESIG/$group/_BV  -I BV -O BV  -out_step 10 -gzo  -o tmp/WIGGLESIG/$group/sig
  endif
  \rm  tmp/WIGGLESIG/$group/_BV

endif

exit $err
