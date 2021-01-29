#!/bin/tcsh
set run=$1
foreach lane (`cat Fastc/$run/laneList`)
  set ff=Fastc/$lane.fastc.gz
  mv $ff $ff.old.gz
  dna2dna -i $ff.old.gz -I fastc -O fastc -UMI 33 -o Fastc/$lane -gzo -gzi -count
end
