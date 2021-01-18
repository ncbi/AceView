#!bin/tcsh -ef

set lane=$1

# hack to back process the best hit files
# 

if ((-e  tmp/COUNT/$lane.hits.gz || -e  tmp/COUNT/$lane.hits.gz.mrnaOrder.gz) && ! -e tmp/COUNT/$lane.mrnaOrder.done) then

  if (! -e tmp/COUNT/$lane.hits.gz.mrnaOrder.gz) then
    mv tmp/COUNT/$lane.hits.gz tmp/COUNT/$lane.hits.gz.mrnaOrder.gz
  endif
  if (-e tmp/COUNT/$lane.too_many_hits.gz ) mv tmp/COUNT/$lane.too_many_hits.gz tmp/COUNT/$lane.too_many_hits.gz.ok
  gunzip -c  tmp/COUNT/$lane.hits.gz.mrnaOrder.gz | bin/bestali -o tmp/COUNT/$lane -exportBest -gzo
  \mv tmp/COUNT/$lane.too_many_hits.gz.ok tmp/COUNT/$lane.too_many_hits.gz
  touch tmp/COUNT/$lane.mrnaOrder.done
  \rm  tmp/COUNT/$lane.hits.gz.mrnaOrder.gz
endif

exit 0

# gunzip -c tmp/COUNT/$lane.hits.gz | bin/bestali -filter $species -maxHit 10 -exportVenn -o tmp/COUNT/$lane

# 
foreach lane (`cat MetaDB/$MAGIC/LaneList`)
  gunzip -c tmp/PHITS_*/$lane.hits.gz | sort -k 1,1 -k 2,2nr | bin/bestali -filter $species -maxHit $maxHit -countBest -exportBest -exportVenn -exportSuffix -errorProfile -aliProfile -exportMito -gzo -o tmp/COUNT/$lane
end
