#!bin/tcsh -f

set run=$1
set lane=$2

foreach target ($targets)
  source scripts/target2target_class.txt
  echo " bin/bestali -i tmp/COUNT/$lane.hits.gz -target_class $target_class --run $lane -o  tmp/COUNT/$lane.seqc.$target"
  bin/bestali -i tmp/COUNT/$lane.hits.gz -target_class $target_class  -run $lane.$target -o  tmp/COUNT/$lane.seqc.$target --seqc
end

touch tmp/COUNT/$lane.seqc_by_target.done
