#!bin/tcsh -f

set lane=$1

bin/wiggle -ventilate -I BHIT -O BHIT -i tmp/PHITS_genome/$lane.hits.gz -o tmp/PHITS_genome/$lane -gzo

exit 0
