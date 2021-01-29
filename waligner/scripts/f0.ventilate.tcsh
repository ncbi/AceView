#!bin/tcsh -f

set lane=$1
# 2019-05-15 : surprisingly i notice that this code never worked in the MiniTest system
# may be i do not need it
exit 0
bin/wiggle -ventilate -I BHIT -O BHIT -i tmp/PHITS_genome/$lane.hits.gz -o tmp/PHITS_genome/$lane -gzo

exit 0
