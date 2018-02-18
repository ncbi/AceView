#!bin/tcsh -f

set lane=$1

bin/dnastitch -small -i Fastc/$lane.fastc.gz -o tmp/ShortFragments/$lane

exit 0
