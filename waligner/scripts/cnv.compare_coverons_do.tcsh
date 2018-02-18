#!bin/tcsh -f

set spongeFile=$1
set fuu=$2
set out=$3

bin/geneelements -sponge 10 -spongeFile $spongeFile -wiggle $fuu > $out

exit 0
