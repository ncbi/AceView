#!bin/tcsh -f 

set run=$1

set mcn=""
if ($?maxChromNameLength) then
  set mcn="-maxChromNameLength $maxChromNameLength"
endif

bin/geneelements -ventilate -run $run $mcn

exit 0
