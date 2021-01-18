#!bin/tcsh -f 

set run=$1

set mcn=""
# if you open this 'hacking' parameter here. then add a line in scripts/LIMIT to set it to zero as a default
if ($?maxChromNameLengthZZZ == 1) then
  set mcn="-maxChromNameLength $maxChromNameLength"
endif

echo "bin/geneelements -ventilate -run $run $mcn"
      bin/geneelements -ventilate -run $run $mcn

exit 0
