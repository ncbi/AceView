#!bin/tcsh -f

set out=$1
set inFileList=$2

bin/bestali  -geneSupport2ace  -gzo -o $out -inFileList $inFileList

exit 0
