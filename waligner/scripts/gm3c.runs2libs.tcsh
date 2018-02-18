#! bin/tcsh -f
set inFileList=$1
set out=$2

bin/bestali  -aceRuns2Lib  -gzo -o $out -inFileList $inFileList

