#!bin/tcsh -ef
set target=$1
set manip=$2
set tissue=$3
set lane=$4


cat tmp/COUNT/$manip/$tissue/f.14.$target.list | gawk -F '\t' '/chondria/{p=$1;i=length(p);zz=0;if(substr(p,n-1,1)=="/")print}' >  tmp/Transplicing/$tissue.$lane.pairs

