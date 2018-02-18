#!bin/tcsh -ef

set target=$1
set manip=$2
set tissue=$3
set lane=$4

    cat tmp/UHITS_$target/$manip.$tissue.$lane.$target.u.hits | gawk -F '\t' '{ng=$1;p=$2;score=$3;nn[p]++;if(nn[p]==1)gg[ng]++;}END{for(i=1;i<=10;i++)printf("%d\t",gg[i]);printf("\n");}' > tmp/MULTIPLICITY/LANES/$manip.$tissue.$lane.$target.mult


