#!bin/tcsh -ef

set manip=$1
set tissue=$2
set lane=$3

cat Fasta/$manip/$tissue.$lane.*fa* | gawk '/^[@>]/{print substr($1,2);}' >! tmp/NOL3/totom2f.$manip.$tissue.$lane
echo ZZZZZ >> tmp/NOL3/totom2f.$manip.$tissue.$lane
cat tmp/NOL3/$manip.$tissue.$lane.list | gawk -F '\t' '{print $2}' | sort -u  >> tmp/NOL3/totom2f.$manip.$tissue.$lane
cat   tmp/NOL3/totom2f.$manip.$tissue.$lane | gawk -F '\t' '/^ZZZZZ/{z=1;next;}{if(z==0){p[$1]=1;next;}if(p[$1]==1)p[$1]=2;}END{for(k in p) if(p[k]==2)print k}' >! tmp/NOL3/totom2f.$manip.$tissue.$lane.list
echo ZZZZZ >>  tmp/NOL3/totom2f.$manip.$tissue.$lane.list
cat  Fasta/$manip/$tissue.$lane.*fa* >>   tmp/NOL3/totom2f.$manip.$tissue.$lane.list
cat  tmp/NOL3/totom2f.$manip.$tissue.$lane.list  | gawk -F '\t' '/^ZZZZZ/{z=1;next;}/^[>@]/{s=substr($1,2);ok=p[s];if(ok==1)printf(">%s\n",s);next;}{if(z==0){p[$1]=1;next;}if(ok==1)print;ok=0;next;}' >!  tmp/NOL3/$manip.$tissue.$lane.fasta
\rm  tmp/NOL3/totom2f.$manip.$tissue.$lane.list tmp/NOL3/totom2f.$manip.$tissue.$lane

