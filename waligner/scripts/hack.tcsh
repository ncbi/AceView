#!/bin/tcsh -f

set lane=$1
zcat tmp/COUNT/$lane.best_score.gz ZZZZZ.gz  tmp/BACTERIA/$lane.hits.gz | gawk -F '\t'  '/^ZZZZZ/{zz++;next;}{if(zz<1){s[$1]=$2;next;}}{if($2 > 0+s[$1])print;}'| gzip  > tmp/BACTERIA/$lane.clean_hits.gz

