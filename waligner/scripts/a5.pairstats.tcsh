#!bin/tcsh -f

set run=$1

## run the statistics
scripts/c2.alistats.tcsh $run COUNT2

## construct the pair length profile

cat tmp/COUNT2/$run/*.pairStats | gawk '/^#/{next;}/^$/{next}{t=$1;i=t2i[t];if(i<1){imax++;i=imax;t2i[t]=i;i2t[i]=t}nam[i]=$1;}/^Orphan/{if($2 == "Any")nn[i]+=$3;next;}{nn[i]+=$2;}END{for(i=1;i<=imax;i++) {if(nn[i]>0) printf("%s\t%s\t%s\n", r,nam[i],nn[i]) ;}}' r=$run > tmp/COUNT2/$run/runPairStats.txt
# exclude the last point, which is the remnant of the distribution
cat tmp/COUNT2/$run/*.pairStats | gawk '/^#Number/{for(i=2;i<NF;i++)n[i]+=$i;if(NF > imax)imax=NF;}END{printf("Insert_length_in_pairs\t%s",r);for(i=2;i<=imax;i++)printf("\t%d",n[i]);printf("\n");}' r=$run >  tmp/COUNT2/$run/runPairHisto.txt

