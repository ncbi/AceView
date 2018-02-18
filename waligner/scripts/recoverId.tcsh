#!bin/tcsh -ef

set run=$1

set mytmp=$TMPDIR/aceview.recoverId.$run.$$
echo $mytmp

mkdir $mytmp
mkdir $mytmp/ALIa
mkdir $mytmp/ALIa/$run
mkdir $mytmp/ALI
mkdir $mytmp/ALI/$run

phasea4:
 echo -n "Phase a4: $run recover the original identifiers "
 date

set nTarget=0

#######################################################################################
## ventilate

foreach lane (`gawk '{i=index($1,"/");if(run == substr($1,1,i-1))print $1;}' run=$run MetaDB/$MAGIC/LaneList`)
  gunzip -c tmp/COUNT/$lane.hits.gz | bin/bestali -ventilate 16 -id_file Fastc/$lane.id.gz -o $mytmp/ALIa/$lane
end

## regroup
foreach ii (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
  cat  $mytmp/ALIa/*.$ii.hits | sort -k 1,1 -k 2,2nr -k 8,8 | gzip -1 >  $mytmp/ALI/$run/g.$ii.hits.gz
end

#######################################################################################
## cleanup

foreach target ($targets)
  \rm $mytmp/PHITS_$target/*.hits.gz

  if (! -d tmp/ALI) mkdir tmp/ALI
  if (! -d tmp/ALI/$run) mkdir tmp/ALI/$run

  mv $mytmp/ALI/$run/* tmp/ALI/$run
end

touch tmp/ALI/$run/recoverId.done

# \rm -rf $mytmp

echo -n "done "
date

exit 0


##############################################################
