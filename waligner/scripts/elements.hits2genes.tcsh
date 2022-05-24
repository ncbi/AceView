#!bin/tcsh -ef

set GM=$1
set target=$2
set run=$3
set COUNT=$4
set lane=$5
set stranded=$6
set SUBSAMPLING=$7

set uGeneSupport=""
if ($8 != "") set uGeneSupport=" -unique_gene_support $8"

set test=2
if ($target == mito || $target == SpikeIn) set test=0
if ($GM == GENE) set test=1
if ($GM == MRNA) set test=1
if ($GM == MRNAH) set test=1
if ($test == 2) then
  echo "elements.hits2gene first param=$GM, should be GENE or MRNA"
  exit 1
endif
echo hello1
set type=g

set ss=""
if ($stranded == 1) set ss='-stranded'
if ($stranded == -1) set ss='-antistranded'

source scripts/target2target_class.txt


 #  -pureNsStrand 1 => reject the reads aligning once on + strand of a gene, once on - strand of another gene
 #                2 => accept the special case of 2 genes on opposite strands, those may be remapped in the non stranded genomic wiggle
set toto=tmp/$COUNT/$lane.hits.gz
set remap=""
if ($target == mito || $target == chloro || $target == SpikeIn) then
  set toto=tmp/$COUNT/$lane.mito.gz
  if ($GM != GENE) continue
else
  if (-e tmp/METADATA/mrnaRemap.gz) then
    set remap=" -mrnaRemap tmp/METADATA/mrnaRemap.gz "
  endif
endif

set gm2='-geneSupport -intergenicSupport'
if ($GM == MRNA) set gm2='-mrnaSupport'
if ($GM == MRNAH) set gm2='-mrnaSupport -hierarchic'

set st="-strategy $Strategy"

set pair=""
foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
  if ($run == $run2) then
    set pair="-pair 500"
  endif
end

set splitM=""
if ($GM == GENE && -e TARGET/GENES/$species.$target.split_mrnas.txt) then
  set splitM="-split_mRNAs TARGET/GENES/$species.$target.split_mrnas.txt"
endif

set isIlm=`cat  MetaDB/$MAGIC/runs.ace | gawk '/^Run /{gsub(/\"/,"",$2);ok=0;if($2==run)ok=1;}/^Illumina/{if(ok==1)okk=1;}END{print 0+okk;}' run=$run`
if ($isIlm == 1) then
  set maxWigErr=3
  set maxWigErrRate=3
else
  set maxWigErr=0
  set maxWigErrRate=0
endif

   set m8kb=""
   if (-e tmp/METADATA/$target.selected8kbTranscriptList.txt && ($GM == MRNA || $GM == MRNAH)) set m8kb="-selected8kbList tmp/METADATA/$target.selected8kbTranscriptList.txt"
   set m5kb=""
   if (-e tmp/METADATA/$target.selected5kbTranscriptList.txt && ($GM == MRNA || $GM == MRNAH)) set m5kb="-selected5kbList tmp/METADATA/$target.selected5kbTranscriptList.txt"

if (-e $toto) then
        echo "bin/bestali -i $toto $gm2 -run $run -target_class $target_class $splitM -pureNsStrand 1 -maxErr $maxWigErr -maxErrRate $maxWigErrRate $pair $ss $st $uGeneSupport $m5kb $m8kb -gzo -o  tmp/GENELANES/$lane.$target.$GM $remap -subsampling $SUBSAMPLING"
              bin/bestali -i $toto $gm2 -run $run -target_class $target_class $splitM -pureNsStrand 1 -maxErr $maxWigErr -maxErrRate $maxWigErrRate $pair $ss $st $uGeneSupport $m5kb $m8kb -gzo -o  tmp/GENELANES/$lane.$target.$GM $remap  -subsampling $SUBSAMPLING
endif

exit 0

foreach run (`cat MetaDB/$MAGIC/RunList`)
  source scripts/mkDir  COUNT2 $run
  foreach lane (`cat Fastc/$run/LaneList`)
     scripts/submit  tmp/COUNT2/$lane "bin/bestali -i tmp/COUNT/$lane.hits.gz -o  tmp/COUNT2/$lane -pair 300  -strategy $Strategy -exportBest -gzo"
  end
end

grep Status tmp/COUNT2/*/*.err | grep -v '= 0'
foreach lane (`cat MetaDB/$MAGIC/LaneList`)
  mv   tmp/COUNT2/$lane.hits.gz tmp/COUNT/$lane.hits.gz
end




