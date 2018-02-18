#!bin/tcsh -f

set dd=$1
set ss=$2
 set srr=`echo $ss | gawk -F '___' '{print $1}'`
 set ff=`echo $ss | gawk -F '___' '{print $2}'`
 set p=`echo $ss | gawk -F '___' '{print $3}'`

echo "$srr $ff $p"
pushd $dd

  if ($p == Paired_end) then
    echo "fastq-dump --split-files --fasta 0 --gzip $srr"
    fastq-dump --split-files --fasta 0 --gzip $srr
  else
    echo "fastq-dump --fasta 0 --gzip $srr"
    fastq-dump --fasta 0 --gzip $srr
  endif

 
popd


 
