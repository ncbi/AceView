#!bin/tcsh -f

# do not use -e, you block if some file is missing

set rg=$1
set run=$2
set ln=$3

set mytmp=$TMPDIR/aceview.s5.countEdited.$$
mkdir $mytmp

touch  $mytmp/tata1 $mytmp/tata2
echo "ZZZZZ\nzzzzz" >  $mytmp/zzzzz

if ($rg == run) then

  gunzip -c tmp/NEWHITS_snp/$run/*.[mw].hits.gz | gawk '/^#/{next}{print}' | sed -e 's/gi|251831106|ref|NC_012920.1|/mito/g'      | sort > $mytmp/tata1
  gunzip -c  tmp/NEWHITS_snp/$run/*.w.multi.hits.gz  | gawk '/^#/{next}{print}' | sed -e 's/gi|251831106|ref|NC_012920.1|/mito/g' | sort > $mytmp/tata2
endif

cat  $mytmp/tata1  $mytmp/zzzzz  $mytmp/tata2 |  bin/snp -countEdited $ln -minAliPerCent 90 -run $run -o  tmp/NEWHITS_snp/$run/s590.multi 


\rm -rf  $mytmp


