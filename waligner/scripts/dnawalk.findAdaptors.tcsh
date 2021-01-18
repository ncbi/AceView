#!bin/tcsh -f

set run=$1
set read=$2

set r2=""
if ($read == 2) set r2='-read2'

echo " gunzip -c Fastc/$run/f*.fastc.gz  | bin/dnawalk -prefix -walk -wordLength 11 -vary 30 -run $run $r2 -o tmp/Adaptors/$run/r.$read"
gunzip -c Fastc/$run/f*.fastc.gz  | bin/dnawalk -prefix -walk -wordLength 11 -vary 30 -run $run $r2 -o tmp/Adaptors/$run/r.$read

exit 0

echo "# " > tmp/$MAGIC.prefix
foreach run (`cat MetaDB/$MAGIC/RunList`)
  set ff=tmp/Profiles/$run/f2.readsBeforeAlignment.1.txt
  if (! -e $ff) continue
  cat $ff | gawk -F '\t' '/^#/{next;}{k++;if(NF==11)for(i=7;i<=10;i++){if($i>80)printf("%s\t%d\t%s\t%s\n",run,k,substr("atgc",i-6,1),$i);}}' run=$run >>tmp/$MAGIC.prefix
end
