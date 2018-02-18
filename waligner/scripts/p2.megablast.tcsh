#!bin/tcsh -f

set run=$1
set lane=$2

gunzip -c tmp/Unaligned/$lane.fastc.gz | gawk '/^>/{split($1,aa,"#");if(aa[2]>10){ok=1;print;}next;}{if(ok==1)print;ok=0;}' > tmp/Unaligned/$lane.unmapped.fastc.short
gunzip -c tmp/Unaligned/$lane.fastc.gz | gawk '/^>/{split($1,aa,"#");if(aa[2]>5 && aa[2]<=10){ok=1;print;}next;}{if(ok==1)print;ok=0;}' | head -1000 >> tmp/Unaligned/$lane.unmapped.fastc.short
gunzip -c tmp/Unaligned/$lane.fastc.gz | gawk '/^>/{n++;split($1,aa,"#");if(n%10==0 && aa[2]<5 ){ok=1;print;}next;}{if(ok==1)print;ok=0;}' | head -1000 >> tmp/Unaligned/$lane.unmapped.fastc.short


if (! -e tmp/Unaligned/$lane.megablast.out) then
  megablast -d nt -i  tmp/Unaligned/$lane.unmapped.fastc.short -m 9  -o  tmp/Unaligned/$lane.megablast.out -D 0 -W 32 -f T -p 100
endif

cat tmp/Unaligned/$lane.megablast.out | sed -e "s/'/ /g" | gawk -F '|' '{if ($1==" gi"){gi=$2;split($5,aa," ");r=aa[2];zz[gi,r]++;if(zz[gi,r]==1)nn[gi]++;}}END{for (k in nn)printf("%s\t%d\n",k,nn[k]);}' | sort -u | sort -k 2nr >  tmp/Unaligned/$lane.gi.list

exit 0

