#!bin/tcsh -ef
set manip=$1
set tissue=$2

set manip2=$manip
if($manip == any) set manip2='*'

if (! -e  tmp/mito_wiggle/pA/$manip.$tissue.double) then
         gunzip -c  tmp/PHITS_mito/$manip2/$tissue/*.hits.gz  | gawk -F '\t' '{if($5 ~ /p[AT]/){i=index($1,"#");n=0+substr($1,i+1);s=$7;a=$8;ns[a]++;if($5 == "pT"){if(s=="Forward"){nrs[a]++;nrt[a]+=n;}else {nfs[a]++;nft[a]+=n;}} if($5 == "pA"){if(s=="Forward"){nfs[a]++;nft[a]+=n;}else {nrs[a]++;nrt[a]+=n;}}}}END{for (a in ns)printf("%d\t%d\t%d\t%d\t%d\n",a,nft[a],nrt[a],nfs[a],nrs[a]);}' | sort -k 1n >  tmp/mito_wiggle/pA/$manip.$tissue.double
endif

          cat tmp/mito_wiggle/pA/$manip.$tissue.double | gawk '{if($2>0)printf("%d\t%d\n",$1,$2);}' | bin/wiggle -I BV -O BV -trackName "$manip.$tissue.f" -out_step 1  > tmp/mito_wiggle/pA/$manip.$tissue.tag.f.bv
          cat tmp/mito_wiggle/pA/$manip.$tissue.double | gawk '{if($3>0)printf("%d\t%d\n",$1,$3);}' | bin/wiggle -I BV -O BV -trackName "$manip.$tissue.r" -out_step 1   > tmp/mito_wiggle/pA/$manip.$tissue.tag.r.bv

          cat tmp/mito_wiggle/pA/$manip.$tissue.double | gawk '{if($4>0)printf("%d\t%d\n",$1,$4);}' | bin/wiggle -I BV -O BV -trackName "$manip.$tissue.f" -out_step 1  > tmp/mito_wiggle/pA/$manip.$tissue.seq.f.bv
          cat tmp/mito_wiggle/pA/$manip.$tissue.double | gawk '{if($5>0)printf("%d\t%d\n",$1,$5);}' | bin/wiggle -I BV -O BV -trackName "$manip.$tissue.r" -out_step 1   > tmp/mito_wiggle/pA/$manip.$tissue.seq.r.bv

