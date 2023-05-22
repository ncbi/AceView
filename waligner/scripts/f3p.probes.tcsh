#!bin/tcsh -f
set chrom=$1
setenv ici `pwd`

set ff=$MAGIC.probes.hits
if (! -e $ff.gz) then
  goto done
endif

set dd=tmp/X.$MAGIC/XH$chrom
if (! -d $dd) then
  echo "missing directory $dd"
  goto done
endif

if (! -d $dd/database) then
  echo "missing directory $dd/database"
  goto done
endif

set probes_map=$dd/f3p.probes.map.ace.gz
if (-e $ff) then
  goto done
endif

if (! -e $dd/f3p.cosmid.map.gz) then
  $ici/bin/tacembly $dd << EOF
    query find sequence Is_cosmid && IntMap == $chrom
    select -o $dd/f3p.cosmid.map  s,m,a1,a2 from s in @, m in s->IntMap, a1 in m[1], a2 in m[2]
    quit
EOF
  gzip $dd/f3p.cosmid.map
endif
 


zcat $dd/f3p.cosmid.map.gz ZZZZZ.gz $ff.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){nn++;c[nn]=$1;chr=$2;c1[nn]=$3;c2[nn]=$4;next;} if(chr!=$11)next;p=$1;p1=$12;p2=$13; for (n=n;n>0 && c1[n] > p1;n--){};for (n=n;n<=nn && c2[n] < p1;n++){};printf("%s\t%s\t%d\t%d\t%d\t%d\n",c[n],p,c1[n],c2[n],p1,p2);}' > $dd/f3p.cosmid2probes.txt
cat $dd/f3p.cosmid2probes.txt | sort -k 1,1 -k 5,5n |  gawk -F '\t' '{c=$1;if(c!=oldc)printf("\nSequence %s\n",c);oldc=c;printf("Probe_hit \"%s\" %d %d\n",$2,$5-$3+1,$6-$3+1);}END{printf("\n");}'  > $dd/f3p.cosmid2probes.ace
cat $dd/f3p.cosmid2probes.txt | sort -k 1,1 -k 5,5n |  gawk -F '\t' '{p=$2;p1=$5;p2=$6;printf("Probe %s\nIntMap %s %d %d\nCapture\n\n", p,chr,p1,p2);}' chr=$chrom  > $dd/f3p.cosmid2probes.ace

gzip -f $dd/f3p.cosmid2probes.ace

echo "pparse $dd/f3p.cosmid2probes.ace.gz"  | bin/tacembly $dd -no_prompt


done:
 echo done
