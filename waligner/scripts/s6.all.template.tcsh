#!bin/tcsh -f
set selectReject=$1
set n0=$2
set n1=$3

cat tmp/KL.$MAGIC/templateList | gawk '{n++;if(n > n0 && n <= n0 + n1) print}' n0=$n0 n1=$n1 > tmp/KL.$MAGIC/s6/s6.$selectReject.$n0

foreach template (`cat  tmp/KL.$MAGIC/s6/s6.$selectReject.$n0`)
  if (-d  tmp/KL.$MAGIC/s6/$template) continue
  if (! -d tmp/KL.$MAGIC/s6/$template) mkdir tmp/KL.$MAGIC/s6/$template
  if (! -e tmp/KL.$MAGIC/s6/$template/any.$selectReject.txt) then
    echo -n "$template "
    date
    ## clean up
    if (-e tmp/KL.$MAGIC/s6/$template/$selectReject.1) \rm tmp/KL.$MAGIC/s6/$template/$selectReject.1

    ## grab the data for each group
    foreach group (`cat MetaDB/$MAGIC/GroupSnpList MetaDB/$MAGIC/RunSnpList`)
      if (-d tmp/KL.$MAGIC/$group) then
        gunzip -c  tmp/KL.$MAGIC/tri.$selectReject.gz ZZZZZ.gz tmp/KL.$MAGIC/$group/$template.gz | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz=1;next;}/Incompatible_strands/{a=$12/($12+$13+1);b=$14/($14+$15+1);if(a < b - .5 || a > b + .5){if(rj == "rejected")print ; next;}}{if(zz<1){ok[$1]=1;next;} z = $1 ":" $2 ":" $3  ; if(ok[z]==1) print }' rj=$selectReject >> tmp/KL.$MAGIC/s6/$template/$selectReject.1
      endif
    end

    ## sort and cumulate per template
    cat   tmp/KL.$MAGIC/s6/$template/$selectReject.1 | sort -k 1,1 -k 2,2n -k 3,3 -k 6,6 >  tmp/KL.$MAGIC/s6/$template/$selectReject.txt

    ## cleanup
    \rm   tmp/KL.$MAGIC/s6/$template/$selectReject.1 

    ## compute the global count
    cat tmp/KL.$MAGIC/s6/$template/$selectReject.txt | bin/snp -merge -run $MAGIC > tmp/KL.$MAGIC/s6/$template/any.$selectReject.txt

  endif
end

