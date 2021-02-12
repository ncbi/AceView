#!bin/tcsh -f

set run=$1
set justMito=$2
set poolDummy=$3
set solidDummy=$4
set collect=$5
set zone=$6

set uu=u

set minSnpFrequency2=$minSnpFrequency
if (-e tmp/SNP_BRS/$MAGIC.minFrequency.txt) then
  # reset minSnpFrequency to at least twice the percent error rate, the error_profile is given in error per million 
  set z=`cat tmp/SNP_BRS/$MAGIC.minFrequency.txt | gawk '{if ($1 == run) z=($4 + 0)/5000 ; if (mf<z)mf=z;}END{print int(mf);}' mf=$minSnpFrequency run=$run`
  set minSnpFrequency2=$z
endif

echo "scripts/snp.collect.tcsh $1 $2 $3 $4 $5 $6 minSnpCover=$minSnpCover minSnpCount=$minSnpCount minFrequency=$minSnpFrequency2"
if (0 &&  -e  tmp/SNP_BRS/$run/$zone.BRS.$uu.gz && ! -e  tmp/SNP/$run/$zone.plus.$uu.BF.gz && ! -e  tmp/SNP/$run/$zone.minus.$uu.BF.gz ) then
# this code is superseded by a direct export in snp.c in case $collect == count
      echo -n "wiggle strand plus start : " 
      date
      gunzip -c tmp/SNP_BRS/$run/$zone.BRS.$uu.gz | gawk -F '\t' '/^nAli/{next}/^nBpAli/{next}/^ERR/{next}{if($1 != "~")gene=$1;x=$2;if(length($3)>1)next;if($4 != strand)next;printf("%s\t%d\t%d\t%d\n",gene,x-1,x,$6);}' strand="+" | bin/wiggle -I BG -O BF -out_step 10 -gzo -o tmp/SNP/$run/$zone.plus.$uu

      echo -n "wiggle strand minus start : "
      date
      gunzip -c tmp/SNP_BRS/$run/$zone.BRS.$uu.gz | gawk -F '\t' '/^nAli/{next}/^nBpAli/{next}/^ERR/{next}{if($1 != "~")gene=$1;x=$2;if(length($3)>1)next;if($4 != strand)next;printf("%s\t%d\t%d\t%d\n",gene,x-1,x,$6);}' strand="-" | bin/wiggle -I BG -O BF -out_step 10 -gzo -o tmp/SNP/$run/$zone.minus.$uu
      echo -n "wiggle done : "
      date
endif

# avoid analysing MetaDB/*List in submitted scripts
set solid=""
if ($solidDummy == solid) set solid="-solid"

set pool=""
if ($poolDummy == pool) set pool="-pool"

# echo the command, then parse the BRS table and export a list of snp counts in target coordinates detailling both strands
# the program collects the snps not wild type at risk 1/100, function snpNotLow() 
# and all the counts at the positions listed in -snp_list even if in the current run we see a wild type

    set vdb=""
    if ($collect == count && -e  tmp/SNP_LIST/$zone/$MAGIC.Variant.list) set vdb="-snp_list   tmp/SNP_LIST/$zone/$MAGIC.Variant.list"
    if ($collect == count && -e  ../Global.Variant.list) set vdb="-snp_list  ../Global.Variant.list"
   
####  Main work
echo -n "$collect main work : "
date

    set out=tmp/SNP/$run/$zone.$collect.$uu.snp
    if ($collect == count) set out=tmp/SNP/$run/$MAGIC.$zone.$collect.$uu.snp

    set select8kb=""
    set target=av
    if ($collect == detect && -e tmp/METADATA/$target.selected8kbTranscriptList.txt) then
      set select8kb="-selected8kbList  tmp/METADATA/$target.selected8kbTranscriptList.txt"
    endif
    if (! -e $out.gz) then
  
      # in the new method, we add the BRS of the runs into the groups, so we can use the high thresholds at the detect stage
      set mins=" -minCover $minSnpCover -minMutant $minSnpCount -minFrequency $minSnpFrequency2"
      echo " bin/snp -BRS_$collect $solid $mins -run $run $pool $vdb -strategy $Strategy  $select8kb -i  tmp/SNP_BRS/$run/$zone.BRS.$uu.gz  -o $out -gzo" 
             bin/snp -BRS_$collect $solid $mins -run $run $pool $vdb -strategy $Strategy  $select8kb -i  tmp/SNP_BRS/$run/$zone.BRS.$uu.gz  -o $out -gzo 

    endif

######

    if ($collect == detect) then
      gunzip -c tmp/SNP/$run/$zone.$collect.$uu.snp.gz |  gawk '/^#/{next}/Incompatible_strands/{next;}{gsub(/>/,"2",$3);c=$9;m=$10;if (c>=minCover && m>=minCount && 100*m >= minFreq*c) printf("%s:%s_%s\n",$1,$2,$3);}' minCover=$minSnpCover minCount=$minSnpCount minFreq=$minSnpFrequency2 | sort -u >   tmp/SNP_LIST/$zone/Variant.$run.$uu.list
    endif

stats:

echo -n "$collect stats : "
date

set zone1="$zone"
if ($collect == count) set zone1="$MAGIC.$zone"

    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz  |  gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}' run=$run minC=$minSnpCover minF=$minSnpFrequency2 | sort -k 5n >   tmp/SNP/$run/$zone1.$collect.$uu.stats
    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz  |  gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}'  run=$run minC=$minSnpCover minF=95 | sort -k 5n >>   tmp/SNP/$run/$zone1.$collect.$uu.stats


    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz | gawk -F '\t' '/Incompatible_strands/{next;}{type=$3;if($10 + $11 > 0){n[type]++;m[type]+=$10;w[type]+=$11;}}END{for(t in n)printf("%s\t%s\t%d\t%d\t%d\n",t,run,n[t],m[t],w[t]);}' run=$run | sort >  tmp/SNP/$run/$zone1.$collect.types.$uu.txt
  
  if ($collect == count) then
    gunzip -c  tmp/SNP/$run/$MAGIC.*.$collect.$uu.snp.gz | bin/histo -snp -o  tmp/SNP/$run/$zone1.snp_frequency
  endif

if ($justMito == 1) then
  touch tmp/SNP/$run/$zone1.snp.$collect.justmito.done
else
  touch tmp/SNP/$run/$zone1.snp.$collect.done
endif

echo -n "$collect done : "
date
echo done

exit 0


##########################################################
## verif scripts

foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/zoneg.1.detect.u.snp.gz) then
    echo tmp/SNP/$run
  endif
end


foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/zoneg.1.$collect.u.snp.gz) then
    echo tmp/SNP/$run
  else
    touch tmp/SNP/$run/snp.$collect.done
  endif
end

foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/s1.err) echo tmp/SNP/$run
  grep Exit tmp/SNP/$run/s1.err
end

# grep Exit tmp/SNP/*/s1.err | grep -v '= 0' | gawk '{split($1,aa,"/");printf("\\rm -rf tmp/SNP*/%s\n",aa[3]);}' > _rm

exit 0

# hack to look for saturation
\rm _toto
foreach ii (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
  cat B$ii/Global.B$ii.Variant.list >> _toto
end
cat _toto | sort -u > Global.B0_15.Variant.list

# by hank i run the scripts 15 times stopping at say 15 .. 8 7 .. 0
\rm _toto
foreach ii (0 1   )
  cat B$ii/Global.B$ii.Variant.list >> _toto
end
cat _toto | sort -u > Global.B0_1.Variant.list

