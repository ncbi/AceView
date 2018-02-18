#!bin/tcsh -ef

set group=$1
set zone=$2
set collect=$3
set phase=$4
set uu=u

set zoneType=""
if ($Strategy != RNA_seq)   set zoneType="zoneg."
if (-e tmp/SNPG/$group/$zoneType$zone.$collect.u.snp.gz) goto stats
if (-e  tmp/SNPG/$group/zoneList.$zone) \rm tmp/SNPG/$group/zoneList.$zone

foreach run (`cat  tmp/SNPG/$group/RunList`)
  if ($zone == mito) then
    echo tmp/SNP/$run/$zone.$collect.u.snp.gz >>  tmp/SNPG/$group/zoneList.$zone
  else
    if ($Strategy == RNA_seq) then
      echo tmp/SNP/$run/$zone.$collect.u.snp.gz  >>  tmp/SNPG/$group/zoneList.$zone
    else
      if ( -e  tmp/SNP_ZONERUN/$run/chom2zone.$zone) then
        foreach zone (`cat tmp/SNP_ZONERUN/$run/chom2zone.$zone`)
          echo   tmp/SNP/$run/zoneg.$zone.$collect.u.snp.gz >>  tmp/SNPG/$group/zoneList.$zone
        end
      endif    
    endif
  endif
end

# toto1: accumulate the data, toto2: a hack to remove the N in the snpnets

if (-e tmp/SNPG/$group/toto1.$zone) \rm  tmp/SNPG/$group/toto1.$zone

if (0) then
########### this horrible hack is supposed to fix an anomaly with the snipnets in case they contain N
# however we cannot afford to do a sort of the whole file
# the snipnet fix should be incorporated in the bin/snp -merge 
 foreach ff (`cat  tmp/SNPG/$group/zoneList.$zone`)
  if (-e $ff) then
    echo $ff
    gunzip -c $ff >>  tmp/SNPG/$group/toto1.$zone
  endif
end
# ~~~~ is after z in ascii order i.e. value 126, this way we get the last gene

# cat  tmp/SNPG/$group/toto1.$zone tmp/SNPG/zzzz | sort -k 1,1 -k 2,2n -k 3,3   | gawk -F '\t' '/^#/{nt[$0]++;if(nt[$0]<2)print;next;}{if(NF > nfMax)nfMax=NF;z= $1 "\t" $2 "\t" $3 ; if(z != oldz){ if(index(bestP bestN,"NN")==0){for(n=1;n<=0+nn;n++){printf("%s\t%s\t%s",oldz,bestP,bestS);for(i=6;i<=nfMax;i++)printf("\t%s",dd[n,i]);printf("\n");}}nn=0;bestP=$4;bestS=$5;oldz=z;pHasN=split(bestP,aa,"N")-1;sHasN=split(bestS,aa,"N")-1;}nn++;for(i=6;i<=nfMax;i++)dd[nn,i]=$i;if(pHasN>0){pHasN2=split($4,aa,"N")-1;if(pHasN2<pHasN){pHasN=pHasN2;bestP=$4;}}if(sHasN>0){sHasN2=split($5,aa,"N")-1;if(sHasN2<sHasN){sHasN=sHasN2;bestS=$5;}}}' | bin/snp -merge -run $group -minCover $minSnpCover | gzip >  tmp/SNPG/$group/$zoneType$zone.$collect.u.snp.gz
##########end of hottible hack
endif

set nmax=`wc -l tmp/SNPG/$group/zoneList.$zone`
   set mins="  -minCover  $minSnpCover -minMutant $minSnpCount -minFrequency $minSnpFrequency"
   if ($collect == count) set mins=" -minCover  $minSnpCover"
   gunzip -c  `cat  tmp/SNPG/$group/zoneList.$zone`  | bin/snp -merge -run $group $mins | gzip >  tmp/SNPG/$group/$zoneType$zone.$collect.u.snp.gz

 touch tmp/SNPG/$group/$phase.$zone.done1

# accumulate the wiggles
if ($collect == count) then 
  set err=0
  if (-e tmp/SNPG/$group/_BV) \rm tmp/SNPG/$group/_BV
  foreach run  (`cat  tmp/SNPG/$group/RunList`)
    if (-e  tmp/SNP/$run/$zone.$collect.BV.gz) then
      gunzip -c   tmp/SNP/$run/$zone.$collect.BV.gz >>  tmp/SNPG/$group/$zone.$collect._BV
    else
      set err=1
    endif
  end
  if (-e tmp/WIGGLESIG/$group/$zone.$collect._BV) then
    if ($err == 0) then
       bin/wiggle -i  tmp/SNPG/$group/$zone.$collect._BV  -I BV -O BV  -out_step 10 -gzo  -o tmp/WIGGLESIG/$group/$zone.$collect
    endif
    \rm  tmp/WIGGLESIG/$group/$zone.$collect._BV
  endif
endif

stats:
   if ($collect == detect) then
      gunzip -c tmp/SNPG/$group/$zone.detect.$uu.snp.gz | gawk '/^#/{next}/Incompatible_strands/{next;}{gsub(/>/,"2",$3);c=$9;m=$10;if (c>=minCover && m>=minCount && 100*m >= minFreq*c) printf("%s:%s_%s\n",$1,$2,$3);}' minCover=$minSnpCover minCount=$minSnpCount minFreq=$minSnpFrequency | sort -u >  tmp/SNP_LIST/$zone/Variant.$group.$uu.list 

    endif

    gunzip -c tmp/SNPG/$group/$zoneType$zone.$collect.u.snp.gz | gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}'  run=$group minC=$minSnpCover minF=$minSnpFrequency | sort -k 5n >  tmp/SNPG/$group/$zoneType$zone.$collect.u.stats
    gunzip -c tmp/SNPG/$group/$zoneType$zone.$collect.u.snp.gz | gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}'  run=$group minF=95 minC=$minSnpCover | sort -k 5n >>  tmp/SNPG/$group/$zoneType$zone.$collect.u.stats

touch tmp/SNPG/$group/$phase.$zone.done2

touch tmp/SNPG/$group/$phase.$zone.done

exit 0
