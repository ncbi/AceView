#!bin/tcsh -f

set run=$1
set groupLevel=$2

if ($groupLevel == 0) then

  if  (! -e tmp/INTRONRUNS/$run/$run.u.list) then
       foreach sublib (`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($1 == run)print $2;}' run=$run | sort  -u`)
          if ( -e tmp/INTRONRUNS/$sublib/$sublib.u.list) then
             cat tmp/INTRONRUNS/$sublib/$sublib.u.list >> tmp/INTRONRUNS/$run/$run.u.list
          endif
       end
  endif

  if (! -e  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz) then
      if (-e  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts) \rm  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts
      foreach sublib (`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($1 == run)print $2;}' run=$run | sort  -u`)
           if ( -e  tmp/INTRONRUNS/$sublib/$sublib.u.intronSupport.counts.gz) then
              gunzip -c  tmp/INTRONRUNS/$sublib/$sublib.u.intronSupport.counts.gz  >>  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts
           endif
      end
      cat tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 ; nn[z] += $4;}END {for (z in nn) printf ("%s\t%d\n", z, nn[z]) ;}' | gzip >  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz
      \rm tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts 
  endif

else

  if (! -e  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz) then
    if (-e  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts) \rm  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts
    foreach r (`cat MetaDB/$MAGIC/g2r | gawk -F '\t' '{if($1 == g)print $2 ;}' g=$run | sort -u`)
      ls -ls  tmp/INTRONRUNS/$r/$r.u.intronSupport.counts.gz
      if ( -e  tmp/INTRONRUNS/$r/$r.u.intronSupport.counts.gz) then
              gunzip -c  tmp/INTRONRUNS/$r/$r.u.intronSupport.counts.gz  >>  tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts
           endif
      end
      cat tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 ; nn[z] += $4;}END {for (z in nn) printf ("%s\t%d\n", z, nn[z]) ;}' | gzip >  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz
      \rm tmp/INTRONRUNS/$run/$run.u.intronSupport.precounts 
  endif

endif

foreach target ($Etargets)
  if (! -e tmp/INTRONRUNS/$run/known_introns.$target.ace) then

          gunzip -c tmp/METADATA/$target.fr.introns.gz ZZZZZ.gz  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1,$2,$3]=1;next;}if (ok[$1,$2,$3]==1) {n++;nn+=$4;}}END{printf("Ali %s\nKnown_introns %s %d introns_supported_by %d reads\n\n", run, target, n,nn);}' run=$run target=$target  >  tmp/INTRONRUNS/$run/known_introns.$target.ace
  endif

end

exit 0
