#!bin/tcsh -f

echo -n "f2 tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt start : "
date

if (! -e tmp/X.$MAGIC/EHITS/f2.done) then
  if (-e tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt1) \rm tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt1
    foreach run (`cat MetaDB/$MAGIC/RunsList`)
      echo $run
      set ff=tmp/INTRONRUNS/$run/$run.u.doubleIntronSupport.ace.gz
      if (-e $ff) then
        gunzip -c $ff  | gawk '{gsub(/\"/,"",$0);}/^DoubleIntron/{g=$2;}/^Run_U/{gg[g]+= $6;uu[g]+= $6;next;}/^Run_nU/{gg[g]+= $6;nu[g]+= $6;next;}END{for(g in gg)if(gg[g]>0)printf("%s\t%d\t%d\t%d\n",g,gg[g],uu[g],nu[g]);}' >>   tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt1
      else
        echo "missing file $ff"
      endif
    end 
 
  cat tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt1 | gawk '{n2[$1]+=$2;n3[$1]+=$3;n4[$1]+=$4;}END{for(k in n2)printf("%s\t%d\t%d\t%d\n",k,n2[k],n3[k],n4[k]);}' | gzip > tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt.gz
  \rm   tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt1 
endif

touch  tmp/X.$MAGIC/EHITS/f2.done

exit 0
