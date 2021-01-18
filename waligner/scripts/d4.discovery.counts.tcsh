#!bin/tcsh -f

set group=$1
set all_counts=$2

set toto=tmp/introns/$group/d4.discovery.counts
echo $toto
echo "# $2 $3 $4 date " > $toto
date >> $toto

  echo  "d4:group:$group "
  if (-e $toto.1) \rm $toto.1
  touch $toto.1
# 50 20 10 5 3 2 1 0 -1
  foreach ii (50 20 10 5 3 2 1 0 -1)
    if ($ii == -1) \rm $toto.1
    touch $toto.1
    echo -n "$group"  >> $toto
    if ($ii > 1) echo -n "\t$ii duos"  >> $toto
    if ($ii == 1) echo -n "\t$ii duo"  >> $toto
    if ($ii == 0) echo  -n "\tuno+duo"  >> $toto
    if ($ii == -1) echo  -n "\tuno"  >> $toto
    if (-e $toto.c) \rm $toto.c

    if ($ii >= 0 && -e tmp/INTRONRUNS/$group/$group.u.intronSupport.counts.gz) then
      gunzip -c  tmp/INTRONRUNS/$group/$group.u.intronSupport.counts.gz | gawk -F '\t' '{if($4 >= ii) printf("%s\t%s\t%s\n",$1,$2,$3);}' ii=$ii >> $toto.1
    endif

      foreach chrom ($chromSetAll)
        foreach jj (50 20 10 5 3 2 1)
          if($ii != $jj) continue
          if($ii != 0 && $jj < $ii) continue
          if ($group != ANY) then
            foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if($1 == g)print $2}END{print g;}' g=$group | sort -u`)
              if (-e tmp/OR/$run2/$chrom/$run2.$jj.0.txt) then
                cat tmp/OR/$run2/$chrom/$run2.$jj.*.txt | gawk -F '\t' '/^INTRON/{if($9 ~ /g[ct]_ag/)printf("%s\t%s\t%s\n",$3,$4,$7);}'  >> $toto.1
              endif
            end
          else
            foreach run2 (`cat MetaDB/$MAGIC/RunList`)
              if (-e tmp/OR/$chrom/$run2.$jj.0.txt) then
                cat tmp/OR/$run2/$chrom/$run2.$jj.*.txt | gawk -F '\t' '/^INTRON/{if($9 ~/g[ct]_ag/)printf("%s\t%s\t%s\n",$3,$4,$7);}'  >> $toto.1
              endif
            end
          endif
        end
      end
      if ($ii <= 0) then
        if ($group != ANY) then
          foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if($1 == g)print $2}END{print g;}' g=$group | sort -u`)
            if (-e  tmp/OR/$run2/d1.$run2.de_uno.txt.gz) then
              gunzip -c tmp/OR/$run2/d1.$run2.de_uno.txt.gz | gawk -F '\t' '{if($9 ~ /g[ct]_ag/){a2=$5;b1=$7;if(a2<b1)da=1;else da=-1;a2+=da;b1-=da;z=chrom "\t" a2 "\t" b1 ; if (nn[z]<1){nn[z]+=$11;printf("%s\t%09d\t%09d\n",$3,a2,b1);}}}'  >> $toto.1
            endif
          end
        else
          foreach run2 (`cat MetaDB/$MAGIC/RunList`)
            if (-e  tmp/OR/$run2/d1.$run2.de_uno.txt.gz) then
              gunzip -c tmp/OR/$run2/d1.$run2.de_uno.txt.gz | gawk -F '\t' '{if($9~/g[ct]_ag/){a2=$5;b1=$7;if(a2<b1)da=1;else da=-1;a2+=da;b1-=da;z=chrom "\t" a2 "\t" b1 ; if (nn[z]<1){nn[z]+=$11;printf("%s\t%09d\t%09d\n",$3,a2,b1);}}}'  >> $toto.1
            endif
          end
        endif
      endif

    if (-e $toto.c) \rm $toto.c
    touch  $toto.c
    if (-e $toto.1) then
      cat  $toto.1 | sed -e 's/CHROMOSOME_//' |  sort -u  >> $toto.c
    endif
    set ni_count=`wc $toto.c | gawk '{printf ("%d", $1);}'`
    echo -n "\t$ni_count" >> $toto
    foreach target ($Etargets)
      set target_count=`echo $all_counts | gawk '{i=index($1,t);split(substr($1,i+1),aa,"-");split(aa[2],bb,"_");print bb[1];}' t=$target`
      if (-e tmp/METADATA/$target.f.introns) then
        echo "$ii $target $target_count"
        cat   tmp/METADATA/$target.[fr].introns | sed -e 's/CHROMOSOME_//' > $toto.$target
        echo ZZZZZ >> $toto.$target
        cat  $toto.c  >> $toto.$target
        cat $toto.$target | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz==0){u=$1 0+$2 0+$3;uu[u]=1;next;}u=$1 0+$2 0+$3;if(uu[u]==1) print; uu[u]=0;}' | wc | gawk '{ni2=ni;c2=cDNA;if(ni2==0)ni2=1;if(c2==0)c2=1;printf ("\t%d\t%d\t%.2f\t%.2f", $1,ni-$1,100.0*$1/ni2,100.0*$1/c2);}' cDNA=$target_count ni=$ni_count  >> $toto
      endif
    end
    echo >> $toto
  end

\rm $toto.*


