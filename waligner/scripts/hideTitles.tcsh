#!bin/tcsh -f

echo Hiding private titles before expotation


cat << EOF > hideTitles.awk
/^#Sorting_title/{next;}
/^#Sample/{next;}
/^#Title/{printf("%s",\$1);for(i=2;i<=NF;i++){z=\$i;j=index(\$i,"/");if(j>0)z=substr(\$i,1,j-1);printf("\t%s",z);}printf("\n");next;}
{print}
EOF

if (! -d RESULTS/Expression/Clean) mkdir RESULTS/Expression/Clean
foreach uu (u nu)
  if ($uu == u) set uu2=unique
  if ($uu == nu) set uu2=quasi_unique
  foreach type (reads_aligned_per_gene expression_index)
    foreach GM (GENE)
      foreach target (av)
        set traget2=$target
        if ($target == av) set target2=AceView
        set f=RESULTS/Expression/$uu2/$target/$MAGIC.$target2.$GM.$uu.$type.txt.gz   
        if (-e $f) then
          echo $f
          gunzip -c $f | gawk -F '\t' -f hideTitles.awk | gzip > RESULTS/Expression/Clean/$MAGIC.$target2.$GM.$uu.$type.txt.gz  
        endif
      end
    end
  end
end



