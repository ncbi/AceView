#!bin/tcsh -f

# Merge all compare files of the project, for gene, mRNA and micro-array
foreach target (av AGLuK)
 foreach GM (GENE MRNAH MA)
  set toto=RESULTS/Expression/unique/$target/$MAGIC.$target.$GM.compared.fused.txt
  echo "## Fused table of all DEG lists in $target project $MAGIC " > $toto
  echo "## Genes up and down in each comparison are listed in succession." >> $toto
  echo "## Details of the histograms of all comparison are available in the files" >> $toto

  foreach cc (`cat MetaDB/$MAGIC/ccc_pair.list | gawk -F ',' '{print $1}'`)
     foreach ff (`ls RESULTS/Expression/unique/$target/$MAGIC.*.$GM.u.$cc.compare.txt`)
       cat $ff  | head -30 | gawk '/^##/{if(index($0, "file")>0)print}' >> $toto
     end
  end
  echo >> $toto
  echo >> $toto
  foreach cc (`cat MetaDB/$MAGIC/ccc_pair.list | gawk -F ',' '{print $1}'`)
    foreach ff (`ls RESULTS/Expression/unique/$target/$MAGIC.*.$GM.u.$cc.compare.txt`)
      cat $ff  | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 >> $toto
    end
  end
  set n=`cat $toto | gawk -F '\t' '/^#/{next;}{if(NF > 5)n++;}END{print n}'
  if ($n < 40)\rm $toto
  \rm $toto.1
 end
end

## create a fused table
set toto=RESULTS/Expression/unique/$MAGIC.compare.fused.txt
echo > $toto.1
foreach target (av AGLuK)
  foreach GM (GENE MRNAH MA)
    foreach cc (`cat MetaDB/$MAGIC/ccc_pair.list | gawk -F ',' '{print $1}'`)
      foreach ff (`ls RESULTS/Expression/unique/$target/$MAGIC.*.$GM.u.$cc.compare.txt`)
        cat $ff  | gawk -F '\t' '/^#/{next;}{if($13 == "" || $14 == "")next;printf ("%s\t",GM);print;}' GM=$GM >> $toto.1
      end
    end
  end
end

cat $toto.1 | gawk -F '\t' -f scripts/ma4.compare.MOA.awk




