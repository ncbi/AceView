#!bin/tcsh -f

set support2ace=$1
set run=$2
set target=$3
set GM=$4
set gm=$5
set uu=$6
set phase=$7

if (! -e tmp/GENERUNS/$run/$run.$target.$GM.$uu."$gm"Support.ace.gz) then   
  if (-e tmp/GENERUNS/$MAGIC.$target.$GM.$uu.ace) \rm  tmp/GENERUNS/$MAGIC.$target.$GM.$uu.ace
  echo "  bin/bestali $support2ace -gzo -o tmp/GENERUNS/$run/$run.$target.$GM.$uu -inFileList tmp/GENERUNS/$run/$run.$target.$GM.$uu.list -run $run"
          bin/bestali $support2ace -gzo -o tmp/GENERUNS/$run/$run.$target.$GM.$uu -inFileList tmp/GENERUNS/$run/$run.$target.$GM.$uu.list -run $run
endif

foreach kb (8kb 5kb)

  set tutu=tmp/GENERUNS/$run/$run.coverage_of_"$kb"_transcripts
  if (-e $tutu.txt.gz) gunzip -f $tutu.txt.gz  
  if ($phase == m2bH && ! -e $tutu.txt) then
    set ok=0
    echo "... constructing tmp/GENERUNS/$run/$run.3pHisto.$kb.txt"
    echo -n "## $run\tfile $tutu.txt\t" > $tutu.txt
    date	>> $tutu.txt
    echo ' ' > tmp/GENERUNS/$run/toto.$$
    foreach target (av $Etargets)
      if (-e tmp/METADATA/$target.selected"$kb"TranscriptList.txt) then
        cat  tmp/METADATA/$target.selected"$kb"TranscriptList.txt > tmp/GENERUNS/$run/toto1.$$
        foreach ff (`cat  tmp/GENERUNS/$run/$run.$target.$GM.3pHisto.$kb.list`)
          set ok=1
          cat $ff >>  tmp/GENERUNS/$run/toto.$$
        end
      endif
    end
    set nMax=`cat tmp/GENERUNS/$run/toto1.$$ | gawk '/^#/{next;}{n++}END{print n}'`
    if ($ok > 0) then
      cat  tmp/GENERUNS/$run/toto.$$ | gawk -F '\t' -f scripts/m2bH.3pHisto.awk nMax=$nMax kb=$kb run=$run > tmp/GENERUNS/$run/toto2.$$
      cat  tmp/GENERUNS/$run/toto2.$$ | head -12 | gawk '/^##/{print}' >>  $tutu.txt
      cat  tmp/GENERUNS/$run/toto2.$$ | head -12 | gawk '/^##/{next}/^#/{print}' >> $tutu.txt
      cat  tmp/GENERUNS/$run/toto2.$$ | gawk '/^#/{next}{print}' | scripts/tab_sort -k 4nr >> $tutu.txt
    else
      \rm  $tutu.txt
    endif
   # gzip $tutu.txt
   \rm tmp/GENERUNS/$run/toto.$$ tmp/GENERUNS/$run/toto1.$$ tmp/GENERUNS/$run/toto2.$$
  endif

end
