#!bin/tcsh -f
set target=$1
set parity=$2
set out=$3

set roc=$out/_roc.$parity
set km=$out/_km.$parity

  if (-e $roc) \rm $roc
  echo -n "# " > $roc
  date >> $roc
  cat $roc > $km

foreach kk (1 2 3 4)

  echo "A continuous classifier is contructed for 503 runs (513 minus 10 runs misclassified for sex), using several methods" >> $roc
  echo "The 503 runs are then ordered according to this classifier" >> $roc
  echo "The cumulative count of the known samples belonging to the ordered list is evaluated, yielding a ROC curve" >> $roc

  if ($kk == 1) echo "ROC curve for favorable and unfavorable outcome in the 96u + 180f known samples, ordered by their continuous classifier" >> $roc
  if ($kk == 2) echo "ROC curve for event-free survival, ordered by their continuous classifier" >> $roc
  if ($kk == 3) echo "ROC curve for overall survival, ordered by their continuous classifier" >> $roc
 # 3 4 5 6 7 8
  if ($kk == 1) set title="favorable versus unfavorable (extremes)"
  if ($kk == 2) set title="event free survival"
  if ($kk == 3) set title="overal survival"
  if ($kk == 4) set title="sex"

  foreach ii (4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
    set kk1=`echo $kk | gawk '{print $1+25}'`
    set mm=`cat $out/$MAGIC.outcome.$parity.txt | gawk -F '\t' '/^#/{printf("%s", $ii)}' ii=$ii`
    if ($ii < 12 || ($ii < 23 && $ii > 14)) then
      echo -n "$target\t" >> $roc
      if ($kk == 1) echo -n "Outcome u/f\t" >> $roc
      if ($kk == 2) echo -n "Outcome EF\t" >> $roc
      if ($kk == 3) echo -n "Outcome OS\t" >> $roc
      if ($kk == 4) echo -n "Outcome Sex\t" >> $roc
       echo -n "Classifier: $mm\t" >> $roc
       cat $out/$MAGIC.outcome.$parity.txt |   grep -v Z | sort -k $ii'n' | gawk -F '\t' -f scripts/NB.rocplot.awk  title="$title" kk=$kk1 clf="$mm" ii=$ii >> $roc
    endif
    if (0) then
      echo -n "$target\t" >> $roc
      if ($kk == 1) echo -n "Outcome u/f\t" >> $roc
      if ($kk == 2) echo -n "Outcome EF\t" >> $roc
      if ($kk == 3) echo -n "Outcome OS\t" >> $roc
      if ($kk == 4) echo -n "Outcome Sex\t" >> $roc
      echo -n "Classifier: $mm\t" >> $roc
       cat $out/$MAGIC.outcome.$parity.txt |  grep OSINTERNAL|  grep -v Z | sort -k $ii'n' | gawk -F '\t' -f scripts/NB.rocplot.awk  title="$title" kk=$kk1 clf="$mm" ii=$ii >> $roc
    endif
  if ($parity == 1 && (($kk == 2 && ($ii == 11 || $ii == 12 || $ii == 22 || $ii == 23)) || ($kk == 3 && ($ii == 13 || $ii == 14 || $ii == 24 || $ii == 25)))) then

    echo -n "$target\t" >> $km
    if ($kk == 1) echo -n "Outcome u/f\t" >> $km
    if ($kk == 2) echo -n "Outcome EF\t" >> $km
    if ($kk == 3) echo -n "Outcome OS\t" >> $km
    if ($kk == 4) echo -n "Outcome Sex\t" >> $roc
       set titleKm="among high risk patients" 
      echo -n "Classifier: $mm\t" >> $km
      cat $out/$MAGIC.outcome.$parity.txt | grep _HR |  sort -k $ii'n' | gawk -F '\t' -f scripts/NB.rocplot.awk  title="$title" kk=$kk1  ii=$ii | grep uf >> $km 
    echo -n "$target\t" >> $km
    if ($kk == 1) echo -n "Outcome u/f\t" >> $km
    if ($kk == 2) echo -n "Outcome EF\t" >> $km
    if ($kk == 3) echo -n "Outcome OS\t" >> $km
    if ($kk == 4) echo -n "Outcome Sex\t" >> $roc
       echo -n "Classifier: $mm\t" >> $km
      cat $out/$MAGIC.outcome.$parity.txt | grep _HR |  sort -k $ii'n' | gawk -F '\t' -f scripts/NB.rocplot.awk  title="$title" kk=$kk1  ii=$ii | grep MCC >> $km
      echo -n "$target\tClassifier HR\t" >> $km
      cat $out/$MAGIC.outcome.$parity.txt | grep _HR | sort -k $ii'n' | gawk -F '\t' -f scripts/NB.KaplanMeier.awk  title="$title" JJ=5 title="$titleKm" >> $km
    endif
  end
end

if ($parity == 1) then
  grep uf $km >> $roc
endif
