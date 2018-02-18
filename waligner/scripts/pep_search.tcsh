#!bin/tcsh -f

setenv ff $1
set pepRef=$2

set nn=`cat $ff | gawk '{n++;}END{print n}'`
set ii=0
set titi7=$3.$ff

touch $titi7
\rm  $titi7
while ($ii < $nn)
  @ ii = $ii + 1
  cat $ff | gawk '{n++;if (n==ii)printf("%s\t%s\t%s",$1,$2,$3);}' ii=$ii >> $titi7
  set pep=`cat $ff | gawk '{n++;if (n==ii)printf("%s",$1);}' ii=$ii`
  cat  $pepRef | gawk /$pep/'{printf("\t%s\t%d",$1,index($2,pep));}' pep=$pep >> $titi7
  echo >> $titi7
end
