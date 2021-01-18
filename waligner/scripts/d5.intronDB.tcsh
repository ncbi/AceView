#!bin/tcsh -f

set phase=$1
set project=$2
echo "phase=$phase  project=$project"

if ($phase == cumul) then
  bin/tacembly MetaDB <<EOF
    bql -o MetaDB/$MAGIC/GroupIntronList select g from p in  ?project where p == "$MAGIC", g in p->run where g#intron
    bql -o tmp/OR/d5.$MAGIC.minS select g,m from p in  ?project where p == "$MAGIC", g in p->run where g#intron, a in g->ali, m in a->Candidate_introns[18]
EOF


set toto = tmp/OR/d5.$MAGIC.de_uno
echo ' ' > $toto.txt
foreach run (`cat MetaDB/$MAGIC/GroupIntronList`)
  set ff=tmp/OR/$run/d4.de_uno.txt.gz
  if (! -e $ff) set ff=tmp/OR/$run/d5.de_uno.txt.gz
  if (-e  $ff) then
    set minS=`cat  tmp/OR/d5.$MAGIC.minS | gawk -F '\t' 'BEGIN{n=1;}{if($1==run)n=$2+0;}END{print n}' run=$run`
    zcat $ff | gawk -F '\t' '{a1=$2+0;a2=$3+0;ii=$1 "__" a1 "_" a2 ; if($4>=minS)printf ("%s\t%s\t%s\t%s\t%d\n",ii,$5,"any",run,$4);}' run=$run minS=$minS>> $toto.txt
  endif
end
cat $toto.txt | sort -V | gzip > $toto.txt.gz
\rm  $toto.txt 
zcat $toto.txt.gz  | gawk -F '\t' '{ii=$1; if($5+0<1)next;if(ii!=old){if (n>0)printf("RNA_seq %d\n",n);n=0;printf("\nIntron \"%s\"\n",ii);}old=ii;printf("de_uno %s %d\n",$4,$5);n+=$5;}END{if (n>0)printf("RNA_seq %d\n",n);printf("\n");}' | gzip > $toto.ace.gz


echo "pparse $toto.ace.gz" | bin/tacembly GeneIndexDB -no_prompt

bin/tacembly GeneIndexDB << EOF
  pparse   tmp/OR/d5.$MAGIC.de_uno.ace.gz
  bql -o tmp/OR/d5.$MAGIC.all_introns select ?Intron
  save
  quit
EOF

cat  tmp/OR/d5.$MAGIC.all_introns | gawk '{split($1,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){b1=a1;b2=a2;}else{b1=a2;b2=a1;}printf("Intron %s\nIntMap %s %d %d\nLength %d\n\n",$1,chrom,a1,a2,b2-b1+1);}' | gzip >  tmp/OR/d5.$MAGIC.all_introns.intmap.ace.gz

echo "pparse tmp/OR/d5.$MAGIC.all_introns.intmap.ace.gz" | bin/tacembly GeneIndexDB -no_prompt

# donor acceptor
bin/tacembly GeneIndexDB << EOF
  find intron
  show -a -f tmp/OR/d5.$MAGIC.any_intron.ace
  quit
EOF

cat tmp/OR/d5.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nDonor %s__%d_%s\n",chrom, a1,s);printf("IntMap %s %d %d\n",chrom, a1-ds,a1);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("gt\n");next;}/^gc_ag/{printf("gc\n");next;}/^ct_ac/{printf("ct\n");next;}/^at_ac/{printf("at\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/OR/d5.$MAGIC.any_donor.ace

cat tmp/OR/d5.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nAcceptor %s__%d_%s\n",chrom, a2,s);printf("IntMap %s %d %d\n",chrom, a2,a2+ds);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("ag\n");next;}/^gc_ag/{printf("ag\n");next;}/^ct_ac/{printf("ac\n");next;}/^at_ac/{printf("ac\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/OR/d5.$MAGIC.any_acceptor.ace

bin/tacembly GeneIndexDB << EOF
  pparse tmp/OR/d5.$MAGIC.any_donor.ace
  pparse tmp/OR/d5.$MAGIC.any_acceptor.ace
  save
  quit
EOF
bin/tacembly GeneIndexDB << EOF
  bql -o tmp/OR/d5.$MAGIC.donors select d,t from d in ?Donor, t in d#gt
  bql -o tmp/OR/d5.$MAGIC.acceptors select a,t from a in ?Acceptor, t in a#ag
  quit
EOF
cat tmp/OR/d5.$MAGIC.donors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/OR/d5.$MAGIC.donors.next_gt.ace
cat tmp/OR/d5.$MAGIC.donors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/OR/d5.$MAGIC.donors.next_gt.ace
cat tmp/OR/d5.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/OR/d5.$MAGIC.donors.next_gt.ace
cat tmp/OR/d5.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/OR/d5.$MAGIC.donors.next_gt.ace

cat tmp/OR/d5.$MAGIC.acceptors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/OR/d5.$MAGIC.acceptors.next_ag.ace
cat tmp/OR/d5.$MAGIC.acceptors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/OR/d5.$MAGIC.acceptors.next_ag.ace
cat tmp/OR/d5.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/OR/d5.$MAGIC.acceptors.next_ag.ace
cat tmp/OR/d5.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/OR/d5.$MAGIC.acceptors.next_ag.ace

if (-e  tmp/OR/d5.$MAGIC.captured_genes.txt) \rm  tmp/OR/d5.$MAGIC.captured_genes.txt

bin/tacembly GeneIndexDB << EOF
  read-models
  parse tmp/OR/d5.$MAGIC.donors.next_gt.ace
  parse tmp/OR/d5.$MAGIC.acceptors.next_ag.ace 
  save
  quit
EOF

### CAPTURE
if (-e tmp/METADATA/$MAGIC.captured_genes.ace) then 
  bin/tacembly GeneIndexDB << EOF
    parse  tmp/METADATA/$MAGIC.captured_genes.ace
    save
    bql -o  tmp/OR/d5.$MAGIC.captured_introns.txt select ii,c from g in ?Gene, c in g->capture where c, ii in g->Intron  
    quit
EOF

  cat tmp/OR/d5.$MAGIC.captured_introns.txt | gawk -F '\t' '{if($1 != old)printf ("\nIntron %s\n",$1);old=$1;printf("Capture %s\n",$2);}END{printf("\n");}' > tmp/OR/d5.$MAGIC.captured_introns.ace
  bin/tacembly GeneIndexDB << EOF
    read-models
    parse tmp/OR/d5.$MAGIC.captured_introns.ace
    save
    quit
EOF
endif
### END CAPTURE


if (-e  RESULTS/Expression/AceFiles/$project.introns.INTRON.u.ace.gz2) then
  set toto=GeneIndexDB/$project.intron_confirmation.ace
  echo $toto
  if (! -e $toto) then
    if (-e  RESULTS/Expression/AceFiles/$project.introns.INTRON.u.ace.gz2) then
      zcat RESULTS/Expression/AceFiles/$project.introns.INTRON.u.ace.gz | gawk '/^Intron/{z=$0}/_SumOfAllReadsInProject/{if ($4 >0) printf("%s\nRun_U %s %s %s %s %s %s\n\n",z,magic,$3,$4,$5,$6,$7);}' magic=$project > $toto
      echo "pparse $toto" | bin/tacembly GeneIndexDB -no_prompt
    else
      echo "Missing file  RESULTS/Expression/AceFiles/$project.introns.INTRON.u.ace.gz"
    endif
  endif
  goto phaseLoop
endif

touch GeneIndexDB/d5.cumul.done

phaseLoop:
  echo d5.intronDB phase $phase  done
