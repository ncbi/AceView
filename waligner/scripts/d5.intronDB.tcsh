#!bin/tcsh -f

set phase=$1
set project=$2
set chrom=$3

echo "phase=$phase  project=$project"

if ($phase == chromCumul) then
  set toto = tmp/introns/d5.$MAGIC.de_uno.$chrom
  echo ' ' > $toto.1
  echo ' ' > $toto.txt
  foreach run (`cat MetaDB/$MAGIC/GroupIntronList`)
    set ff=tmp/OR/$run/d4.de_uno.txt.gz
    if (! -e $ff) set ff=tmp/OR/$run/d5.de_uno.txt.gz
    if (-e  $ff) then
      set minS=`cat  tmp/introns/d5.$MAGIC.minS | gawk -F '\t' 'BEGIN{n=1;}{if($1==run)n=$2+0;}END{print n}' run=$run`
      zcat $ff | gawk -F '\t' "/^$chrom\t/"'{a1=$2+0;a2=$3+0;ii=$1 "__" a1 "_" a2 ; if($7=="known" || ($5 == "gt_ag" && $4>=minS))printf ("%s\n",ii);}' run=$run minS=$minS chrom=$chrom >> $toto.1
    endif
  end
  cat $toto.1 | cut -f 1 | sort -u | gzip > $toto.2.gz
  foreach run (`cat MetaDB/$MAGIC/GroupIntronList`)
    set ff=tmp/OR/$run/d4.de_uno.txt.gz
    if (! -e $ff) set ff=tmp/OR/$run/d5.de_uno.txt.gz
    if (-e  $ff) then
      zcat $toto.2.gz ZZZZZ.gz $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}'"/^$chrom\t/"'{a1=$2+0;a2=$3+0;ii=$1 "__" a1 "_" a2 ; if(ok[ii]==1)printf ("%s\t%s\t%s\t%s\t%d\n",ii,$5,"any",run,$4);}' run=$run minS=$minS chrom=$chrom >> $toto.txt
    endif
  end

  cat $toto.txt | sort -V | gzip > $toto.txt.gz
  \rm  $toto.txt  $toto.1 $toto.2.gz
  zcat $toto.txt.gz  | gawk -F '\t' '{ii=$1; if($5+0<1)next;if(ii!=old){if (n>0)printf("RNA_seq %d\n",n);n=0;printf("\nIntron \"%s\"\n",ii);}old=ii;printf("de_uno %s %d\n",$4,$5);n+=$5;}END{if (n>0)printf("RNA_seq %d\n",n);printf("\n");}' | gzip > $toto.ace.gz
  
  goto phaseLoop
endif

if ($phase == parseCumul) then
  if (-e  tmp/introns/_r.d5) \rm  tmp/introns/_r.d5
  foreach chrom ($chromSetAll)
    echo "pparse tmp/introns/d5.$MAGIC.de_uno.$chrom.ace.gz" >> tmp/introns/_r.d5
  end
  echo 'save' >> tmp/introns/_r.d5
  echo 'quit' >> tmp/introns/_r.d5
  
  cat  tmp/introns/_r.d5 | bin/tacembly GeneIndexDB -noprompt

  bin/tacembly GeneIndexDB << EOF
    bql -o tmp/introns/d5.all_introns select ii from ii in ?Intron where (!ii#length || ! ii#IntMap)
    quit
EOF

  cat  tmp/introns/d5.all_introns | gawk '{split($1,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){b1=a1;b2=a2;}else{b1=a2;b2=a1;}printf("Intron %s\nIntMap %s %d %d\nLength %d\n\n",$1,chrom,a1,a2,b2-b1+1);}' | gzip >  tmp/introns/d5.all_introns.intmap.ace.gz

  echo "pparse tmp/introns/d5.all_introns.intmap.ace.gz" | bin/tacembly GeneIndexDB -no_prompt

  touch  tmp/introns/d5.$MAGIC.parse.done

  if (! -e GeneIndexDB/d5.parseGenome.done) then
    echo "pparse TARGET/Targets/$species.genome.fasta.gz" | bin/tacembly GeneIndexDB -no_prompt
    touch GeneIndexDB/d5.parseGenome.done
  endif

  bin/tacembly GeneIndexDB <<EOF
    query find intron IS 22_* &&  ! type
    spush
    bql -o GeneIndexDB/d5.intron_feet.R.txt select i,m,s,x,y ,f1,f2 from i in @,m in i->intmap, s in ?Sequence where s==m,x in m[1],y in m[2] where y<x, f1 in DNA(s,x,x-1), f2 in DNA(s,y+1,y)
    sxor
    bql -o GeneIndexDB/d5.intron_feet.F.txt select i,m,s,x,y ,f1,f2 from i in @,m in i->intmap, s in ?Sequence where s==m,x in m[1],y in m[2] where x<y, f1 in DNA(s,x,x+1), f2 in DNA(s,y-1,y)
    quit
EOF

EOF
  goto phaseLoop
endif

if ($phase == donorAcceptor) then
# donor acceptor
  bin/tacembly GeneIndexDB << EOF
    find intron
    show -a -f tmp/introns/d5.$MAGIC.any_intron.ace
    quit
EOF

  cat tmp/introns/d5.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nDonor %s__%d_%s\n",chrom, a1,s);printf("IntMap %s %d %d\n",chrom, a1-ds,a1);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("gt\n");next;}/^gc_ag/{printf("gc\n");next;}/^ct_ac/{printf("ct\n");next;}/^at_ac/{printf("at\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/introns/d5.$MAGIC.any_donor.ace

  cat tmp/introns/d5.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nAcceptor %s__%d_%s\n",chrom, a2,s);printf("IntMap %s %d %d\n",chrom, a2,a2+ds);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("ag\n");next;}/^gc_ag/{printf("ag\n");next;}/^ct_ac/{printf("ac\n");next;}/^at_ac/{printf("ac\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/introns/d5.$MAGIC.any_acceptor.ace

  bin/tacembly GeneIndexDB << EOF
    pparse tmp/introns/d5.$MAGIC.any_donor.ace
    pparse tmp/introns/d5.$MAGIC.any_acceptor.ace
    save
    bql -o tmp/introns/d5.$MAGIC.donors select d,t from d in ?Donor, t in d#gt
    bql -o tmp/introns/d5.$MAGIC.acceptors select a,t from a in ?Acceptor, t in a#ag
    quit
EOF
  cat tmp/introns/d5.$MAGIC.donors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/introns/d5.$MAGIC.donors.next_gt.ace
  cat tmp/introns/d5.$MAGIC.donors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/d5.$MAGIC.donors.next_gt.ace
  cat tmp/introns/d5.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/d5.$MAGIC.donors.next_gt.ace
  cat tmp/introns/d5.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/introns/d5.$MAGIC.donors.next_gt.ace

  cat tmp/introns/d5.$MAGIC.acceptors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/introns/d5.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/d5.$MAGIC.acceptors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/d5.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/d5.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/d5.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/d5.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/introns/d5.$MAGIC.acceptors.next_ag.ace

  if (-e  tmp/introns/d5.$MAGIC.captured_genes.txt) \rm  tmp/introns/d5.$MAGIC.captured_genes.txt

  bin/tacembly GeneIndexDB << EOF
    read-models
    parse tmp/introns/d5.$MAGIC.donors.next_gt.ace
    parse tmp/introns/d5.$MAGIC.acceptors.next_ag.ace 
    save
    quit
EOF
  touch tmp/introns/d5.$MAGIC.donorAcceptor.done
  goto phaseLoop
endif

### CAPTURE
if ($phase == capture) then
  if (-e tmp/METADATA/$MAGIC.av.captured_genes.ace) then 
    bin/tacembly GeneIndexDB << EOF
      parse  tmp/METADATA/$MAGIC.av.captured_genes.ace
      save
      bql -o  tmp/introns/d5.$MAGIC.captured_introns.txt select ii,c from g in ?Gene, c in g->capture where c, ii in g->Intron  
      quit
EOF

    cat tmp/introns/d5.$MAGIC.captured_introns.txt | gawk -F '\t' '{if($1 != old)printf ("\nIntron %s\n",$1);old=$1;printf("Capture %s\n",$2);}END{printf("\n");}' > tmp/introns/d5.$MAGIC.captured_introns.ace
    bin/tacembly GeneIndexDB << EOF
      read-models
      parse tmp/introns/d5.$MAGIC.captured_introns.ace
      save
      quit
EOF
  endif
  ### END CAPTURE
  touch tmp/introns/d5.$MAGIC.capture.done
  goto phaseLoop
endif


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
