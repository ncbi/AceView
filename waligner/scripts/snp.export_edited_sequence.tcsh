#!bin/tcsh -ef

set zone=$1

# 4, 3, 20 are also the default values in the code

set minSnpCount=$2
set minSnpCover=$3
set minSnpFrequency=$4

  set n=0
  if (-e tmp/SNP/$zone.snp.list1) \rm tmp/SNP/$zone.snp.list1

  foreach run (`cat MetaDB/$MAGIC/RunList`)
    if (-e  tmp/SNP/$run/$zone.detect.u.txt.gz) then
      gunzip -c tmp/SNP/$run/$zone.detect.u.txt.gz | gawk '/^#/{next}/Incompatible_strands/{next;}{if($9>=100*minCover && $10>=100*minMutant && $8>=minFrequency)printf("%s:%s%s\n",$1,$2,$3);}' minCover=$minSnpCover minMutant=$minSnpCount  minFrequency=$minSnpFrequency  >> tmp/SNP/$zone.u.snp.list1
    endif
    set n=`cat Fastc/$run/Max_probe_length | gawk '{if(x<1)x=n;if ($1>200)$1=1 ;if($1>n)x=$1;}END{print x}' n=$n`
  end 

# tmp/VariantDB.zoneg.1/Variant.compare
  cat TARGET/Targets/$species.*.mask.txt ZZZZZ  tmp/SNP/$zone.u.snp.list1 | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){n++;if(bad1[$1]>0){bad2[$1]=n;}else{bad1[$1]=n;bad2[$1]=n;}gg[n]=$1;a1[n]=$2;a2[n]=$3;next;} split($1,aa,"[:A-Z]");g=aa[1];x=aa[2]; n1=bad1[g];n2=bad2[g];if(n1>0){for(i=n1;i<=n2;i++)if(gg[i]==g && x>=a1[i] && x<=a2[i])next;}print}' >  tmp/SNP/$zone.u.snp.list2

  cat tmp/SNP/$zone.u.snp.list2  | sed -e 's/>/2/' |  sort -u > tmp/SNP/$zone.u.snp.list
  \rm tmp/SNP/$zone.u.snp.list[12]

  set toto=tmp/SNP_ZONE/$zone.fasta.gz
  if (! -e $toto) set toto=TARGET/Targets/$species.$zone.fasta.gz
  if (-e $toto) then
    echo "bin/snp -editSequence $n -fasta $toto -snp_list tmp/SNP/$zone.u.snp.list -newIntrons tmp/Transloc/$MAGIC.transloc.txt   -gzo -o tmp/SNP/$zone.edited_sequence -run $zone"
    bin/snp -editSequence $n -fasta $toto -snp_list tmp/SNP/$zone.u.snp.list -newIntrons tmp/Transloc/$MAGIC.transloc.txt   -gzo -o tmp/SNP/$zone.edited_sequence -run $zone
  endif

exit

