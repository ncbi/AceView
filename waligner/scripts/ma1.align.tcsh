#!bin/tcsh -f

set pp=$1
set isAffy=`echo $pp | gawk '{i=index($1,"Affy");print i;}'`
  if (! -e  MicroArray/Oligos/Fasta/$pp.count) then
    bin/dna2dna -i MicroArray/Oligos/Fasta/$pp.fasta.gz -I fasta -O count -o MicroArray/Oligos/Fasta/$pp
  endif

# we coud also use target EBI, midway between av and RefSeq
  foreach target (av RefSeq genome)
    if (! -e TARGET/Targets/$species.$target.fasta.gz) continue
    source scripts/target2target_class.txt
    set target2=$target
    if ($species == rn && $target == RefSeq) set target2=RefSeq_NX_MR_2013_05_27
    if (! -e  MicroArray/Hits/$pp.$target2.hits.gz2) then
      set strand=""
      if ($target != genome && $pp =~ AGL* ) set strand="-strand"

      set mm="-minAli 35"
      if ($pp =~ Affy* ) set mm="-errRateMax 10"   # Affymetrix
      echo "clipalign -i MicroArray/Oligos/Fasta/$pp.fasta.gz $mm -seedLength 16 -seedOffset 1 -seedShift 5 -maxHit 30 -t TARGET/Targets/$species.$target2.fasta.gz -o MicroArray/Hits/$pp.$target -gzo -target_class $target_class $strand"
      bin/clipalign -i MicroArray/Oligos/Fasta/$pp.fasta.gz $mm -seedLength 16 -seedOffset 1 -seedShift 5 -maxHit 30 -t TARGET/Targets/$species.$target2.fasta.gz -o MicroArray/Hits/$pp.$target -gzo -target_class $target_class $strand
    endif
  end

  gunzip -c MicroArray/Hits/$pp.*.hits.gz | bin/bestali -exportBest -exportVenn -gzo -o MicroArray/Hits/best.$pp.EBI
   
# dedicated hack 2012_06_21
if (0) then
    bin/clipalign -i MicroArray/Oligos/Fasta/$pp.fasta.gz -errRateMax 5 -seedLength 16 -seedOffset 1 -seedShift 5 -maxHit 10 -t TARGET/Targets/RefSeq.pg.2012.fasta.gz -o MicroArray/Hits/$pp.RefSeq.pg.2012 -gzo -target_class pg -strand
     bin/clipalign -i MicroArray/Oligos/Fasta/$pp.fasta.gz  -errRateMax 5 -seedLength 16 -seedOffset 1 -seedShift 5 -maxHit 10 -t TARGET/Targets/rn.magic_june21.fasta.gz -o MicroArray/Hits/$pp.magic.2012_06_21 -gzo -target_class magic -strand


# locate the zones assayed by Affy, the probeset hits inside 100kb and above 30bp and is touching at most 9 sites
# 2013_01_13, we add the conditions nn[z] > 7  and $5<50, to remove single probe hits
# 2013_05_24, we remove the limit at 8 probes but we export nn[z]

 gunzip -c best.Affy.Rat230_2.hits.gz | gawk '/^#/{next}/genome/{if($10>=10)next;i=index($1,"_at");printf("%s\t%s\t%d\t%d\n", substr($1,1,i+2),$11,$12,$13);}' | sort -k 1,1 -k 2,2 -k 3,3n | gawk '{a1=$3;a2=$4;s="+";if(a1>a2){a0=a1;a1=a2;a2=a0;s="-";}z=$1 "\t" s $2;if(! aa1[z] || aa1[z] > a1)aa1[z]=a1; if(! aa2[z] || aa2[z] < a2)aa2[z]=a2;nn[z]++;}END{for (z in aa1) if(nn[z]>0) printf("%s\t%d\t%d\t%d\t%d\n",z,aa1[z],aa2[z], aa2[z] - aa1[z]+1,nn[z]);}' | sort | gawk '{n++;p=$1;z=$0;z2p[z]=p;nn[p]++;}END{for (z in z2p)printf("%d\t%s\n",nn[z2p[z]],z);}' | sort -k 2,2 -k 7,7nr | gawk -F '\t' '{if($2==p && $7<n)next;p=$2;n=$7;print}' | sort -k 1n > $pp.probeset2chrom.txt

# 2013_05_24, limit now to probeset hitting < 10 chromosomes over a region < 10M
cat $pp.probeset2chrom.txt | gawk '{if($1<10 && $6 > 20 && $6 < 10000000)print}' |  sort | cut -f 2,3,4,5,6,7 | gawk '{n++;p=$1;z=$0;z2p[z]=p;nn[p]++;mm[p]=$7;}END{for (z in z2p)printf("%d\t%s\t%d\n",nn[z2p[z]],z,mm[z2p[z]]);}' >  $pp.probeset2chrom.compact.txt


\rm $pp.probeset2chrom.compact.txt.gz
gzip  $pp.probeset2chrom.compact.txt
# count the probes actually part of each compact zone
gunzip -c  $pp.probeset2chrom.compact.txt.gz ZZZZZ.gz  best.Affy.Rat230_2.hits.gz  | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz<1){ps=$2;s=substr($3,1,1);c=substr($3,2);a1=$4;a2=$5;np[ps]++;n=np[ps];aac[ps,n]=c;aas[ps,n]=s;aa1[ps,n]=a1;aa2[ps,n]=a2;next;}}/genome/{p=$1;i=index($1,"_at");ps=substr($1,1,i+2);if(0)print ps,np[ps] ;c=$11;a1=$12;a2=$13;s="+";if(a1>a2){a0=a1;a1=a2;a2=a0;s="-";}for(n=1;n<=np[ps];n++){if(0)print a1,aa1[ps,n];if(c==aac[ps,n] && s==aas[ps,n] && a1>=aa1[ps,n] && a2 <= aa2[ps,n]){okn[ps,n]++;if(index(ok[ps,n],$1 ";")<1)oknu[ps,n]++;ok[ps,n] = ok[ps,n] $1 ";" ;}if(0)print ok[ps,n];}}END{for (ps in np) for(n = 1 ; n <= np[ps];n++){a1=aa1[ps,n];a2=aa2[ps,n];s=aas[ps,n];gsub(/>/,"",ok[ps,n]);gsub(ps,"",ok[ps,n]);printf("Affy_%s\t%s\t%d\t%d\t%d\t%d\t%s\n",ps,aac[ps,n],a1,a2,oknu[ps,n],okn[ps,n],ok[ps,n]);}}' | sort >  $pp.probeset2chrom.compact.count_probes.txt


# construct the coordinates of the probes in the probeset
gunzip -c  $pp.probeset2chrom.compact.txt.gz ZZZZZ.gz  best.Affy.Rat230_2.hits.gz  | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz<1){ps=$2;s=substr($3,1,1);c=substr($3,2);a1=$4;a2=$5;np[ps]++;n=np[ps];aac[ps,n]=c;aas[ps,n]=s;aa1[ps,n]=a1;aa2[ps,n]=a2;next;} p=$1;i=index($1,"_at");ps=substr($1,1,i+2);if(0)print ps,np[ps] ;c=$11;a1=$12;a2=$13;s="+";if(a1>a2){a0=a1;a1=a2;a2=a0;s="-";}for(n=1;n<=np[ps];n++){if(0)print a1,aa1[ps,n];if(c==aac[ps,n] && s==aas[ps,n] && a1>=aa1[ps,n] && a2 <= aa2[ps,n])ok[ps,n]=ok[ps,n] ":" a1 "/" a2;if(0)print ok[ps,n];}}END{for (ps in np) for(n = 1 ; n <= np[ps];n++){n1=split(ok[ps,n],aa,":");if(0)print n1,ok[ps,n];if(n1>4)kk[ps]++;}for (ps in np){k=0; for(n = 1 ; n <= np[ps];n++){n1=split(ok[ps,n],aa,":");if(n1>4){ k++; su="";if(kk[ps]>1)su="#"k"/"kk[ps];a1=aa1[ps,n];a2=aa2[ps,n];if(aas[ps,n]=="-"){a0=a1;a1=a2;a2=a0;}printf("Sequence %s\tSubsequence \"Affy_%s%s\" %d %d\n\n",aac[ps,n],ps,su,a1,a2);ii=split(ok[ps,n],aa,":");for(i=2;i<=ii;i++){split(aa[i],bb,"/");if(aas[ps,n]=="+"){b1=bb[1]-a1+1;b2=bb[2]-a1+1;}else{b1=a1-bb[2]+1;b2=a1-bb[1]+1;}printf("Exon\tAffy_%s%s\t%09d\t%09d\n",ps,su,b1,b2);} }}}}' >  $pp.probeset2chrom.compact.counts.txt1

# transform into  Affy__pg*
cat $pp.probeset2chrom.compact.counts.txt1 | grep Sequence | sort | gawk '/^Sequence/{s=$2;if(s!=old)printf("\nSequence %s\n",s);old=s;printf("Subsequence %s %d %d\n",$4,$5,$6);;}END{printf("\n");}' > $pp.probeset2chrom.compact.counts.ace
cat $pp.probeset2chrom.compact.counts.ace | gawk '/^Sequence/{seq=$2;next;}/^Subseq/{p=$2;gsub(/\"/,"",p);gsub(/Affy_/,"",p);i=index(p,"#");if(i>1)p=substr(p,1,i-1);printf("Probeset %s\nIntMap %s %d %d \n\n",p,seq,$3,$4);}' > $pp.probeset2intmap.ace
cat $pp.probeset2chrom.compact.counts.txt1 | grep Exon | sort | gawk '/^Exon/{s=$2;if(s!=old)printf("\nSequence \"%s\"\nIs_predicted_gene\nMethod Affy\n",s);old=s;printf("Source_exons %d %d\n",0+$3,0+$4);}END{printf("\n");}' >> $pp.probeset2chrom.compact.counts.ace

# there are 8111 TAQ, 8080 in best.TAQ

   foreach target (av RefSeq)
     if ($isAffy >= 1) then
     gunzip -c best.Affy.Rat230_2.hits.gz | gawk '/^#/{next}{print $1}' | sort -u | gawk '{i=index($1,"_at");print substr($1,1,i);}' | sort | gawk '{n[$1]++;}END{for(k in n)if (n[k]>=8)print k}' | sort > well_mapped_probeSet.txt

     gunzip -c MicroArray/Hits/best.$pp.hits.gz | gawk -F '\t' '/^#/{next;}{if (index($8,target)<3)next; if($9=="-")$9=$11;p=$1;i=index(p,"_at");ps=substr(p,1,i+2);printf("%s\t%s\t%s\n",ps,p,$9);}' target=$target | sort -u > MicroArray/Hits/$pp.$target.probe2gene.txt
        cat MicroArray/Hits/$pp.$target.probe2gene.txt | gawk '{nn[$1 "\t" $3]++}END{for(k in nn)printf ("%s\t%d\n",k,nn[k]);}' | sort -u >  MicroArray/Hits/$pp.$target.probeset2gene.txt 
 
    gunzip -c MicroArray/Hits/$pp.$target.hits.gz | gawk -F '\t' '/^#/{next;}{if (index($8,target)<3)next; if(1 ||  $9=="-")$9=$11;p=$1;i=index(p,"_at");ps=substr(p,1,i+2);printf("%s\t%s\t%s\n",ps,p,$9);}' target=$target | sort -u > MicroArray/Hits/$pp.$target.probe2mrna.txt
       cat MicroArray/Hits/$pp.$target.probe2mrna.txt | gawk '{nn[$1 "\t" $3]++}END{for(k in nn)printf ("%s\t%d\n",k,nn[k]);}' | sort -u >  MicroArray/Hits/$pp.$target.probeset2mrna.txt 


       cat MicroArray/Hits/$pp.$target.probe2gene.txt | gawk '{nn[$2]++;g[$2]=g[$2] ":" $3;}END{for(k in nn){u[nn[k]]++;}for (k in u)print k,u[k];}' | sort -k 1n > MicroArray/Hits/$pp.$target.probe2gene.multiplicity.txt 
       cat MicroArray/Hits/$pp.$target.probeset2gene.txt | gawk '{nn[$1]++;g[$1]=g[$1] ":" $3;}END{for(k in nn){u[nn[k]]++;}for (k in u)print k,u[k];}' | sort -k 1n > MicroArray/Hits/$pp.$target.probeset2gene.multiplicity.txt 

       cat  MicroArray/Hits/$pp.$target.probeset2gene.txt | gawk -F '\t' '{n[$3]++}END{for(k in n)printf ("%s\t%d\n",k,n[k]);}' | sort -k 1n > MicroArray/Hits/$pp.$target.probeset2gene.mapped_probe_per_probeset.txt

# we require 8 probes per probeset matching the gene
       cat MicroArray/Hits/$pp.$target.probeset2gene.txt | gawk '{nn[$1]++;g[$1]=g[$1] ":" $2 "\t" $3 ;}END{for(k in nn){if (nn[k]==1)printf("%s\t%s\n",k,substr(g[k],2));if(nn[k]==2){split(g[k],aa,":");g1=aa[2];g2=aa[3];if (g1<"ZZ" && g2 >= "a")printf("%s\t%s\n",k,g1);if (g2<"ZZ" && g1 >= "a")printf("%s\t%s\n",k,g2);}}}' | gawk -F '\t' '{if($3 >= 8)print}' > MicroArray/Hits/$pp.$target.probeset2gene.unique.txt
       cat MicroArray/Hits/$pp.$target.probeset2gene.txt | gawk '{if(substr($2,1,3)!="NM_" && substr($2,1,3)!="XM_" && substr($2,1,3)!="NR_"&& substr($2,1,3)!="XR_")nn[$1]++;nn[$1]=0+nn[$1];g[$1]=g[$1] ":" $2 "\t" $3 ;}END{for(k in nn){if (nn[k]<=5){n1=split(g[k],aa,":");for(i=2;i<=n1;i++)printf("%s\t%s\n",k,aa[i]);}}}' | gawk -F '\t' '{if($3 >= 8)print}' > MicroArray/Hits/$pp.$target.probeset2gene.quasi_unique.txt


       cat MicroArray/Hits/$pp.$target.probeset2gene.unique.txt | gawk '{if($2 == "-")next; nm=$2;nm1=substr(nm,1,3); if(nm1== "NM_" || nm1== "NR_" || nm1== "XM_" || nm1== "XR_") next; printf ("Probeset %s\nGene %s %d\n\n",$1,$2,$3);}' >  MicroArray/Hits/$pp.$target.probeset2gene.unique.ace
       cat MicroArray/Hits/$pp.$target.probeset2gene.quasi_unique.txt | gawk '{if($2 == "-")next;nm=$2;nm1=substr(nm,1,3); if(nm1== "NM_" || nm1== "NR_" || nm1== "XM_" || nm1== "XR_") next;  printf ("Probeset %s\nGene %s %d \n\n",$1,$2,$3);}' >  MicroArray/Hits/$pp.$target.probeset2gene.quasi_unique.ace
        cat MicroArray/Hits/$pp.$target.probeset2gene.quasi_unique.txt | gawk '{nm=$2;nm1=substr(nm,1,3); if(nm1== "NM_" || nm1== "NR_" || nm1== "XM_" || nm1== "XR_") {z=$2;i=index(z,".");printf ("Probeset %s\nNM_id %s %d\n\n",$1,substr(z,1,i-1),$3);}}' >  MicroArray/Hits/$pp.$target.probeset2nmid.quasi_unique.ace
     else
       gunzip -c MicroArray/Hits/best.$pp.hits.gz | gawk -F '\t' '/^#/{next;}{if ($8 == tgcl)printf("%s\t%s\n",$1,$9);}'  tgcl=$target_class  | sort -u | sed -e 's/>//' > MicroArray/Hits/$pp.$target.probe2gene.txt

# accept case -2: 2 genes en antisense et garder le premier par ordre alphabetique unix, donc celui en majuscules
#  sort | gawk '{if($1==old)next;old=$1;print;}'
       gunzip -c MicroArray/Hits/best.$pp.hits.gz | gawk -F '\t' '/^#/{next;}{if ($8 == tgcl && $10==1)printf("%s\t%s\n",$1 ,$9);}'   tgcl=$target_class | sort -u | sed -e 's/>//' | sort > MicroArray/Hits/$pp.$target.probe2gene.unique.txt
 
       gunzip -c MicroArray/Hits/best.$pp.hits.gz | gawk -F '\t' '/^#/{next;}{if($8==tgcl && $10 ==1)printf("%s\t%s\n",$1,$11);}' tgcl=$target_class | sort -u | sed -e 's/>//' | sort > MicroArray/Hits/$pp.$target.probe2mrna.unique.txt

if (0) then
  # comptage des probes stranded or not
  set ss=ns
  set ss=s

  if ($ss == s)  set toto=RESULTS/Probe_fate.stranded.txt
  if ($ss == ns)  set toto=RESULTS/Probe_fate.non_stranded.txt
  echo -n "# " > $toto
  date >> $toto

    set ff=MicroArray/Hits/best.AGLuK.hits.gz 

    echo -n "Probes uniquely mapping in an AceView Gene without mismatch in sense, accepted " >> $toto
    gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 == 60 && $12<$13) print $1;}' | sort -u | wc | gawk '{printf("\t%d\n",$1);}'  >> $toto

    if ($ss == ns) then
      echo -n "Probes uniquely mapping in an AceView Gene without mismatch in antisense, accepted " >> $toto
      gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 == 60 && $12>$13) print $1;}' | sort -u | wc | gawk '{printf("\t%d\n",$1);}'  >> $toto
    endif

    echo -n "Probes uniquely mapping in an AceView Gene over their entire length with mismatch, accepted" >> $toto
    gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 < 60 && $5 == 60) print $1;}' | sort -u | wc | gawk '{printf("\t%d\n",$1);}'  >> $toto
    echo -n "Probes uniquely mapping in an AceView Gene over at least 35 bases, accepted" >> $toto
    gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 < 60 && $5 < 60) print $1;}' | sort -u | wc | gawk '{printf("\t%d\n",$1);}'  >> $toto

    if ($ss == ns) then
      echo -n "Probes mapping at a single locus but touching 2 genes antisense to each other, rejected" >> $toto
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av" && $10 == -2) print $1;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto
    endif

    echo -n "Probes mapping on transcripts but not on genome, hance accross an exon juction" >> $toto
     gunzip -c $ff | gawk -F '\t' '{ if($8 == "Z_genome")gg[$1]=1;if($8 == "ET_av" && $10 == 1 && $2 == 60 && $12<$13) ok[$1]=1;}END{for(k in ok)if(gg[k]<1)nn++;}END{printf("\t\t\t%d\n",nn);}'


    gunzip -c $ff | gawk -F '\t' '{ if($8 == "Z_genome")gg[$1]=1;if($8 == "ET_av" && $10 == 1 && $2 <= 60 && $12<$13) ok[$1]=1;}END{for(k in ok)if(gg[k]<1)ns++;}END{print ns}'
4077

    if (0 && $ss == s) then   
      # this case does not occur, since in standed mapping the probe is atributed to a single gene
      echo -n "Probes mapping at a single locus but touching 2 genes antisense to each other, accepted" >> $toto
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av" && $10 == -2) print $1;}' | sort -u | wc | gawk '{printf("\t%d\n",$1);}'  >> $toto
    endif

    echo -n "Probes mapping in several distinct genes, rejected" >> $toto
    gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av" && $10 > 1) print $1;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto
    echo -n "Probes mapping uniquely on the genome, outside genes, without mismatch, rejected " >> $toto
    gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2==60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto
    echo -n "Probes mapping uniquely on the genome, outside genes, over their entire length, with mismatch, rejected " >> $toto
    gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<60 && $5 == 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto
    echo -n "Probes mapping uniquely on the genome, outside genes, at least 35 bases, with mismatch, rejected " >> $toto
    gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<60 && $5 < 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto
    echo -n "Probes mapping in 2 to 9 positions on the genome, outside any gene, rejected " >> $toto
    gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 > 1)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p ;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  >> $toto

    cat $toto | gawk -F '\t' '/accepted/{na+=$2;}/rejected/{nr+=$3}{nn+=$2+$3;}END{printf("Total\t%d\t%d\nUnmapped %d\tTotal %d\n",na,nr,43291-nn,43291);}' > $toto.1
    cat $toto.1 >> $toto
    \rm $toto.1

    mv $toto $toto.1 
    cat $toto.1 | gawk -F '\t' '{printf("%s\t%s\t%s\n",$2,$3,$1);}' > $toto
    \rm $toto.1
endif

  if (0) then
  # hand debugging
      echo -n "Probes mapping uniquely on the genome, outside genes, without mismatch, rejected "
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2==60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  
      echo -n "Probes mapping uniquely on the genome, outside genes, over their entire length, with mismatch, rejected "
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<60 && $5 == 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}' 
      echo -n "Probes mapping uniquely on the genome, outside genes, at least 35 bases, with mismatch, rejected " 
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<60 && $5 < 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'  

      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av" || $8 == "KT_RefSeq") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<=60 && $5 <= 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'
      gunzip -c $ff | gawk -F '\t' '{  if($8 == "ET_av") av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<=60 && $5 <= 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'
       gunzip -c MicroArray/Hits/best.$pp.EBI.hits.gz  | gawk -F '\t' '{  if($8 == "ET_av" || $8 == "KT_RefSeq" || ($10==1 && $8 == "MT_EBI")) av[$1]=1; if($8 == "Z_genome" && $10 == 1 && $2<=60 && $5 <= 60)gg[$1]=1;}END{for(p in gg)if(av[p]<1)print p;}' | sort -u | wc | gawk '{printf("\t\t%d\n",$1);}'

  endif


       cat MicroArray/Hits/$pp.$target.probe2gene.txt | gawk '{nn[$1]++;g[$1]=g[$1] ":" $2;}END{for(k in nn){u[nn[k]]++;}for (k in u)print k,u[k];}' | sort -k 1n > MicroArray/Hits/$pp.$target.probe2gene.multiplicity.txt 
     
       cat MicroArray/Hits/$pp.$target.probe2gene.txt | gawk '{nn[$1]++;g[$1]=g[$1] ":" $2;}END{for(k in nn){if (nn[k]==1)printf("%s\t%s\n",k,substr(g[k],2));if(nn[k]==2){split(g[k],aa,":");g1=aa[2];g2=aa[3];if (g1<"ZZ" && g2 >= "a")printf("%s\t%s\n",k,g1);if (g2<"ZZ" && g1 >= "a")printf("%s\t%s\n",k,g2);}}}' > MicroArray/Hits/$pp.$target.probe2gene.unique.txt
       cat MicroArray/Hits/$pp.$target.probe2gene.txt | gawk '{if(substr($2,1,3)!="NM_" && substr($2,1,3)!="XM_" && substr($2,1,3)!="NR_"&& substr($2,1,3)!="XR_")nn[$1]++;nn[$1]=0+nn[$1];g[$1]=g[$1] ":" $2;}END{for(k in nn){if (nn[k]<=5){n1=split(g[k],aa,":");for(i=2;i<=n1;i++)printf("%s\t%s\n",k,aa[i]);}}}' > MicroArray/Hits/$pp.$target.probe2gene.quasi_unique.txt

       cat MicroArray/Hits/$pp.$target.probe2gene.unique.txt | gawk '{printf ("Gene %s\nMicroArray %s\n\n",$2,$1);}' >  MicroArray/Hits/$pp.$target.probe2gene.unique.ace
       cat MicroArray/Hits/$pp.$target.probe2gene.quasi_unique.txt | gawk '{printf ("Gene %s\nMicroArray %s\n\n",$2,$1);}' >  MicroArray/Hits/$pp.$target.probe2gene.quasi_unique.ace
 
 # map the yet unassigned probes to genes via the extent of the geneboxes, respecting the strand
       cat MicroArray/Hits/$pp.av.probe2gene.quasi_unique.txt ZZZZZ TARGET/GENES/av.gene2intMap.txt ZZZZZ MicroArray/Hits/$pp.probe2chrom.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{if(zz<2){s="+";if(0+$3>0+$4){s="-";x=$3;$3=$4;$4=x;}printf("%s\t%s\t%d\t%d\tGENE\t%s\n",$2,s,$3,$4,$1); next;}}{if(ok[$1]==1)next;gsub(/-/,":",$2);split($2,aa,":");s="+";if(0+aa[2]>0+aa[3]){s="-";x=aa[2];aa[2]=aa[3];aa[3]=x;}printf("%s\t%s\t%d\t%d\tPROBE\t%s\n",substr(aa[1],4),s,aa[2],aa[3],$1);}' | sort -k 1,1 -k 2,2 -k 3,3n -k 4,4n |  gawk -F '\t' '{if(c!= $1 || s != $2){c=$1;s=$2;g1=0;g2=0;p1=0;p2=0;probe=0;gene=0;}}{if($5 == "GENE"){g1=$3;if(g1>g2)gene=$6;g2=$4;if(g1<p2)printf("%s\t%s\n",probe,gene);next;}probe=$6;p1=$3;p2=$4;if(p1<g2)printf("%s\t%s\n",probe,gene);}' | sort | gawk -F '\t' '{if(index(gg[$1],$2)<1)gg[$1]=gg[$1] "," $2;}END{for (g in gg)printf("%s\t%s\n",g,substr(gg[g],2));}'  > MicroArray/Hits/$pp.$target.probe2gene.quasi_unique.supplementary.txt

     endif

   end

if (0) then
  # how many probesete do not map in av but map to the genome
  cat $pp.probeset2chrom.txt | cut -f 2 | sort -u > ps.g
  cat $pp.av.probeset2gene.quasi_unique.txt |  cut -f 1 | sort -u > ps.av
  cat $pp.RefSeq.probeset2gene.quasi_unique.txt |  cut -f 1 | sort -u > ps.rs

 cat ps.av ZZZZZ ps.rs | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.rs_noav
 cat ps.g  ps.av ZZZZZ ps.rs | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.rs_nog
 cat ps.av ZZZZZ ps.g | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.g_noav
 cat ps.rs ZZZZZ ps.g_noav | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])print $1;}' | sort -u > ps.g_noav_rs
 cat ps.rs ZZZZZ ps.g_noav | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.g_noav_nors
 cat ps.rs ZZZZZ ps.av | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])print $1;}' | sort -u > ps.avrs
 cat ps.rs ZZZZZ ps.av | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.av_nors
 cat ps.rs ZZZZZ ps.g_noav | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])print $1;}' | sort -u > ps.g_noav_rs
 cat ps.rs ZZZZZ ps.g_noav | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.g_noav_nors
 cat ps.g ps.rs ps.av ZZZZZ ps.all | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1])next;print $1;}' | sort -u > ps.un

 # verif we must get 31099

 cat ps.un ps.rs_nog ps.g_noav_rs ps.g_noav_nors ps.av_nors ps.avrs | sort -u | wc
 cat ps.un | gawk '{printf("%s\tUnmapped\n",$1);}' > ps.types
 cat ps.rs_nog | gawk '{printf("%s\tRefSeq_not_on_Rn4\n",$1);}' >> ps.types
 cat ps.g_noav_rs | gawk '{printf("%s\tRefSeq_mapped_noAceView\n",$1);}' >> ps.types
 cat ps.g_noav_nors  | gawk '{printf("%s\tCluster_on_Rn4\n",$1);}' >> ps.types
 cat ps.av_nors | gawk '{printf("%s\tAceView_noRefSeq\n",$1);}' >> ps.types
 cat ps.avrs | gawk '{printf("%s\tAceView_RefSeq\n",$1);}' >> ps.types

# create the mapping table in the MicroArray/HITS dir
# we create 4 tables ps.types, ps.map.g ps.map.av ps.map.rs with 2 columns

  cat $pp.probeset2chrom.compact.count_probes.txt    | gawk -F '\t' '{if($5<3)next;m=$2 ":" $3 "-" $4 "(" $5 ")"; ps=$1;gsub(/^Affy_/,"",ps);if(mm[ps]) mm[ps]=mm[ps] ";" ;  mm[ps]=mm[ps] m  ;}END{for (k in mm)printf("%s\t%s\n",k,mm[k]);}' >  ps.map.g
  cat $pp.av.probeset2gene.quasi_unique.txt   | gawk -F '\t' '{m=$2 "(" $3 ")" ; if(mm[$1]) mm[$1]=mm[$1] ";" ;  mm[$1]=mm[$1] m ;}END{for (k in mm)printf("%s\t%s\n",k,mm[k]);}' >  ps.map.av
  cat ps.rs ZZZZZ $pp.RefSeq.probeset2mrna.txt   | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{if(ok[$1]<1)next;m=$2 "(" $3 ")" ; if(mm[$1]) mm[$1]=mm[$1] ";" ;  mm[$1]=mm[$1] m ;}END{for (k in mm)printf("%s\t%s\n",k,mm[k]);}' | sed -e 's/X__//g' >  ps.map.rs

  cat ps.types ZZZZZ ps.map.g ZZZZZ ps.map.av ZZZZZ ps.map.rs | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<4){mmm[$1]=1;mm[0+zz,$1]=$2;next;}}END{printf("#Probeset\tType\tRn4 genome\tAceView 2008\tRefSeq 2012\n");for(k in mmm){printf("%s",k);for(i=0;i<4;i++){u=mm[i,k];if(length(u)<1)u="-";printf("\t%s",u);}printf("\n");}}' | sort > ps.mapping_type.txt

# join the previous table with the AceView gene coordinates and sort by genomic position (disgusting code but nice output
set toto=Affy_AceView.correspondance_and_mapping.txt
echo "# " > $toto
date >>  $toto
echo "#Location on genome\tAffymetrix probeset\tType\tProbeset position on Rn4 genome (number of mapped probes in this region)\tAceView gene 2008 (number of probes mapped to this gene)\tgeneId\tAceView gene position and orientation\tRefSeq-2012 transcript(number of probes mapping to this transcript)\tRefSeq Entrez gene [GeneId]" >>  $toto
cat  ../../TARGET/GENES/av.gene.ace ZZZZZ  ps.mapping_type.txt | gawk '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1)gsub(/\"/,"",$0);}/GeneId/{gid[g]=$2;next;}/^Gene /{g=$2;next;}/^IntMap/{ok[g]=1;map[g]= $2 ":" $3 "-" $4;if(0)print " UUUU " g " map " map[g];next;}{if(zz<1)next;}{printf("%s\t%s\t%s",$1,$2,$3);ng=split($4,ggg,";");if(0)print " ZZ " $4 " xxx " ng ; m="";id="";for(i=1;i<=ng;i++){gg=ggg[i]; if(0)print "YYYY " gg; if(length(gg)>0){split(gg,gga,"(");g=gga[1];ok[g]=2;m1=map[g];if(length(m1)>1){if(length(m)>1)m=m ";";m=m m1}id1=gid[g];if(id1){if(i>1)id=id ";";id=id id1}}}printf("\t%s\t%s\t%s",$4,id,m);printf("\t%s\n",$5);}END{for(g in ok)if(ok[g]==1)printf("-\tAceView_noAffy\t-\t%s\t%s\t%s\t-\n",g,gid[g],map[g]);}' | gawk -F '\t' '{m=$3;if(m=="-")m=$6;split(m,mm,";");split(mm[1],cc,":");split(cc[2],aa,"-");a1=aa[1];a2=aa[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}c="zz";if(cc[1])c=cc[1];if(c=="X")c="XX";if(c=="Y")c="XY";if(length(c)==1)c = "0" c;printf("%s\t%09d\t%09d\t",c,a1,a2);print;}' | sort | gawk -F '\t' '{c=$1;gsub(/^0/,"",c);gsub(/^X/,"",c);gsub(/zz/,"-",c);a1=0+$2;a2=0+$3;if(c=="-" && substr($5,1,7)=="Cluster")$5="Affy_not_well_mapped";if(substr($5,1,7)=="Cluster" && a2-a1>100000)$5="Affy_overspread_on_Rn4";split($6,vv,"(");split(vv[2],vv2,")");np=0+vv2[1];if(np<8 && ( substr($5,1,7)=="Cluster" || $5 == "Affy_over_spread_on_Rn4"))$5="Affy_partly_mapped_on_Rn4";printf("%s:%d-%d",c,a1,a2);if(index($5,"Affy")<1)$5="Affy_" substr($5,1); for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >  $toto.1

# add the rn5 nm->gid correspondance
cat ../../TARGET/MRNAS/rn5.RefSeq2GeneId2Gene.txt ZZZZZ $toto.1 | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{print;}{if(zz<1){nm2gid[$1]=$3;nm2gene[$1]=$4;next;}}{nm8="";n8=split($8,gg8,";");for(i=1;i<=n8;i++){split(gg8[i],ggg8,".");gid8=nm2gid[ggg8[1]];if(gid8 && index(nm8,"[" gid8 "]")<1){if(length(nm8)>0 )nm8=nm8 ";";nm8=nm8 nm2gene[ggg8[1]] "[" gid8 "]" ;}} printf("%s",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);if(nm8=="")nm8="-";printf("\t%s\n",nm8);}' >>  $toto

\rm $toto.1


# \cp $toto ~/rat2012/DeepLeming/RESULTS
# count the categories in the table demultiplexing the ; and clearing the (count-probes)
cat $toto | grep AceView_ | gawk -F '\t' '{if($5=="-")next;gsub(/;/,"\n",$5);print $5;}' | gawk '{split($1,aa,"(");print aa[1];}' | sort -u | wc 


endif

if ($pp == qPCR) then
# complete panic for Biogazelle qPCR, only 8k probes hit the gens, the other 8k hit the genome
# we try to associate them to the genes via their cordinates


endif


# Are the AGLuK probes stranded or not stranded
# check this question by compute the correlation coeff of + and - probes with RNA-seq
if (0) then

# list the probes uniquely aligning over their whole length on + or - strand
set pp=AGLuK
set ff=MicroArray/Hits/best.AGLuK.hits.gz

gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 == 60 && $12<$13) print $1 "\t" $9 ;}' | sort -u | sed -e 's/>//' > MicroArray/Hits/AGLuK.forward.list
gunzip -c $ff | gawk -F '\t' '{ if($8 == "ET_av" && $10 == 1 && $2 == 60 && $12>$13) print $1 "\t" $9 ;}' | sort -u | sed -e 's/>//' > MicroArray/Hits/AGLuK.reverse.list

foreach fr (forward reverse)
# grab the expression index of these probes in all runs
  cat MicroArray/Hits/AGLuK.$fr.list ZZZZZ OTHER_PIPELINES/AGLuK.ace | gawk  '/^ZZZZZ/{zz++;next;}{if(zz<1){p2g[$1]=$2;next;}}/^Gene/{ok=0;if(p2g[$2]){p=$2;ok=1}next;}/^Run_U/{if(ok<1)next;printf("%s:%s:%s\t%f\n",p,p2g[p],$2,$3);}' > MicroArray/Hits/AGLuK.$fr.ma.idx

# grab the expression index of the corresponding genes in all runs
  gzip  MicroArray/Hits/AGLuK.$fr.list
gunzip -c  MicroArray/Hits/AGLuK.$fr.list.gz ZZZZZ.gz RESULTS/Expression/unique/av/NB.AceView.GENE.u.ace.gz | gawk  '/^ZZZZZ/{zz++;next;}{if(zz<1){g2p[$2]=$1;next;}}/^Gene/{gsub(/\"/,"",$2);ok=0;if(g2p[$2]){g=$2;ok=1}next;}/^Run_U/{if(ok<1)next;printf("%s:%s:%s\t%f\n",g2p[g],g,$2,$3);}'  > MicroArray/Hits/AGLuK.$fr.av.idx2
  gunzip  MicroArray/Hits/AGLuK.$fr.list.gz

end

foreach fr (forward reverse)
  echo -n "$fr "
  cat   MicroArray/Hits/AGLuK.$fr.ma.idx  > tutu
  cat ZZZZZ >> tutu
  cat   MicroArray/Hits/AGLuK.$fr.av.idx2  >> tutu
  cat  tutu | gawk '/^ZZZZZ/{zz++;next;}{pp[$1]=1;p2x[zz+1,$1]=$2;next;}END{for(p in pp){a=p2x[1,p];b=p2x[2,p];if(a>0 && b>0){N++;X+=a;X2+=a*a;Y+=b;Y2+=b*b;XY+=a*b;}}x =X2/N - X*X/(N * N) ;y = Y2/N - Y*Y/(N * N) ; w = XY/N - X*Y/(N * N) ;w = w / sqrt (x*y) ;printf("N=%d %f\n",N,w);}'
end

# compute run by run the correlation AGLuK/av of the variation of the gene index, centering each gene around its mean
\rm titi7
foreach fr (forward)
    echo $fr
    cat   MicroArray/Hits/AGLuK.$fr.ma.idx  > tutu7
    cat ZZZZZ >> tutu7
    # cat   MicroArray/Hits/AGLuK.$fr.av.idx  | sed -e 's/Rhs955/RhsXXX/' | sed -e 's/Rhs1200/Rhs955/' | sed -e 's/RhsXXX/Rhs1200/' | sed -e 's/Rhs1210/RhsXXX/' | sed -e 's/Rhs1124/Rhs1210/' | sed -e 's/RhsXXX/Rhs1124/' >> tutu7 
    cat   MicroArray/Hits/AGLuK.$fr.av.idx >> tutu7 
    cat tutu7 | ~/GO/bin/ma_av_correl $fr >> titi7
end
cat titi7 | sort -k 4nr > titi8
echo -n "# RESULTS/$MAGIC.ma.av.correl.stranded.new.txt " > titi9
date >> titi9
echo "# Strand\tRun\tNb of genes\tCorrelation coefficient of the centered index for these genes in AceView and $pp\tTitle\tRunId\tSample\tMicroarray" >> titi9
cat r2t2s2m.txt ZZZZZ titi8 | gawk -F '\t' '/^#/{next;}/ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){r2t[$1]=$2;r2rid[$1]=$3;r2s[$1]=$4;r2m[$1]=$5;next;}}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,r2t[$2],r2rid[$2],r2s[$2],r2m[$2]);}' >> titi9
\cp titi9 RESULTS/$MAGIC.ma.av.correl.stranded.new.txt

### 2014_05_05  : correlation Agilent RNA-seq for some specific runs
foreach run (Rhs800 Rhs1000 Rhs1040 Rhs1130 Rhs1170 Rhs1180 Rhs1200)
  cat  ZZZZZ MicroArray/Hits/AGLuK.forward.ma.idx ZZZZZ MicroArray/Hits/AGLuK.forward.av.idx | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{split($1,aa,":");if(aa[3]!=run)next;g=aa[2];gg[g]=1;gx[g,zz]=$2;}END{for(g in gg)printf("%s\t%.2f\t%.2f\n",g,gx[g,1],gx[g,2]);}' run=$run > MicroArray/Hits/AGLuK.forward.av.correl.$run.txt &
end


R

 a1200=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1200.txt")
 mm1200=as.matrix(a1200[a1200[,2]+a1200[,3]>10,2:3])
 mm1200=as.matrix(,2:3])
 cor(mm1200)

 a1200=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1200.txt")
 mm1200=as.matrix(a1200[a1200[,2]+a1200[,3]>10,2:3])
 mm1200=as.matrix(,2:3])
 cor(mm1200)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs800.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1000.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1040.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1130.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1170.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1180.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 aa=read.table("MicroArray/Hits/AGLuK.forward.av.correl.Rhs1200.txt")
 mm=as.matrix(aa[,2:3])
 cor(mm)

 ii=(aa[,2]<aa[,3]+1) & (aa[,2]>aa[,3]-1)
 cor(aa[ii,2:3])
endif
