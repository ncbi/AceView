#!bin/tcsh -f

if ($1 == CG) goto phaseCG
if ($1 == iter0) goto phaseiter0
if ($1 == iter1) goto phaseiter1
if ($1 == iter2) goto phaseiter2
if ($1 == iter2bis) goto phaseiter2bis
if ($1 == iter3) goto phaseiter3
if ($1 == iter4) goto phaseiter4
if ($1 == iter5) goto phaseiter5
if ($1 == chi2) goto phasechi2
goto phaseLoop

# FATIGUE

Data from columbia, suggested bu Lipman on Feb 6, 2013

we downloaded 4480 fastq.gz file by sftp

===
Feb 18, verif du download
 sftp | tee tyty
 
ls -ls on the remote machine 
mget *
ls -ls locally
# compare
 ls -ls *.fastq.gz | gawk '{printf("%s\t%s\n",$6,$10);}' | sort > titi1
cat tyty | gawk '/fastq.gz/{split($9,aa,"/");printf("%s\t%s\n",$5,aa[3]);}' | sort > titi2

 cat titi1 titi2 | sort -u | wc
diff titi1 titi2
# perfect

###########
## Archive in ARCHIVE/
## Archive 4480 fastq.gz file for 200 patients

# Feb 18, convert to fastq

# list of files
ls *.R1.fastq.gz | sed -e 's/\.R[12].fastq.gz//' > r.list
cp r/list ..
cd ..

mkdir Raw
foreach ff (`cat r.list`)
  scripts/submit Raw/$ff "$bin/dna2dna -i1 DATA/$ff.R1.fastq.gz -i2 DATA/$ff.R2.fastq.gz -I fastq -O raw | sort | gzip > Raw/$ff.raw.gz"
end

# create a single file per person, we have 2200 paired-end raw files, for 200 persons, 120M read per person

cat r.list | sort | gawk '{split($1,aa,"_");p=aa[1]; if(p != old) n++; old=p; printf("Run Rhs_%s\nTitle %s\nPaired_end\nRNA\nFile raw Fatigue.2013_01/%s.raw.gz \n\n",p,p,$1);}' | head

# FATIGUE

Data from columbia, suggested bu Lipman on Feb 6, 2013

we downloaded 4480 fastq.gz file by sftp

===
Feb 18, verif du download
 sftp | tee tyty
 
ls -ls on the remote machine 
mget *
ls -ls locally
# compare
 ls -ls *.fastq.gz | gawk '{printf("%s\t%s\n",$6,$10);}' | sort > titi1
cat tyty | gawk '/fastq.gz/{split($9,aa,"/");printf("%s\t%s\n",$5,aa[3]);}' | sort > titi2

 cat titi1 titi2 | sort -u | wc
diff titi1 titi2
# perfect

###########
## Archive in ARCHIVE/
## Archive 4480 fastq.gz file for 200 patients

# Feb 18, convert to fastq

# list of files
ls *.R1.fastq.gz | sed -e 's/\.R[12].fastq.gz//' > r.list
cp r/list ..
cd ..

mkdir Raw
foreach ff (`cat r.list`)
  scripts/submit Raw/$ff "$bin/dna2dna -i1 DATA/$ff.R1.fastq.gz -i2 DATA/$ff.R2.fastq.gz -I fastq -O raw | sort | gzip > Raw/$ff.raw.gz"
end

# create a single file per person, we have 2200 paired-end raw files, for 200 persons, 120M read per person

cat r.list | sort | gawk '{split($1,aa,"_");p=aa[1]; if(p != old) n++; old=p; printf("Run Rhs_%s\nTitle %s\nPaired_end\nRNA\nFile raw Fatigue.2013_01/%s.raw.gz \n\n",p,p,$1);}' | head

###############################################################################################
###########
## create a bunch of test sets

set s=NoLib3
set s=OvB7
set s=each53
set s=fatig
set s=HRD
set s=hr80
set s=Fb1
set s=Aspi
set s=Clopi
set s=Lipid
set s=Salt

\rm  True_groups.ace VV_groups.ace Rdm_groups.ace
foreach p (109 127 181 211 227 277 283 349 421 443 569 599 601 641 673 727 751 811 919 929)
  echo 1 | gawk '{z="FC" s "."  p ;printf("Run %s.1\n-D Runs\n\nRun %s.2\n-D Runs\n\nRun %s.3\n-D Runs\n\nRun %s.4\n-D Runs\n\n",z,z,z,z,z,z);}' p=$p s=$s >> True_groups.ace
  echo 1 | gawk '{z="VV" s "."  p ;printf("Run %s.1\n-D Runs\n\nRun %s.2\n-D Runs\n\nRun %s.3\n-D Runs\n\nRun %s.4\n-D Runs\n\n",z,z,z,z,z,z);}' p=$p s=$s >> VV_groups.ace

  echo 1 | gawk '{z="FC" s "." p ;printf("Run %s.1\nCompare_to %s.2 %s.1 %s.2 %s.3 %s.4\n\n",z,z,z,z,z,z);}' p=$p s=$s >> True_groups.ace
  echo 1 | gawk '{z="FC" s "." p ;printf("Run %s.3\nCompare_to %s.4 %s.3 %s.4 %s.1 %s.2\n\n",z,z,z,z,z,z);}' p=$p s=$s >> True_groups.ace
  echo 1 | gawk '{z="VV" s "." p ;printf("Run %s.1\nCompare_to %s.3 %s.1 %s.3 %s.2 %s.4\n\n",z,z,z,z,z,z);}' p=$p s=$s >> VV_groups.ace
  echo 1 | gawk '{z="VV" s "." p ;printf("Run %s.2\nCompare_to %s.4 %s.2 %s.4 %s.1 %s.3\n\n",z,z,z,z,z,z);}' p=$p s=$s >> VV_groups.ace

  cat A.list ZZZZZ B.list | gawk '/^ZZZZZ/{zz++}/^Run/{ng[zz+1]++;rr[zz+1,ng[zz+1]]=$2;}END{for (zz=1;zz<=2;zz++){n1=ng[zz];n=int(613*sqrt(p))%n1;for(i=0;i<n1;i++){z=int(2*i/n1);n+=p;n=n%n1;printf("Run FC%s.%d.%d\t//n1=%d n=%d\nRuns %s\nProject %s\n\n",s,p,2*z+zz,n1, n,rr[zz,n],proj);}}}' p=$p  proj=$MAGIC s=$s >> True_groups.ace
  cat A.list ZZZZZ B.list | gawk '/^ZZZZZ/{zz++}/^Run/{ng[zz+1]++;rr[zz+1,ng[zz+1]]=$2;}END{for (zz=1;zz<=2;zz++){n1=ng[zz];n=int(613*sqrt(p))%n1;for(i=0;i<n1;i++){z=int(2*i/n1);n+=p;n=n%n1;printf("Run VV%s.%d.%d\t//n1=%d n=%d\nRuns %s\nProject %s\n\n",s,p,2*z+zz,n1, n,rr[zz,n],proj);}}}' p=$p  proj=$MAGIC s=$s >> VV_groups.ace
 
end

if (0) then
 
set s=Fav ; set NN=91
set s=HRDS ; set NN=94
set s=MNA ; set NN=92
set s=1 ; set NN=117
set s=4S ; set NN=47
set s=4 ; set NN=116
set s=F10G ; set NN=100
set s=XS ; set NN=24
set s=AB_BGI ; set NN=80
set s=AMB ; set NN=80

\rm  VV_groups.ace 
# 109 127 181 211 227
  foreach p ( 569 599 601 641 673 727 751 811 919 929)

    echo 1 | gawk '{z="VV" s "."  p ;printf("Run %s.1\n-D Runs\n\nRun %s.2\n-D Runs\n\nRun %s.3\n-D Runs\n\nRun %s.4\n-D Runs\n\n",z,z,z,z,z,z);}' p=$p s=$s >> VV_groups.ace

  echo 1 | gawk '{z="VV" s "." p ;printf("Run %s.1\nCompare_to %s.3\n\n",z,z);}' p=$p s=$s >> VV_groups.ace

    cat A.list ZZZZZ B.list | gawk '/^ZZZZZ/{zz++}/^Run/{ng[zz+1]++;rr[zz+1,ng[zz+1]]=$2;}END{nk=ng[1]+ng[2];k=int(31417*sqrt(ng[1]*ng[2]))%nk;for (zz=1;zz<=2;zz++){n1=ng[zz];n1a=int(n1*ng[1]/nk);kk=0;n=int(613*sqrt(p))%n1;for(i=0;i<n1;i++){n+=p;n=n%n1;k+=p;k=k%nk;z=0;if(kk>=n1a) z=1;kk++;printf("Run VV%s.%d.%d\t//n1=%d n=%d\nRuns %s\nProject %s\n\n",s,p,2*z+1, n1,n,rr[zz,n+1],proj);}}}' p=$p  proj=$MAGIC s=$s >> VV_groups.ace

  end

endif

echo toto | gawk 'END{n=79;n1=80;;for(i=0;i<n1;i++){n+=p;n=n%n1;print n}' p=$p

 # cat A.list  B.list | gawk '/^ZZZZZ/{zz++}/^Run/{ng[zz+1]++;rr[zz+1,ng[zz+1]]=$2;}END{for (zz=1;zz<=1;zz++){n1=ng[zz];n=int(613*sqrt(p))%n1;for(i=0;i<n1;i++){z=int(4*i/n1);n+=p;n=n%n1;printf("Run Rdm%s.%d.%d\t//n1=%d n=%d\nRuns %s\nProject %s\n\n",s,p,z+zz,n1, n,rr[zz,n],proj);}}}' p=$p proj=$MAGIC  s=$s >> Rdm_groups.ace


#############################################################################################
########### FREQUENCY ANALYSIS OF GENES , removing baddies

####################################
# grab the files with AUC2 > 68 correct code
setenv filterLimit 62
setenv filterLimit 68
setenv filterLimit 70
setenv filterLimit 72
setenv filterLimit 78
setenv filterLimit2 78
setenv limits 3
setenv filterByGene 1
setenv limits "30 60 80 100 120 140 160 180 200"
setenv limits "10 20 30 40 50 60 70 80 90 100"
setenv limits "10 15 20 25 30 35 40 45 50 55 60"
setenv limits "10 12 14 16 18 20 22 24"
setenv limits "10 12 14 16 18 20 22 24 26 28 30"
setenv limits "5 8 10 12 14 16 18 20"
setenv limits "2 3 4 5 6 7 8 9 10"
phaseiter0:

foreach type (FC VV)

  set toto=RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.$type.txt
  echo > $toto.1
  set nstrata=0
  foreach uu (u)
    set qu=""
    if ($uu == nu) set qu="quasi_"
# RefSeq av EBI seqc av RefSeq seqc 
    foreach target (av RefSeq seqc snp)
      if ($target != av && $target != snp && ! -d tmp/PHITS_$target) continue
      if ($target == snp && ! -d tmp/SNPH) continue
      foreach pm (RefSeq av Plus Minus)
        # GENE MRNAH SNP
        foreach typ (GENE MRNAH SNP)
          foreach beta (0  )
            echo "$target\t$uu\t$pm\t$typ\t$beta  RESULTS/Expression/"$qu"unique/$target/$MAGIC.*.$typ.$uu.$type*.beta.$beta.txt"
            cat  RESULTS/Expression/$qu'unique'/$target/$MAGIC.*.$typ.$uu.$type*.beta.$beta.txt | gawk '/# AUC/{ntype=0;ok=0;if(1 && $3<1 && $7>98)next;if($7>filterLimit)ok=1;next;}/^Genes:/{ntype++;type="ZERO";if(ntype==1)type="Minus";if(ntype==2) type="Plus";ok1=0;if(type==pm)ok1=1;next;}/^$/{ok1=0;next;}'"/$pm genes/"'{next;if(ok==1)for(i=5 ; i<=NF;i++){n[$i]++;}next;}{if(ok+ok1==2)n[$1]+=$2;}END{for(k in n)printf("%s\t%d\t%s\t%s\t%s\t%s\t%d\n",k,n[k],target,uu,pm,typ,beta);}' target=$target beta=$beta uu=$uu pm=$pm typ=$typ  filterLimit=$filterLimit >> $toto.1
            if ($pm == Plus) then
              @ nstrata = $nstrata + `cat  RESULTS/Expression/$qu'unique'/$target/$MAGIC.*.$typ.$uu.$type*.beta.$beta.txt | gawk '/# AUC/{if($7>filterLimit)ok++;next;}END{print ok}' filterLimit=$filterLimit `
            endif
          end
        end
      end
    end
  end
 
  set toto=RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.$type.txt
  echo -n "# " > $toto
  date >> $toto
  echo "# $MAGIC using $nstrata strata" >> $toto
  echo "# Gene\tInstance\tPM\tType\tGene Total\tsubtotal\tRefSeq u no iter\tRefSeq u iter\tRefSeq nu no iter\tRefSeq nu iter\tAceView u no iter\tAceView u iter\tAceView nu no iter\tAceView nu iter" >> $toto

  cat ZZZZZ $toto.1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;bad[g]="Baddy";gene=g;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}}dubious[gene]="Dubious";next;}}{g=$1;if(substr(g,1,2)=="FC" && length(g)==7)next;if(0+g>0)next;gg[g]+=$2;tt[$6]=1;ggt[g,$5,$6]+=$2;n[g,$3,$4,$5,$6,0+$7]=$2;gene=g;if(substr(gene,1,2)=="S_"){i=index(gene,".");if(i>3)gene=substr(gene,1,i-1);gene=substr(gene,3);}if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0 && substr(gene,i-2,1)==".")gene=substr(gene,1,i-3);if(i>0 && substr(gene,i-3,1)==".")gene=substr(gene,1,i-4);}}ngene[gene]+=$2;g2gene[g]=gene;}END{for(pmi=0;pmi<2;pmi++){pm="Plus";if(pmi==1)pm="Minus";for(g in gg)for (t in tt)if(ggt[g,pm,t]>0){printf("%s\t%s\t%s\t%s\t%d\t%d",g2gene[g],g,pm,t,ngene[g2gene[g]],ggt[g,pm,t]);for (ira=0;ira<2;ira++){ra="RefSeq";if(ira>0)ra="av";{for(uui=0;uui<2;uui++){uu="u";if(uui>0)uu="nu"; for(iter=0;iter<2;iter++)printf("\t%d",n[g,ra,uu,pm,t,iter]);}}}printf("\n");}}}' | gawk '{if($1=="used" || $1 == "runs")next;print;}' | sort -k 5nr  >> $toto

# list the best element per gene

  set toto2=`echo $toto | sed -e 's/\.txt/.singleInstancePerGene.txt/'`
  cat $toto | gawk '/^#/{print}' > $toto2
# sort genes col 1 by number of occurence of instance col 6 preferring MRNAH to GENE col 4
  if ($filterByGene == 0) then
    cat $toto | sort  -k 6,6nr -k 4,4r | gawk -F '\t' '{g=$1;if(ok[g]<1)print ; ok[g]=1;}' >> $toto2
  else
    cat $toto | sort  -k 5,5nr -k 6,6nr -k 4,4r | gawk -F '\t' '{g=$1;if(ok[g]<1)print ; ok[g]=1;}' >> $toto2
  endif


  foreach limit  (2 3 4 5 6 8  10 12 14 16 18 20 22 24 26 28 30 50)
    cat $toto2 | gawk -F '\t' '{if($5 > limit)ng++;if($6>100*limit)nt++;}END{print limit "\t" ng "\t" nt ;}' limit=$limit | tee -a $toto2
  end

  echo $toto $toto2
end 

touch RESULTS/baddy_batchies.txt 

set limitb=10
cat ZZZZZ RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.FC.singleInstancePerGene.txt ZZZZZ RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.VV.singleInstancePerGene.txt  ZZZZZ  RESULTS/baddy_batchies.txt | gawk  '/^ZZZZZ/{zz++;next;}{if ($5<limit && zz==1)next;if (2*$5<limit && zz==2)next;zz1=zz;g=$2;if(zz==3){zz1=2;g=$1;}nn[g]=1;nz[g,zz1]=1;next;}END{for(g in nn)if(nz[g,1] + nz[g,2]==2)print g}' limit=$limitb | sort >  RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.limit_$limitb.FC_killed_by_batchies_andVV.txt

goto phaseLoop

##########################
# given the list of genes/mrna/snp selected in iter0, construct in SIG_* the index.ace file for these objects

phaseiter1:
set toto=RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.FC.singleInstancePerGene.txt

set limit0=0

if (0) then
  cat RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.FC.singleInstancePerGene.txt | grep SNP | gawk -F '\t' '{if($5>5)print $2;}' > SIG_TR.$MAGIC.$limit/$MAGIC.SNP.list
  cat SIG_TR.$MAGIC.$limit/$MAGIC.SNP.list ZZZZZ tmp/SNPH/$MAGIC_ALL.snp.cover.10.filtered  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){split($1,aa,":");gsub(/2/,">",aa[3]);g=aa[1] ":" aa[2] ":" aa[3] ;s[g]=1; next;}g=$1 ":" $2 ":" $3 ; if(s[g]==1)print}' >  SIG_TR.$MAGIC.$limit/$MAGIC.SNP.snp
endif

# $limits
foreach limit ($limits)
  \rm -rf SIG_TR.$MAGIC.$limit 
  mkdir SIG_TR.$MAGIC.$limit  
  if ($limit0 > 0) continue
  if ($limit0 == 0) set limit0=$limit

  touch RESULTS/baddy_batchies.txt

  cat ZZZZZ RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.FC.singleInstancePerGene.txt ZZZZZ RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.VV.singleInstancePerGene.txt  ZZZZZ  RESULTS/baddy_batchies.txt | gawk  '/^ZZZZZ/{zz++;next;}{if ($5<100*limit && zz==1)next;if (2*$5<100*limit && zz==2)next;zz1=zz;g=$2;if(zz==3){zz1=2;g=$1;}nn[g]=1;nz[g,zz1]=1;next;}END{for(g in nn)if(nz[g,1] > 0 &&  nz[g,2]>0)print g}' limit=$limit | sort >   SIG_TR.$MAGIC.$limit/baddies

  foreach pm (Plus Minus)
    cat SIG_TR.$MAGIC.$limit/baddies ZZZZZ $toto | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){bad[$1]=1;next;}}{if(bad[$2]==1)next;j=6 - filterByGene ; if($3 == PM && ok[$1]<1 && $j >= limit){ok[$1]=1; printf("%s\n",$2);}}' PM=$pm limit=$limit filterByGene=$filterByGene >  SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm
  end

  if (! $?MAGIC_all) then
    echo "please setenv MAGIC_all to the project containing all reads (NBS Fatigue etc)"
    goto phaseLoop
  endif

# hack for predicting $MAGIC given a hand selected list of genes and MRNAS
# accelerate the code by searching the smaller (limit>3) lists in the already exported larger (limit==3) list
  # GENE MRNAH  RefSeq av EBI seqc
  foreach GM (GENE MRNAH SNP)
    foreach target (av RefSeq seqc snp)
      if ($target != snp && ! -d tmp/PHITS_$target) continue
      if ($GM == SNP && $target != snp) continue
      if ($GM != SNP && $target == snp) continue

      foreach pm (Plus Minus) 
        if ($limit == $limit0) then
          if (-e SIG_TR.$MAGIC.$limit/$MAGIC.all.sig_genes.ace) continue
          if (-e SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.$target.$GM.ace) continue
          gzip SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm
          if (-e RESULTS/Expression/unique/RefSeq/$MAGIC_all.RefSeq.GENE.u.ace.gz) then
            echo "using final $pm.$target.$GM"
            gunzip -c SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.gz ZZZZZ.gz RESULTS/Expression/unique/$target/$MAGIC_all.*.$GM.u.ace.gz | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}/^Gene/{v=0;gsub(/\"/,"",$2);v=ok[$2];}/^Tran/{v=0;gsub(/\"/,"",$2);v=ok[$2];$1="Gene";}{if(v==1)print}' > SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.$target.$GM.ace
          else
             echo "NOT using final"
             gunzip -c SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.gz ZZZZZ.gz RESULTS/Expression/unique/$target/$MAGIC.*.$GM.u.ace.gz | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}/^Gene/{v=0;gsub(/\"/,"",$2);v=ok[$2];}/^Tran/{v=0;gsub(/\"/,"",$2);v=ok[$2];$1="Gene";}{if(v==1)print}' > SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.$target.$GM.ace
          endif
          gunzip SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.gz
        else
          cat SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm ZZZZZ  SIG_TR.$MAGIC.$limit0/$MAGIC.sig_genes.$pm.$target.$GM.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}/^Gene/{v=0;gsub(/\"/,"",$2);v=ok[$2];}/^Tran/{v=0;gsub(/\"/,"",$2);v=ok[$2];$1="Gene";}{if(v==1)print}' > SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes.$pm.$target.$GM.ace
        endif
      end
    end
  end

  cat  SIG_TR.$MAGIC.$limit/$MAGIC.sig_genes*.ace > SIG_TR.$MAGIC.$limit/$MAGIC.all.sig_genes.ace
end

goto phaseLoop

##############  run the prediction program on the mixed (gene+mrna+snp) selected object in directories SIG_* at limit==$limit0
phaseiter2:
echo "Run again on about 1000 best mix of genes and transcripts, still not using the excluded runs"
 touch  toto.reject
 if (-e RESULTS/$MAGIC.candidate_genes.filtered) \cp $filterLimit.limit_20.FC_killed_by_batchies_andVV.txt toto.reject


scripts/geneindex.tcsh snp4 000
set limit0=0
foreach limit ($limits)
  if ($limit0 > 0) continue
  if ($limit0 == 0) set limit0=$limit
  if (! -d SIG_TR.$MAGIC.$limit) mkdir SIG_TR.$MAGIC.$limit
  set info=tmp/GENEINDEX/$MAGIC.snp.info.ace
  set runlist=MetaDB/$MAGIC/RunListSorted
  set dg=SIG_TR.$MAGIC.$limit0/$MAGIC.all.sig_genes.ace
 
  echo "bin/geneindex -deepGene $dg -u -runList $runlist -runAce $info  -o  SIG_TR.$MAGIC.$limit/$MAGIC.allOut.sig_genes -gzo -pA -keepIndex  -compare -iterate -rejectMarker toto.reject"
  bin/geneindex -deepGene $dg -u -runList $runlist -runAce $info  -o  SIG_TR.$MAGIC.$limit/$MAGIC.allOut.sig_genes -gzo -pA -keepIndex  -compare -iterate -rejectMarker toto.reject

end

goto phaseLoop

###########  reselect the best objects 

set limit0=0
setenv filterLimit2 74
setenv limits "3 5 8 10 12 14 16 18 20"        
setenv limits "3 4 5 6 7 8 9 10 12 14"

setenv filterLimit2 66
setenv limits "2 3 4 5 6 7 8 9 10 12 14"
setenv limits "5 8 10 12 14 16 18 20"

phaseiter3:
set limit0=0

foreach limit ($limits)
  if ($limit0 > 0) continue 
  if (! -d SIG_TR.$MAGIC.$limit) continue
  if ($limit0 == 0) set limit0=$limit
end

foreach type (FC VV)
  set toto=RESULTS/$MAGIC.second_round_genes.filtered$filterLimit2.$type.txt
  echo > $toto.1
  set nstrata=0
      foreach pm (Plus Minus)
          foreach beta (0  )
            cat  SIG_TR.$MAGIC.$limit0/$MAGIC.allOut.sig_genes.$type*.beta.$beta.txt | gawk '/# AUC/{ntype=0;ok=0;if(1 && $3<1 && $7>98)next;if($7>filterLimit)ok=1;next;}/^Genes:/{ntype++;type="ZERO";if(ntype==1)type="Minus";if(ntype==2) type="Plus";ok1=0;if(type==pm)ok1=1;next;}/^$/{ok1=0;next;}'"/$pm genes/"'{next;if(ok==1)for(i=5 ; i<=NF;i++){n[$i]++;}next;}{if(ok+ok1==2)n[$1]+=$2;}END{for(k in n)printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\n",k,n[k],target,uu,pm,typ,beta);}' target=mix beta=$beta uu=uu pm=$pm   filterLimit=$filterLimit2 >> $toto.1
            if ($pm == Plus) then
              @ nstrata = $nstrata + `cat SIG_TR.$MAGIC.$limit0/$MAGIC.allOut.sig_genes.$type*.beta.$beta.txt  | gawk '/# AUC/{if($7>filterLimit)ok++;next;}END{print ok}' filterLimit=$filterLimit2 `
            endif
          end
      end

  echo -n "# " > $toto
  date >> $toto
  echo "# $MAGIC using $nstrata strata" >> $toto
  echo "# Gene\tInstance\tPM\tType\tGene Total\tsubtotal\tRefSeq u no iter\tRefSeq u iter\tRefSeq nu no iter\tRefSeq nu iter\tAceView u no iter\tAceView u iter\tAceView nu no iter\tAceView nu iter" >> $toto

  cat ZZZZZ $toto.1 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;bad[g]="Baddy";gene=g;if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0)gene=substr(gene,1,i-3);}}dubious[gene]="Dubious";next;}}{g=$1;if(substr(g,1,2)=="FC" && length(g)==7)next;if(0+g>0)next;gg[g]+=$2;tt[$6]=1;ggt[g,$5,$6]+=$2;n[g,$3,$4,$5,$6,0+$7]=$2;gene=g;if(substr(gene,1,2)=="S_"){i=index(gene,".");if(i>3)gene=substr(gene,1,i-1);gene=substr(gene,3);}if(substr(gene,1,3_)=="X__"){gene=substr(gene,4);}else{if(substr(gene,1,1)=="_"){split(gene,aa,"_");gene=aa[2];}else{i=index(gene,"Aug");if(i>0 && substr(gene,i-2,1)==".")gene=substr(gene,1,i-3);if(i>0 && substr(gene,i-3,1)==".")gene=substr(gene,1,i-4);}}ngene[gene]+=$2;g2gene[g]=gene;}END{for(pmi=0;pmi<2;pmi++){pm="Plus";if(pmi==1)pm="Minus";for(g in gg)for (t in tt)if(ggt[g,pm,t]>0){printf("%s\t%s\t%s\t%s\t%d\t%d",g2gene[g],g,pm,t,ngene[g2gene[g]],ggt[g,pm,t]);for (ira=0;ira<2;ira++){ra="RefSeq";if(ira>0)ra="av";{for(uui=0;uui<2;uui++){uu="u";if(uui>0)uu="nu"; for(iter=0;iter<2;iter++)printf("\t%d",n[g,ra,uu,pm,t,iter]);}}}printf("\n");}}}' | gawk '{if($1=="used" || $1 == "runs")next;print;}' | sort -k 5nr  >> $toto


# list the best element per gene

  set toto2=`echo $toto | sed -e 's/\.txt/.singleInstancePerGene.txt/'`
  cat $toto | gawk '/^#/{print}' > $toto2
# sort genes col 1 by number of occurence of instance col 6 preferring MRNAH to GENE col 4
  cat $toto | sort  -k 5,5nr  | gawk -F '\t' '{g=$1;if(ok[g]<1)print ; ok[g]=1;}' >> $toto2


  foreach limit  ($limits)
    cat $toto2 | gawk -F '\t' '{if($5 > limit)ng++;if($6>100*limit)nt++;}END{print limit "\t" ng "\t" nt ;}' limit=$limit | tee -a $toto2
  end

  echo $toto $toto2

end
goto phaseLoop

phaseiter4:
set limit0=0
foreach limit ($limits)
  if ($limit0 > 0) continue
  if (! -d SIG_TR.$MAGIC.$limit) continue
  if ($limit0 == 0) set limit0=$limit
end

# setenv limits "20 25 30 34 37"
# setenv limits "27 30 38 41 57"
setenv MAGIC1 $MAGIC
setenv MAGIC $MAGIC1'final'
scripts/geneindex.tcsh snp4 000
setenv MAGIC $MAGIC1

foreach limit ($limits)
  if (! -d SIG_TR2.$MAGIC.$limit) mkdir SIG_TR2.$MAGIC.$limit
  cat RESULTS/$MAGIC.second_round_genes.filtered$filterLimit2.FC.singleInstancePerGene.txt | gawk -F '\t' '{if($5>=100*limit)print $2;}' limit=$limit >  SIG_TR2.$MAGIC.$limit/geneList.txt
  set info=tmp/GENEINDEX/$MAGIC.snp.info.ace
  if (-e tmp/GENEINDEX/$MAGIC'final'.snp.info.ace) set info=tmp/GENEINDEX/$MAGIC'final'.snp.info.ace 
  set runlist=MetaDB/$MAGIC/RunListSorted
  if (-e tmp/GENEINDEX/$MAGIC'final'.RunListSorted)  set runlist=tmp/GENEINDEX/$MAGIC'final'.RunListSorted
  echo "bin/geneindex -deepGene SIG_TR.$MAGIC.$limit0/$MAGIC.all.sig_genes.ace -u -runList $runlist -runAce $info  -o  SIG_TR2.$MAGIC.$limit/$MAGIC.second_round -gzo -pA -keepIndex  -compare -iterate -rejectMarker toto.reject -selectMarker SIG_TR2.$MAGIC.$limit/geneList.txt"
  bin/geneindex -deepGene SIG_TR.$MAGIC.$limit0/$MAGIC.all.sig_genes.ace -u -runList $runlist -runAce $info  -o  SIG_TR2.$MAGIC.$limit/$MAGIC.second_round -gzo -pA -keepIndex  -compare -iterate -rejectMarker toto.reject -selectMarker SIG_TR2.$MAGIC.$limit/geneList.txt
end

goto phaseLoop

phaseiter5:

date > RESULTS/SIG_TR2.$MAGIC.AUC2.txt

foreach limit ($limits)
  pushd  SIG_TR2.$MAGIC.$limit
  date > AUC2.txt
  foreach ff (`ls *.histo.*.txt`)
    echo $limit | gawk '{gsub("$MAGIC.allOut.sig_genes.","",f);gsub(/...histo.[01].txt/,"",f);printf("%s\t%s\t",$1,f);}' f=$ff >> AUC2.txt
    grep AUC2 $ff >> AUC2.txt
  end
  cat AUC2.txt >> ../RESULTS/SIG_TR2.$MAGIC.AUC2.txt
  popd
end

if ($filterByGene == 1) then
  set nam=by_Gene_score 
else
  set nam=by_Transcript_score 
endif
echo $nam
set toto7=RESULTS/SIG_TR2.$MAGIC.AUC2.summary.filtered$filterLimit2.$nam.txt
echo -n "# " > $toto7
date >> $toto7
echo -n "Limit\tPlus genes\tMinus genes\tAny" >> $toto7
echo  "\tfiltered $filterLimit\t$nam" >> $toto7

foreach limit  ($limits)
  set np=`wc SIG_TR2.$MAGIC.$limit/$MAGIC.sig_genes.Plus | gawk '{print $1}'`
  set nm=`wc SIG_TR2.$MAGIC.$limit/$MAGIC.sig_genes.Minus | gawk '{print $1}'`
  set nn=`echo "$np $nm" | gawk '{print $1+$2}'`
  echo "$limit\t$np\t$nm\t$nn" >> $toto7
end

foreach iter (iter0 iter1)
  cat RESULTS/SIG_TR2.$MAGIC.AUC2.txt | gawk -F '\t' '{if($11!=iter)next;g=$2;n=$1;ns[n]=1;x=$4;y=$8;gg[g]=1;gnx[g,n]=x;gny[g,n]=y;if(index(g,"FC")>0){nfc[n]++;nfcx[n]+=x;nfcy[n]+=y;}if(index(g,"VV")>0){nvv[n]++;nvvx[n]+=x;nvvy[n]+=y;}}END{printf("0\tIteration");for(n=1;n<=300;n++)if(ns[n]==1)printf("\tTraining using X genes occuring in %d strates",n);for(n=1;n<=300;n++)if(ns[n]==1)printf("\tTest using X genes occuring in %d strates",n);printf("\nAverage strata\t%s",iter);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%.2f",nfcx[n]/nfc[n]);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%.2f",nfcy[n]/nfc[n]);printf("\nAverage transverse control\t%s", iter);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%.2f",nvvx[n]/nvv[n]);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%.2f",nvvy[n]/nvv[n]);for(g in gg){printf("\n%s\t%s",g,iter);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%s",gnx[g,n]);for(n=1;n<=300;n++)if(ns[n]==1)printf("\t%s",gny[g,n]);}printf("\n");}' iter=$iter | sort >>  $toto7
end

foreach limit ($limits)
  if (! -d  RESULTS/SIG_TR2.$MAGIC) mkdir RESULTS/SIG_TR2.$MAGIC
  mkdir RESULTS/SIG_TR2.$MAGIC/filtered$filterLimit.filtered2$filterLimit2.$limit
  \cp SIG_TR2.$MAGIC.$limit/*cluded*.beta.?.txt RESULTS/SIG_TR2.$MAGIC/filtered$filterLimit.filtered2$filterLimit2.$limit
  if (! -d SIG_TR2.$MAGIC.$limit) continue
  \cp SIG_TR2.$MAGIC.$limit/baddies RESULTS/SIG_TR2.$MAGIC/baddies.limit$limit.filtered$filterLimit.$nam.txt

  pushd  SIG_TR2.$MAGIC.$limit
    foreach ff (`ls $MAGIC.allOut.sig_genes.*cluded* $MAGIC.allOut.sig_genes.*227.3* $MAGIC.allOut.sig_genes.*Chronic* `)
      \cp $ff  ../RESULTS/SIG_TR2.$MAGIC/$ff.limit$limit.filtered$filterLimit.$nam
    end
  popd
end

echo  $toto7

goto phaseLoop

phaseiter40:
  
set bestLimit=15
echo "Select again the most frequently used objects in the best SIG_TR"

set toto=RESULTS/SIG_TR.$MAGIC.$bestLimit/candidate_genes.iterated.txt
echo > $toto.1

    foreach pm (Plus Minus)
        foreach beta (0 1 )
          cat SIG_TR.$MAGIC.$bestLimit/$MAGIC.allOut.sig_genes.FC*.beta.$beta.txt | gawk '/# AUC/{ok=0;if($7>filterLimit)ok=1;next;}'"/$pm genes/"'{if(ok==1)for(i=5 ; i<=NF;i++){n[$i]++;}}END{for(k in n)printf("%s\t%d\t%s\t%d\n",k,n[k],pm,beta);}' beta=$beta pm=$pm filterLimit=$filterLimit >> $toto.1
        end
    end

echo -n "# " > $toto
date >> $toto
echo "# $MAGIC" >> $toto
echo "# Gene\tInstance\tPM\tType\tGene Total\tsubtotal\tRefSeq u no iter\tRefSeq u iter\tRefSeq nu no iter\tRefSeq nu iter\tAceView u no iter\tAceView u iter\tAceView nu no iter\tAceView nu iter" >> $toto

cat $toto.1 | gawk -F '\t' '{g=$1;if(0+g>0)next;if(g=="runs")next;if(g=="used")next;n=$2;pm=$3;beta=$4;pms[pm]=1;nn[g]+=n;ngpm[g,pm]+=n;nnn[g,pm,beta]+=n;}END{for(pm in pms){for(g in nn)if(ngpm[g,pm]>0){printf("%s\t%s\t%d\t%d\t%d\t%d\n",pm,g,nn[g],ngpm[g,pm],nnn[g,pm,0],nnn[g,pm,1]);}}}' > $toto.2

set limit=10
cat $toto.2 | sort -k 4,4n  | grep Plus | tail -$limit > $toto.3
cat $toto.2 | sort -k 4,4n  | grep Minus | tail -$limit >> $toto.3

mkdir SIG2_TR.$MAGIC.$limit
cat $toto.3 | cut -f 2 | sort -u > SIG2_TR.$MAGIC.$limit/gene.list 
cat SIG2_TR.$MAGIC.$limit/gene.list ZZZZZ  SIG_TR.$MAGIC.3/$MAGIC.all.sig_genes.ace  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}/^Gene/{v=0;gsub(/\"/,"",$2);v=ok[$2];}{if(v==1)print}'  > SIG2_TR.$MAGIC.$limit/$MAGIC.all.sig_genes.ace
 
 bin/geneindex -deepGene SIG2_TR.$MAGIC.$limit/$MAGIC.all.sig_genes.ace -u -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.snp.info.ace  -o  SIG2_TR.$MAGIC.$limit/$MAGIC.allOut.sig_genes -gzo -pA -keepIndex  -compare -iterate 


# list the best element per gene

set toto2=`echo $toto | sed -e 's/\.txt/.singleInstancePerGene.txt/'`
cat $toto | gawk '/^#/{print}' > $toto2
cat $toto | sort -k 6nr | gawk -F '\t' '{g=$1;if(ok[g]<1)print ; ok[g]=1;}' >> $toto2


goto phaseLoop
ls SIG_TR.$MAGIC.*/Fatigue.* > tata
cat tata | gawk '{f=$1;gsub(/\/Fatigue/,"/FatigTH15d5",f);print "mv " $1 " " f;}' > _m
source _m

## prepare a list of mrna for wich we will construct the mRNA wiggles
cat RESULTS/$MAGIC.candidate_genes.filtered$filterLimit.txt | gawk -F '\t' '{if($4!="MRNAH" || $6<15)next;m=$2;c="ET_av";if(substr(m,1,1)=="_")c="KT_RefSeq";print c "\t" m;}' | sort -u > mRNA2Wiggle.txt

#######################################
# export a selection of histograms accross the whole population
set toto=RESULTS/RefSeq.candidates.histo.txt 
echo -n "# " > $toto
date >> $toto
foreach ff ( `ls RESULTS/Expression/*unique/RefSeq/Differential_genes/Gr1_2.RefSeq.GENE.*.1.txt`)
  head -3  $ff | gawk '{print}' >> $toto
  cat RESULTS/RefSeqGenes.txt ZZZZZ $ff | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){ok["X__" $1]=1;next;}if(ok[$1])print}' >>  $toto
end
cat $toto | cut -f 1,2,3,4,5,6,7,8

cat MetaDB/Fatigue/runs.ace  | gawk '/^Runs/{next;}/^Run/{r=$2;next;}/^Sample/{printf("%s\t%s\n",r,$2);}' | sed -e 's/\"//g' | gzip > MetaDB/$MAGIC/run2sample.ace.gz

set toto=RESULTS/RefSeq.candidates.index.txt 
echo -n "# " > $toto
date >> $toto
\rm $toto.1
foreach ff ( `ls RESULTS/Expression/*unique/RefSeq/Gr1_2.RefSeq.GENE.*.ace.gz`)
  gunzip -c  MetaDB/$MAGIC/run2sample.ace.gz ZZZZZ.gz RESULTS/RefSeqGenes.txt.gz ZZZZZ.gz $ff | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){r2s[$1]=$2;next;}}{if(zz<2){ok["X__" $1]=1;next;}}{gsub(/\"/,"",$0);}/^Gene/{g=$2;gok=ok[g];next;}/Run_nU/{if($11=="NAw" || gok<1)next;r=$2;rr[r]=1;gsub("NA/","",$3);ngr[g,r]=$3;}END{printf("# Gene");for (ir = 3991; ir<=4190;ir++) {r="Rhs"ir;printf("\t%s",r);}printf("\n# Sample");for (ir = 3991; ir<=4190;ir++) {r="Rhs"ir;printf("\t%s",r2s[r]);}for(g in ok)if(ok[g]==1){printf("\n%s",g);for (ir = 3991; ir<=4190;ir++){r="Rhs"ir;printf("\t%s",ngr[g,r]);}}printf("\n");}' >>  $toto.1
end
head -1 $toto.1 >> $toto
tail -n +2 $toto.1 | sort >> $toto
\rm  $toto.1
cat $toto | cut -f 1,2,3,4,5,6,7,8

#######################################
# hack to archive the CH data

foreach ii (20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38)

cat toto | grep CHARGE_RNASeq$ii | gawk '{n=split($2,aa,"/");printf("mv %s CHARGE_RNASeq%d\n",aa[n],ii);}' ii=$ii >> _m
end
head _m
source _m
 rsync -avh --whole-file ~/aut/tmp/DATA/CHARGE_RNASeq 2013_06.Human.CH &

qusage 10 
foreach dd (`ls -d B? B??`)
  pushd $dd
  setenv MAGIC $dd
  MAGIC s10
  popd
end

 gunzip -c RESULTS/Expression/unique/AceView/Fatigue.av.MRNAH.nu.ace.gz  | gawk '/^Gene/{print;next;}/ChronicFatigue/{print;next;}' | head -60


g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Tran/{g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Transcript/{g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Group_/{if($2=="ChronicFatigueAndControls_199"){z0=$3;s0=$10;}if($2=="ChronicFatigueGroup1"){z1=$3;s1=$10;}if($2=="ChronicFatigueGroup2"){z2=$3;s2=$10;if(z0 >= 5)printf("%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\n",g,type,z0,z1,z2,z2-z1,s0,s1,s2);}}'

cat RESULTS/29PredictorGenes.txt | grep Aug | gawk '/_/{printf("ET_av\t%s\n",$1);}' > titi
cat RESULTS/29PredictorGenes.txt | grep -v Aug | gawk '/_/{printf("KT_RefSeq\t%s\n",$1);}' >> titi
cat titi TARGET/MRNAS/mRNA2Wiggle.txt | sort -u > mRNA2Wiggle.txt



cat TARGET/MRNAS/hs.NM2gene.txt ZZZZZ RESULTS/2084ForIngenuity.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$4;nm[g]=nm[g] "," $2;gid[g]=gid[g] "," $3;next;}}{gsub (/\r/,"",$NF);g=$1;printf("%s\t%s\t%s",g,substr(gid[g],2),substr(nm[g],2));for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > RESULTS/2084ForIngenuity.gid.nm.txt


====

# rattrapage des fasta de LOC39
foreach run (`cat MetaDB/$MAGIC/RunList | head -9`)
  foreach lane (`cat Fastc/$run/LaneList`)
    echo $lane
    if (! -e tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc) then
      gunzip -c tmp/MRNAWIGGLELANE/$lane.remapped.*LOC399942*.?.minerr0.hits.gz | gawk '/^#/{next}{print substr($1,1,length($1)-1);}' | sort -u >  tmp/MRNAWIGGLELANE/$lane.LOC399942.list
      bin/dna2dna -select tmp/MRNAWIGGLELANE/$lane.LOC399942.list -i Fastc/$lane.fastc.gz -I fastc -O fastc -o tmp/MRNAWIGGLELANE/$lane.LOC399942
      mv tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc1 
      cat tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc1 | gawk '/^>/{printf(">%s_%s\n",run,substr($1,2));next;}{print}' run=$run >  tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc
    endif
    cat tmp/MRNAWIGGLELANE/$lane.LOC399942.fastc | gawk '/^>/{seq=$1;next;}{i=index($1,"><");printf("%s.f\n%s\n%s.r\n%s\n",seq,substr($1,1,i-1),seq,substr($1,i+2));}' >  tmp/MRNAWIGGLELANE/$lane.LOC399942.fasta
  end
end

##################################################################
#### Order all patients

\rm toto
foreach ff (`ls RESULTS/Expression/*unique/*/$MAGIC.*.FC*.*...beta.0.txt`)
  cat $ff | gawk -F '\t' '/^Run:/{print}' | scripts/transpose | gawk '/^Rhs/{n++;print n,$1;}' >> toto
end

# plot the order of all patients
set toto1=RESULTS/$MAGIC.patient_ordering_in_multiple_strates.av.MRNAH.txt
echo -n "# " > $toto1
date >> $toto1
cat << EOF >> $toto1
#  $MAGIC.patient_ordering_in_multiple_strates.txt av MRNAH
# In each strate, a gene signature is constructed by randomly selecting half the controls and half the patients
# then, using this signature all patients and controls are ordered
# The present table reports the centered average position of each person
# In a perfect situation all members of the one group would have a negative score, all members of the other group, a positive score
EOF
echo "# Run\tSample\tGroup\tNumber of strates\tScore\tScore of group 1 patients\tSore of group2 patients" >> $toto1
cat MetaDB/$MAGIC/samples.ace ZZZZZ GR12.ace ZZZZZ toto | gawk '/^ZZZZZ/{zz++;next;}{if(zz==2){if(r2g[$2]==1)nnn++;n[$2]++;nn[r2g[$2],$2]+=$1;next;}}{gsub(/\"/,"",$2);}/^Runs/{r2g[$2]=gr;next;}/^Run/{if(zz==1){gr++;ggr[gr]=$2;}else r2s[$2]=sample;next;}/^Sample/{sample=$2;next;}END{zero=nnn;for (r in n){printf("%s\t%s\t%s\t%d\t%d",r,r2s[r],ggr[r2g[r]],n[r],nn[1,r]+nn[2,r]-zero);if(r2g[r]==1)printf("\t%d\t\n",(nn[1,r]-zero)/n[r]);else printf("\t\t%.2f\n",(nn[2,r]-zero)/n[r]);}}' | sort -k 5n  >> $toto1



cat $toto1 | gawk -F '\t' '{printf("%d\t%d\n",100+$6,100+$7);}' | bin/histo -plain -o toto99 -plot -columns 6,7

##################################################################
####

phaseLoop:
 echo done
exit 0

##################################################################
##################################################################

#  MAGIC g4 m4H s4a ; ~/Fatigue/README CG ; 
# setenv MAGIC NBS80_20final ; MAGIC ii1 ; scripts/geneindex.tcsh snp4 000 ; 
# setenv MAGIC NBS80_20 ; ~/Fatigue/README iter1 ; ~/Fatigue/README iter2 ; ~/Fatigue/README iter3


# catch all reads in LOC100996385 and in TMSB10.aAug10
#  wm1H wait wm2H wait wm3H

#  gunzip -c /home/mieg/Fatigue/tmp/MRNAWIGGLELANEH/$lane.remapped.ET_av:TMSB10.aAug10.r.minerr0.hits.gz ZZZZZ 
#  gunzip -c /home/mieg/Fatigue/tmp/MRNAWIGGLELANEH/Rhs3991/f.4.remapped.KT_RefSeq:_LOC100996385_1.f.minerr0.hits.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}'
# take the intersect
# this gives a set of coordinates in  TMSB10.aAug10
# take all reads in that area of type < (they align on plus strand)
# extract their fasta file
# assemble using saucisse on the 27 bp or using the acembly code

###########

# can we optimize the concordance of the 2 lists
cat RESULTS/NBS.candidate_genes.filtered$filterLimit.txt > toto1
cat RESULTS/NBS80_20.candidate_genes.filtered$filterLimit.txt > toto2

cat toto1 | sort -k 5nr -k 6nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $5>=30)print ; ok[g]=1;}' | cut -f 1 > toto1u
cat toto2 | sort -k 5nr -k 6nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $5>=43)print ; ok[g]=1;}' | cut -f 1 > toto2u
wc toto[12]u
cat toto[12]u | sort -u | wc 

cat toto1 | sort -k 6nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $6>=15)print ; ok[g]=1;}'  | cut -f 1 > toto1u
cat toto2 | sort -k 6nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $6>=17)print ; ok[g]=1;}'  | cut -f 1 > toto2u
wc toto[12]u
cat toto[12]u | sort -u | wc 

cat toto1 |  gawk -F '\t' '{g=$1;print g "\t" $5+$6;}' | sort -k 2nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $2>=45)print ; ok[g]=1;}'  | cut -f 1 > toto1u
cat toto2 |  gawk -F '\t' '{g=$1;print g "\t" $5+$6;}' | sort -k 2nr | gawk -F '\t' '{g=$1;if(ok[g]<1 && $2>=60)print ; ok[g]=1;}'  | cut -f 1 > toto2u
wc toto[12]u
cat toto[12]u | sort -u | wc 

#######
## snp use  F50_50.mysnp  tmp/SNPH/Fatigue.snp.cover.10.filtered

# bin/geneindex -deepSNP  $MAGIC.mysnp -u -mask tmp/GENEINDEX/Fatigue.snp.u.mask  -runList MetaDB/Fatigue/RunListSorted -runAce tmp/GENEINDEX/Fatigue.snp.info.ace  -o  tmp/GENEINDEX/Results/Fatigue.snp..SNP.u -gzo -pA -method snp -export aitv -exportDiffGenes  -compare  -iterate

bin/geneindex -deepSNP  SIG_TR.$MAGIC.5/$MAGIC.SNP.snp -u   -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.snp.info.ace  -o  SIG_TR.$MAGIC.5/$MAGIC.SNP -gzo -pA -method snp            -export aitv   -exportDiffGenes  -compare  -iterate
 
\rm toto
foreach pm (Plus Minus)
  foreach beta (0 1 )
    cat  SIG_TR.$MAGIC.5/$MAGIC.SNP.FC*.beta.$beta.txt | gawk '/# AUC/{ok=0;if($7>filterLimit)ok=1;next;}'"/$pm genes/"'{if(ok==1)for(i=5 ; i<=NF;i++){n[$i]++;}}END{for(k in n)printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\n",k,n[k],target,uu,pm,typ,beta);}' target=snp beta=$beta uu=u pm=$pm typ=SNP  filterLimit=62 >> toto
  end
end
cat toto | sort -k 2nr |  gawk '/:/{if($2>10)print;}' >  SIG_TR.$MAGIC.$limit/$MAGIC.SNP.list1
cat SIG_TR.$MAGIC.$limit/$MAGIC.SNP.list1 ZZZZZ SIG_TR.$MAGIC.$limit/$MAGIC.SNP.snp  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){split($1,aa,":");gsub(/2/,">",aa[3]);g=aa[1] ":" aa[2] ":" aa[3] ;s[g]=1; next;}g=$1 ":" $2 ":" $3 ; if(s[g]==1)print}' >  SIG_TR.$MAGIC.$limit/$MAGIC.SNP.snp1

bin/geneindex -deepSNP  SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1 -u   -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.snp.info.ace  -o  SIG_TR.$MAGIC.5/$MAGIC.SNP1 -gzo -pA -method snp            -export aitv   -exportDiffGenes  -compare  -iterate

cat SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1 ZZZZZ SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1 |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{g=$1 ":" $2 ":" $3 ;if(zz<1){n[g]++;if($19 =="Both_strands" && $22==$23)good[g]++;next;}if(3*good[g]>n[g])print;}' > SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1.compatible

cat SIG_TR.$MAGIC.5/$MAGIC.SNP.snp ZZZZZ SIG_TR.$MAGIC.5/$MAGIC.SNP.snp |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{g=$1 ":" $2 ":" $3 ;if(zz<1){n[g]++;if($19 =="Both_strands" && $22==$23)good[g]++;next;}if(3*good[g]>n[g])print;}' > SIG_TR.$MAGIC.5/$MAGIC.SNP.snp.filtered

bin/geneindex -deepSNP  SIG_TR.$MAGIC.5/$MAGIC.SNP.snp.filtered -u   -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.snp.info.ace  -o  SIG_TR.$MAGIC.5/$MAGIC.SNP_compatible -gzo -pA -method snp            -export aitv   -exportDiffGenes  -compare  -iterate

 
\rm toto
foreach pm (Plus Minus)
  foreach beta (0 1 )
    cat  SIG_TR.$MAGIC.5/$MAGIC.SNP_compatible.FC*.beta.$beta.txt | gawk '/# AUC/{ok=0;if($7>filterLimit)ok=1;next;}'"/$pm genes/"'{if(ok==1)for(i=5 ; i<=NF;i++){n[$i]++;}}END{for(k in n)printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\n",k,n[k],target,uu,pm,typ,beta);}' target=snp beta=$beta uu=u pm=$pm typ=SNP  filterLimit=62 >> toto
  end
end
cat toto | sort -k 2nr |  gawk '/:/{if($2>10)print;}' >  SIG_TR.$MAGIC.$limit/$MAGIC.SNP_compatible.list1
cat SIG_TR.$MAGIC.$limit/$MAGIC.SNP_compatible.list1 ZZZZZ SIG_TR.$MAGIC.$limit/$MAGIC.SNP.snp  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){split($1,aa,":");gsub(/2/,">",aa[3]);g=aa[1] ":" aa[2] ":" aa[3] ;s[g]=1; next;}g=$1 ":" $2 ":" $3 ; if(s[g]==1)print}' >  SIG_TR.$MAGIC.$limit/$MAGIC.SNP_compatible.snp1

bin/geneindex -deepSNP  SIG_TR.$MAGIC.5/$MAGIC.SNP_compatible.snp1 -u   -runList MetaDB/$MAGIC/RunListSorted -runAce tmp/GENEINDEX/$MAGIC.snp.info.ace  -o  SIG_TR.$MAGIC.5/$MAGIC.SNP1_compatible  -gzo -pA -method snp            -export aitv   -exportDiffGenes  -compare  -iterate

mkdir  RESULTS/SIG_SNP
\cp SIG_TR.$MAGIC.5/$MAGIC.SNP.snp RESULTS/SIG_SNP/$MAGIC.SNP.snp.txt
\cp SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1 RESULTS/SIG_SNP/$MAGIC.SNP.snp1.txt
\cp  SIG_TR.$MAGIC.5/$MAGIC.SNP1*cluded* RESULTS/SIG_SNP
\cp  SIG_TR.$MAGIC.5/$MAGIC.SNP.snp1 RESULTS/SIG_SNP/$MAGIC.SNP.snp1.txt
\cp SIG_TR.$MAGIC.5/$MAGIC*Chronic*  RESULTS/SIG_SNP

## AUC2 du test = 72%

#####################
#####################
## hack to eliminate the genes that go in the batch effect
ls RESULTS/Expression/unique/AceView/Fatigue.av.MRNAH.nu | grep LibBatch1Group1_39_LibBatch2G
roup1_58...beta.0.txt

echo " " > titi
foreach gg (LibBatch1Group1_39_LibBatch2Group1_58 LibBatch1Group2_34_LibBatch2Group2_20 LibBatch1Group2_34_LibBatch3Group2_45 LibBatch2Group2_20_LibBatch3Group2_45)
  foreach ff (`ls RESULTS/Expression/unique/*/Fatigue.*.[GM]*.nu.$gg...beta.0.txt`)
    cat $ff  |  gawk '/^Genes:/{ok++;next;}{if(ok>0)printf ("%s\t%s\n",$1,$2);}' >> titi
  end
end
cat titi | gawk '{s[$1]+=$2;}END{for (g in s) if (s[g]>=200)print g"\t" s[g];}' | sort -k 2nr >  RESULTS/baddy_batchies.txt
# use as:  geneindex -reject baddy_batchies.txt

cat RESULTS/baddy_batchies.txt ZZZZZ RESULTS/Expression/*unique/*/Fatigue.*ChronicFatigueGroup1_ChronicFatigueGroup2...histo.0.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){bad[$1]=1;next;}}{if(bad[$1]==1 || bad[$1]==2){bad[$1]++;print}}' > RESULTS/baddy_batchies.histo.txt

##################################################################
#### in geneindex.tcsh, if the gene split per zone is available (file tmp/METADATA/$target.GENE.ns.sponge )
#### then we export the zone index in  tmp/GENEINDEX/Results/NB.av.GENE.u.geneClusters.ace.gz
#### we can then run the predictor program using thse zone-index as

# gunzip -c tmp/GENEINDEX/Results/HRmini.av.MRNAH.u.geneClusters.ace.gz | gawk '/^Run_U/{$3=14 + 3* ($3-10);}{print}' >  tmp/GENEINDEX/Results/HRmini.av.MRNAH.u.geneClusters.dilated.ace

# bin/geneindex -keepIndex -deepGene tmp/GENEINDEX/Results/HRmini.av.MRNAH.u.geneClusters.dilated.ace -u -runList MetaDB/HRmini/RunListSorted -runAce tmp/GENEINDEX/HRmini.snp.SNP.info.ace -o toto89 -gzo -pA -target_class any -exportDiffGenes -compare


##################################################################
####

phaseLoop:
 echo done
exit 0

##################################################################
##########################
## vahan

cat tmp/SNPH/*/KTwins.snp.sorted.filtered | grep KKP |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(1 || g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==200)printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}' > RESULTS/KTwins.comparison.min2.txt

cat tmp/SNPH/mito/KTwins.snp.sorted.filtered | grep KKP |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==200)printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}' > RESULTS/KTwins.comparison.compatible.mito.min2.txt

cat tmp/SNPH/*/KTwins.snp.sorted | grep KKP |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==200)printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}' > RESULTS/KTwins.comparison.min2.txt


# limit to indels
cat tmp/SNPH/*/KTwins.snp.sorted.filtered2 | grep KKP |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;if(substr($3,2,1)==">")next;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==200)printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}' >  RESULTS/KTwins.comparison.compatible_indels.minCoveron2.txt2

# limit to transitions
cat tmp/SNPH/*/KTwins.snp.sorted.filtered2 | grep KKP |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;if($3 == "A>G" || $3 == "T>C" || $3 == "G>A" || $3 == "C>T") ; else next;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==200)printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}' >  RESULTS/KTwins.comparison.compatible_transitions.minCoveron2.txt2

# get the score in all the runs
cat tmp/SNPH/*/KTwins.snp.sorted | grep KKPGP9 |  gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g;srg[snp,run]=g;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;printf("\n\n%s %s", r1, r2);for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==0 && srg[s,r2]==2 || srg[s,r1]==2 && srg[s,r2]==0 )printf("%s %s\n",s,r1);if(srg[s,r2]==200)printf("%s %s\n",s,r2);for(g2 in gs)printf("\t%s",g2);for(g1 in gs){printf("\n%s",g1);for(g2 in gs)printf("\t%d",nn[r1,r2,g1,g2]);}printf("\n");}}'

cat tmp/SNPH/*/KTwins.snp.sorted.filtered | grep KKPGP9 | gawk -F '\t' '{snp=$1 "\t" $2 "\t" $3;if($3 == "A>G" || $3 == "T>C" || $3 == "G>A" || $3 == "C>T") ; else next;run=$6;g=$7;snps[snp]=1;runs[run]=1;if(g==0 || g==1 || g==2)gs[g]=g+10;srg[snp,run]=g+10;}END{for(r1 in runs)for(r2 in runs){if(r1 >= r2)continue ; if(length(r2)<3)continue ;for(s in snps)nn[r1,r2,srg[s,r1],srg[s,r2]]++;for(s in snps)if(srg[s,r1]==10 && srg[s,r2]==12 || srg[s,r1]==12 && srg[s,r2]==10 )printf("chrom\t%s\t%s:%d %s:%d\n",s,r1,srg[s,r1],r2,srg[s,r2]);}}'  >> homo_hetero

KKPGP90 KKPGP91  chrom 1        103530884       G>A :: KKPGP90:10 KKPGP91:12

gzip homo_hetero
zcat  homo_hetero.gz ZZZZZ.gz tmp/SNP/Ghs*/zoneg.*.count.u.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$2 ":" $3 ":" $4]=1;next; }}{if(ok[$1 ":" $2 ":" $3] == 1)print}' > bigsnp2
cat bigsnp2 | cut -f 1,2,3,6,7,8,9,10,11 | sort -k 1,1 -k 2,2 -k 3,3 -k 4,4 | gawk -F '\t' '{r=$4;if(r<"Ghs558")next;if(r<"Ghs565")u="##";else u="--";if(old != "" && $2!=old){printf("\n");old="";}if($7<5)next;old=$2;printf("%s::%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3, $4,u,$5,$6,$7,$8,$9);}' > bigsnp3
cat bigsnp3   | gawk  -F '\t' '{s=$1;if($3=="##")k[s]=1;if($3=="--")m[s]=1;vv[s]=vv[s] "\n" $0; next;}END{for(s in k)if(m[s]==1)print vv[s];}'




cat tmp/SNPH/*/KTwins.snp.sorted.filtered | grep KKPGP | cut -f 9,10 > toto4




