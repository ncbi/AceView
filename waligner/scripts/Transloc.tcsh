##=============================
##=============================
##=============================
## gene pairs in NB

# collect all gene-gene links
echo " collect all gene-gene links in all runs"
if (! -d  tmp/GeneLinks)  mkdir tmp/GeneLinks

if (! -e tmp/GeneLinks/$MAGIC.links_all.support) then
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    if (! -e tmp/GeneLinks/$run.geneLinks) then
      echo "    collecting $run"
      gunzip -c  tmp/COUNT/$run/*.geneLinks.gz | gawk -F '\t' '{if($1<$2)nn[$1 "\t" $2]+=$3;else nn[$2 "\t" $1]+=$3;}END{for(k in nn)printf("%s\t%s\t%d\n",run,k,nn[k]);}' run=$run | sort -k 4nr > tmp/GeneLinks/$run.geneLinks
    endif
    cat  tmp/GeneLinks/$run.geneLinks >> tmp/GeneLinks/$MAGIC.links_all
  end

# add by hand a few cases
  echo "Rhs811\tSGTA\tTCF3\t1000" >>  tmp/GeneLinks/$MAGIC.links_all
endif

# global counts accross all runs
echo "collate per gene-pair accross all runs"
if (! -e  tmp/GeneLinks/$MAGIC.links_all.support) then
  cat  tmp/GeneLinks/$MAGIC.links_all.support | sort -k 1,1 > tmp/GeneLinks/$MAGIC.links_all.sorted
  mv tmp/GeneLinks/$MAGIC.links_all.sorted  tmp/GeneLinks/$MAGIC.links_all
  cat tmp/GeneLinks/$MAGIC.links_all |  gawk -F '\t' '{if($4>=10)nn[$2 "\t" $3]+=$4;}END{for(k in nn)printf("Any\t%s\t%d\n",k,nn[k]);}'  | sort -k 4nr > tmp/GeneLinks/$MAGIC.links_all.global

  # re-collate in all runs just the links with at least 10 supports, with an enumeratio of the sample names -> 10565 cases 
  cat  MetaDB/$MAGIC/_q2 ZZZZZ tmp/GeneLinks/$MAGIC.links_all.global tmp/GeneLinks/$MAGIC.links_all |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){r2nb[$1]=$7;next;}}{k=$2 "\t" $3;}/^Any/{ok[k]=1;next;}{if(ok[k]<1)next;else{if($4>=10){nn10[k]+=$4;nnP10[k]++;}nn1[k]+=$4;nnP[k]++;nb=$1;if(length(r2nb[$1])>2)nb=r2nb[$1];nn2r[k] = nn2r[k] "," nb "(" $4 ")";}}END{for(k in nn1)printf("Any\t%s\t%d\t%d\t%d\t%d\t%s\n",k,nn1[k],nn10[k],nnP[k],nnP10[k],nn2r[k]);}'  | sort -k 4nr > tmp/GeneLinks/$MAGIC.links_all.support
  \rm tmp/GeneLinks/$MAGIC.links_all tmp/GeneLinks/$MAGIC.links_all.global
endif

echo -n "   collated N gene-pair candidates \tN="
wc -l tmp/GeneLinks/$MAGIC.links_all.support

echo -n "Selecting the significant pairs \tN="
# nice doccument presenting these supporting samples, limited to the AceVew genes (excluding EBI and RefSeq -> 5734 cases 

if (! -e  tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt) then
  echo -n "# " > tmp/GeneLinks/$MAGIC.GeneLinks.txt
  date >> tmp/GeneLinks/$MAGIC.GeneLinks.txt
  echo "# Type\tChromosome\tChromosome\tDistance\tGene1\tGene2\tLocus1\tLocus2\tSupport\tSupport in samples with at least 10\tNb sample with support\tNb sample with support > 10\tDetails" >> tmp/GeneLinks/$MAGIC.GeneLinks.txt
cat   tmp/METADATA/av.GENE.info.ace ZZZZZ tmp/GeneLinks/$MAGIC.links_all.support | gawk '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}/^Gene /{g=$2;gsub(/\"/,"",g);next;}/^IntMap/{gsub(/\"/,"",$0);isg2m[g]=1;g2m[g]= $2 ":" $3 "_" $4 ;chrom[g]=$2;a1[g]=$3;if($3<$4)st[g]=1;else st[g]=-1;next;}{g1=$2;g2=$3;g1a=g1;g2a=g2;gsub(/[0-9]/,"",g1a);gsub(/[0-9]/,"",g2a);gsub(/P$/,"",g1a); gsub(/P$/,"",g2a);gsub("zzzzz","",g1a); gsub("zzzzz","",g2a);if(zz<1)next;if(isg2m[g1]<1)next;type="0";da=a1[g1]-a1[g2];if(da<0)da=-da; if(substr(g1a,1,3)!="LOC" && (index(g1a,g2a)>0 || index(g2a,g1a)>0))type="4_famille"; else if(chrom[g1]==chrom[g2] && st[g1] == st[g2] && da<200000)type="5_extension";else if(chrom[g1]==chrom[g2] && st[g1] != st[g2])type="2_inversion";if(chrom[g1]!=chrom[g2]){da="0";if(type=="0")type="1_links2chroms"}if(type=="0")type="0_deletion" ;printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n",type,chrom[$2],chrom[$3],da,$2,$3,g2m[$2],g2m[$3],$4,$5,$6,$7,$8);}' | sort >> tmp/GeneLinks/$MAGIC.GeneLinks.txt

  bin/rsync tmp/GeneLinks/$MAGIC.GeneLinks.txt RESULTS

# select the mUnder1, with high support in not too many individuals
  set nR=`cat MetaDB/$MAGIC/RunList | wc -l`
  cat tmp/GeneLinks/$MAGIC.GeneLinks.txt | gawk -F '\t' '/^#/{print;next;}{n=split($13,aa,",");suma=0;sumb=0;n1=nR;best="";nx=0;n11=0;mx=0;for(i=1;i<=n;i++){split(aa[i],bb,"(");split(bb[2],cc,")");x=cc[1];nx++;if(x<10){suma+=cc[1];if(x>=1)n11++;}else {n1--;best=best "," aa[i]}if(x>mx)mx=x;sumb+=x;}if(n1<1 || 4*nX > nR)next;if(mx<n11+8)next;ma=suma/n1;if(ma<.25 && 2*n1>nR){ printf("%s",$1);for(i=2;i<=NF-1;i++)printf("\t%s",$i);printf("\t%s",best);printf("\t%s\n",$NF);}}' nR=$nR > tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt
  bin/rsync tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt RESULTS
endif
wc -l  tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt

# look for recurence of runs
cat tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt | gawk -F '\t' '/^#/{print;next;}{n=split($13,aa,",");for(i=1;i<=n;i++){split(aa[i],bb,"(");split(bb[2],cc,")");x=cc[1];nx++;if(x>=10)rr[bb[1]]++;}}END{for(r in rr) printf("%s\t%d\n",r,rr[r]);}' | sort -k 2n | tail -100 > tmp/GeneLinks/$MAGIC.GeneLinks.big
cat MetaDB/$MAGIC/_q2 ZZZZZ tmp/GeneLinks/$MAGIC.GeneLinks.big | gawk '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){t2r[$7]=$1;next;}printf("Run %s\n", t2r[$1]);}' > tmp/GeneLinks/$MAGIC.GeneLinks.big.list
cat tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt | gawk -F '\t' '/^#/{print;next;}{n=split($13,aa,",");for(i=1;i<=n;i++){split(aa[i],bb,"(");split(bb[2],cc,")");x=cc[1];nx++;if(x>=10)rr[bb[1]]=rr[bb[1]] "," $5 ":" $6;}}END{for(r in rr)printf("%s\t%s\n",r,rr[r]);}' >  tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1togene

cat MetaDB/$MAGIC/_q2 ZZZZZ tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1togene | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){t2r[$2]=$1;next;}r=$1;if(t2r[r])r=t2r[r];printf("Run %s // %s\t",r,$1);n=split($2,aa,",");for(i=2;i<=n;i++){gsub(/:/," ",aa[i]);printf("Links_2_genes %s\t",aa[i]);}printf("\t\n");}' | sort -k 3 | sed -e 's/\t/\n/g' > tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1togene.ace

# set nR=`cat MetaDB/$MAGIC/RunList | wc -l`
# cat tmp/GeneLinks/$MAGIC.GeneLinks.txt | grep 5605067 |gawk -F '\t' '/^#/{print;next;}{n=split($13,aa,",");suma=0;sumb=0;n1=nR;for(i=1;i<=n;i++){split(aa[i],bb,"(");split(bb[2],cc,")");x=cc[1];if(x<10){suma+=cc[1];n1--;}sumb+=x;}m=sum/n1;print m,n1,suma,sumb,n,suma/n1,sumb/nR;}' nR=$nR



# cat tmp/GeneLinks/$MAGIC.GeneLinks.txt | gawk -F '\t' '{n=split($13,aa,",");for(i=1;i<=n;i++){split(aa[i],bb,"(");split(bb[2],cc,")");if(cc[1]>0)printf("%d\n", cc[1]);};}'  | bin/histo -plain -o rear_support_per_person -plot

# Find the exact splice events by realigning all genes
# take a pair of gene , raboute les mrna.a, exporte les pairs qui mappent sur g1.a et g2.a
# realigne en mode saut pour trouver la vraie deletion

# construction de la gene-pair
if (! -d tmp/TranslocGenePair)  mkdir  tmp/TranslocGenePair
cat tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt | gawk -F '\t' '{g1=$5 "::" $7;g2=$6 "::" $8;n=split($13,aa,",");for(i=2;i<=n;i++)printf("%s\t%s\t%s\n",g1,g2,aa[i]);}' >  tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1a
cat  MetaDB/$MAGIC/_q2 ZZZZZ tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1a | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){nb2r[$7]=$1;next;}}{split($3,aa,"(");printf("%s\t%s\t%s\n",$1,$2,nb2r[aa[1]]);}' > tmp/TranslocGenePair/g2g2r.txt

cat tmp/TranslocGenePair/g2g2r.txt | sort | gawk '{split($1,aa,"::");split($2,bb,"::");ok[aa[1],bb[1]]++;if(ok[aa[1],bb[1]]>1)next;printf("if(! -d  tmp/TranslocGenePair/GenePair.%s.%s) then\n  mkdir  tmp/TranslocGenePair/GenePair.%s.%s\n",aa[1],bb[1],aa[1],bb[1]);printf(" scripts/submit  tmp/TranslocGenePair/GenePair.%s.%s/prepare \" scripts/Transloc.GenomePair.tcsh %s %s %s\"\nendif\n",aa[1],bb[1],aa[1],bb[1],"ZZZZZ");}' > tmp/TranslocGenePair/g2g2r.go
cat tmp/TranslocGenePair/g2g2r.txt | sort | gawk '{split($1,aa,"::");split($2,bb,"::");run=$3;ok[aa[1],bb[1],run]++;if(ok[aa[1],bb[1],run]>1)next;printf("if( -d  tmp/TranslocGenePair/GenePair.%s.%s) then\n",aa[1],bb[1]); printf("  scripts/submit  tmp/TranslocGenePair/GenePair.%s.%s/cumul \"scripts/Transloc.GenomePair.tcsh %s %s %s\"\nendif\n",aa[1],bb[1],aa[1],bb[1],run);}' > tmp/TranslocGenePair/g2g2r.go2

goto phaseLoop
# create the directories and the acedb databases, wait, run the alignments on the farm
source  tmp/TranslocGenePair/g2g2r.go
scripts/submit wait
source  tmp/TranslocGenePair/g2g2r.go2

scripts/submit wait

# analysis of results
# reject cases where the same pair appears as F1.F2 and F2.F1
cat tmp_old/TranslocGenePair/Gene*/*/all.genome.introns | gawk '{r=$1;t=$2;g=$5;if(g2r[g,r]<1){g2r[g,r]=1;gg2r[g]=gg2r[g] "," r ;}if(g2t[g,t]<1){ gg[g]=gg[g] "," t  ;}g2t[g,t]+=$13;ggn[g] += $13;}END{for(g in gg){split(g,aa,"."); printf("%s", ggn[g] "\t" aa[2] "\t" aa[3] "\t" g  "\t" substr(gg2r[g],2)  ) ; n=split(gg[g],aa, ",");for(i=2;i<=n;i++)printf("\t%s\t%d",aa[i],g2t[g,aa[i]]);printf("\n");}}' | sort -k 1n > RESULTS/TranslocGenePairs.txt


cat MetaDB/$MAGIC/_q2 ZZZZZ tmp/GeneLinks/$MAGIC.GeneLinks.mUnder1.txt | gawk -F '\t'  '/^ZZZZZ/{zz++;next;}{gsub(/\"/,"",$0);}{if(zz<1){nb2r[$7]=$1;next;}}{g1=$5;g2=$6;z[g1,g2]=$0;n=split($13,aa,",");for(i=2;i<=n;i++){split(aa[2],bb,"(");printf("tmp/TranslocGenePair/GenePair.%s.%s/%s\t%s\t%s\t%s\t",g1,g2,nb2r[bb[1]],g1,g2,nb2r[bb[1]]);print;}}'  >   tmp/TranslocGenePair/g2g2r.1

head -1 tmp/TranslocGenePair/g2g2r.1 > tmp/TranslocGenePair/g2g2r.2
foreach dd (`cat tmp/TranslocGenePair/g2g2r.1 | gawk '/^#/{next;}{print $1;}'`)
  if (-e $dd/summary.txt) then
    foreach ii (1 2 5 6 7 8 9 10)
      cat tmp/TranslocGenePair/g2g2r.1 |  gawk -F '\t' '{if($1==dd){for(i=2;i<=NF;i++)printf("\t%s",$i);}}' dd=$dd  >>tmp/TranslocGenePair/g2g2r.2
      cat $dd/all.genome.introns | tail -1 | gawk -F '\t' '{printf("\t%s\t%s\t%d",$2,$11,$13);}' >> tmp/TranslocGenePair/g2g2r.2
      echo $ii | gawk '{tt="GENE\tDNA\tCHROMOSOME\tCOORD\tWIGGLE\tEXON_1\tEXON_2\tLINK";split(tt,aa,"\t");printf("\t%s",aa[$1]);}' >> tmp/TranslocGenePair/g2g2r.2
      cat $dd/summary.txt | cut -f $ii | scripts/transpose >> tmp/TranslocGenePair/g2g2r.2
    end
  endif
end

echo -n "# " >  RESULTS/TranslocGenePair.summary.single_individual.txt
date >>  RESULTS/TranslocGenePair.summary.single_individual.txt
echo "\tGene1\tGene2\Run\tChromosome1\tChromosome2\tDistance\tGene1\tGene2\tLocus1\tLocus2\tSupport in any run\tSupport in samples with at least 10\tNumber of samples with some support\tNumber of samples with support >= 10\tRuns with candidate rearrangments (number of supporting "  >> RESULTS/TranslocGenePair.summary.single_individual.txt
cat  tmp/TranslocGenePair/g2g2r.2 | gawk -F '\t' '{if($16==1)print;}' >> RESULTS/TranslocGenePair.summary.single_individual.txt





phaseLoop:
  echo done
