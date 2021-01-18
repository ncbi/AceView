#!bin/tcsh -f

#############################################
## generic metadata file, should be replaced by tmp/METADATA

#  export from 37lm5 the genes with pastille_regulation_antisens
  # cat toto.as.list ZZZZZ TARGET/GENES/av.geneTypes.txt | gawk '/ZZZZZ/{zz++;next;}{if(zz<1){g=$2;gsub(/\"/,"",g);as[g]=1;next;}}{t=$2;gsub("conserved","ancient",t);gsub("novel","recent",t);g=$1;if(index(t,"spliced")>0){if(as[g]==1)t=t "_antisense";else t=t "_free";}printf("%s\t%s\n",g,t);}' > TARGET/GENES/av.geneTypesAntiSense.txt
 set toto=TARGET/GENES/av.metadata.txt

  gunzip -c TARGET/Targets/hs.av.fasta.gz | gawk '/^>/{split($1,aa,"|");print aa[3];}' | $tab_sort -u > $toto.1
  cat ZZZZZ TARGET/GENES/av.gene2intMap.txt >> $toto.1
  cat ZZZZZ TARGET/GENES/av.gene2geneid.txt >> $toto.1
  # cat ZZZZZ TARGET/GENES/av.geneTypesAntiSense.txt >> $toto.1
  cat ZZZZZ TARGET/GENES/av.gene2title.txt | grep -v 'cloud gene' >> $toto.1
  cat ZZZZZ   >> $toto.1
  gunzip -c TARGET/MRNAS/hs.av.transcript2gene.txt.gz  | sed -e 's/\"//g' >> $toto.1
  cat ZZZZZ MicroArray/Hits/AGLuK.av.probe2gene.unique.txt   >> $toto.1
  
  echo -n "# " > $toto
  date >> $toto
  echo "# Metadada for the 55836 AceView-2011 genes (a subset of the AceView 2010 genes) used in the SEQC project." >> $toto
  printf "# The corresponding fasta file is at ftp-SEQC/SEQC_Reference_Targets/human.AceView.2010.selected.v4.fasta.gz\n" >> $toto
  printf "# Gene\tChromosome\tGeneId\tType\tTitle\tTranscripts\tUniquely mapped Agilent probes\n" >> $toto
  cat $toto.1 |  gawk -F '\t'  '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){gg[$1]=1;next;}}{if(zz<2){g2chrom[$1]="chr" $2 ":" $3 "-" $4;next;}}{if(zz<3){g2gid[$1]=g2gid[$1] ";" $2;next;}}{if(zz<4){g2typ[$1]=$2;next;}}{if(zz<5){g2title[$1]=$2;next;}}{if(zz<6){g2mrna[$2]=g2mrna[$2] ";" $1;next;}}{if(zz<7){g2probe[$2]=g2probe[$2] ";" $1;next;}}END{for(g in gg){typ=g2typ[g];if(typ=="")typ="single_exon";printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",g,g2chrom[g],substr(g2gid[g],2),typ,g2title[g],substr(g2mrna[g],2),substr(g2probe[g],2));}}' | $tab_sort >> $toto

\rm $toto.1


# comparison between the DEG for the groups listed in MOA.list between the expression-index measured by methods Target[1:2]
set Target1=$1
set Target2=$2
set GM0=$3
set MOA_genes=MOA.$Target1.$Target2.$GM0

if (! -d tmp/MOA) mkdir tmp/MOA
if (! -d tmp/$MOA_genes) mkdir tmp/$MOA_genes
if (! -d RESULTS/$MOA_genes) mkdir RESULTS/$MOA_genes
\rm  RESULTS/$MOA_genes/*

set qunique=unique

tbly MetaDB <<EOF
  query find Experiment project == $MAGIC // NB_DEG
  bql -a -o tmp/MOA/experiments.$GM0.txt  select e,r1,r2 from e in @, r1 in e->DEG_high_run, r2 in e->DEG_low_run
EOF

# ATTENTION, use capitalization of Experiments names so each pair is ordered from good to bad phenotype
cat tmp/MOA/experiments.$GM0.txt | gawk -F '\t' '{gsub(/\"/,"",$0);s=1;if($2>$3){z=$2;$2=$3;$3=z;s=-1;}printf("%s\t%s\t%s\t%d\n",$1,$2,$3,s);}' | sort -k 2,2 -k 3,3 -k 1,1 | gawk -F '\t' 'BEGIN{u=-1;}{u=-u;s=$6;if(u!=s){z=$2;$2=$3;$3=z;}printf("%s\t%s\t%s\t%d\n",$1,$2,$3,s);}' > tmp/MOA/experiments.sorted.$GM0
cat tmp/MOA/experiments.sorted.$GM0 | gawk -F '\t' '{printf("%s_%s\n",$2,$3);printf("%s_%s\n",$3,$2);}' > tmp/MOA/MOA.list

###############################################
## grab the most significant on/off genes transcripts and probes

set toto=RESULTS/MOA.close_to_ON_OFF_elements.txt
echo -n "#Date " > $toto
date >> $toto
echo "# This table identifies the top differential genes with low or very low expression in one group, aiming to identify ON/OFF genes or transcripts" >> $toto 
printf "#Compared groups\tDetection method: GENE RNA-seq genes, MRNAH RNA-seq transcripts, MA microarray\tDifferential score (range 0 to 200)\tAverage expression index of low group\tAverage expression index of high group\tlog2 of Fold_change\tGene or transcript\tMicroarray probe\n" >> $toto

foreach moa (`cat tmp/MOA/MOA.list`)
    set ff1=$moa
echo "--- $sign $moa $ff1"
  foreach gm (GENE MRNAH )
    set target=av
    if ($gm == MA) set target=AGLuK
    set ff=RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target.$gm.u.$ff1.diffGenes.0.txt 
    if (! -e $ff) then
       echo "FATAL ERROR $moa $gm"
       ls -ls $ff
       continue
    endif
echo "$sign $ff1"

    # isG==2 for GENE , 0 for MRNAH because stupid table has a different format
    set isG=0
    if ($gm == MA) set isG=1
    if ($gm == GENE) set isG=2
    if ($MAGIC == Pain && $gm == GENE) set isG=1

    set p2g=MicroArray/Hits/AGLuK.av.probe2gene.unique.txt 
    cat $p2g ZZZZZ $ff  | gawk -F '\t' '/^#/{next}/^ZZZZZ/{zz++;next;}{if(zz<1){p2g[$1]=$2;next;}}{g=$1;score[g]=$2;fc=$(10+isG);if(fc<0){fc=-fc;low[g]=$(11+isG);}else high[g]=$(11+isG);FC[g]=fc;}END{for(g in score)if(FC[g]>1 && low[g]<8.5 && score[g]>=70){g1=g;g2="";if(gm=="MA"){g1=p2g[g];g2=g;}printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",moa,gm,score[g],low[g],high[g],FC[g],g1,g2);}}' moa=$ff1 gm=$gm isG=$isG   | sort -u | sort -k 3nr  >> $toto
  end
end
mv $toto $toto.1
cat $toto.1 | sort -k 1,1 -k 3nr -k 4n > $toto
cat $toto |  gawk -F '\t' '/^#/{print;next;}{if($2 == "GENE" && $3>190 && $6 < 6.5)print}' > $toto.GENE.txt

###############


foreach target ($Target1 $Target2)
  if (-e tmp/$MOA_genes/tutu.$target) \rm tmp/$MOA_genes/tutu.$target
  echo XXXXX >> tmp/$MOA_genes/tutu.$target
  set sign= -1
  \rm DEG.$target.list
  set GM=MA
  if ($target == av || $target == RefSeq || $target == EBI) set GM=$GM0
  foreach ff (`cat tmp/MOA/MOA.list`)
    @ sign =  - $sign
    if ($sign == 1) then
      set ff1=$ff
      cat tmp/MOA/experiments.$GM0.txt  | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\t%s\t%s\t%s\n",$2,$3,$4,$5,$6);}' | grep $ff > tmp/$MOA_genes/tutu7
    else
      cat tmp/MOA/experiments.$GM0.txt  | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\t%s\t%s\t%s\n",$3,$2,$4,$5,$6);}' | grep $ff > tmp/$MOA_genes/tutu7
      set ff1=`cat tmp/MOA/experiments.$GM0.txt  | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s\t%s\t%s_%s\n",$2,$3,$3,$2);}' | grep $ff |  gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\n",$1,$2);}'`
    endif
    if ($target == av) then
      set limit=`cat tmp/$MOA_genes/tutu7 | cut -f 2`
    else
      set limit=`cat tmp/$MOA_genes/tutu7 | cut -f 3`
    endif
    echo $target,$ff,$sign,$limit
    set dd="RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target.$GM.u"
    echo  $dd.$ff.diffGenes.1.txt >> DEG.$target.list
    # echo $target,$sign,$ff1,$dd
    # echo $dd.$ff1.diffGenes.1.txt
    set p2g=MicroArray/Hits/AGLuK.av.probe2gene.unique.txt 
    if ($GM0 == MRNAH) set p2g=MicroArray/Hits/AGLuK.av.probe2mrna.unique.txt 
    echo " " > RESULTS/$MOA_genes/hs.av.transcript2gene.txt
    if ($GM0 == MRNAH)  gunzip -c TARGET/MRNAS/hs.av.transcript2gene.txt.gz | sed -e 's/\"//g' >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
    set GMS=GENES
    if ($GM0 == MRNAH)  set GMS=MRNAS

    # add the meta data and reorder the AGL lines so that they are keyed on the gene
    if (0) then 
      cat  TARGET/GENES/av.gene2geneid.txt ZZZZZ $p2g ZZZZZ TARGET/GENES/av.geneTypesAntiSense.txt ZZZZZ TARGET/GENES/av.gene2title.txt ZZZZZ TARGET/$GMS/av.gene2intMap.txt ZZZZZ RESULTS/$MOA_genes/hs.av.transcript2gene.txt ZZZZZ $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gid[$1]=gid[$1] ";" $2 ; next;}}{if(zz==1){affy[$2]=affy[$2] ";" $1; next;}}{if(zz==2){gType[$1]=gType[$1] ";" $2; next;}}{if(zz==3){gTitle[$1]=$2; next;}}{if(zz==4){gChrom[$1]= "chr" $2 ":" $3 "-" $4; next;}}{if(zz==5){t2g[$1]=$2; next;}}{if(substr($1,1,5)=="UKv4_")affy[$1]=";"$1;t1=substr(gType[$1],2);t2="single_exon";if(t1)t2=t1;printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,ff,$2,substr(gid[$1],2),substr(affy[$1],2),t2,gTitle[$1] gTitle[t2g[$1]],gChrom[$1],$11);}' ff=$ff | gawk -F '\t' '/^#/{print;next;}/^Genes/{print;next;}{if($3>limit)print}' limit=$limit sign=$sign >> tmp/$MOA_genes/tutu.$target
    endif

   # same, better metadata
   cat TARGET/GENES/av.metadata.txt ZZZZZ RESULTS/$MOA_genes/hs.av.transcript2gene.txt ZZZZZ $p2g ZZZZZ $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g = $1 ; gChrom[g]=$2 ; gid[g]=$3 ;gType[g]=$4 ;gTitle[g]=$5;next;}}{if(zz==1){t2g[$1]=$2; next;}}{if(zz==2){affy[$2]=affy[$2] ";" $1; next;}}{if(substr($1,1,5)=="UKv4_")affy[$1]=";"$1;t2=gType[$1];printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,ff,$2,gid[$1],substr(affy[$1],2),t2,gTitle[$1] gTitle[t2g[$1]],gChrom[$1],$11);}' ff=$ff | gawk -F '\t' '/^#/{print;next;}/^Genes/{print;next;}{if($3>limit)print}' limit=$limit sign=$sign >> tmp/$MOA_genes/tutu.$target


  end
end

if (0) then

cat TARGET/GENES/av.metadata.txt ZZZZZ RESULTS/$MOA_genes/hs.av.transcript2gene.txt ZZZZZ $p2g ZZZZZ $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g = $1 ; gChrom[g]=$2 ; gid[g]=$3 ;gType[g]=$4 ;gTitle[g]=$5;next;}}{if(zz==1){t2g[$1]=$2; next;}}{if(zz==2){affy[$2]=affy[$2] ";" $1; next;}}{if(substr($1,1,5)=="UKv4_")affy[$1]=";"$1;t1=substr(gType[$1],2);t2="single_exon";if(t1)t2=t1;printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,ff,$2,substr(gid[$1],2),substr(affy[$1],2),t2,gTitle[$1] gTitle[t2g[$1]],gChrom[$1],$11);}' ff=$ff | gawk -F '\t' '/^#/{print;next;}/^Genes/{print;next;}{if($3>limit)print}' limit=$limit sign=$sign >> tmp/$MOA_genes/tutu3.$target

endif


# TARGET/GENES/av.metadata.txt
#
foreach target ($Target1 $Target2)
  cat tmp/$MOA_genes/tutu.$target | gawk -F '\t' '/^XXXXX/{sign=-1;next;}/^# Param/{next;}/^# AUC/{next;}/^Genes: overexpressed in/{sign=-sign;if (0)print "TTTTT",sign ; next;}{g=$1;f=$2;ff[f]=1;jj=f2jj[f];if(jj<1){jjMax++;jj=jjMax;f2jj[f]=jj;jj2f[jj]=f;}sco=$3;if(0 && index(g,"Rat230")==0 && index(g,":")>0)next;gg[g]=1;score[f,g]=sign*sco;if(0)print "TTTT",g,score[f,g];if(1 || substr(target,1,2) == "av"){if($4)gid[g]= ";" $4;if($5)affy[g]=$5;if($6)gType[g]=$6;if($7)title[g]=$7;if($8)chrom[g]=$8;}if($9){nV++; if(nV %2 == 1)valeur[g]=valeur[g] "(" $9;else valeur[g]=valeur[g] ":" $9 ")";}}END{printf("0\t0\t0\tGene\tNCBI GeneId\tMicroarray probe(s) mapped to the gene\tGene description\tType\tMap chromosome strand from to on build 37\tValues");for(jj=1;jj<=jjMax;jj++) {f=jj2f[jj];printf("\t%s",f);}for (g in gg){z=0;z2=0;nnp=0;nnn=0;col=1;for(jj=1;jj<=jjMax;jj++) {f=jj2f[jj];col=2*col;z1=score[f,g];if(z1<0){z1=-z1;if(z1>0)nnn++;}else{if(z1>0)nnp++};z+=z1;if(z1>0)z2+=col;}nnnn=0;if(nnn>0 && nnp > 0)nnnn=100000;if(nnn>1)nnn=2;if(nnp>2)nnp=2;printf("\n%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s",nnnn+1000+0*(nnn+nnp),100*z2+nnn,z,g,substr(gid[g],2),affy[g],title[g],gType[g],chrom[g],valeur[g]);for(jj=1;jj<=jjMax;jj++) {f=jj2f[jj];printf("\t%s",score[f,g]);}}printf("\n");}'  target=$target | sort -k 1,1n -k 2,2nr -k 3nr | gawk -F '\t' '{if(old > 0 && $2 != old)printf("\n");old=$2;print;}' > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$target.txt
end


set nMoa=`cat tmp/MOA/MOA.list | wc -l`

# gene types :: obsolete, we have a better count by gene-type further down
# set toto=RESULTS/$MOA_genes/$MOA_genes.gene_types_new_novel.txt
# echo -n '# ' > $toto
# date >> $toto

#    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$Target1.txt | gawk -F '\t' '{g=$4;t=$8;v=$9;nnn[t]++;n[t,0]++;if($5)known[t,0]++;else novel[t,0]++;for(ii=1;ii<=nMoa/2;ii++){i=10+ii;z=$i;if(1){n[t,ii]++;if($5)known[t,ii]++;else novel[t,ii]++;}}}END{printf("Limit\tType");for(ii=0;ii<=nMoa/2;ii++)printf("\t\tAceView\tKnown\tNovel");printf("\n");for(k in nnn)if(length(k)>3){printf("%d\t%s",1, k);for(ii=0;ii<=nMoa/2;ii++)printf("\t\t%d\t%d\t%d",n[k,ii],known[k,ii],novel[k,ii]); printf("\n");}}' nMoa=$nMoa >> $toto

set toto=RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.txt 
echo -n '# ' > $toto.1
date >>  $toto.1

 cat  $p2g  ZZZZZ TARGET/GENES/av.metadata.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$Target2.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$Target1.txt  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){af2gene[$1]=$2;gene2af[$2]=gene2af[$2] ";" $1;next;}}{if(zz==1){g=$1;gid[g]=$3 ;gMap[$1]=$2;gType[$1]=$4; gTitle[$1]=$5; next;}}{ if(zz==2){if(0)next;af=$4;iiAf++;ii2af[iiAf]=af;afscore[af]=$1;for(i=2;i<=3;i++)afscore[af]=afscore[af] "\t" $i;affy[af]="";for(i=10;i<=NF;i++)affy[af]=affy[af] "\t" $i; next;}}{ if(NF<5)next;if(NF>nf){nf=NF;}g=$4;if(substr(g,1,5)=="UKv4_")gene2af[g]=";"g;na=split(gene2af[g],aa,";");printf("%s",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);for(iii=2;iii<=7;iii++){afok[aa[iii]]=1;if(affy[aa[iii]]){printf("\t##");if(affy[aa[iii]]){printf("\t%s",aa[iii]);printf("\t%s",affy[aa[iii]]);}else{for(i=0;i<nf-7;i++)printf("\t");}}}printf("\t###\n");}END{for(g in gene2af){  na=split(gene2af[g],aa,";"); }for(ii=1;ii<=iiAf;ii++){af=ii2af[ii];if(afok[af]<1 && affy[af]){g=af2gene[af];t1="single_exon";t2=gType[g];if(t2)t1=t2;g2=af;if(g)g2=g;  i3Max=split(gene2af[g],aa3,";");if(i3Max<1){i3Max=2;aa3[2]=af;}if(0 && g=="AKT2")print "YYYY",g,gene2af[g],i3Max;for(i3=2;i3<=i3Max;i3++){af3=aa3[i3]; afok[af3]=1;if(0 && g=="AKT2")print "ZZZ",g,gene2af[g],i3Max,i3,af3,affy[af3];if(affy[af3]=="")continue;printf("%s\t%s\t%s\t%s\t%s\t%s\t%s",afscore[af3],g2,gid[g2],af3,gTitle[g2],t1,gMap[af2gene[af3]]);for(i=0;i<nf-9;i++)printf("\t");printf("\t##\t%s\t%s\t###\n",af3,affy[af3]);}}}}' |  gawk -F '\t' '{if(old > 0 && $2 != old)printf("\n");old=$2;print;}' | gawk -F '\t' '{printf("%s",$1);for(i=2;i<=6;i++)printf("\t%s",$i); for(iii=0;iii<6;iii++){i=12+nMoa/2+iii*(nMoa/2+4);printf("\t%s",$i);}for(i=7;i<=NF;i++)printf("\t%s",$i);printf("\n");}'  nMoa=$nMoa |  sort -k 4,4 >>  $toto.1


# BUG on a rate la fusion des lignes du gene AKT2 (et pleins de cas similaires)
# cat TARGET/GENES/av.gene2geneid.txt ZZZZZ $p2g   ZZZZZ TARGET/GENES/av.geneTypesAntiSense.txt ZZZZZ TARGET/GENES/av.gene2title.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$Target2.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.$Target1.txt  | gawk -F '\t' '/^ZZZZZ/{print}/AKT2/{print}' > tata

# cat  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.txt | grep AKT2 > tata

if (1) then
  # we need to fuse the probes mapped to the same gene , in MRNA case they are eliminated above
  cat  $toto.1 | sort -k 4,4 >   $toto.1s
  cat $toto.1s ZZZZZ $toto.1s |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{g=$4;if(zz<1){ng[g]++;next;}if(ng[g]<2 && $1>0)print;}'  | sort -k 4,4 > $toto.2
  cat $toto.1s ZZZZZ $toto.1s |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{g=$4;if(zz<1){ng[g]++;next;}if(ng[g]>=2 && $1>0)print;}' | sort -k 4,4 > $toto.3

  echo "1\t1\t1\tZZZZ"  >> $toto.3

  cat $toto.3 | grep IGHV_ > tata
  cat $toto.3 | grep SLC8A1.eAug10 > tata
  echo "1\t1\t1\tZZZZ" >> tata

  if ($GM0 == GENE)  then
    cat  $toto.3  | gawk -F '\t' '{g=$4;if(nng>0 && g!=oldg){printf ("%d\t%d\t%d\t%s\t%d\t%s",nn1,nn2,nn3/nng,oldg,gid,substr(ppp,2));if(nng>6)nng=6;for(i=1;i<=nng;i++)printf("\t%s",pp[i]);for(i=nng+1;i<=6;i++)printf("\t");for(i=13;i<21;i++)printf("\t%s",uu[i]);for(jj=1;jj<=nng;jj++)for(i=1;i<=8;i++)printf("\t%s",vv[i,jj]);printf("\t###\n");nn3=0;nng=0;ppp="";}oldg=g;nng++;nn1=$1;nn2=$2;nn3+=$3;gid=$5;ppp=ppp ";" $6 ; pp[nng]=$6; for(i=1;i<21;i++)uu[i]=$i;for(i=1;i<=8;i++)vv[i,nng]=$(i+20);}' > $toto.4
  else
    cat  $toto.3 | gawk -F '\t' '{t=$4;if(nng>0 && t!=oldt){printf ("%d\t%d\t%d\t%s\t%s\t%s",nn1,nn2,nn3/nng,oldt,gid,substr(ppp,2));if(nng>6)nng=6;for(i=1;i<=nng;i++)printf("\t%s",pp[i]);for(i=nng+1;i<=6;i++)printf("\t");for(i=13;i<21;i++)printf("\t%s",uu[i]);for(jj=1;jj<=nng;jj++)for(i=1;i<=8;i++)printf("\t%s",vv[i,jj]);printf("\t###\n");nn3=0;nng=0;ppp="";}oldt=t;nng++;nn1=$1;nn2=$2;nn3+=$3;gid=$5;ppp=ppp ";" $6 ; pp[nng]=$6; for(i=1;i<21;i++)uu[i]=$i;for(i=1;i<=8;i++)vv[i,nng]=$(i+20);}' > $toto.4
  endif

  cat $toto.2 $toto.4 > $toto
  \rm $toto.1  $toto.1s  $toto.2  $toto.3 $toto.4
endif

## AGL
set pp=AGLuK
if (! -e  MicroArray/Hits/$pp.probe2chrom.unique.txt) then
  gunzip -c MicroArray/Hits/$pp.genome.hits.gz | gawk -F '\t' '/Z_genome/{gsub(">","",$1);s="+";if($12>$13)s="-";printf("%s\tchr%s:%d-%d\n",$1,$11,$12,$13);}' > MicroArray/Hits/$pp.probe2chrom.txt
  cat  MicroArray/Hits/$pp.probe2chrom.txt ZZZZZ  MicroArray/Hits/$pp.probe2chrom.txt | gawk '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){nn[$1]++;next;}if(nn[$1]==1)print;}' > MicroArray/Hits/$pp.probe2chrom.unique.txt

  gunzip -c MicroArray/Oligos/Fasta/$pp.fasta.gz | gawk '/^>/{printf ("Probe %s\n",substr($1,2));next;}{printf("Motif %s\nLength %d\n\n",$1,length($1));}' > MicroArray/Oligos/original/$pp.motif.ace
  cat MicroArray/Oligos/original/$pp.ln_entropy_TM.txt | gawk '{gsub(/>/,"",$1);printf("Probe %s\nTM %s\n\n",$1,$4);}' > MicroArray/Oligos/original/$pp.ln_entropy_TM.ace

  cat MicroArray/Hits/$pp.probe2chrom.unique.txt | gawk '/^#/{next;}{split($2,aa,":");split(aa[2],bb,"-");gsub(/chr/,"",aa[1]);printf ("Probe %s\nIntMap %s %d %d\n\n",$1,aa[1],bb[1],bb[2]);}' > MicroArray/Hits/$pp.probe2chrom.unique.ace
endif

echo ZZZZZ >  RESULTS/$MOA_genes/hs.av.transcript2gene.txt
if ($GM0 == MRNAH) gunzip -c TARGET/MRNAS/hs.av.transcript2gene.txt.gz | sed -e 's/\"//g' >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
echo ZZZZZ >>  RESULTS/$MOA_genes/hs.av.transcript2gene.txt
if ($GM0 == MRNAH) cat   TARGET/GENES/av.gene2geneid.txt >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
echo ZZZZZ >>  RESULTS/$MOA_genes/hs.av.transcript2gene.txt
if ($GM0 == MRNAH) then
  if (-e  TARGET/GENES/av.geneTypesAntiSense.txt) cat  TARGET/GENES/av.geneTypesAntiSense.txt  >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
  if (-e  TARGET/GENES/av.geneTypesAntiSense.txt.gz) gunzip -c  TARGET/GENES/av.geneTypesAntiSense.txt.gz  >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
endif
echo ZZZZZ >>  RESULTS/$MOA_genes/hs.av.transcript2gene.txt
cat $p2g >>  RESULTS/$MOA_genes/hs.av.transcript2gene.txt
echo ZZZZZ >>  RESULTS/$MOA_genes/hs.av.transcript2gene.txt 
cat TARGET/MRNAS/av.gene2intMap.txt >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt

  set sign=-1
  set nn=0
  foreach moa (`cat tmp/MOA/MOA.list`)
    @ sign = - $sign
    if ($sign > 0) then
      @ nn = $nn + 1

      echo -n '# ' > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
      date >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
      cat MicroArray/Hits/$pp.probe2chrom.txt RESULTS/$MOA_genes/hs.av.transcript2gene.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.txt | gawk -F '\t'  '/^ZZZZZ/{zz++;next;}/^#/{if(0)next;}{if(zz<1){probe2map[$1]=$2;next;}}{if(zz<2){mrna2gene[$1]=$2;next;}}{if(zz<3){gene2gid[$1]=$2;next;}}{if(zz<4){gene2typ[$1]=$2;next;}}{if(zz<5){p2g[$1]=$2;g2p[$2]=$1;next;}}{if(zz<6){m2chrom[$1]="chr" $2 ":" $3 "-" $4;next;}}{line++;ok = 0;z1=1000;z2=0;z3=0;z4=0;z5=0;nm=0;np=0;nma=0;npa=0;nmr=0;npr=0;if(line==1){gsub("_"," ",moa);gsub("MYCN not amplified"," ",moa);z1=0;ok=1;i=16+nn;$i="RNA-seq " moa ;for(iii=1;iii<=6;iii++){i=16+nn+iii*(nMoa/2+4);$i="Probe " iii " " moa;}}if($15);else $15=m2chrom[$4];if($15);else $15=probe2map[$4];i=16+nn;z=$i+0;if(z<0){nm++;nmr++;z=-z;z3++;z4=1;}else if(z>0){np++;npr++;z3++;z4=1;}z2+=z;for(iii=1;iii<=6;iii++){i=16+nn+iii*(nMoa/2+4);z=$i+0;if(z<0){nm++;nma++;z=-z;z3++;z5=1;}else if(z>0){np++;npa++;z3++;z5=1;}z2+=z;}if(nm>0)z1=2000;if(nm*np>0){z1=3000;if(nma*npa>0)z1=4000;} z6=1;if(p2g[$4]||g2p[$4])z6=0;if(nm+np+ok>0){z2=z2/(z3+.001);printf("%d\t%d\t%d\t%s",z1+40-10*z4-20*z5+100*z6,z2,z3,$4);typ=gene2typ[mrna2gene[$4]];if(typ);else typ="single_exon";if(GM=="MRNAH")printf("\t%s\t%s\t%s",mrna2gene[$4],gene2gid[mrna2gene[$4]],typ);for(i=5;i<=16;i++)printf("\t%s",$i);i=16+nn;printf("\t%s",$i);for(iii=1;iii<=6;iii++){i=16+nn+iii*(nMoa/2+4);printf("\t%s",$i);}printf("\n");}}' nMoa=$nMoa nn=$nn  moa=$moa GM=$GM0 | sort -k 1,1n -k 2,2nr  | gawk -F '\t' '{if(old > 0 && $1 != old)printf("\n");old=$1;print;}' >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
      wc RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt 
    endif
  end

cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.*.txt | cut -f 4 | gawk '/^UK/{print $1}' | sort -u  > RESULTS/$MOA_genes/orphan_probes.list

## verify the thresholdsadd metadata to the fused table

foreach moa (`cat tmp/MOA/MOA.list`)
  echo -n " $moa plus "
  set n=17
  if ($GM0 == MRNAH) set n=20
  cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | cut -f $n | gawk '{if($1<0){$1=-$1;if($1>0)print}}' | sort -u | sort -k 1n | head -1
  echo -n " $moa plus "
  cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | cut -f $n | gawk '{if($1+0>0)print}' | sort -u | sort -k 1n | head -1
end


## add metadata to the fused table


if ($GM0 == MRNAH) then

  set toto=RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.txt 
  set toto1=RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.with_metadata.txt 

  cat MicroArray/Hits/$pp.probe2chrom.txt >  $toto.9
  echo ZZZZZ >>  $toto.9
  gunzip -c TARGET/MRNAS/hs.av.transcript2gene.txt.gz | sed -e 's/\"//g' >> $toto.9
  echo ZZZZZ >>  $toto.9 
  cat $p2g >>   $toto.9  
  echo ZZZZZ >>    $toto.9
  cat TARGET/MRNAS/av.gene2intMap.txt >>   $toto.9  
  echo ZZZZZ >>    $toto.9
  cat TARGET/GENES/av.metadata.txt >>  $toto.9

  cat $toto.9 ZZZZZ $toto |  gawk -F '\t'  '/^ZZZZZ/{zz++;next;}/^#/{if(NF>4){printf("# Differential feature\tMap\tGene\tNCBI GeneId\tGene type");for(i=5;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}}{if(zz<1){probe2map[$1]=$2;next;}}{if(zz<2){mrna2gene[$1]=$2;next;}}{if(zz<3){p2g[$1]=$2;g2p[$2]=$1;next;}}{if(zz<4){m2chrom[$1]="chr" $2 ":" $3 "-" $4;next;}}{if(zz<5){g=$1;gene2gid[g]=$3;gene2typ[g]=$4;gene2title[g]=$5;next;}}{t=$4;g=mrna2gene[t];if(length(g)<2)g=p2g[t];if(length(g)<2)g=t;map=m2chrom[t];$14="";$15="";$16=t;$13=gene2title[g];if(length(map)<2)map=m2chrom[t];if(length(map)<2)map=probe2map[t];printf("%s\t%s\t%s\t%s\t%s",t,map,g,gene2gid[g],gene2typ[g]);for(i=6;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $toto1.1
  echo -n '# ' > $toto1
  date >>  $toto1
  cat $toto1.1  | gawk '/Parameters/{next;}/^#/{print;}' | sort -u >> $toto1
  cat $toto1.1  | gawk '/^#/{next;}/^UKv4_/{next;}{print;}' | sort  >> $toto1
  cat $toto1.1  | gawk '/^#/{next;}/^UKv4_/{print;}' | sort  >> $toto1

endif

echo "establish RMP"

# establish type RMP : R:RnaSeq  0:absent, 1:up, 2:down, 3:both,  M:microArray 0:absent, 1:up, 2:down, 3:both, P:has probe o:absent, 1:exists
foreach moa (`cat tmp/MOA/MOA.list | sort -u`)  
  cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | sort -k 5,5 -k 4,4 > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.6
  # add the metadata

  if ($GM0 == MRNAH) then
    mv RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.6 RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.5
    cat $toto.9 ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.5 |  gawk -F '\t'  '/^ZZZZZ/{zz++;next;}/^#/{if(zz==5)print;next;}{if(zz<1){probe2map[$1]=$2;next;}}{if(zz<2){mrna2gene[$1]=$2;next;}}{if(zz<3){p2g[$1]=$2;g2p[$2]=$1;next;}}{if(zz<4){m2chrom[$1]="chr" $2 ":" $3 "-" $4;next;}}{if(zz<5){g=$1;gene2gid[g]=$3;gene2typ[g]=$4;gene2title[g]=$5;next;}}{t=$4;g=mrna2gene[t];if(length(g)<2)g=p2g[t];if(length(g)<2)g=t;$6=gene2gid[g];$7=gene2typ[g];$16=gene2title[g];map=m2chrom[t];if(length(map)<2)map=probe2map[t];$17="";$18=map;$19=t;printf("%s",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.6 
  endif

  cat  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.6 ZZZZZ RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.6 | gawk -F '\t' '/^#/{if(zz==1)print;next;}/^ZZZZZ/{zz++;next;}{g=$4;if(GM=="MRNAH")g=$5;if(g!=oldg){Ru=0;Rd=0;Mu=0;Md=0;P=0;}oldg=g;if(zz<1){u=$17;if(GM=="MRNAH")u=$20;if(u+0>0)Ru=1;if(u+0<0)Rd=2;for(jj=0;jj<=6;jj++){k=18+jj;if(GM=="MRNAH")k=21+jj;if(k<=NF){u=$k;if(u+0>0)Mu=1;if(u+0<0)Md=2;}}u=$6;if(GM=="MRNAH")u=$9;if(length(u)>4)P=1;score[g]=100*(Ru+Rd)+10*(Mu+Md)+P;next;}printf("%03d",score[g]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' GM=$GM0 | sort -k 1,1n -k 2,2nr -k 5,5 -k 4,4 > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
  \rm  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.[56]
end

# final sorting by type

  foreach moa (`cat tmp/MOA/MOA.list | sort -u`)  
    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | gawk -F '\t' '{if(substr($4,1,5)=="UKv4_")next;print;}' > titi.$$
    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | gawk -F '\t' '{if(substr($4,1,5)=="UKv4_")print;}' > titi.1.$$

    cat titi.$$ | gawk -F '\t' '/^#/{print}' >  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n '# ' > RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    date >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "# Sorting column, 3 digits: hundreds [0:not differential, 1:up, 2:down"  >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == MRNAH) echo -n ", 3:complex some transcripts up and some down" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "] in RNA-seq, tens [0,1,2,3] in Agilent, units [0: no probe, 1:at least one probe]" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tAverage of the absolute value of the differential score in all measurments, RNA or Agilent">> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tNumber of measures">> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == MRNAH) echo -n "\tTranscript" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tGene\tNCBI GeneId" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == MRNAH) echo -n "\tGene type\t" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == MRNAH) echo -n "\tMicroarray probe(s) mapped to the transcript" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == GENE)  echo -n "\tMicroarray probe(s) mapped to the gene" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tAgilent probe 1\tAgilent probe 2\tAgilent probe 3" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tAgilent probe 4\tAgilent probe 5\tAgilent probe 6" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tGene description\tGene type\tMap chromosome from to on build 37" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == GENE)  echo -n "\tGene" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    if ($GM0 == MRNAH)  echo -n "\tTranscript" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo -n "\tRNA-seq $moa" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    foreach i (1 2 3 4 5 6)
      echo -n "\tProbe $i $moa" >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    end
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==111)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==11)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==101)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==100)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt

    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==221)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
     echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
   cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==21)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
     echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
   cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==201)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
     echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
   cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==200)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt

    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==331)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==301 || $1==311 || $1==321)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==300)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt

    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==31 || $1==131 || $1==231)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    cat titi.$$ | gawk -F '\t' '/^#/{next;}{if($1==121 || $1==211)print}' | sort -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
     cat titi.1.$$ | gawk -F '\t' '/^#/{next;}{printf("%d",1000+$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' | sort -k 1,1n -k 2,2nr -k 4,4 >> RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
    echo >>  RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt
 

    \rm  titi.$$  titi.1.$$
  end


# count the associated genes: 1000 = up, 2000 = down, 3000 = contradictions, 
# for each Up-Up, Up agilent missed in rna, up agilent not measired in rna, up rna-seq missed in agilent, up in rna-seq not assessed in agilent

if (-e RESULTS/$MOA_genes/_tutuG) \rm RESULTS/$MOA_genes/_tutuG
if (-e RESULTS/$MOA_genes/_tutuT) \rm RESULTS/$MOA_genes/_tutuT
  set moa2=""
  foreach moa (`cat tmp/MOA/MOA.list`)  
    if ($moa == $moa2) continue
    set moa2=$moa
    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | gawk -F '\t' '/^#/{next;}{g=$4;if(GM=="MRNAH")g=$5;s=$1;if(s+0<1)next;g2s[g]=s;ns1[s]++;if(ns1[s]==1){k++;kk[k]=s;}npr=0;probes=$6;if(GM=="MRNAH")probes=$9;if(probes)npr=split(probes,aa,";");if(npr>gpr[g])gpr[g]=npr;nprv=0;i0=18;if(GM=="MRNAH")i0=21;for(i=i0;i<NF;i++)if($i+0!=0)nprv++;if(nprv>gprv[g])gprv[g]=nprv;}END{for(g in g2s){x=g2s[g];ns[x]++;nspr[x]+=gpr[g];nsprv[x]+=gprv[g];}for(i=1;i<=k;i++){s=kk[i];printf("%s\t%d\t%d\t%d\t%d\n",moa,s,ns[s],nspr[s],nsprv[s]);for(g in g2s)if(g2s[g]>=3000)printf("#\t%s\t%s\t%d\n",moa,g,g2s[g]);}}' moa=$moa GM=$GM0  >> RESULTS/$MOA_genes/_tutuG
    if($GM0 == MRNAH) then
      cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt | gawk -F '\t' '/^#/{next;}{g=$4;                    s=$1;if(s+0<1)next;g2s[g]=s;ns1[s]++;if(ns1[s]==1){k++;kk[k]=s;}npr=0;probes=$9;if(probes)npr=split(probes,aa,";");if(npr>gpr[g])gpr[g]=npr;nprv=0;i0=21;for(i=i0;i<NF;i++)if($i+0!=0)nprv++;if(nprv>gprv[g])gprv[g]=nprv;}END{for(g in g2s){x=g2s[g];ns[x]++;nspr[x]+=gpr[g];nsprv[x]+=gprv[g];}for(i=1;i<=k;i++){s=kk[i];printf("%s\t%d\t%d\t%d\t%d\n",moa,s,ns[s],nspr[s],nsprv[s]);for(g in g2s)if(g2s[g]>=3000)printf("#\t%s\t%s\t%d\n",moa,g,g2s[g]);}}' moa=$moa GM=$GM0  >> RESULTS/$MOA_genes/_tutuT
    endif
  end


cat RESULTS/$MOA_genes/_tutuG | gawk '/^#/{next}{m=$1;s=$2;n=$3;npr=$4;nprv=$5;mm[m]=1;ns[s]++;if(ns[s]==1){k++;kk[k]=s;}ss[s]+=n;sm[s,m]=n;smpr[s,m]=npr;smprv[s,m]=nprv;}END{for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%03d",s);}  for(m in mm){printf("\n%s",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",sm[s,m]);}}for(m in mm){printf("\n%s.obs",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smprv[s,m]);}}for(m in mm){printf("\n%s.possible",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smpr[s,m]);}}for(m in mm){printf("\n%s.missed",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smpr[s,m]-smprv[s,m]);}}printf("\n");}' | scripts/transpose  > RESULTS/$MOA_genes/transcripts_stats_per_gene.txt

  if($GM0 == MRNAH) then
    cat RESULTS/$MOA_genes/_tutuT | gawk '/^#/{next}{m=$1;s=$2;n=$3;npr=$4;nprv=$5;mm[m]=1;ns[s]++;if(ns[s]==1){k++;kk[k]=s;}ss[s]+=n;sm[s,m]=n;smpr[s,m]=npr;smprv[s,m]=nprv;}END{for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%03d",s);}  for(m in mm){printf("\n%s",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",sm[s,m]);}}for(m in mm){printf("\n%s.obs",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smprv[s,m]);}}for(m in mm){printf("\n%s.possible",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smpr[s,m]);}}for(m in mm){printf("\n%s.missed",m);for(i=1;i<=k;i++){s=kk[i];if(ss[s]>0)printf("\t%s",smpr[s,m]-smprv[s,m]);}}printf("\n");}' | scripts/transpose  > RESULTS/$MOA_genes/transcripts_stats_per_transcript.txt
  endif

foreach gt (gene transcript)
  if ($gt == transcript && $GM0 == GENE) continue
  \cp  RESULTS/$MOA_genes/transcripts_stats_per_$gt.txt  RESULTS/$MOA_genes/transcripts_stats_per_$gt.sum.txt
  cat  RESULTS/$MOA_genes/transcripts_stats_per_$gt.txt | gawk -F '\t' '{if(GM=="MRNAH" && gt="gene"){g=$5;genes[g]++;if(genes[g]>1)next;}s=$1;if(s==100 || s==101 || s==111 || s==11)z="1 Gene up";else if (s==200 || s==201 || s==221 || s==21)z="2 Gene down";else if(s==331)z="9 Complex RNA and agilent";else if(s==31 || s==131 || s==231)z="10 Complex agilent";else if(s==300 || s==301 || s==311 || s==321)z="11 Complex RNA";else if(s>1000)z="1000 Orphan probe not mapped or not mapped to a single gene";else z="12 Conflict";for(i=2;i<=NF;i++){u[i,z]+=$i;if(s<1000)uu[i]+=$i;u[i,0+s]=$i;}uuu[z]=1;nf=NF;}END{for(i=2;i<=nf;i++){u[i,"5 RNA and AGL"]=u[i,111]+u[i,221];u[i,"6 AGL not RNA"]=u[i,11]+u[i,21];u[i,"7 RNA not AGL"]=u[i,101]+u[i,201];u[i,"8 RNA discovery"]=u[i,100]+u[i,200];}uuu["5 RNA and AGL"]=1;uuu["6 AGL not RNA"];uuu["7 RNA not AGL"]=1;uuu["8 RNA discovery"]=1;for(z in uuu){printf("\n%s",z);for(i=2;i<=nf;i++)printf("\t%d",u[i,z]);}printf("\n98 Total");for(i=2;i<=nf;i++)printf("\t%d",uu[i]);printf("\n");printf("\n99 Total measurable both in Agilent and RNA-seq");for(i=2;i<=nf;i++)printf("\t%d",uu[i]-u[i,"8 RNA discovery"]);printf("\n");}' GM=$GM0 gt=$gt | sort -k 1,1n >> RESULTS/$MOA_genes/transcripts_stats_per_$gt.sum.txt
end

# finally we clipout the UK only

  foreach moa (`cat tmp/MOA/MOA.list | sort -u`)  
    mv RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.9
    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.9 | gawk -F '\t' '/^#/{print;next;}/^toto/{next;}/^Sorting/{print;next;}{if($1+0<1)next;if (index($1,"3")>0)next;if (substr($1,1,2)=="21" || substr($1,1,2)=="12")next;if($1>1000)next;if(old+0>0 && $1!=old)printf("\n");old=$1;if ($1>0)print}' >  RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt
    echo toto | gawk '{for(i=1;i<=500;i++)printf("\n");}' >> RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt
    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.9 | gawk -F '\t' '/^#/{next;}/^Sorting/{next;}/^toto/{next;}{if($1+0<1)next;if (substr($1,1,2)=="21" || substr($1,1,2)=="12")next;if (index($1,"3")==0)next;if($1>1000)next;old=$1;print}' | sort -k 5,5 -k 2,2nr  -k 4,4 -k 1,1 >>  RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt
    echo toto | gawk '{for(i=1;i<=500;i++)printf("\n");}' >> RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt
   cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.9 | gawk -F '\t' '/^#/{next;}/^Sorting/{next;}/^toto/{next;}{if($1+0<1)next;if (substr($1,1,2)=="21" || substr($1,1,2)=="12");else next;if (index($1,"3")>0)next;if($1>1000)next;old=$1;print}' | sort -k 5,5 -k 2,2nr  -k 4,4 -k 1,1 >>  RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt
    \rm RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.$moa.txt.9
  end

endif

# to get the noMA counts edit if(zz==4){if(0)next;...} to if(1)next;  in the construction of the FDRlimit.fused file
# mv _tutu _tutu.noMA
# mv  RESULTS/$MOA_genes/transcripts_stats_per_gene.txt  RESULTS/$MOA_genes/transcripts_stats_per_gene.noMA.txt

\rm RESULTS/$MOA_genes/_tutu[TG]

# RESULTS/MOA.av.AGLuK.MRNAH/MOA.av.AGLuK.MRNAH.FDRlimit.fused.*


###########################################################
## count the differential probes per type

set toto=RESULTS/$MOA_genes/differential_probes_per_gene_category.txt
 set sign=-1
  set nn=0
  foreach moa (`cat tmp/MOA/MOA.list`)
    @ sign = - $sign
    if ($sign > 0) then
      @ nn = $nn + 1
      cat MicroArray/Hits/AGLuK.av.probe2gene.unique.txt ZZZZZ RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt ZZZZZ  RESULTS/MOA.av.AGLuK.MRNAH/MOA.av.AGLuK.MRNAH.FDRlimit.AGLuK.txt  | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}/Sorting/{next;}{if(zz<1){p2g[$1]=$2 ; next;}}{if(zz<2){g=$4;if(GM=="MRNAH")g=$5;types[$1]=1;g2typ[g]=$1;next;}}{p=$4;s=$(10+nn);g=p2g[p];ss=0;if(s+0>0)ss=1;if(s+0<0)ss=2;typ=g2typ[g];if(ss>0 && typ)nTyp[typ,ss]++;}END{printf("# Number of probes per given gene type\n#Experiment\tType\tPlus probes\tMinus probes");for(typ in types)printf("\n%s\t%s\t%d\t%d",moa,typ,nTyp[typ,1],nTyp[typ,2]);printf("\n");}' nn=$nn GM=$GM0 moa=$moa >> $toto

    endif
  end

set toto=RESULTS/$MOA_genes/All_differential_genes_seen_via_$GM0.txt
\rm $toto.1
foreach moa (`cat tmp/MOA/MOA.list`)
  cat RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt | gawk -F '\t' '/^#/{next;}/Sorting/{next;}{g=$4;if(GM=="MRNAH")g=$5;s=$1;if(s==100 || s==101 || s==111 || s==11)z="1 Gene up";else if (s==200 || s==201 || s==221 || s==21)z="2 Gene down";else if(s==331)z="9 Complex RNA and agilent";else if(s==31 || s==131 || s==231)z="10 Complex agilent";else if(s==300 || s==301 || s==311 || s==321)z="11 Complex RNA";else if(s>1000)z="1000 Orphan probe not mapped or not mapped to a single gene";else z="12 Conflict";if(s<1000)printf("%s\t%s\t%s\t%s\n",GM,moa,z,g);}' GM=$GM0 moa=$moa > $toto.1
end
echo -n "# " > $toto
date >> $toto
echo "# Each gene is counted only once in each experiment, the counts are therefore additive" >> $toto
echo "Experiment\tType\tGene name" >> $toto
cat $toto.1 | sort > $toto
cat $toto.1 | gawk '{printf("%s\tAny\tany\t%s\n",GM,$4);}'  GM=$GM0  | sort >> $toto
\rm $toto.1


###########################################################
## gene types statistics

echo -n "# " > RESULTS/$MOA_genes/Venn.MOA_gene_types.txt
date >> RESULTS/$MOA_genes/Venn.MOA_gene_types.txt
echo "# MOA\tType\tRNA-seq with no geneId no probe\tRNA-seq geneId no probe\tRNA-seq with probe and no geneid\tRNA_seq with probe and geneid\tAceView and RefSeq\tAveView only\tGene with AGL probe\tNo Agilent probe\tAny" >> RESULTS/$MOA_genes/Venn.MOA_gene_types.txt
echo toto > toto.$$
foreach moa (`cat tmp/MOA/MOA.list | sort -u`)
  if (! -e   RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt) continue
   cat RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt >> toto.$$
   cat RESULTS/$MOA_genes/$MOA_genes.fused.$moa.txt  | gawk -F '\t' '/^#/{next;}{g=$4;gid=$5;ma=$6;t=$14;if(GM=="MRNAH"){g=$5;genes[g]++;if(genes[g]>1)next;gid=$6;t=$7;ma=$9;}z=0;if(length(g)<1)next;if(gid>0)z=1;tt[t]=1;if(length(ma)>1)g2ma[g]=2;if(g==ma)next;gg[g]=1;gt[g,t]=z+1;}END{for(g in gg)for(t in tt)if(gt[g,t]>0)nn[t,gt[g,t]+g2ma[g]]++;for(t in tt)for(i=1;i<=4;i++)nnn[i]+=nn[t,i];for(t in tt)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",moa,t,nn[t,1],nn[t,2],nn[t,3],nn[t,4],nn[t,2]+nn[t,4],nn[t,1]+nn[t,3],nn[t,3]+nn[t,4],nn[t,1]+nn[t,2],nn[t,1]+nn[t,2]+nn[t,3]+nn[t,4]);printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",moa,"Any",nnn[1],nnn[2],nnn[3],nnn[4],nnn[2]+nnn[4],nnn[1]+nnn[3],nnn[3]+nnn[4],nnn[1]+nnn[2],nnn[1]+nnn[2]+nnn[3]+nnn[4]);}' moa=$moa GM=$GM0 | sort  >> RESULTS/$MOA_genes/Venn.MOA_gene_types.txt
end
cat toto.$$   | gawk -F '\t' '/^#/{next;}{g=$4;gid=$5;ma=$6;t=$14;if(GM=="MRNAH"){g=$5;gid=$6;t=$7;ma=$9;}z=0;if(length(g)<1)next;if(gid>0)z=1;tt[t]=1;if(length(ma)>1)g2ma[g]=2;if(g==ma)next;gg[g]=1;gt[g,t]=z+1;}END{for(g in gg)for(t in tt)if(gt[g,t]>0)nn[t,gt[g,t]+g2ma[g]]++;for(t in tt)for(i=1;i<=4;i++)nnn[i]+=nn[t,i];for(t in tt)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",moa,t,nn[t,1],nn[t,2],nn[t,3],nn[t,4],nn[t,2]+nn[t,4],nn[t,1]+nn[t,3],nn[t,3]+nn[t,4],nn[t,1]+nn[t,2],nn[t,1]+nn[t,2]+nn[t,3]+nn[t,4]);printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",moa,"Any",nnn[1],nnn[2],nnn[3],nnn[4],nnn[2]+nnn[4],nnn[1]+nnn[3],nnn[3]+nnn[4],nnn[1]+nnn[2],nnn[1]+nnn[2]+nnn[3]+nnn[4]);}' moa=Union GM=$GM0 | sort  >> RESULTS/$MOA_genes/Venn.MOA_gene_types.txt
\rm  toto.$$

gunzip -c TARGET/Targets/hs.av.fasta.gz | gawk '/^>/{split($1,aa,"|");print aa[3];}' | sort -u > titi
cat titi ZZZZZ TARGET/GENES/av.geneTypesAntiSense.txt  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{if(ok[$1]==1)n[$2]++}END{for(k in n)printf("%s\t%d\n",k,n[k]);}' | sort  >> RESULTS/$MOA_genes/Venn.MOA_gene_types.txt


########################
## histogram of index and of delta in the 4 Affy and av

foreach target ($Target1 $Target2)
  if (-e tmp/$MOA_genes/tutu.$target) \rm tmp/$MOA_genes/tutu.$target
  echo XXXXX >> tmp/$MOA_genes/tutu.$target
  set sign=-1
  \rm DEG.$target.list
  set GM=MA
  if ($target == av || $target == RefSeq || $target == EBI) set GM=$GM0
  foreach moa (`cat tmp/MOA/MOA.list`)
    @ sign = - $sign
    if ($sign == 1) then
      set ff1=$moa
      cat tmp/MOA/experiments.$GM0.txt | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\t%s\t%s\t%s\n",$2,$3,$4,$5,$6);}' | grep $moa > tmp/$MOA_genes/tutu7
    else
      cat tmp/MOA/experiments.$GM0.txt | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\t%s\t%s\t%s\n",$3,$2,$4,$5,$6);}' | grep $moa > tmp/$MOA_genes/tutu7
      set ff1=`cat tmp/MOA/experiments.txt | gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s\t%s\t%s_%s\n",$2,$3,$3,$2);}' | grep $moa |  gawk -F '\t' '{gsub(/\"/,"",$0);printf("%s_%s\n",$1,$2);}'`
    endif
    if ($target == av) then
      set limit=`cat tmp/$MOA_genes/tutu7 | cut -f 3`
    else
      set limit=`cat tmp/$MOA_genes/tutu7 | cut -f 3`
    endif
    echo $target,$moa,$sign,$limit
    set dd="RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.$target.$GM.u"
    echo  $dd.$moa.diffGenes.1.txt >> DEG.$target.list
    # echo $target,$sign,$ff1,$dd
    # echo $dd.$ff1.diffGenes.1.txt
    set p2g=MicroArray/Hits/AGLuK.av.probe2gene.unique.txt 
    if ($GM0 == MRNAH) set p2g=MicroArray/Hits/AGLuK.av.probe2mrna.unique.txt 
    echo " " > RESULTS/$MOA_genes/hs.av.transcript2gene.txt
    if ($GM0 == MRNAH)  gunzip -c TARGET/MRNAS/hs.av.transcript2gene.txt.gz | sed -e 's/\"//g' >> RESULTS/$MOA_genes/hs.av.transcript2gene.txt
    set GMS=GENES
    if ($GM0 == MRNAH)  set GMS=MRNAS
    if ($target == AGLuK) set GMS=AGL
    # add the meta data and reorder the AGL lines so that they are keyed on the gene

    cat  TARGET/GENES/av.gene2geneid.txt ZZZZZ  $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{print;next;}{if(zz<1){gid[$1]=gid[$1] ";" $2 ; next;}}{if($2<limit)next;delta=$11;if(GMS=="GENES"){delta=$12;}if(delta<=0)next;nn[delta]++;}END{for(d in nn)printf("Delta\t%s\t%s\t%.1f\t%d\n",moa,GMS,d,nn[d])}' limit=$limit sign=$sign  moa=$moa GMS=$GMS | sort -k 3n >> RESULTS/PrevalenceProfiles/index_profile
    cat  TARGET/GENES/av.gene2geneid.txt ZZZZZ  $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{print;next;}{if(zz<1){gid[$1]=gid[$1] ";" $2 ; next;}}{if($2<limit)next;delta=$11;idx=int(2*$12+.49)/2;if(GMS=="GENES"){delta=$12;idx=int(2*$13+.49)/2;}if(delta<=0)next;nn[idx]++;}END{for(d in nn)printf("Average of High group\t%s\t%s\t%.1f\t%d\n",moa,GMS,d,nn[d])}' limit=$limit sign=$sign  moa=$moa GMS=$GMS | sort -k 3n >> RESULTS/PrevalenceProfiles/index_profile
    cat  TARGET/GENES/av.gene2geneid.txt ZZZZZ  $dd.$ff1.diffGenes.1.txt | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{print;next;}{if(zz<1){gid[$1]=gid[$1] ";" $2 ; next;}}{if($2<limit)next;delta=$11;idx=int(2*$12+.49)/2;if(GMS=="GENES"){delta=$12;idx=int(2*$13+.49)/2;}if(delta=0)next;nn[idx]++;}END{for(d in nn)printf("Average of Low group\t%s\t%s\t%.1f\t%d\n",moa,GMS,d,nn[d])}' limit=$limit sign=$sign  moa=$moa GMS=$GMS | sort -k 3n >> RESULTS/PrevalenceProfiles/index_profile


 
 cat $dd.$ff1.diffGenes.1.txt | gawk -F '\t' '{print}' | cut -f 1,2,3 | head -8



  end
end

cat  RESULTS/PrevalenceProfiles/index_profile | grep Delta | gawk -F '\t' '/^#/{next}{dd[$1]=1;mm[$2]=1;gm[$3]=1;nn[$4]=1;z[$1,$2,$3,$4]=$5;}END{for(d in dd) for (m in mm) for(g in gm)printf("\t%s %s %s",d,m,g);for(n in nn){if(n+0>0)printf("\n%.1f",n);for(d in dd) for (m in mm) for(g in gm)printf("\t%d",z[d,m,g,n]);}printf("\n");}' | sort -k 1n >  RESULTS/PrevalenceProfiles/index_profile.histo.txt
cat  RESULTS/PrevalenceProfiles/index_profile | grep -v Delta | gawk -F '\t' '/^#/{next}{dd[$1]=1;mm[$2]=1;gm[$3]=1;nn[$4]=1;z[$1,$2,$3,$4]=$5;}END{for(d in dd) for (m in mm) for(g in gm)printf("\t%s %s %s",d,m,g);for(n in nn){if(n+0>0)printf("\n%.1f",n);for(d in dd) for (m in mm) for(g in gm)printf("\t%d",z[d,m,g,n]);}printf("\n");}' | sort -k 1n >>  RESULTS/PrevalenceProfiles/index_profile.histo.txt


###################################################################
# Triple delight


set toto=RESULTS/MOA.Gene.MRNAH.AGLuK.txt

cat RESULTS/MOA.av.AGLuK.GENE/MOA.av.AGLuK.GENE.FDRlimit.fused.txt | gawk -F '\t' '/^#/{next;}/^\(\+\)eQC-42/{next;}{g=$4;gid[g]=$5 "\t" $15 "\t" $14 "\t" $13 ;for(i=1;i<=nMoa/2;i++)z[g,i]=$(i+16);}END{for(g in gid) {printf("GENE\t%s\t%s\t%s",g,g,gid[g]);for(i=1;i<=nMoa/2;i++)printf("\t%s",z[g,i]);printf("\n");}}' nMoa=$nMoa >  $toto.1

cat RESULTS/MOA.av.AGLuK.MRNAH/MOA.av.AGLuK.MRNAH.FDRlimit.fused.with_metadata.txt | gawk -F '\t' '/^#/{next;}{t=$1;g=$3;gid2[g]=$4 "\t" $2 "\t" $5 "\t" $13 ;ng2t[g]++;g2t[g,ng2t[g]]=t;for(i=1;i<=nMoa/2;i++)z[t,i]=$(i+16);}END{for(g in ng2t)for(k=1;k<=ng2t[g];k++){t=g2t[g,k];printf("MRNAH\t%s\t%s\t%s",g,t,gid2[g]);for(i=1;i<=nMoa/2;i++)printf("\t%s",z[t,i]);printf("\n");}}' nMoa=$nMoa >>  $toto.1

cat RESULTS/MOA.av.AGLuK.GENE/MOA.av.AGLuK.GENE.FDRlimit.fused.txt |  gawk -F '\t' '/^#/{next;}{p=$22;g=$4;gid3[g]=$5 "\t" $15 "\t" $14 "\t" $13 ;for(kk=0;kk<6;kk++){p=$(22+(4+nMoa/2)*kk);if(p){ng2p[g]++;g2p[g,ng2p[g]]=p;for(i=1;i<=nMoa/2;i++)z[p,i]=$(i+24+(4+nMoa/2)*kk);}}}END{for(g in ng2p)for(k=1;k<=ng2p[g];k++){p=g2p[g,k];printf("AGLuK\t%s\t%s\t%s",g,p,gid3[g]);for(i=1;i<=nMoa/2;i++)printf("\t%s",z[p,i]);printf("\n");}}' nMoa=$nMoa  >>  $toto.1

echo -n '# ' > $toto
date >> $toto

echo -n "#n-score\tTop-score\t1\t4\t4s\tMNA\tGene\tGeneId\tMap\tType\tTitle" >> $toto
foreach moa (Stg1 Stg4 Stg4s MNA)
  foreach gm (Gene MRNAH AGLuK)
     echo -n "\t$moa $gm plus" >> $toto
     echo -n "\t$moa $gm minus" >> $toto
  end
end
foreach moa (Stg1 Stg4 Stg4s MNA)
  foreach gm (Gene MRNAH AGLuK)
     echo -n "\t$moa $gm plus" >> $toto
     echo -n "\t$moa $gm minus" >> $toto
  end
end
echo >> $toto


# just count the number of + or - via the 3 methods
cat $toto.1 | gawk -F '\t' '{gma=substr($1,1,1);g=$2;if(substr(g,1,5)=="UKv4_")next;if(gid[g]=="" || (index($5,"chr")>0))gid[g]=$4 "\t" $5 "\t" $6 "\t" $7 ;for(i=1;i<=nMoa/2;i++){z=$(7+i);if(z+0<0){if(gi[g,i]<2)gi[g,i] += 2;ggg[g]++;gnc[g,gma,i]++;if(-z>gn[g,gma,i])gn[g,gma,i]=-z;if(-z>ts[g])ts[g]=-z;}if(z+0>0){ggg[g]++;if(gi[g,i]%2==0)gi[g,i] += 1;gpc[g,gma,i]++;if(z>gp[g,gma,i])gp[g,gma,i]=z;;if(z>ts[g])ts[g]=z;}}}END{for(g in gid) {if(ggg[g]<1)continue ;k=0;for(i=1;i<5;i++){if(gi[g,i]==3)k+=1000;if(gi[g,i]>=1)k+=100;}printf("%d\t%.1f\t%d\t%d\t%d\t%d\t%s\t%s",k+0*(gi[g,1]+gi[g,2]+gi[g,3]+gi[g,4]),ts[g],gi[g,1],gi[g,2],gi[g,3],gi[g,4],g,gid[g]);for(i=1;i<=nMoa/2;i++)for(igma=1;igma<=3;igma++){gma=substr("GMA",igma,1);printf("\t+%d\t-%d",gp[g,gma,i],gn[g,gma,i]);}for(i=1;i<=nMoa/2;i++)for(igma=1;igma<=3;igma++){gma=substr("GMA",igma,1);printf("\t+%d\t-%d",gpc[g,gma,i],gnc[g,gma,i]);}printf("\n");}}'  nMoa=$nMoa | sort  -k 1nr -k 3n -k 4n -k 5n -k 6n -k 2nr -k 7 | gawk -F '\t' '{a=$3 $4 $5 $6;if(a != old) printf("\n");old=a;print}'  >> $toto


if (0) then

  set sign=-1
  set nn=0
  if (-e $toto.1) \rm $toto.1
  foreach moa (`cat tmp/MOA/MOA.list`)
    @ sign = - $sign
    if ($sign < 0) continue
    @ nn = $nn + 1
    echo -n "$moa"  >> $toto.1

    # export gene type moa score on a line, then later create the nice looking table

    # self should be 1 when comparing AGL_at to AGL_microarray
    set self=1
    set slf=0

    cat RESULTS/$MOA_genes/$MOA_genes.FDRlimit.fused.txt | gawk -F '\t' '{line++;ok = 0;if(line<2 || NF <17)next;g=$4;if($6 && (self || $6 != $4))rnaAffy[g]=1;i=16+nn;z=$i+0;if(z!=0){ok=1;rna[g]=1;}for(iii=1;iii<=6;iii++){i=16+nn+iii*(nMoa/2+4);z=$i+0;if(z!=0){ok+=2;affy[g]=1;next;}}}END{for(g in rna){nR++;if(rnaAffy[g]>0)nR2++;}for(g in affy){if(rna[g]>0){nRA++;if(rnaAffy[g]>0)nRA2++;}else {nA++;if(rnaAffy[g]>0)nA2++;}}nR=nR-nRA;nR2=nR2-nRA2;printf("\t\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f",nR+1*nRA,nA+1*nRA,nR+nRA+nA,nRA,100.0*nRA/(.01+nR+nRA),100.0*nRA/(.01+nA+nRA),100*nR/(.01+nA));printf("\t\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f",nR2+1*nRA2,nA2+1*nRA2,nR2+nRA2+nA2,nRA2,100.0*nRA2/(.01+nR2+nRA2),100.0*nRA2/(.01+nA2+nRA2),100*nR2/(.01+nA2));}'  nMoa=$nMoa nn=$nn moa=$moa self=$self >> $toto.1

    echo >> $toto.1
  end
  cat $toto.1 >> $toto
  cat $toto.1 | gawk -F '\t' '{if (NF > nf)nf=NF;for(i=1;i<=NF;i++)s[i]+=$i;}END{printf("Sum");for(i=2;i<=nf-3;i++)printf("\t%d",s[i]);printf("\t\t\t\n");}' >> $toto
  echo >> $toto

cat $toto  | gawk -F '\t' '/^Sum/{if (NF > nf)nf=NF;for(i=1;i<=NF;i++)s[i]+=$i;}END{printf("SuperSum");for(i=2;i<=nf-3;i++)printf("\t%d",s[i]);printf("\t\t\t\n");}' > $toto.1
cat $toto.1 >> $toto

endif

\cp $toto RESULTS/$MOA_genes 

########################

####################################################################################
### Triple delight av, AGLuK, RNA_at_AGLuK


set dd=RESULTS/DEG_Lists_exported
set ff=$dd/Integrated_list.txt

set dd=RESULTS
set ff=$dd/MOA.Gene.MRNAH.AGLuK.txt

cat TARGET/GENES/av.metadata.txt ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;g2map[g]=$2;g2typ[g]=$4;g2title[g]=$5;g2gid[g]=$3;n=split($6,aa,";");for(i=1;i<=n;i++)t2g[aa[i]]=g;g2probes[g]=$7;n=split($7,aa,";");for(i=1;i<=n;i++)p2g[aa[i]]=g;g2np[g]=n;next;}}{g=$7;printf("%s",$1);for(i=2;i<=7;i++)printf("\t%s",$i);printf("\t%s",g2gid[g]);if(length($7)>0)printf("\t%d\t%s\t%s\t%s",0+g2np[g],g2map[g],g2typ[g],g2title[g]);else printf("\t\t\t\t");for(i=12;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $dd/Integrated_list_with_probes.txt

cat  TARGET/GENES/CancerCensus.txt ZZZZZ $dd/Integrated_list_with_probes.txt |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){c=$1 "::" ;if($3) c="S:" $3 ;if($4) c="G:" $4 ;if($3 && $4) c="S:" $3 "; G:" $4 ;if (length(c)>2){if($1)g2cancer[$1]=c;if($1)gid2cancer[$1]=c;}next;}}{if(length($7)>1 && $11=="")$11="single_exon";printf("%s",$1);for(i=2;i<=12;i++)printf("\t%s",$i);c="";n=split($8,aa,";");for(i=1;i<=n;i++){c1=gid2cancer[aa[i]];if(length(c1>1))c=c1;}c1=g2cancer[$7];if(length(c)==0 && length(c1)>1)c=c1;printf("\t%s",c);for(i=13;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $dd/Integrated_list_with_cancer.txt

# call type 31 if 3 in 3 in RNA-seq, 32 if 3 in Agilent, 33 if 3 in both, 4 if RNA-seq contradicts Agilent
cat $dd/Integrated_list_with_cancer.txt   |  gawk -F '\t' '/^#/{print;next;}{for(ii=0;ii<4;ii++){dx=13+6*ii;r=0;m=0;if($(dx+1)+$(dx+3)+0>0)r+=1;if($(dx+2)+$(dx+4)+0<0)r+=2;if($(dx+5)>0)m+=1;if($(dx+6)<0)m+=2;u=0;if(r==3){u=31;if(m==3)u=33;}else if(m==3){u=32;if(r==3)u=33;}else if(m*r==0)u=m+r;else if (m==r)u=m;else u=4;if($1)$(3+ii)=u;} printf("%s",$1); for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $dd/Integrated_list_with_cancer.4.txt

cat $dd/Integrated_list_with_cancer.4.txt  | gawk -F '\t' '/^#/{print;next;}{for(kk=0;kk<=4;kk++){if(kk==0)t=$3 "\t" $4 "\t" $5 "\t" $6;if(kk==1){if($3<1)continue;t=$3 "\t" "-" "\t" "-" "\t" "-";} if(kk==2){if($4<1)continue;t="-" "\t" $4 "\t" "-" "\t" "-";} if(kk==3){if($5<1)continue;t="-" "\t" "-" "\t" $5 "\t" "-"; }if(kk==4){if($6<1)continue;t="-" "\t" "-" "\t" "-" "\t" $6;} tt[t]=1;ng[t]++;$8=0;n=$9;if(n+0>6)$9=6;for(i=7;i<=NF;i++){nt[t,i]+=$i;if($i+0!=0)nt2i[t,i]++;}if(NF>nf)nf=NF;}}END{for(t in tt){x=2;if(index(t,"-")>0)x=1;printf("%d\t\t%s\t%d\t%d\t%d\t\t\t\t%d\t%d",x,t,ng[t],0+nt[t,8],0+nt[t,9],0+nt[t,14],0+nt[t,15]);for(i=16;i<=nf;i++)printf("\t%d in %d genes",0+nt[t,i],0+nt2i[t,i]);printf("\n");}}' | sort -k 1,1n -k 3,3n -k 4,4n -k 5,5n -k 6,6n > $dd/Integrated_stats.txt


cat TARGET/GENES/av.metadata.txt ZZZZZ toto23 |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g=$1;g2typ[g]=$4;g2gid[g]=$3;n=split($6,aa,";");for(i=1;i<=n;i++)t2g[aa[i]]=g;g2probes[g]=$7;n=split($7,aa,";");for(i=1;i<=n;i++)p2g[aa[i]]=g;g2np[g]=n;next;}}{g=$1;if(g2typ[g])print g,g2typ[g];}'

cat $dd/Integrated_list_with_cancer.4.txt | gawk -F '\t' '/^#/{print;next;}' > $dd/Integrated_stats.gene_type.txt
cat $dd/Integrated_list_with_cancer.4.txt | gawk -F '\t' '/^#/{next;}{if($11=="")next;for(kkk=0;kkk<2;kkk++)for(kk=0;kk<=4;kk++){if(kkk==1)$11="Any";if(kk==0)t=$11 "\t" "=" "\t" "-" "\t" "-" "\t" "-" ;if(kk==1){if($3<1)continue;t=$11 "\t" $3 "\t" "-" "\t" "-" "\t" "-" ;} if(kk==2){if($4<1)continue;t=$11 "\t" "-" "\t" $4 "\t" "-" "\t" "-" ; }if(kk==3){if($5<1)continue;t=$11 "\t" "-" "\t" "-" "\t" $5 "\t" "-" ;} if(kk==4){if($6<1)continue;t=$11 "\t" "-" "\t" "-" "\t" "-" "\t" $6 ;}  tt[t]=1;ng[t]++;if($13)nc[t]++;if(kk==0 && kkk==0)$8=split($8,aa,";");n=$9;if(n+0>6)$9=6;for(i=7;i<=NF;i++){nt[t,i]+=$i;if($i && $i!=0)nt2i[t,i]++;}}if(NF>nf)nf=NF;}END{for(t in tt){if(index(t,"0")>0)continue;x=4;if(index(t,"=")>0)x=3;if(index(t,"Any")>0)x-=2;printf("%d\t%s\t%d\t%d in %d genes\t%d in %d genes\t\t\t\t%d",x,t,ng[t],0+nt[t,8],0+nt2i[t,8],0+nt[t,9],0+nt2i[t,9],nc[t]);for(i=14;i<=nf;i++)printf("\t%d in %d genes",0+nt[t,i],0+nt2i[t,i]);printf("\n");}}' | sort -k 1,1n -k 3,3n -k 4,4n -k 5,5n -k 6n,6n -k 2,2 >> $dd/Integrated_stats.gene_type.txt

###########################################################
### Venn diagram of Consistency between the GENE MRNAH AGLuK
### new method 2014_03_05, directly from the integrated table

set toto=RESULTS/Venn.Integrated_list.DEG_withDETorDEP.txt
echo -n "# " > $toto
date >> $toto
 cat RESULTS/Integrated_list_with_cancer.4.txt | gawk -F '\t' -f scripts/ma4.counts.awk avOnly=0 >> $toto

set toto=RESULTS/Venn.Integrated_list.DEG_withDETorDEP.withGid.txt
echo -n "# " > $toto
date >> $toto
 cat RESULTS/Integrated_list_with_cancer.4.txt | gawk -F '\t' -f scripts/ma4.counts.awk avOnly=1 >> $toto

set toto=RESULTS/Venn.Integrated_list.DEG_withDETorDEP.avOnly.txt
echo -n "# " > $toto
date >> $toto
 cat RESULTS/Integrated_list_with_cancer.4.txt | gawk -F '\t' -f scripts/ma4.counts.awk avOnly=2 >> $toto

###########################################################
### Venn diagram of Consistency between the GENE MRNAH AGLuK



set dd=RESULTS

# call type 31 if 3 in 3 in RNA-seq, 32 if 3 in Agilent, 33 if 3 in both, 4 if RNA-seq contradicts Agilent
echo -n "# " >  $dd/Integrated_Venn.txt
date >>  $dd/Integrated_Venn.txt
cat $dd/Integrated_list_with_cancer.txt   |  gawk -F '\t' '/^#/{next;}{if($11=="")next;for(ii=0;ii<4;ii++){dx=13+6*ii;g=0;t=0;a=0;if($(dx+1)>0)g+=1;if($(dx+2)<0)g+=2;if($(dx+3)>0)t+=1;if($(dx+4)<0)t+=2;if($(dx+5)>0)a+=1;if($(dx+6)<0)a+=2;nn[ii,g,t,a]++; }}END{printf("\tGene\tTranscript\tAgilent\tGT\tTA\tAG\tGTA");split("Stage 1;Stage 4;Stage 4s;MNA",aa,";");for(ii=0;ii<4;ii++)for(k=1;k<=2;k++){printf("\n%s Type %d",aa[ii+1],k);printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",nn[ii,k,0,0],nn[ii,0,k,0],nn[ii,0,0,k],nn[ii,k,k,0],nn[ii,0,k,k],nn[ii,k,0,k],nn[ii,k,k,k],nn[ii,k,0,0]+nn[ii,0,k,0]+nn[ii,0,0,k]+nn[ii,k,k,0]+nn[ii,0,k,k]+nn[ii,k,0,k]+nn[ii,k,k,k]);}printf("\n");}' >> $dd/Integrated_Venn.txt
cat $dd/Integrated_list_with_cancer.txt   |  gawk -F '\t' '/^#/{next;}{if($11=="")next;if($9+0<1)next;for(ii=0;ii<4;ii++){dx=13+6*ii;g=0;t=0;a=0;if($(dx+1)>0)g+=1;if($(dx+2)<0)g+=2;if($(dx+3)>0)t+=1;if($(dx+4)<0)t+=2;if($(dx+5)>0)a+=1;if($(dx+6)<0)a+=2;nn[ii,g,t,a]++; }}END{printf("\tGene\tTranscript\tAgilent\tGT\tTA\tAG\tGTA");split("Stage 1;Stage 4;Stage 4s;MNA",aa,";");for(ii=0;ii<4;ii++)for(k=1;k<=2;k++){printf("\n%s Type %d with probe",aa[ii+1],k);printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",nn[ii,k,0,0],nn[ii,0,k,0],nn[ii,0,0,k],nn[ii,k,k,0],nn[ii,0,k,k],nn[ii,k,0,k],nn[ii,k,k,k],nn[ii,k,0,0]+nn[ii,0,k,0]+nn[ii,0,0,k]+nn[ii,k,k,0]+nn[ii,0,k,k]+nn[ii,k,0,k]+nn[ii,k,k,k]);}printf("\n");}' >> $dd/Integrated_Venn.txt

# Top genes with most differential transcripts


set toto=$dd/MOA.mostDiffAlternativeTranscripts.txt
echo -n "#Date " > $toto
date >> $toto
cat $dd/Integrated_list_with_cancer.4.txt | gawk -F '\t' '/^#/{printf("#Class\tBest delta\tRNA-seq delta\tAgilent delta\t"); print ; next;}' >> $toto

if (-e $toto.1) \rm $toto.1

  set sign=-1
  set nn=0

  foreach moa (`cat tmp/MOA/MOA.list`)
    @ sign = - $sign
    if ($sign < 0) continue
    @ nn = $nn + 1

    cat $dd/Integrated_list_with_cancer.4.txt  | gawk -F '\t' '/^#/{ next }{if($1 < 1000) next; t=$(2+nn) ; if(t != 31 && t != 32) next ; g=$7 ; ngid = index($8,";"); if (ngid > 1000) next;  rnaP = $(10 + 6*nn); rnaN = $(11 + 6*nn); aglP = $(12 + 6*nn); aglN = $(13 + 6*nn) ; rnaD = aglD = 0 ;  if (rnaP * rnaN < 0) rnaD = rnaP - rnaN  ; if (aglP * aglN < 0) aglD = rnaP - rnaN  ; if(rnaD<0)rnaD = -rnaD;if(aglD<0)aglD = -aglD; if(0) print g,t,rnaP,rnaN,rnaD; if (rnaD > MN || aglD > MN) { x=rnaD ; if (aglD > x)x=aglD; printf("%s\t%d\t%d\t%d\t",moa,x,rnaD,aglD); print ;}}' nn=$nn moa=$moa MN=60 >> $toto.1


end

cat $toto.1 | sort -k 1,1 -k 2,2nr >> $toto
\rm $toto.1


# grab the histograms
\rm  RESULTS/MOA.mostDiffAlternativeTranscripts.*.histos.txt
cat $toto | gawk -F '\t' '{m=substr($1,2);i=index(m,"Stg");m2=substr(m,i) "_" substr($1,1,i-1);printf("cat RESULTS/Expression/unique/av/Diff_genes/NB.av.MRNAH.u.%s.diffGenes.1.txt   | grep %s >> RESULTS/MOA.mostDiffAlternativeTranscripts.%s.histos.txt\n",$1,$11,$1);}' > titatou2
cat $toto | gawk -F '\t' '{m=substr($1,2);i=index(m,"Stg");m2=substr(m,i) "_" substr($1,1,i-1);printf("cat RESULTS/Expression/unique/av/Diff_genes/NB.av.MRNAH.u.%s.diffGenes.1.txt   | grep %s >> RESULTS/MOA.mostDiffAlternativeTranscripts.%s.histos.txt\n",m2,$11,m2);}' >> titatou2
source titatou2



goto phaseLoop

##############################################
# histos, i put them here for the moment

foreach target ($Etargets)
  set target2=$target
  if($target == av) set target2=AceView
  gunzip -c RESULTS/Expression/unique/$target2/$MAGIC.$target.GENE.u.ace.gz | gawk 'BEGIN{print "0\n400";}/^Run_U/{if($11 != "ok" || 0+$3 < 0)next; n++;z=0+$3;print 5*int(2*z);}' | bin/histo -plot -w 80 -plain -o RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo
  gunzip -c RESULTS/Expression/unique/$target2/$MAGIC.$target.GENE.u.ace.gz | gawk 'BEGIN{print "0\n400";}/^Run_U/{if($11 != "NA" || 0+$3 < 0)next; n++;z=0+$3;print 5*int(2*z);}' | bin/histo -plot -w 80 -plain -o RESULTS/Expression/unique/$MAGIC.$target.just_NA.newindex.histo
end
foreach target ($Etargets)
  set target2=$target
  if($target == av) set target2=AceView
  \cp RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo.txt tutu
  tail -n +3 tutu | gawk -F '\t' '/^400/{next;}{x=(0+$1)/10;n=0+$2;nn+=n;xx+=n*x;}END{printf("Average\t%s\n",xx/nn);}' >> RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo.txt
  \rm RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo.correl.txt
  \rm RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo.ps
  \rm RESULTS/Expression/unique/$MAGIC.$target.no_NA.newindex.histo.tcsh
end


exit 0

#### enddraftVenn  $MOA_genes
## same but we cut by MOA

foreach moa (`cat tmp/MOA/MOA.list`)

  cat tmp/MOA.av_newindex_$th.AGLuK_rescaled80.newindex_10/tutu.av_newindex_$th | grep $moa | gawk -F '\t' '/^#/{next}/Gene/{next;}{if($3>=40)print $1}' | sort -u > tata.$moa.av 
cat tmp/MOA.av_newindex_$th.AGLuK_rescaled80.newindex_10/tutu.AGLuK_rescaled80.newindex_10 |   grep $moa | gawk -F '\t' '/^#/{next;}/Gene/{next;}{if($3>=40)print $1}' | sort -u >  tata.$moa.AGLuK_rescaled80
  cat tmp/MOA.av_newindex_$th.AGLuK_at_newindex_10/tutu.AGLuK_at_newindex_10 |  grep $moa | gawk -F '\t' '/^#/{next}/Gene/{next;}{if(0 && $3>=40)print $1}' | sort -u > tata.$moa.AGLuK_at   

# Venn diagram of probes and gene between av_th12, AGLuK_at_th10 and AGLuK_rescaled80_th10
  wc tata.$moa.*
 
# passons en gene
  cat  MicroArray/Hits/AGLuK.av.probe2gene.unique.txt ZZZZZ tata.$moa.AGLuK_at         | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){p2g[$1]=$2;next;}p=$1;g=p;if(p2g[$1])g=p2g[$1];print g;}'| sort -u > tata.$moa.AGLuK_at.genes
  cat  MicroArray/Hits/AGLuK.av.probe2gene.unique.txt ZZZZZ tata.$moa.AGLuK_rescaled80 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){p2g[$1]=$2;next;}p=$1;g=p;if(p2g[$1])g=p2g[$1];print g;}'| sort -u > tata.$moa.AGLuK_rescaled80.genes
cat  tata.$moa.av | sort -u > tata.$moa.av.genes

  wc tata.$moa.*.genes

#Venn
  cat  tata.$moa.AGLuK_at.genes ZZZZZ  tata.$moa.AGLuK_rescaled80.genes ZZZZZ tata.$moa.av.genes |  gawk -F '\t' 'BEGIN{zz=1;}/^#/{next}/^ZZZZZ/{zz=2*zz;next;}{if($1)gg[$1]+=zz;}END{printf("8 Class\tNumber of genes\tAGLuK_at limit 10\tAGLuK_rescaled_80 limit 10\tAceView genes limit 10");for(g in gg){nnn++;nn[gg[g]]++;} for(i=1;i<=7;i++){printf("\n%d\t%d",i,nn[i]);a="";b="";c="";if(i%2==1)a="*";if(int(i/2)%2==1)b="*";if(int(i/4)%2==1)c="*";printf("\t%s\t%s\t%s",a,b,c);}printf("\n0\t%d\t.\t.\t.\n",nnn);}' | sort -k 1nr > RESULTS/$MOA_genes/Venn/Venn.gene_counts.$moa.txt

#Actual lists
  cat  tata.$moa.AGLuK_at.genes ZZZZZ  tata.$moa.AGLuK_rescaled80.genes ZZZZZ tata.$moa.av.genes |  gawk -F '\t' 'BEGIN{zz=1;}/^#/{next}/^ZZZZZ/{zz=2*zz;next;}{gg[$1]+=zz;}END{printf("# Class\tGenes\tAGLuK_at limit 10\tAGLuK_rescaled_80 limit 10\tAceView genes limit 10");for(g in gg){i=gg[g];printf("\n%d\t%s",i,g);a="";b="";c="";if(i%2==1)a="*";if(int(i/2)%2==1)b="*";if(int(i/4)%2==1)c="*";printf("\t%s\t%s\t%s",a,b,c);}printf("\n");}' | sort -k 1nr > tata.$moa.all_genes

  cat  TARGET/GENES/av.gene2geneid.txt ZZZZZ TARGET/GENES/av.geneTypesAntiSense.txt  ZZZZZ TARGET/GENES/av.gene2title.txt ZZZZZ TARGET/GENES/av.gene2intMap.txt MicroArray/Hits/AGLuK.probe2chrom.txt ZZZZZ tata.$moa.all_genes | gawk -F '\t' 'BEGIN{printf("# Class\tAGLuK_at limit 10\tAGLuK_rescaled_80 limit 10\tAceView genes limit 10\tGene\tMap\tGeneId\tType\tTitle\n");}/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){gid[$1]=gid[$1] ";" $2 ; next;}}{if(zz==1){gType[$1]=gType[$1] ";" $2; next;}}{if(zz==2){gTitle[$1]=$2; next;}}{if(zz==3){if(NF>2)gChrom[$1]= "chr" $2 ":" $3 "-" $4; else gChrom[$1]=$2;next;}}{g=$2;if(g){typ=substr(gType[g],2);if(length(typ)<1){if(substr(g,1,3)=="UKv"){if(gChrom[g])typ="Probe not in gene";else typ="Unmapped probe";}else typ="Single exon";}printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$3,$4,$5,$2,gChrom[g],substr(gid[g],2),typ,gTitle[g]);}}' > tata.$moa.all_genes.txt

  \cp  tata.$moa.all_genes.txt RESULTS/$MOA_genes/Venn/Venn.gene_list.$moa.txt

end

# clean up parasites
\rm RESULTS/$MOA_genes/Venn/Venn.gene_list.Stg1_MYCN_not_amplified_117_Stg_AnyNot1noMNA_256.txt RESULTS/$MOA_genes/Venn/Venn.gene_list.Stg4s_MYCN_not_amplified_48_Stg_AnyBut4SnoMNA_325.txt RESULTS/$MOA_genes/Venn/Venn.gene_list.Stg_AnyBut4noMNA_257_Stg4_MYCN_not_amplified_116.txt RESULTS/$MOA_genes/Venn/Venn.gene_list.Stg_AnyNotMNA_281_StgAny_MYCN_amplified_92.txt 


  
foreach target ($Etargets)
  echo "$target nb_genes_not_NA "
  set target2=$target
  if($target == av) set target2=AceView
  ls -ls RESULTS/Expression/unique/$target2/$MAGIC.$target.GENE.u.ace.gz
continue
  gunzip -c RESULTS/Expression/unique/$target2/$MAGIC.$target.GENE.u.ace.gz | gawk 'BEGIN{if(0)print "0\n400";}/^Gene /{g=$2;next;}/^Run_U/{if($11 != "ok" || 0+$3 < 0)next; gg[g]++;if(gg[g]==1){n++;if(1)printf("Gene\t%s\n",g);}}END{if (0) print n}' | grep -v 'ERCC-' > RESULTS/Expression/unique/$MAGIC.$target.GENE.list
end

set toto=RESULTS/Expression/unique/$MAGIC.transcript.counts.txt
date > $toto
foreach target ($Etargets)
  set target2=$target
  if($target == av) set target2=AceView
  echo -n "$target nb_transcripts_not_NA " >> $toto
  gunzip -c  RESULTS/Expression/unique/$target2/$MAGIC.$target.MRNAH.u.ace.gz  |  gawk '/^Tran /{g=$2;next;}/^Run_U/{if($11 != "ok" || 0+$3 < 0)next; gg[g]++;if(gg[g]==1){n++;if(0) print g,n;}}END{print "\t",n} '>> $toto
end
foreach target ($Etargets)
  gunzip -c TARGET/Targets/$species.$target.fasta.gz | gawk '/^>/{n++}END{printf("Nb transcripts in %s\t%d\n",target,n);}' >> $toto
end

# the program below counts introns supported at least n times in at least one run
# the code  histo of intron supports with origin of target in geneindex.tcsh uses the cumul which is more interesting
# grep fot 'true log' in geneindex.tcsh

# av  223292
# RefSeq nb_transcripts_not_NA Found  36654
# EBI nb_transcripts_not_NA Found  154668

 
goto phaseLoop

########################
## hack Venn gene_types statistics reusing the files exported for Matthias, probably this step is useless, the next step is better
set dd=RESULTS/DEG_Lists_exported
set toto=$dd/stats.txt
echo -n "# " > $toto
date >> $toto

foreach GM (Gene transcripts Agilent)
  set ff=$dd/$GM.txt
  if ($GM == Agilent) set ff=$dd/Gene.txt
  foreach i1 (1 9 17 25)
    set i2=$i1
    if ($GM == Agilent)  then
      @ i1 = $i1 + 1 
      @ i2 = $i1 + 5 
    endif
    echo $GM $i1 $i2

    foreach hasGid (0 1 2)
     cat TARGET/GENES/av.metadata.txt ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){g=$1;g2map[g]=$2;g2typ[g]=$4;g2gid[g]=$3;n=split($6,aa,";");for(i=1;i<=n;i++)t2g[aa[i]]=g;g2probes[g]=$7;n=split($7,aa,";");for(i=1;i<=n;i++)p2g[aa[i]]=g;next;}}{if(zz<1)next;}/^Gene/{title=$i1;next;}/^Transcript/{title=$i1;next;}{for(i=i1;i<=i2;i++){g=$i;if(GM=="transcripts")g=t2g[g];if(GM=="Agilent"){p=g;g=p2g[g];}if(length(g)>0){gg[g]=1;pp[p]=1;}}next;}END{has[0]="";has[1]=" with GeneId"; has[2]=" no GeneId" ;printf("Title\t%s %s\n",title, has[hasGid]);for(g in gg){if(hasGid==0 || (hasGid==1 && g2gid[g]) || (hasGid==2 && g2gid[g]=="")){ng++;nt[g2typ[g]]++;n=split(g2probes[g],aa,";");ng2p[g2typ[g]]+=n;}}for(p in pp){g=p2g[p];if(hasGid==0 || (hasGid==1 && g2gid[g]) || (hasGid==2 && g2gid[g]=="")){npp[g2typ[g]]++;}}for(t in nt){ng1+=nt[t];printf("%s\t%d\t%d\t%d\n",t,nt[t],ng2p[t],npp[t]);}printf("zNo type\t%d\n",ng-ng1);printf("zAny\t%d\n",ng)}' i1=$i1 i2=$i2 GM=$GM hasGid=$hasGid | sort >> $toto
      echo >> $toto
    end
  end
end

foreach GM (Gene transcripts Agilent)
  if ($GM == Agilent) then
    cat $dd/Gene.txt | cut -f 2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,23,26,27,28,29,30 > $dd/$GM.union
  else
    cat $dd/$GM.txt | cut -f 1,9,17,25 > $dd/$GM.union
  endif
end

foreach GM (Gene transcripts Agilent)
  set ff=$dd/$GM.union
  set i1=1
  set i2=100

    foreach hasGid (0 1 2)
      cat TARGET/GENES/av.metadata.txt ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){g=$1;g2typ[g]=$4;g2gid[g]=$3;n=split($6,aa,";");for(i=1;i<=n;i++)t2g[aa[i]]=g;g2probes[g]=$7;n=split($7,aa,";");for(i=1;i<=n;i++)p2g[aa[i]]=g;next;}}{if(zz<1)next;}/^Gene/{title=$i1;next;}/^Transcript/{title=$i1;next;}{for(i=i1;i<=i2;i++){g=$i;if(GM=="transcripts")g=t2g[g];if(GM=="Agilent"){p=g;g=p2g[g];}if(length(g)>0){gg[g]=1;pp[p]=1;}}next;}END{has[0]="";has[1]=" with GeneId"; has[2]=" no GeneId" ;printf("Title\tUnion %s %s\n",GM, has[hasGid]);for(g in gg){if(hasGid==0 || (hasGid==1 && g2gid[g]) || (hasGid==2 && g2gid[g]=="")){ng++;nt[g2typ[g]]++;n=split(g2probes[g],aa,";");ng2p[g2typ[g]]+=n;}}for(p in pp){g=p2g[p];if(hasGid==0 || (hasGid==1 && g2gid[g]) || (hasGid==2 && g2gid[g]=="")){npp[g2typ[g]]++;}}for(t in nt){ng1+=nt[t];printf("%s\t%d\t%d\t%d\n",t,nt[t],ng2p[t],npp[t]);}printf("zNo type\t%d\n",ng-ng1);printf("zAny\t%d\n",ng)}' i1=$i1 i2=$i2 GM=$GM hasGid=$hasGid | sort >> $toto
      echo >> $toto
    end
 
end

cat $toto | sed -e 's/probe 1//' | gawk -F '\t' '/^Title/{ti=$2;nti++;tis[nti]=ti;next;}{t=$1;if(nti==1){ii++;tt[ii]=t;}nt[ti,t]=$2;}END{for(i=1;i<=nti;i++)printf("\t%s",tis[i]);for(j=1;j<=ii;j++){t=tt[j];printf("\n%s",t);for(i=1;i<=nti;i++)printf("\t%d",nt[tis[i],t]);}printf("\n");}'> RESULTS/DEG_Lists_exported/jolie_stats.txt



goto phaseLoop
#################################################################

# construct the metadata file 

if (! -e TARGET/GENES/av.metadata.txt) then
  tbly ~/mm_37_2/ZZ <<EOF
    find gene
    show -a -f TARGET/GENES/toto.ace
    quit
EOF
  echo "Gene toto" >> TARGET/GENES/toto.ace
  cat  TARGET/GENES/toto.ace | gawk '{gsub(/\"/,"",$0);}/^Gene /{if (g)printf("%s\t%s\t%s\t%s\t%s\n",g,substr(gid,2),map,type,title);g=$2;gid="";map="";type="";title="";next}/^Title/{title=substr($0,7);next;}/^GeneId/{gid=gid ";" substr($0,8);next;}/^IntMap/{map=$2 ":" $3 "-" $4;next;}' > TARGET/GENES/av.metadata.txt 
endif

#################################################################
#################################################################
phaseLoop: 
  echo done


