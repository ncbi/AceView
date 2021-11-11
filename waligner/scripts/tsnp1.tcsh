#!bin/tcsh -f

set run=$1
set zone=$2
set MB=$3
set Strategy=$4
set target_fasta=$5
set method=$6
set part=$7
set t1=$8
set t2=$9


set SAM=""
if ($MB == MB) then
  set SAM="-SAM "
  set method=MagicBlast
endif

set nAna=1

set chrom=`echo $zone | sed -e 's/zone[rG]\.//' -e 's/_[a-z]$//'` 

set minSnpFrequency2=$minSnpFrequency
if (-e tmp/TSNP/$MAGIC.minFrequency.txt) then
  # reset minSnpFrequency to at least twice the percent error rate, the error_profile is given in error per million 
  set z=`cat tmp/TSNP/$MAGIC.minFrequency.txt | gawk '{if ($1 == run) z=($4 + 0)/5000 ; if (mf<z)mf=z;}END{print int(mf);}' mf=$minSnpFrequency run=$run`
  set minSnpFrequency2=$z
endif

set trgt=""
if ($Strategy == RNA_seq) then
  set trgt="-target_class ET_av"
else
  set t12=`echo "x $t1 $t2" | gawk '{z="";if(NF==3)z=" -t1 " $2 " -t2 " $3 ;print z}'`
  echo $t12
  set trgt="-target_class Z_genome  -t $chrom $t12"
endif

echo "bin/tricoteur -run $run -laneList tmp/TSNP/$run/LaneList.$part $trgt -target_fasta $target_fasta -method $method -minSnpCover $minSnpCover -minSnpCount $minSnpCount -minSnpFrequency $minSnpFrequency2 -dx 8 -minAliPerCent 70 -o tmp/TSNP/$run/$zone/tsnp1.$MB.$part -nAna $nAna -qc -uno $SAM"
      bin/tricoteur -run $run -laneList tmp/TSNP/$run/LaneList.$part $trgt -target_fasta $target_fasta -method $method -minSnpCover $minSnpCover -minSnpCount $minSnpCount -minSnpFrequency $minSnpFrequency2 -dx 8 -minAliPerCent 70 -o tmp/TSNP/$run/$zone/tsnp1.$MB.$part -nAna $nAna -qc -uno $SAM
touch 

exit 0


######################################################################################
##################   CORONA dedicated scripts june 2020

bin/tricoteur -run COV-20200315-P10-C09-P -laneList tmp/TSNP/COV-20200315-P10-C09-P/LaneList.3 -target_class Z_genome -t NC_045512  -target_fasta tmp/SNP_ZONE/zoneG.NC_045512_a.fasta.gz -method MagicBlast -minSnpCover 100 -minSnpCount 6 -minSnpFrequency 5 -dx 8 -o tmp/TSNP/COV-20200315-P10-C09-P/zoneG.NC_045512_a/tsnpMB1.3 -nAna 1 -SAM
bin/tricoteur -run COV-20200315-P10-C09-P -laneList tmp/TSNP/COV-20200315-P10-C09-P/LaneList.3 -target_class Z_genome -t NC_045512  -target_fasta tmp/SNP_ZONE/zoneG.NC_045512_a.fasta.gz -method MagicBlast -minSnpCover 100 -minSnpCount 6 -minSnpFrequency 5 -dx 8 -o titi -nAna 1 -SAM
bin/tricoteur -run COV-20200315-P10-C09-P -laneList tmp/TSNP/COV-20200315-P10-C09-P/LaneList.3 -target_class Z_genome -t NC_045512  -target_fasta tmp/SNP_ZONE/zoneG.NC_045512_a.fasta.gz -method MagicBlast -minSnpCover 100 -minSnpCount 6 -minSnpFrequency 5 -dx 8 -o titi2 -nAna 1 


# essai: 111 snp vu chez NYC et dans TRicoteur, on veut extraire dans les 2 cas la liste des run ou on le voir nous et eux et faire le Venn

cat RESULTS/SNV_111_NYC-Magic.txt | sed -e 's/\r//g' > _a
# cherche la liste chez nous ou on a vu le snp a > 90%
cat _a  ZZZZZ RESULTS/SNV/MMM.zoneG.NC_045512_a.snp_frequency_table.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}


############################### compare the list of SNPS in NYC
cat variants_June5_vcf_jor_vars.clades.table.txt | cut -f 2 | sort -u > _a
cat MetaDB/MMM/RunsList | sed -e 's/COV-/COVHA-/' -e 's/_/-/g' | sort -u > _b

echo '# 1: runs just in variants_June5_vcf_jor_vars.clades.table.txt, 2: just in Magic analysis, 3: in both' > runs2NYCsnps.txt
cat _a ZZZZZ _b | gawk '/ZZZZZ/{zz++;next;}{n[$1]+=1+zz;}END{for (k in n) printf("%s\t%d\n",k,n[k]);}' | sort -k 2n >> runs2NYCsnps

# project the NYC table to kkep only the runs we analyzed
cat runs2NYCsnps ZZZZZ variants_June5_vcf_jor_vars.clades.table.txt | gawk -F '\t' '/ZZZZZ/{zz++;next;}{if(zz<1){if($2 >= 2)ok[$1]=1;next;}if(ok[$2]+0>=0)print}' >  variants_June5_vcf_jor_vars.clades.table.limited_to_runs_analyzed_in_Magic.txt
# get the snp list
cat runs2NYCsnps ZZZZZ variants_June5_vcf_jor_vars.clades.table.txt | gawk -F '\t' '/ZZZZZ/{zz++;next;}{if(zz<1){if($2 >= 2)ok[$1]=1;next;}if(ok[$2]+0>=0)print}' | cut -f 11 | sort -u  > _c
# retag them as NYC_JOR
# i hand edit thenmes of the deletions
 cat _c | gawk '{split($1,aa,":");split(aa[2],bb,"-");split(aa[3],cc,">");if (length(aa[3])!=3 && length(cc[2])<length(cc[1]))printf("Variant \"NC_045512:%d:Del_%d:%s:%s\"\n",bb[1],length(cc[1])-length(cc[2]),cc[1],cc[2]);}' | sort -V | wc

Variant "NC_045512:683:Del_9:CTAAAGTCAT:C"
Variant "NC_045512:684:Del_9:TAAAGTCATT:T"
Variant "NC_045512:685:Del_9:AAAGTCATTT:A"
Variant "NC_045512:685:Del_9:AAAGTCATTT:A"
Variant "NC_045512:685:Del_9:AAAGTCATTT:A"
Variant "NC_045512:1831:Del_2:TAT:T"
Variant "NC_045512:5535:Del_3:AACA:A"
Variant "NC_045512:5909:Del_7:TATAAATT:T"
Variant "NC_045512:7813:Del_2:TTC:T"
Variant "NC_045512:11270:Del_16:ATGGTTGATACTAGTTT:A"
Variant "NC_045512:11538:Del_2:TTG:T"
Variant "NC_045512:15858:Del_2:TAA:T"
Variant "NC_045512:17521:Del_2:ATG:A"
Variant "NC_045512:20296:Del_67:AAATTAGAAGGCTATGCCTTCGAACATATCGTTTATGGAGATTTTAGTCATAGTCAGTTAGGTGGTTT:A"
Variant "NC_045512:21451:Del_5:GGTCAA:G"
Variant "NC_045512:21452:Del_5:GTCAAA:G"
Variant "NC_045512:21463:Del_2:GAT:G"
Variant "NC_045512:21565:Del:GT:G"
Variant "NC_045512:21859:Del_3:CATA:C"
Variant "NC_045512:21956:Del_4:GAATT:G"
Variant "NC_045512:21974:Del_18:GATCCATTTTTGGGTGTTT:G"
Variant "NC_045512:22288:Del_6:TGCTTTA:T"
Variant "NC_045512:22455:Del_32:AAACAAAGTGTACGTTGAAATCCTTCACTGTAG:A"
Variant "NC_045512:23006:Del_2:GGT:G"
Variant "NC_045512:23552:Del_16:ATACCCATTGGTGCAGG:A"
Variant "NC_045512:24127:Del_13:TAACGGCCTTACTG:T"
Variant "NC_045512:25133:Del_2:AAG:A"
Variant "NC_045512:26071:Del_2:CAT:C"
Variant "NC_045512:26568:Del_30:CTCCTTGAACAATGGAACCTAGTAATAGGTT:C"
Variant "NC_045512:26573:Del_30:TGAACAATGGAACCTAGTAATAGGTTTCCTA:T"
Variant "NC_045512:26622:Del:CT:C"
Variant "NC_045512:26774:Del_56:GGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTT:G"
Variant "NC_045512:26774:Del_57:GGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTTT:G"
Variant "NC_045512:26783:Del_35:TGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTC:T"
Variant "NC_045512:27131:Del_7:CTATAAAT:C"
Variant "NC_045512:27194:Del_27:GACAACAGATGTTTCATCTCGTTGACTT:G"
Variant "NC_045512:27791:Del_3:CTTT:C"
Variant "NC_045512:27899:Del:AT:A"
Variant "NC_045512:29749:Del_10:ACGATCGAGTG:A"

# remainder

######### july 5, compare with the Substitution of the Texas paper

# example: RESULTS/NYC-Texas_SV_SNV_VCF/VCF/NYC_0.02AF_filtered/COVHA-20200403-P1-D03-P_bwamem.bam.lowfreq.vcf

cat RESULTS/NYC-Texas_SV_SNV_VCF/VCF/NYC_0.02AF_filtered/*_bwamem.bam.lowfreq.vcf | gawk '/^NC_/{printf("%s__%s__%s__%s\n",$1,$2,$4,$5);}'  > _a

cat _a | sort -u | gawk -F "__" '{c="NC_045512";x=$2;A=$3;G=$4;printf("Variant \"%s:%d:Sub:%s:%s\"\nIntMap %s %d %d \"Base %s2%s %d is modified\"\nVCF %d %s %s\nTreagen_TX\nTyp %s2%s\n%s2%s\n\n", c,x,A,G,c,x,x+1,A,G,x,x,A,G,A,G,A,G);}' > Treagen_TX.ace 

# Parent_sequence Found_in_genome

# Treagen with counts if > .95 : DP=250;AF=0.072000;SB=33;DP4=142,88,18,0

cat RESULTS/NYC-Texas_SV_SNV_VCF/VCF/NYC_0.02AF_filtered/*_bwamem.bam.lowfreq.vcf | gawk '/^##source/{run=$NF;gsub(/_bwamem.bam/,"",run);}/^NC_/{split($8,aa,";");split(aa[1],bb,"=");c=bb[2];split(aa[4],bb,"=");split(bb[2],cc,",");v=cc[3]+cc[4];r=cc[1]+cc[2];printf("%s__%s__%s__%s__%s__%d__%d__%d\n",$1,$2,$4,$5,run,v,r,c);}'  > _a
cat _a | sort -u | gawk -F "__" '{chr="NC_045512";x=$2;A=$3;G=$4;printf("Variant \"%s:%d:Sub:%s:%s\"\nIntMap %s %d %d \"Base %s2%s %d is modified\"\nVCF %d %s %s\nTreagen_TX\nTyp %s2%s\n%s2%s\n\n",chr,x,A,G,chr,x,x+1,A,G,x,x,A,G,A,G,A,G);}' > Treagen_TX.ace 
cat _a | sort -u | gawk -F "__" '{chr="NC_045512";x=$2;A=$3;G=$4;run=$5;gsub(/-/,"_",run);v=$6;r=$7;c=$8;if(v>=0 && 100*v >= -95*c)printf("Variant \"%s:%d:Sub:%s:%s\"\nTCounts %s %d %d %d  Frequency %.2f\n\n",chr,x,A,G,run,v,r,c,100*v/c);}' | grep -v COVWC  > Treagen_TX.TCounts.preace 
cat r2t.txt ZZZZZ Treagen_TX.TCounts.preace | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){t2r[$2]=$1;next;}}/TCounts/{r=t2r[$2];if(r)$2=r;}{print}' > Treagen_TX.TCounts.ace 


#### Magic and Magic Blast counts if AF  > 95%
cat tmp/TSNP_DB/zoneG.NC_045512_a/MMM.zoneG.NC_045512_a.snp_counts.tsf | sort -V | gawk -F '\t' '{v=$1;run=$2;m=($4+$5)/2;r=($6+$7)/2;c=($8+$9)/2;f=($12+$13)/2;if(c>0){if (v != oldV) {oldV=v;printf("\nVariant %s\n", v) ;}if(c<10)f=-10;printf("VCounts %s %d %d %d Frequency %.2f\n",run,m,r,c,f);}}END{printf("\n");}' > VCounts.ace

tace tmp/TSNP_DB/$zone <<EOF
  query find variant Vcounts
  edit -D Vcounts
  parse  VCounts.ace
  save
  quit
EOF
 
tace tmp/TSNP_DB/$zone <<EOF
  query find variant counts 
  show -a -f toto.ace Counts
  quit
EOF
 
set toto=RESULTS/snps.Magic_MagicBlast_Treagen.frequency.txt
echo -n "### $toto : " > $toto
date >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=5){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]>15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv)if(okv[v]==1){printf("\n%s",v);for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++)printf("\t%.2f",ff[v,ir,t]);}printf("\n");}' | gawk '/^#/{print}' >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=5){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv){okv[v]=0;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=3;t++){z=ff[v,ir,t];if(z>=5)okv[v]=1;}}for (v in vv)if(okv[v]==1){printf("\n%s",v); for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++){z=ff[v,ir,t];if(cc[v,ir,t]<-20)z=-20;printf("\t%.2f",z);}}printf("\n");}' | gawk '/^#/{next;}{print}' | sort -V >> $toto
 
set toto=RESULTS/snps.Magic_MagicBlast_Treagen.counts.txt
echo -n "### $toto : " > $toto
date >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=5){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv)if(okv[v]==1){printf("\n%s",v);for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++)printf("\t%.2f",ff[v,ir,t]);}printf("\n");}' | gawk '/^#/{print}' >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=5){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv){okv[v]=0;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=3;t++){z=ff[v,ir,t];if(z>=5)okv[v]=1;}}for (v in vv)if(okv[v]==1){ok2=1;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){for(t1=1;t1<=3;t1++)for(t2=t1+1;t2<=4;t2++){a=ff[v,ir,t1];b=ff[v,ir,t2];if(a>b+20 || b > a+20)ok2=1;}} if(ok2) {printf("\n%s",v);for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++){z=ff[v,ir,t];if(cc[v,ir,t]<-20)z=-20;printf("\tv:%d,r:%d,c:%d,f:%d",mm[v,ir,t],rr[v,ir,t],cc[v,ir,t],z);}}}printf("\n");}' | gawk '/^#/{next;}{print}' | sort -V >> $toto

set toto=RESULTS/snps.Magic_MagicBlast_Treagen.counts_with_contradictions.txt
echo -n "### $toto : " > $toto
date >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=95){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv){okv[v]=0;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++){z=ff[v,ir,t];if(z>=5)okv[v]=1;}}for (v in vv)if(okv[v]==1){printf("\n%s",v);for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++)printf("\t%.2f",ff[v,ir,t]);}printf("\n");}' | gawk '/^#/{print}' >> $toto
cat MetaDB/MMM/RunListSorted ZZZZZ toto.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){ir++;r2i[$1]=ir;i2r[ir]=$1;irMax=ir;next;}}{gsub(/\"/,"",$0);}/^Variant/{v=$2;vv[v]=1;}{t=0;tr=0;}/^MCounts/{t=1;tr=1;}/^MBCounts/{t=2;tr=2;}/^TCounts/{t=4;tr=4;}/^VCounts/{t=5;tr=8;}{if(t>0){run=$2;ir=r2i[run];mm[v,ir,t]=$3;rr[v,ir,t]=$4;cc[v,ir,t]=$5;ff[v,ir,t]=$7;if($7>=95){ok[ir]= or(ok[ir],tr);if(tr>0&&tr<8)okv[v]=1;}}}END{printf("# Run");for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){r=i2r[ir];printf("\tM:%s\tMB:%s\tT:%s\tV:%s",r,r,r,r);}for (v in vv){okv[v]=0;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=3;t++){z=ff[v,ir,t];if(z>=5)okv[v]=1;}}for (v in vv)if(okv[v]==1){ok2=0;for (ir=1;ir<=irMax;ir++)if(ok[ir]==15){for(t1=1;t1<=3;t1++)for(t2=t1+1;t2<=4;t2++){a=ff[v,ir,t1];b=ff[v,ir,t2];if(a>b+20 || b > a+20)ok2=1;}} if(ok2) {printf("\n%s",v);for (ir=1;ir<=irMax;ir++)if(ok[ir]==15)for(t=1;t<=5;t++){z=ff[v,ir,t];if(cc[v,ir,t]<-20)z=-20;printf("\tv:%d,r:%d,c:%d,f:%d",mm[v,ir,t],rr[v,ir,t],cc[v,ir,t],z);}}}printf("\n");}' | gawk '/^#/{next;}{print}' | sort -V >> $toto

