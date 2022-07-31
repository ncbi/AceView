#!bin/tcsh -f

set phase=$1
set type=GENESPX
if ($3 != "") set type=$3 

set toto=RESULTS/Expression/unique/av/$MAGIC.ordinal_norm

if ($phase == count) goto phaseCount
if ($phase == stats) goto phaseStats
goto phaseLoop


phaseCount:
echo "ordinal_norm.tcsh $phase $type"

cat MetaDB/$MAGIC/runs.ace | gawk 'BEGIN{c="NULL";}{gsub(/\"/,"",$0);}/^Run /{ok=0;run=$2;next;}/^Capture/{printf("RUN\t%s\t%s\n",run,$2);}' > $toto.rc
cat tmp/METADATA/$MAGIC.av.captured_genes.ace | gawk 'gsub(/\"/,"",$0);/^Gene/{g=$2;}/^Capture_touch/{next;}/^Capture/{printf("GENE\t%s\t%s\n",g,$2);}' > $toto.gc

set ff=RESULTS/Expression/AceFiles/$MAGIC.AceView.$type.A1.u.ace.gz
ls -ls $ff
if (-e $ff) then
  zcat $ff | gawk '/^$/{g=0;next;}{gsub(/\"/,"",$0);}/^Gene /{g=$2;next;}{r=0;}/^Run_U/{r=$2;}/^Group_U/{r=$2;}{if(g!=0 && r!=0)printf("%s\t%.2f\t%s\n",r,$3,g);next;}' > $toto.1
endif
cat $toto.1 |  grep -v _SumOfAllReads | sort -k 1,1 -k 2,2nr > $toto.2
\rm $toto.1
phaseStats:

cat $toto.rc $toto.gc $toto.2 | gawk -F '\t' '/^RUN/{r2c[$2]=$3;next;}/^GENE/{gc[$2,$3]=1;next;}{r=$1;g=$3;c=r2c[r];cc=gc[g,c]+0;nn[$1,cc]++;n=nn[$1,cc];zn=","n",";if(index(",100,200,300,500,700,1000,2000,3000,5000,7000,10000,",zn)>0) printf("%s\t%d\t%d\t%s\n",$1,cc,n,$2);}' > $toto.3
cat $toto.3 | sort -k 1,1 -k 2,2n -k 3,3n | gawk -F '\t' '{nr[$1]++;if(nr[$1]==1){ir++;i2r[ir]=$1;}nc[$3]++;if(nc[$3]==1){jc++;j2c[jc]=$3;}z[$1,$2,$3]=$4;}END{printf("#Run");for(j=1;j<=jc;j++)printf("\t%s",j2c[j]);printf("\t\t");for(j=1;j<=jc;j++)printf("\t%s",j2c[j]);for(i=1;i<=ir;i++){printf("\n%s",i2r[i]);for(j=1;j<=jc;j++)printf("\t%.2f",z[i2r[i],0,j2c[j]]);printf("\t\t");for(j=1;j<=jc;j++)printf("\t%.2f",z[i2r[i],1,j2c[j]]);}printf("\n");}' > $toto.txt

ls -ls $toto.*
cat $toto.txt

goto phaseLoop





















phaseLoop:
   echo done
  
