
## analysis of SEQC_main, 2012_11_29

cd IntronHisto
\rm _toto _runs
set ali2=0

foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep A_SEQC |  grep SEQC_ROC | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep A_SEQC |  grep Illumina | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep A_SEQC |  grep SOLiD | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep B_SEQC | grep SEQC_ROC | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep B_SEQC |  grep Illumina | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep B_SEQC |  grep SOLiD | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end

foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep C_SEQC |  grep SEQC_ROC | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep C_SEQC |  grep Illumina | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep C_SEQC |  grep SOLiD | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep D_SEQC | grep SEQC_ROC | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep D_SEQC |  grep Illumina | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end
foreach run (`cat MetaDB/$MAGIC/gtitle.txt | grep D_SEQC |  grep SOLiD | cut -f 1 | sed -e 's/\"//g'`)
  echo $run >> _runs
end

\rm _toto3
set oldrun="xxx"

foreach run (`cat _runs | grep Rhs | head -30000`)
  if ($run == $oldrun) continue
  set oldrun=$run
  if (! -e tmp/OR/$run/d1.$run.de_uno.txt.gz) continue
  gunzip -c tmp/OR/$run/d1.$run.de_uno.txt.gz >> _toto
  set ali=`cat  MetaDB/$MAGIC/runAligned.txt |  gawk -F '\t' '{gsub(/\"/,"",$0);if(zz<1&&$1== run){zz++;printf("%d", $2);last;}}' run=$run`
  set kb=`cat  MetaDB/$MAGIC/runAligned.txt |  gawk -F '\t' '{gsub(/\"/,"",$0);if(zz<1&&$1== run){zz++;printf("%d", $3);last;}}' run=$run`
  set sample=`cat  MetaDB/$MAGIC/gtitle.txt |  gawk -F '\t' '{gsub(/\"/,"",$0);if(zz<1&&$1== run){zz++;printf("%s", $6);last;}}' run=$run`
  echo "RUN\t$run\t$sample\t$ali\t$kb" >> _toto3
end

cat _toto3 | grep _A_ > _toto3A
cat _toto3 | grep _B_ > _toto3B
cat _toto3 | grep _C_ > _toto3C
cat _toto3 | grep _D_ > _toto3D

cat _toto3A | grep BGI > _toto3A.BGI
cat _toto3C | grep BGI > _toto3C.BGI
cat _toto3D | grep BGI > _toto3D.BGI
cat _toto3B | grep BGI > _toto3B.BGI

\rm _toto5
foreach i (1 2 3 4 5 6 7 8)
  foreach ss ( _s _t _u)
    foreach pl (AGR BGI CNL COH MAY NVS LIV NWU PSU SQW )
      foreach ii (_1_ _2_ _3_ _4_ _5_)
        cat _toto3 | grep I_$pl | grep $ii | grep $ss$i >> _toto5
      end 
    end
  end
end
 
 
foreach t (A B C D)
  cat TARGET/MRNAS/introns_RefSeq.txt ZZZZZ  TARGET/MRNAS/introns_av.txt ZZZZZ TARGET/MRNAS/introns_ensembl.txt ZZZZZ TARGET/MRNAS/introns_WEHI.txt ZZZZZ  _toto3$t | gawk -F '\t' -f scripts/d5.intronHisto.awk out=intronList.$t tissue=$t > $toto.$t &
end

foreach t (A B C D)
  cat TARGET/MRNAS/introns_RefSeq.txt ZZZZZ  TARGET/MRNAS/introns_av.txt ZZZZZ TARGET/MRNAS/introns_ensembl.txt ZZZZZ TARGET/MRNAS/introns_WEHI.txt ZZZZZ  _toto3$t.BGI | gawk -F '\t' -f scripts/d5.intronHisto.awk out=intronList.$t.BGI tissue=$t > $toto.$t.BGI &
end

# titration=1 to have d5.intronHisto.awk export the titration columns
set toto="intronHistoMainNB.refseq_encode_magic"
set toto="intronHistoMainNB.refseq_encode_magic_just_gt_ag_sorted"
set toto="intronHistoMainNB.refseq_encode_magic_sorted_new"
#set toto="intronHistoMainNB.non_gt_ag"
#set toto="intronHistoTitrationACDB.fuzzy1000"
#set toto="intronHistoT"
#set toto="intronHistoT.countOnce"
echo -n "# " > $toto.txt
date >> $toto.txt
echo "# variation in the number of introns with depth of the experiments" >> $toto.txt
echo "# We count independent discoveries by reads aligning discontinuously on the gene,unknowingly of previous annotation" >> $toto.txt
cat TARGET/MRNAS/introns_RefSeq.txt ZZZZZ  TARGET/MRNAS/introns_av.txt ZZZZZ TARGET/MRNAS/introns_encode.txt ZZZZZ TARGET/MRNAS/NB_introns_100_supports.txt  ZZZZZ  _toto5NB_sorted | gawk -F '\t' -f scripts/d5.intronHisto.awk out=$toto.list tissue=Any titration=0 donor=0 gtag=0 fuzzy=0 countOnce=0 | tee -a  $toto.txt 

### Venn diagram
# compute the Venn diagram at different thresholds
date >  intronsMagic.BGI.Venn.txt
echo "Magic BGI introns using Magic counts" >>  intronsMagic.BGI.Venn.txt
cat intronList.?.BGI | gawk -F '\t' -f scripts/d5.intronVenn.awk limit=10 >> intronsMagic.BGI.Venn.txt
cat intronList.?.BGI | gawk -F '\t' -f scripts/d5.intronVenn.awk limit=5 >> intronsMagic.BGI.Venn.txt
cat intronList.?.BGI | gawk -F '\t' -f scripts/d5.intronVenn.awk limit=2 >> intronsMagic.BGI.Venn.txt
cat intronList.?.BGI | gawk -F '\t' -f scripts/d5.intronVenn.awk limit=1 >> intronsMagic.BGI.Venn.txt


##########################################
##########################################
## R fits-> python fit, see poisson.py

## fit the numer of reads supporting an intron per kb aligned

set toto="intronHisto.txt"
# nb of reads supporting an intron in de-uno per kb aligned
cat $toto | gawk -F '\t' '{if(0+$4>0)printf("%d\t%d\t%.3f\n",$6,$7,100*$7/$4);}' > r2is.txt 
# number of av introns supported at least once
cat $toto | gawk -F '\t' '{if(0+$4>0)printf("%d\t%d\n",$6,$71);}' > av1.txt 

R
rr=read.table("r2is.txt")
# limit to 4Gb sequence pour eviter le trou (perte du disque NB lors du recalcul)
rr2=rr[rr[,1]<4000000000,1:2])
plot(rr2[,1],rr2[,2])
lines(rr2[,1],1.235*rr2[,1])

## A B samples only rate of reads with introns 14,9% in A, 11.5% in B
 rr2=rr[rr[,1]<1320000000,1:2]
plot(rr2[,1],rr2[,2])
lines(rr2[,1],1.490*rr2[,1])
lines(rr2[,1],1.15*rr2[,1]+180000000)

# very nice display showing that B has longer UTRs
pdf("deUnoPerTbAlignedinA_B.pdf")
> plot(rr2[,1],rr2[,2])
> lines(rr2[,1],1.490*rr2[,1])
> lines(rr2[,1],1.15*rr2[,1]+180000000)
> dev.off()

# number of av introns supported at least once
rr=read.table("av1.txt")
plot(rr)
# limit to library A
rr2=rr[rr[,1]<650000000,1:2]
plot(rr2[,1],rr2[,2])
lines(rr2[,1],295000.0+0.0000*rr2[,1] -rr2[,2])
lines(rr2[,1],9000*600000000/(1+rr2[,1]))

plot(rr2[,1],rr2[,2]
lines(rr2[,1],295000.0+0.0000*rr2[,1] -rr2[,2])
 lines(rr2[,1],12000*600000000/(1+rr2[,1]))

rr=read.table("av1.txt")
rr2=rr[rr[,1]>1800000000 & rr[,1]<4000000000,1:2]
plot(rr2)

