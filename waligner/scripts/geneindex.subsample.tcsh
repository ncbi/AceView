#!bin/tcsh -f

set phase=$1
set uu=$2
set kk=$3

if (x$uu != xu && x$uu != xnu) then
  echo "missing second  parameter  u or nu"
  goto done
endif
  
set kks="10k 20k 50k 100k 200k 500k 1M 2M 5M 10M 20M 50M 100M 0k"

if ($phase == run) then
  if  ($uu == u) set UU=unique
  if  ($uu == nu) set UU=quasi_unique

    bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.av.GENE.$uu.ace -$uu -mask tmp/GENEINDEX/$MAGIC.av.$uu.mask  -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.av.GENE.info.ace  -o  tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.$uu -gzo  -pA -method Gene_AceView     -stableGenes TARGET/Targets/hs.av.stable_genes.txt -referenceGenome GRCh37.p10__NCBI_37_5__ANNOTATION_RELEASE.104__2013_02_01  -target_class ET_av -geneGroup TARGET/GENES/Gene_groups.ace  -exportDiffGenes  -compare -correlation    -htmlSpecies hs   -export abitvz -subsample $kk
    touch RESULTS/Expression.$kk/$UU/av/run.done

  goto done

endif

#############################################################################

if ($phase == compute) then
  if  ($uu == u) set UU=unique
  if  ($uu == nu) set UU=quasi_unique

  foreach kk ($kks)
      if (-e RESULTS/Expression.$kk/$UU/av/run.done) continue
      if (! -d RESULTS/Expression.$kk) mkdir  RESULTS/Expression.$kk
      if (! -d RESULTS/Expression.$kk/$UU) mkdir  RESULTS/Expression.$kk/$UU
      if (-d RESULTS/Expression.$kk/$UU/av) \rm -rf  RESULTS/Expression.$kk/$UU/av
      mkdir  RESULTS/Expression.$kk/$UU/av

    scripts/submit tmp/GENEINDEX/Results/EXpresson.$kk/$MAGIC.$myUU  "scripts/geneindex.subsample.tcsh run $uu $kk"
  end

  goto done
endif

#############################################################################

if ($phase == transfer) then
  foreach uu (u nu)
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique
    foreach kk ($kks)
      if (! -e RESULTS/Expression.$kk/$UU/av/run.done) continue
      if (-e RESULTS/Expression.$kk/$UU/av/transfer.done) continue
      \rm -rf  RESULTS/Expression.$kk/$UU/av
      mkdir  RESULTS/Expression.$kk/$UU/av
      \cp tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.$uu.*.profiles.txt RESULTS/Expression.$kk/$UU/av
      touch RESULTS/Expression.$kk/$UU/av/run.done
      touch RESULTS/Expression.$kk/$UU/av/transfer.done
    end
  end
  goto done
endif

#############################################################################

if ($phase == report) then
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique
  echo "report uu=$uu UU=$UU"
  foreach kk ($kks)
    if (! -e RESULTS/Expression.$kk/$UU/av/transfer.done) continue
    if (-e RESULTS/Expression.$kk/$UU/av/report.done) continue
    scripts/geneindex.subsample.tcsh reportOne $uu $kk 
  end
  goto done
endif

#############################################################################

if ($phase == reportOne) then
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique

  if (-e RESULTS/Expression.$kk/$UU/av/report.done) goto done
  set Expression=Expression.$kk

  set ln=geneBox_length
  set ln=mRNA_length
  set toto=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENE.$ln.DEG.$uu.stats.txt
  echo -n "### $toto : subsampling $kk " > $toto
  date >> $toto

  set cap=A1A2I2I3R1R2
  set totocap=RESULTS/$Expression/$UU/av/AECDB_diff.GENE.$cap.$ln.DEG.$uu.profile.stats.txt
  echo -n "### $toto$cap : subsampling $kk $cap " > $totocap
  date >> $totocap

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB AGLR1_AECDB AGLR2_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB BSPR1_AECDB )
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | head -1 >> $toto
      cat $ff | tail -6 >> $toto
    end

    if (-e RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$cap.$uu.AGLR1_AECDB_Profile.score.genes.profiles.txt) then 
      foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$cap.$uu.$comp'_Profile'.score.genes.profiles.txt`)
        cat $ff | head -1 >> $toto$cap
        cat $ff | tail -6 >> $toto$cap
      end
    endif
  end

  set ln=mRNA_length
  set toto=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENE.$ln.DEG.$uu.heatmap.txt
  echo -n "### $toto :" > $toto
  date >> $toto



# grep index min/max et fold change
  cat RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$uu.RNA_Total_ACB_Profile.score.genes.profiles.txt | gawk -F '\t' '/^#/{next;}{printf("LnMiMxFc\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$41,$42,$40);}' > $toto.0

  date > $toto.1
  date > $toto.2

  foreach comp (AGLR1_AECDB AGLR2_AECDB BSPR1_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB)
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$45,$46);}' ff=$ff | grep AECDB_Profile | sed -e "s/RESULTS\/$Expression\/$UU\/av\/AECDB_diff.$kk.AceView.GENE.$uu.//" -e 's/_AECDB_Profile.score.genes.profiles.txt//' >> $toto.1
    end
  end

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB)
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENE.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$29,$30);}' ff=$ff | grep ACB_Profile | sed -e "s/RESULTS\/$Expression\/$UU\/av\/AECDB_diff.$kk.AceView.GENE.$uu.//" -e 's/_ACB_Profile.score.genes.profiles.txt//' >> $toto.2
    end
  end

  cat $toto.[012] |  gawk -F '\t' -f scripts/deg_capture_heatmap.awk outf=RESULTS/$Expression/$UU/av/AECDB_diff    > $toto.3
  cat $toto.3 | head -8 | gawk '/^#/{print}' >> $toto
  cat $toto.3 | gawk '/^#/{next;}{print}' | sort -k 7,7 -k 1,1  >> $toto
  #\rm $toto.[0123]
  ls -ls $toto.[0123]

  set toto2=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENE.$ln.DEG.$uu.heatmap.A2R2I3.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A2")*index(c,"R2")*index(c,"I3")>0)print;}' >> $toto2


  set toto2=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENE.$ln.DEG.$uu.heatmap.A1A2I2I3R1R2.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A1")*index(c,"A2")*index(c,"I2")*index(c,"I3")*index(c,"R1")*index(c,"R2")>0)print;}' >> $toto2

ls -ls $toto
wc $toto
goto done

  echo $toto.done
  touch RESULTS/Expression.$kk/$UU/av/report.done
  goto done
endif

#############################################################################

goto done

date > toto
foreach kk ($kks)
  cat RESULTS/Expression.$kk/$UU/av/AECDB_diff.$kk.GENE.mRNA_length.DEG.$uu.profile_stats.txt | gawk -F '\t' '/^## Stat/{r=$1;n=0;next;}{n++;}{if(kk=="0k")kk=1000000000;gsub("M","000k",kk);gsub("k","",kk);if(n==2)printf ("%s\t%010d\t%d\n", r,kk, $2)}' kk=$kk >> toto
end
cat toto | sort -k 1,1 -k 2nr

#############################################################################

done:
    echo -n "Done "
    date

#############################################################################
#############################################################################
