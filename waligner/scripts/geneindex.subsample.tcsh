#!bin/tcsh -f

set phase=$1
set kk=$2

set kks="10k 20k 50k 100k 200k 500k 1M 2M 5M 10M 20M 50M 100M 0k"

if ($phase == run) then

  if (0) then
    bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.av.GENE.u.ace -u -mask tmp/GENEINDEX/$MAGIC.av.u.mask  -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.av.GENE.info.ace  -o  tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.u -gzo  -pA -method Gene_AceView     -stableGenes TARGET/Targets/hs.av.stable_genes.txt -referenceGenome GRCh37.p10__NCBI_37_5__ANNOTATION_RELEASE.104__2013_02_01  -target_class ET_av -geneGroup TARGET/GENES/Gene_groups.ace    -htmlSpecies hs   -export abitvz -subsample $kk
  else
    bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.av.GENE.nu.ace -nu -mask tmp/GENEINDEX/$MAGIC.av.nu.mask  -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.av.GENE.info.ace  -o  tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.nu -gzo  -pA -method Gene_AceView     -stableGenes TARGET/Targets/hs.av.stable_genes.txt -referenceGenome GRCh37.p10__NCBI_37_5__ANNOTATION_RELEASE.104__2013_02_01  -target_class ET_av -geneGroup TARGET/GENES/Gene_groups.ace -exportDiffGenes  -compare -correlation    -htmlSpecies hs   -export abitvz -subsample $kk
  endif
  touch RESULTS/Expression.$kk/quasi_unique/av/run.done
  goto done
endif

if (0) then
  bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.av.GENE.u.ace -u -mask tmp/GENEINDEX/$MAGIC.av.nu.mask  -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.av.GENE.info.ace  -o  tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.u -gzo  -pA -method Gene_AceView     -stableGenes TARGET/Targets/hs.av.stable_genes.txt -referenceGenome GRCh37.p10__NCBI_37_5__ANNOTATION_RELEASE.104__2013_02_01  -target_class ET_av -geneGroup TARGET/GENES/Gene_groups.ace    -htmlSpecies hs   -export abitvz -subsample $kk
endif

#############################################################################

if ($phase == compute) then
  foreach kk ($kks)
    if (-e RESULTS/Expression.$kk/quasi_unique/av/run.done) continue
    if (! -d RESULTS/Expression.$kk) mkdir  RESULTS/Expression.$kk
    if (! -d RESULTS/Expression.$kk/quasi_unique) mkdir  RESULTS/Expression.$kk/quasi_unique
    if (-d RESULTS/Expression.$kk/quasi_unique/av) \rm -rf  RESULTS/Expression.$kk/quasi_unique/av
    mkdir  RESULTS/Expression.$kk/quasi_unique/av

    scripts/submit tmp/GENEINDEX/Results/$MAGIC.$kk "scripts/geneindex.subsample.tcsh run $kk"
  end
  goto done
endif

#############################################################################

if ($phase == transfer) then
  foreach kk ($kks)
    if (! -e RESULTS/Expression.$kk/quasi_unique/av/run.done) continue
    if (-e RESULTS/Expression.$kk/quasi_unique/av/transfer.done) continue
    \rm -rf  RESULTS/Expression.$kk/quasi_unique/av
    mkdir  RESULTS/Expression.$kk/quasi_unique/av
    \cp tmp/GENEINDEX/Results/$MAGIC.$kk.AceView.GENE.nu.*.profiles.txt RESULTS/Expression.$kk/quasi_unique/av
    touch RESULTS/Expression.$kk/quasi_unique/av/run.done
    touch RESULTS/Expression.$kk/quasi_unique/av/transfer.done
  end
  goto done
endif

#############################################################################

if ($phase == report) then
  foreach kk ($kks)
    if (! -e RESULTS/Expression.$kk/quasi_unique/av/transfer.done) continue
    if (-e RESULTS/Expression.$kk/quasi_unique/av/report.done) continue
    scripts/geneindex.subsample.tcsh reportOne $kk 
  end
  goto done
endif

#############################################################################

if ($phase == reportOne) then

  if (-e RESULTS/Expression.$kk/quasi_unique/av/report.done) goto done
  set Expression=Expression.$kk

  set ln=geneBox_length
  set ln=mRNA_length
  set toto=RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.GENE.$ln.DEG_nu_profile_stats.txt
  echo -n "### $toto : subsampling $kk " > $toto
  date >> $toto

  set cap=A1A2I2I3R1R2
  set totocap=RESULTS/$Expression/quasi_unique/av/AECDB_diff.GENE.$cap.$ln.DEG_nu_profile_stats.txt
  echo -n "### $toto$cap : subsampling $kk $cap " > $totocap
  date >> $totocap

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB AGLR1_AECDB AGLR2_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB BSPR1_AECDB )
    foreach ff (`ls  RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.nu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | head -1 >> $toto
      cat $ff | tail -6 >> $toto
    end

    if (-e RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.$cap.nu.AGLR1_AECDB_Profile.score.genes.profiles.txt) then 
      foreach ff (`ls  RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.$cap.nu.$comp'_Profile'.score.genes.profiles.txt`)
        cat $ff | head -1 >> $toto$cap
        cat $ff | tail -6 >> $toto$cap
      end
    endif
  end

  if (0) then
    set toto=RESULTS/$Expression/unique/av/AECDB_diff.$kk.GENE.$ln.DEG_u_profile_stats.txt
    echo -n "### $toto :" > $toto
    date >> $toto

    foreach ff (`ls  RESULTS/$Expression/unique/av/AECDB_diff.$kk.AceView.GENE.u.*.score.genes.profiles.txt`)
      cat $ff | head -1 >> $toto
      cat $ff | tail -6 >> $toto
    end
  endif

  set ln=mRNA_length
  set toto=RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.GENE.$ln.DEG_nu_heatmap.txt
  echo -n "### $toto :" > $toto
  date >> $toto



# gep index min/max et fold change
  cat RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.nu.RNA_Total_ACB_Profile.score.genes.profiles.txt | gawk -F '\t' '/^#/{next;}{printf("MiMxFc\t%s\t%s\t%s\t%s\n",$2,$41,$42,$40);}' > $toto.0

  date > $toto.1
  date > $toto.2

  foreach comp (AGLR1_AECDB AGLR2_AECDB BSPR1_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB)
    foreach ff (`ls  RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.nu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$45,$46);}' ff=$ff | grep AECDB_Profile | sed -e "s/RESULTS\/$Expression\/quasi_unique\/av\/AECDB_diff.$kk.AceView.GENE.nu.//" -e 's/_AECDB_Profile.score.genes.profiles.txt//' >> $toto.1
    end
  end

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB)
    foreach ff (`ls  RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.AceView.GENE.nu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$29,$30);}' ff=$ff | grep ACB_Profile | sed -e "s/RESULTS\/$Expression\/quasi_unique\/av\/AECDB_diff.$kk.AceView.GENE.nu.//" -e 's/_ACB_Profile.score.genes.profiles.txt//' >> $toto.2
    end
  end

  cat $toto.[012] |  gawk -F '\t' -f scripts/deg_capture_heatmap.awk outf=RESULTS/$Expression/quasi_unique/av/AECDB_diff   >> $toto
  \rm $toto.[012]

  set toto2=RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.GENE.$ln.DEG_nu_heatmap.A2R2I3.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A2")*index(c,"R2")*index(c,"I3")>0)print;}' >> $toto2


  set toto2=RESULTS/$Expression/quasi_unique/av/AECDB_diff.$kk.GENE.$ln.DEG_nu_heatmap.A1A2I2I3R1R2.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A1")*index(c,"A2")*index(c,"I2")*index(c,"I3")*index(c,"R1")*index(c,"R2")>0)print;}' >> $toto2

  echo $toto
  touch RESULTS/Expression.$kk/quasi_unique/av/report.done
  goto done
endif

#############################################################################

date > toto
foreach kk ($kks)
  cat RESULTS/Expression.$kk/quasi_unique/av/AECDB_diff.$kk.GENE.mRNA_length.DEG_nu_profile_stats.txt | gawk -F '\t' '/^## Stat/{r=$1;n=0;next;}{n++;}{if(kk=="0k")kk=1000000000;gsub("M","000k",kk);gsub("k","",kk);if(n==2)printf ("%s\t%010d\t%d\n", r,kk, $2)}' kk=$kk >> toto
end
cat toto | sort -k 1,1 -k 2nr

#############################################################################

done:
    echo -n "Done "
    date

#############################################################################
#############################################################################