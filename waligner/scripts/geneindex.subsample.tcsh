#!bin/tcsh -f

set phase=$1
set uu=$2
set kk=$3

if (x$uu != xu && x$uu != xnu) then
  echo "missing second  parameter  u or nu"
  goto done
endif
  
set kks="10k 20k 50k 100k 200k 500k 1M 2M 5M 10M 20M 50M 100M 0k"
# set kks="10k"

if ($phase == run) then
  if  ($uu == u) set UU=unique
  if  ($uu == nu) set UU=quasi_unique

    bin/geneindex -deepGene tmp/GENEINDEX/$MAGIC.av.GENESP.$uu.ace -$uu -mask tmp/GENEINDEX/$MAGIC.av.$uu.mask  -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.av.GENESP.info.ace  -o RESULTS/Expression.$kk/$UU/av/$MAGIC.$kk.AceView.GENESP.$uu -gzo  -pA -method Gene_AceView     -stableGenes TARGET/Targets/hs.av.stable_genes.txt -referenceGenome GRCh37.p10__NCBI_37_5__ANNOTATION_RELEASE.104__2013_02_01  -target_class ET_av -geneGroup TARGET/GENES/Gene_groups.ace  -exportDiffGenes  -compare -correlation    -htmlSpecies hs   -export abitvz -subsample $kk
    touch RESULTS/Expression.$kk/$UU/av/runSP.done

  goto done

endif

#############################################################################

if ($phase == compute) then
  if  ($uu == u) set UU=unique
  if  ($uu == nu) set UU=quasi_unique

  foreach kk ($kks)
      if (-e RESULTS/Expression.$kk/$UU/av/runSP.done) continue
      if (! -d RESULTS/Expression.$kk) mkdir  RESULTS/Expression.$kk
      if (! -d RESULTS/Expression.$kk/$UU) mkdir  RESULTS/Expression.$kk/$UU
      if (-d RESULTS/Expression.$kk/$UU/av) \rm -rf  RESULTS/Expression.$kk/$UU/av
      mkdir  RESULTS/Expression.$kk/$UU/av

    scripts/submit RESULTS/Expression.$kk/$UU/av/$MAGIC.compute.$myUU  "scripts/geneindex.subsample.tcsh run $uu $kk"
  end

  goto done
endif

#############################################################################

if ($phase == transfer) then
  goto done
  foreach uu (u nu)
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique
    foreach kk ($kks)
      if (! -e RESULTS/Expression.$kk/$UU/av/runSP.done) continue
      if (-e RESULTS/Expression.$kk/$UU/av/transferSP.done) continue
      # \rm -rf  RESULTS/Expression.$kk/$UU/av
      # mv tmp/GENEINDEX/Results/Expression.$kk RESULTS/Expression.$kk/$UU/av
      # touch RESULTS/Expression.$kk/$UU/av/runSP.done
      touch RESULTS/Expression.$kk/$UU/av/transferSP.done
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
    if (! -e RESULTS/Expression.$kk/$UU/av/runSP.done) continue
    if (-e RESULTS/Expression.$kk/$UU/av/reportSP.done) continue
    scripts/geneindex.subsample.tcsh reportOne $uu $kk 
  end
  goto done
endif

#############################################################################

if ($phase == reportOne) then
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique

  set totog=RESULTS/Expression.0k/unique/av/gene_truth.list.0k
  touch $totog

  if (-e RESULTS/Expression.$kk/$UU/av/report.done2) goto done
  set Expression=Expression.$kk

  set ln=geneBox_length
  set ln=mRNA_length
  set toto=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENESP.$ln.DEG.$uu.stats.txt
  echo -n "### $toto : subsampling $kk " > $toto
  date >> $toto

  set cap=A1A2I2I3R1R2
  set totocap=RESULTS/$Expression/$UU/av/AECDB_diff.GENESP.$cap.$ln.DEG.$uu.profile.stats.txt
  echo -n "### $toto$cap : subsampling $kk $cap " > $totocap
  date >> $totocap

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB AGLR1_AECDB AGLR2_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB BSPR1_AECDB )
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENESP.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | head -1 >> $toto
      cat $ff | tail -6 >> $toto
    end

    if (-e RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENESP.$cap.$uu.AGLR1_AECDB_Profile.score.genes.profiles.txt) then 
      foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENESP.$cap.$uu.$comp'_Profile'.score.genes.profiles.txt`)
        cat $ff | head -1 >> $toto$cap
        cat $ff | tail -6 >> $toto$cap
      end
    endif
  end

  set ln=mRNA_length
  set toto=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENESP.$ln.DEG.$uu.heatmap.txt
  echo -n "### $toto :" > $toto
  date >> $toto

  set totoCL=RESULTS/Expression/$UU/av/AECDB_diff.CL.GENESP.DEG.$uu.heatmap.txt
  echo -n "### $totoCL :" > $totoCL
  date >> $totoCL



# grep index min/max et fold change
  cat RESULTS/Expression.0k/$UU/av/AECDB_diff.0k.AceView.GENESP.$uu.RNA_Total_ACB_Profile.score.genes.profiles.txt | gawk -F '\t' '/^#/{next;}{printf("LnMiMxFc\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$41,$42,$40);}' > $toto.0

  date > $toto.1
  date > $toto.2
  date > $toto.3

  foreach comp (AGLR1_AECDB AGLR2_AECDB BSPR1_AECDB ROCR1_AECDB ROCR2_AECDB ILMR1_AECDB ILMR2_AECDB ILMR2_lowQ_AECDB ILMR3_AECDB)
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENESP.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '/^#/{next;}{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$45,$46);}' ff=$ff | grep AECDB_Profile | sed -e "s/RESULTS\/$Expression\/$UU\/av\/AECDB_diff.$kk.AceView.GENESP.$uu.//" -e 's/_AECDB_Profile.score.genes.profiles.txt//' >> $toto.1
    end
  end

  foreach comp (RNA_Total_ACB RNA_PolyA_ACB Nanopore.titr_AGLR2_ACB PacBio2.titr.ccs3_AGLR2_ACB Nanopore.titr_ROCR3_ACB PacBio2.titr.ccs3_ROCR3_ACB)
    foreach ff (`ls  RESULTS/$Expression/$UU/av/AECDB_diff.$kk.AceView.GENESP.$uu.$comp'_Profile'.score.genes.profiles.txt`)
      cat $ff | gawk -F '\t' '/^#/{next;}{printf("%s\t%s\t%s\t%s\t%s\t%s\n",ff,$2,$3,$5,$29,$30);}' ff=$ff | grep ACB_Profile | sed -e "s/RESULTS\/$Expression\/$UU\/av\/AECDB_diff.$kk.AceView.GENESP.$uu.//" -e 's/_ACB_Profile.score.genes.profiles.txt//' >> $toto.2
    end
  end

  if ($kk == 0k) then
    \rm $toto.3
    foreach cl (CL1-Brain-B_priv-2sA1 CL10-Testis-B_priv-2sA1 CL2-Breast-B_priv-2sA1 CL3-Cervix-B_priv-2sA1 CL4-Liver-B_priv-2sA1 CL5-Lipo-B_priv-2sA1 CL6-Blym-B_priv-2sA1 CL7-Tlym-B_priv-2sA1 CL8-Macr-B_priv-2sA1 CL9-Skin-B_priv-2sA1 SumOfCL1toCL10-B_2grps-2_A1  A-UHR-B_priv_4.2sA1 A-UHR-B_4lR3 CL1-Brain-B_4lR3 CL2-Breast-B_4lR3 CL3-Cervix-B_4lR3 CL4-Liver-B_4lR3 CL5-Lipo-B_4lR3 CL6-Blym-B_4lR3 CL7-Tlym-B_4lR3 CL8-Macr-B_4lR3 CL9-Skin-B_4lR3 CL10-Testis-B_4lR3 )
      set ff=RESULTS/Expression/unique/av/CL.AceView.GENESP.u.$cl.score.genes.profiles.txt
      set n1=`cat $ff | head -12 | transpose | grep -n 'Sum of the differential scores of the even columns' | gawk -F : '{print $1}'`
      set n2=`cat $ff | head -12 | transpose | grep -n 'Sum of the differential scores of the odd columns' | gawk -F : '{print $1}'`
      cat $ff | gawk -F '\t' '/^#/{next;}{printf("%s\t%s\t%s\t%s\t%s\t%s\n",cl,$2,$3,$5,$n1,$n2);}' cl=$cl n1=$n1 n2=$n2 >> $toto.3 
    end
    cat $totog $toto.[0123] |  gawk -F '\t' -f scripts/deg_capture_heatmap.awk CL=1 outf=RESULTS/$Expression/$UU/av/AECDB_diff    > $toto.4
    cat $toto.4 | head -8 | gawk '/^#/{print}' >> $totoCL
    cat $toto.4 | gawk '/^#/{next;}{print}' | sort -k 7,7r -k 1,1  | grep -v non-DEG | grep -v Weak >> $totoCL
    cat $toto.4 | gawk '/^#/{next;}{print}' | sort -k 7,7r -k 1,1  | grep Weak >> $totoCL
    cat $toto.4 | gawk '/^#/{next;}{print}' | sort -k 7,7r -k 1,1  | grep non-DEG  >> $totoCL
  endif

  cat $totog $toto.[012] |  gawk -F '\t' -f scripts/deg_capture_heatmap.awk outf=RESULTS/$Expression/$UU/av/AECDB_diff    > $toto.5
  cat $toto.5 | head -8 | gawk '/^#/{print}' >> $toto
  cat $toto.5 | gawk '/^#/{next;}{print}' | sort -k 7,7r -k 1,1  >> $toto

  # \rm $toto.[012345]
  ls -ls $toto.[012345]

  set toto2=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENESP.$ln.DEG.$uu.heatmap.A2R2I3.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A2")*index(c,"R2")*index(c,"I3")>0)print;}' >> $toto2


  set toto2=RESULTS/$Expression/$UU/av/AECDB_diff.$kk.GENESP.$ln.DEG.$uu.heatmap.A1A2I2I3R1R2.txt
  echo -n "### $toto2 :" > $toto2
  date >> $toto2

  cat $toto | gawk '/^###/{next;}/^#/{print;next;}{c=$3;if(index(c,"A1")*index(c,"A2")*index(c,"I2")*index(c,"I3")*index(c,"R1")*index(c,"R2")>0)print;}' >> $toto2

ls -ls $toto
wc $toto
goto done
 

  echo $toto.done
  touch RESULTS/Expression.$kk/$UU/av/reportSP.done
  goto done
endif

#############################################################################

if ($phase == cumul) then
    if  ($uu == u) set UU=unique
    if  ($uu == nu) set UU=quasi_unique

  # establish the truth list
  set totog=RESULTS/Expression.0k/unique/av/gene_truth.list
    if (! -e $totog.ZZ) then
      cat RESULTS/Expression.0k/unique/av/AECDB_diff.0k.GENESP.mRNA_length.DEG.u.heatmap.txt | gawk -F '\t' '{printf("Truth\t%s\t%s\n",$1,$7);}' > $totog.0k
    endif

    set toto=RESULTS/Expression.0k/deg_truth_per_depth.txt
    echo $toto > toto.1
    echo -n "### $toto : " > $toto
    date >> $toto

    set toto2=RESULTS/Expression.0k/deg_truth_per_depth_per_platform.txt
    echo -n "### $toto2 : " > $toto2
    date >> $toto2
    cat $toto | gawk -F '\t' '/^##/{next;}/^#/{print;}' >> $toto2


    set k=0
    foreach kk ($kks)
      set k=`echo $k | gawk '{print $1+1;}'`
      set kkD=$kk
      if ($kk == 0k) set kkD="Entire_dataset" 
      cat RESULTS/Expression.$kk/$UU/av/AECDB_diff.deg_truth.txt | gawk -F '\t' '/^##/{next;}/^#/{for(i=2;i<= NF;i++){tt[i]=$i;}nf=NF;jj=0;next;}{jj++;for(i=2;i<=nf;i++)printf("%d\t%d\t%d\t%s\t%s\t%s\t%d\n",jj,i,k,tt[i],$1,kk,0+$i);}' k=$k kk=$kkD >> toto.1
    end
    cat toto.1 | sort -k 1,1n -k 2,2n -k 3,3n > toto.2
    cat toto.2 | gawk -F '\t' '{i=$1;j=$2;k=$3;a[i]=$5;b[j]=$4;kk[k]=$6;z[i,j,k]=$7;if(i+0>iMax)iMax=i;if(j+0>jMax)jMax=j;if(k+0>kMax)kMax=k;}END{for (i=1;i<=iMax;i++){printf("\n\n##%s\tRun",a[i]);for(k2=1;k2<=kMax;k2++)printf("\t%s",kk[k2]);for (j=1;j<=jMax;j++){if (length(a[i])*length(b[j])>0){ printf("\n%s\t%s",a[i],b[j]);for(k=1;k<=kMax;k++)printf("\t%d",0+z[i,j,k]);}}}printf("\n");}' >> $toto

    cat toto.2 | gawk -F '\t' '{i=$1;j=$2;k=$3;a[i]=$5;b[j]=$4;kk[k]=$6;z[i,j,k]=$7;if(i+0>iMax)iMax=i;if(j+0>jMax)jMax=j;if(k+0>kMax)kMax=k;}END{for (j=1;j<=jMax;j++){printf("\n\n##%s\tRun",b[j]);for(k2=1;k2<=kMax;k2++)printf("\t%s",kk[k2]);for (i=1;i<=iMax;i++){if (length(a[i])*length(b[j])>0){ printf("\n%s\t%s",a[i],b[j]);for(k=1;k<=kMax;k++)printf("\t%d",0+z[i,j,k]);}}}printf("\n");}' >> $toto2

  goto done

   # print
   $2 vide ->jette
   $2 Cumul print
   $2  Traget no total polya 
      trueA trueB newtrieA newtrueB amissed Bmissed InternalIn   Opposite inconsUndec 
   $2  non tragtte
      trueA trueB newtrieA newtrueB InternalIn   Opposite inconsUndec  weak

  goto done
endif

#############################################################################

done:
    echo -n "Done "
    date

#############################################################################
#############################################################################
