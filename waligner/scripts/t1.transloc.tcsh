#!bin/tcsh -f

set phase=$1
set run=$2

if ($phase == t1) goto phaseT1
if ($phase == t2) goto phaseT2
if ($phase == t3) goto phaseT3
if ($phase == t4) goto phaseT4

exit (0)

##### Phase 1 collect the translocations runs
phaseT1:

  foreach target ($Etargets)
    if ($target == magic) continue
    # collate the gene fusions
    source scripts/target2target_class.txt
    if (! -e  tmp/Transloc/$run/t1.transloc.$target.txt.gz && -e tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt) then
      echo ' ' >  tmp/Transloc/$run/t1.transloc.$target.txt
      foreach lane (`cat Fastc/$run/LaneList`)
        echo "t1 $lane"
        if (-e tmp/COUNT/$lane.hits.gz) then
           echo "bin/tricoteur  -hitFile tmp/COUNT/$lane.hits.gz -run $run -lane $lane -target_class  $target_class -geneFusion -o tmp/Transloc/$lane -geneMap tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt"
                 bin/tricoteur  -hitFile tmp/COUNT/$lane.hits.gz -run $run -lane $lane -target_class  $target_class -geneFusion -o tmp/Transloc/$lane -geneMap tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt
          sleep 1
          ls -ls tmp/Transloc/$lane.geneFusion.txt
          cat tmp/Transloc/$lane.geneFusion.txt >>  tmp/Transloc/$run/t1.transloc.$target.txt
	  \rm  tmp/Transloc/$lane.geneFusion.txt
        endif
      end
      gzip tmp/Transloc/$run/t1.transloc.$target.txt
    endif
# count
    zcat  tmp/Transloc/$run/t1.transloc.$target.txt.gz | gawk -F '\t' '/^#/{next;}{z = $3 ;} /PAIR/{nn[z]++;np[z]++;next;}/READ/{nn[z]++;nr[z]++;}END{for(k in nn)if(nn[k]>0)printf("%s\t%d\t%d\n",k,nr[k],np[k]);}' | sort -k 2nr > tmp/Transloc/$run/t1.transloc.$target.count
  end
  
  touch tmp/Transloc/$run/t1.transloc.done

goto phaseLoop

#######################################################################################
## phase t2:   cumul the translocations per runs

phaseT2:

  echo -n "Phase t2 cumul the translocations in each runs, group and project "
  date

  if (-e MetaDB/$MAGIC/r2sublib) then
    foreach run (`cat MetaDB/$MAGIC/r2sublib | cut -f 1 | sort -u`)
      foreach target ($Etargets)
        if (! -e tmp/Transloc/$run/t1.transloc.$target.count) then 
          echo ' ' > tmp/Transloc/$run/_t
          foreach run2 (`cat MetaDB/$MAGIC/r2sublib | gawk '{if($1==run) print $2;}' run=$run | sort -u`)
            cat  tmp/Transloc/$run2/t1.transloc.$target.count >>  tmp/Transloc/$run/_t
          end
          cat tmp/Transloc/$run/_t | gawk -F '\t' '{n1[$1] += $2 ; n2[$1] += $3;}END{for (k in n1) printf("%s\t%d\t%d\n",k,n1[k],n2[k]);}' | sort -k 2nr > tmp/Transloc/$run/t1.transloc.$target.count
          \rm tmp/Transloc/$run/_t
        endif
      end
    end
  endif

  foreach level (1 2 3 4 5 6 7 8 9)
    foreach group (`cat MetaDB/$MAGIC/g2r | gawk -F '\t' '{if ($3 == level)print $1;}' level=$level | sort -u`)
      foreach target ($Etargets)
        if (! -e tmp/Transloc/$group/t1.transloc.$target.count) then 
          echo ' ' > tmp/Transloc/$group/_t
          foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if($1==group) print $2;}' group=$group | sort -u`)
            if (-e tmp/Transloc/$run2/t1.transloc.$target.count) then
            cat  tmp/Transloc/$run2/t1.transloc.$target.count >>  tmp/Transloc/$group/_t
            endif
          end
          cat tmp/Transloc/$group/_t | gawk -F '\t' '{n1[$1] += $2 ; n2[$1] += $3;}END{for (k in n1) printf("%s\t%d\t%d\n",k,n1[k],n2[k]);}' | sort -k 2nr > tmp/Transloc/$group/t1.transloc.$target.count
          \rm tmp/Transloc/$group/_t
        endif
      end
    end
  end

## Create a project table
set minFusion=10
foreach target ($Etargets)
  if ($target == magic) continue
  set toto=tmp/Transloc/t3.$MAGIC.$target.txt
  if (-e $toto.5) continue
  echo $target > tmp/Transloc/t3.$MAGIC.a
  foreach run (`cat MetaDB/$MAGIC/GroupListSorted`)
    set ff=tmp/Transloc/$run/t1.transloc.$target.count
    if (! -e $ff) continue
    cat $ff | gawk -F '\t' '{if($2 >= minF){printf("%s\t",run);print;}}' minF=$minFusion run=$run >> tmp/Transloc/t3.$MAGIC.a
  end
  echo ZZZZZ >>  tmp/Transloc/t3.$MAGIC.a
  foreach run (`cat MetaDB/$MAGIC/RunsListSorted`)
    set ff=tmp/Transloc/$run/t1.transloc.$target.count
    if (! -e $ff) continue
    cat $ff | gawk -F '\t' '{if($2 >= minF){printf("%s\t",run);print;}}' minF=$minFusion run=$run >> tmp/Transloc/t3.$MAGIC.a
  end
  # cat  tmp/Transloc/t3.$MAGIC.a | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz+0<1){if($3>=minC)gg[$1]=1;next;}}' minC=10

  cat  tmp/Transloc/t3.$MAGIC.a ZZZZZ MetaDB/$MAGIC/Run2Title.txt   | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz==2){gsub(/\"/,"",$0);r=$1;if(r2i[r]>0)title[r]=$3;next;}}{if($3+$4<minC)next;}{r=$1;i=r2i[r]+0;if(i==0){iMax++;i=iMax;r2i[r]=i;i2r[i]=r;}k=$2;if(zz == 1){nn[k]+=$3+$4;nn1[k]+=$3;nn2[k]+=$4;}n1[k,i]+=$3;n2[k,i]+=$4;nnn2+=$4;}END{printf("# Run\tAny\tSupporting reads");for(i=1;i<=iMax;i++)printf("\t%s",i2r[i]);if(nnn2>0){printf("\t\t# Run\tAny\tAny2");for(i=1;i<=iMax;i++)printf("\t%s",i2r[i]);}printf("\n# Title\tAny\tAny1");for(i=1;i<=iMax;i++){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}if(nnn2>0){printf("\t\t# Title\tAny\tSupporting read-pairs");for(i=1;i<=iMax;i++){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}}for (k in nn){ printf("\n%s\t%d\t%d",k,nn[k],nn1[k]);for(i=1;i<=iMax;i++)printf("\t%d",n1[k,i]);if(nnn2>0){printf("\t\t%s\t%d\t%d",k,nn[k],nn2[k],nn2[k]);for(i=1;i<=iMax;i++)printf("\t%d",n2[k,i]);}}}END{printf("\n");}' minC=5 > tmp/Transloc/t3.$MAGIC.b
  echo -n "### File $toto : " > $toto
  date >> $toto
  echo "### Table of candidate gene fusions in project $MAGIC. Left table: single read  support, right table: pair support" >> $toto
  cat tmp/Transloc/t3.$MAGIC.b | gawk '/^#/{print}' >> $toto

  if (! -e tmp/Transloc/$target.gene2intmap.txt) then
    cat tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt | gawk -F '\t' '{gene=$5;split($2,aa,":");chrom[gene]=aa[1];split(aa[2],bb,"-");a1=bb[1];a2=bb[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}if(0+aa1[gene]==0){aa1[gene]=a1;aa2[gene]=a2;}if(a1<aa1[gene])aa1[gene]=a1;if(a2>aa2[gene])aa2[gene]=a2;}END{for (g in chrom)printf("%s\t%s\t%d\t%d\n",g,chrom[g],aa1[g],aa2[g]);}' | sort > tmp/Transloc/$target.gene2intmap.txt
  endif
  cat  tmp/Transloc/$target.gene2intmap.txt ZZZZZ tmp/Transloc/t3.$MAGIC.b | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz<1){g=$1;gc[g]=$2;g1[g]=$3;g2[g]=$4;next;}}{if($2<5)next;split($1,aa,"__");ga=aa[1];gb=substr(aa[2],1,length(aa[2])-2);if(gc[ga] == gc[gb]){c1=g1[ga];if(c1<g1[gb])c1=g1[gb];c2=g2[ga];if(c2>g2[gb])c2=g2[gb]; if (c1<c2)next;dc=g1[gb]-g2[ga];if(dc>0 && dc < 100000)next;dc=g1[ga]-g2[gb];if(dc>0 && dc < 100000)next;}printf("%s(%s:%d-%d__%s:%d-%d)",$1,gc[ga],g1[ga],g2[ga],gc[gb],g1[gb],g2[gb]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' | sort -k 2nr >> $toto


  \cp $toto RESULTS/Transloc
end

goto phaseLoop

#####

## confirmation des transloc nanopore par illumina ?

# les coordonnees exactes sont donnes dans le fichier
#  tmp/Transloc/$run/t1.transloc.$target.txt.gz


zcat  tmp/Transloc/TotR5-100_S8/t1.transloc.av.txt.gz  | head
#### Gene fusion candidates	2020-02-21_12:18:30	file=tmp/Transloc/TotR5-100_S8/f2.1.geneFusion.txt
#Run	Read	Fusion	Gene_A	Gene_B	Chrom_A	Chrom_B	Distance	Type	x1 A	x2 A	mRNA_A	from	to	x1 B	x2 B	mRNA_B	from	to	Distinct_supports	Support	Gene_A supports	Gene_B supports	score A	score B	Ali A	Ali B	c1 A	c2 A	c1 B	c2 B
TotR5-100_S8/f2.1	n.223121#1	ABCC4__skerdorbo++	ABCC4	skerdorbo	13	3	0	READ	1	74	ABCC4.aAug10	2595	2668	73	131	skerdorbo.cAug10	266	324	1	1	132	2	0	0	74	59	1	74	73	131
TotR5-100_S8/f2.1	n.223121#1	ABCC4__skerdorbo++	ABCC4	skerdorbo	13	3	0	READ	151	123	ABCC4.aAug10	2640	2668	124	66	skerdorbo.cAug10	266	324	1	1	132	2	0	0	29	59	123	151	66	124

# le fichier tmp/Transloc/t3.crn.av.txt donne la liste des gene fusions interessantes (couver > $minC dans au moins 1 run)
head tmp/Transloc/t3.crn.av.txt
### File tmp/Transloc/t3.crn.av.txt : Thu Feb 20 17:36:37 EST 2020
### Table of candidate gene fusions in project crn. Left table: single read  support, right table: pair support
# Run	Any	Supporting reads	coco	C1_S1	C2_S2	D1_S3	E1_S4	E2_S5	F1_S6	F2_S7	S8	S9	BL10_200_B2_F1_S7	QuartetD5_500_B3_C1_S9	BL10-200-2-BL10-50-1-TotR5-200_S5	D5-25-2-TotR10-50_S7	TotR-100-TotR8-100-BL10-100-2_S8		# Run	Any	Any2	coco	C1_S1	C2_S2	D1_S3	E1_S4	E2_S5	F1_S6	F2_S7	S8	S9	BL10_200_B2_F1_S7	QuartetD5_500_B3_C1_S9	BL10-200-2-BL10-50-1-TotR5-200_S5	D5-25-2-TotR10-50_S7	TotR-100-TotR8-100-BL10-100-2_S8
# Title	Any	Any1	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL		# Title	Any	Supporting read-pairsNULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL
RN7SL2andRN7SL3__LOC723809++(14:50320341-50329627__7:104535070-104567092)	111	111	151	0	0	0	35	38	14	13	0	0	11	0	RN7SL2andRN7SL3__LOC723809++	111	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0

### idee generale
### on lit /t3.$MAGIC.$target.txt, puis  tmp/Transloc/TotR5-100_S8/t1.transloc.av.txt.gz et on compte pour chaque run le
### soutient aux cordonnes exactes au format tsv

### use variant_caller.c : -geneFusion geneFusionFileName
setenv MAGIC crn
set target=av

foreach run (`cat MetaDB/$MAGIC/RunsList`)
  bin/variant_caller -run $run -t $target -geneFusionFile tmp/Transloc/t3.$MAGIC.$target.txt -geneFusionPositions   tmp/Transloc/$run/t1.transloc.$target.txt.gz -o   tmp/Transloc/$run/t4.transloc.$target
end

# debug the fabrication of t1

bin/tricoteur -hitFile toto.hits -run xx -target_class ET_av -geneFusion -geneMap tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt -o tatou


foreach run (COV-20200315-P10-A12-P  COV-20200316-P12-D04-P)
  scripts/submit tmp/Transloc/$run/tricot2 "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 200 -dx 12  -o tmp/Transloc/$run/tricot2"
end
foreach run (`cat MetaDB/M73/RunsList`)
  scripts/submit tmp/Transloc/$run/tricot2 "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 200 -dx 12  -o tmp/Transloc/$run/tricot2"
end
variant_caller  itarget_class Z_genome  -run COV-20200314-P9-B05-P -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 3000 -minSnpCount 150 -dx 12 -t1 25000 -t2 25800


foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricot10.bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricot10 "variant_caller -target_class Z_genome  -run $run -maxLanes 10 -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 200 -minSnpCount 20 -dx 8 -o tmp/Transloc/$run/tricot10"
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricot5.bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricot5 "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 200 -minSnpCount 20 -dx 8 -o tmp/Transloc/$run/tricot5"
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotD.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotD "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 10 -dx 8 -t1 500 -t2 1000 -o tmp/Transloc/$run/tricotD "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotE.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotE "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 20 -dx 12 -min_bridge 15000 -o tmp/Transloc/$run/tricotE "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotG.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotG "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 20 -dx 12  -o tmp/Transloc/$run/tricotG "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotH.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotH "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 10 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotH "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotI.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotI "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 10 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotI "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotJ.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotJ "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 10 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotJ "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotK.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotK "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 100 -minSnpCount 10 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotK "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotL.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotL "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotL "
  endif
end

# added a count in bifurcate
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotM.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotM "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotM "
  endif
end
# added a proba in bifurcate, but removed the count that trimmed iMax in bifurcate
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotN.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotP "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotP "
  endif
end
# in biffurcate, cMax < 128 , also fixed an error on short del 9 bases in tctCompressSegs
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotN.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotQ "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotQ "
  endif
end
# in biffurcate, prune iiStep if arrayMax(hits) > 1000
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotN.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotR "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotR "
  endif
end
foreach run (`cat MetaDB/MM/RunsList | grep COV`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotN.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotS "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 999 -o tmp/Transloc/$run/tricotS "
  endif
end
# 8 threads
foreach run (`cat MetaDB/MMM/RunsList | grep wist`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotN.Bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotT "variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotT "
  endif
end

goto phaseLoop

#######################################################################################
## phase t2:   cumul the translocations per runs
phaseT2Cov:

  echo -n "Phase T2 cumul the translocations in each runs, group and project "
  date

set chrom=mutated_cov_May7
set chrom=NC_045512
foreach run (`cat MetaDB/$MAGIC/RunsList`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotY.bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotY "bin/variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/$species.genome.fasta.gz -t $chrom -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -o tmp/Transloc/$run/tricotE "
  endif
end

set chrom=mutated_cov_May7
set chrom=NC_045512
foreach run (`cat MetaDB/$MAGIC/RunsList`)
  if (! -d  tmp/Transloc/$run) mkdir  tmp/Transloc/$run
  if (! -e tmp/Transloc/$run/tricotZ.bridges.ace) then
    scripts/submit tmp/Transloc/$run/tricotZ "bin/variant_caller -target_class Z_genome  -run $run -target_fasta TARGET/Targets/$species.genome.fasta.gz -t $chrom -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8  -maxLanes 10 -subSampling 10 -o tmp/Transloc/$run/tricotZ "
  endif
end

goto phaseLoop


goto phaseLoop
##############################################


cat <<EOF > tmp/Transloc/t2.d.awk
  /:Del_/ { 
    v = \$1 ; r = \$2 ; m = \$4 ; w1 = \$8 ; w2 = \$9 ;
    z  = -20 ; if (m+w1+w2 >= 20) z = 100.0 * m / (m+w1+w2) ;
    z1 = -20 ; if (m+w1 >= 20) z1 = 100.0 * m / (m+w1) ;
    z2 = -20 ; if (m+w2 >= 20) z2 = 100.0 * m / (m+w2) ;
    printf ("%s\t%s\t7\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n",v,r,z,m,w1,z1,m,w2,z2);
} 
EOF

cat tmp/Transloc/t2.$MAGIC.val.tsf | gawk -F '\t' '{v=$1;r=$2;m=$4;w=$5;m1=$6;m2=$7;w1=$8;w2=$9;if(m+w==0){m=(m1+m2)/2;w=(w1+w2)/2;}w3=m+w;if(w3==0)w3=1;z=100.0*m/w3;if(w3<20)z=-20;w3=m1+w1;z1=-20;if(w3>=20)z1=100.0*m1/w3;z2=-20;m3=m2+w2;if(w3>=20)z2=100.0*m2/w3;printf("%s\t%.2f\t%s\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d\t%d==\n",v,z,r,m,w,z1,m1,w1,z2,m2,w2);}' | sort -k 1,1 -k 2,2nr -k 3,3 > tmp/Transloc/t2.$MAGIC.val.tsf.sorted





cat tmp/Transloc/t2.$MAGIC.val.tsf.sorted | gawk -F '\t' '{if($1 != old){old=$1;printf("\n\nVariant \"%s\"\n-D Validation",$1);}r=$3;z=$2;m=$4;w=$5;w1=$6;w2=$7;printf("\nValidation \"%s\" %.2f %d variant %d reference",r,z,m,w);if(w1+w2>0)printf(" %d refD %d refA",w1,w2);}END{printf("\n\n"); }' > tmp/Transloc/t2.$MAGIC.val.ace
 

if (0) then

An interesting set of long deletions joining the motif 70-75:ACGAAC to each of its 8 repeats was observed in a majority of the samples with good virus coverage. They are supported by  thousands of reads representing up to 40% of  the reads ignoring the breakpoints. They effectively join the 5' end of the virus, known to  contain the transcription initiation site (reference) to the beginning of the majority of the functional annotated open reading frames : M, ORF3a,  (see the table). The notable exception are ORF1ab, ORF7b and ORF10, for which there there are no deletions and no ACGAAC site. The second C of the ACGAAC is immediatly in front of the A of the ATG methionine codon of ORF7a and ORF8, 1 base upstream in S, 2 base upstream in ORF3a and E, 8 bases in N and 44 and 155 in M and ORF6. This raises the possibility of a functonal extension of the translated proteins upstream of the current annotation of M and ORF6. The 15 amino acids encoded upstream of M (NIILVFLFGTLILA) are significantly homologous to xxx. Because of our very deep coverage, our measurments are more precise but similar to reference yyy.   


endif

# Sedlazeck
 cat ../RESULTS/Sedlazeck_merged_svtyper_dist100bp_max1kbp.vcf | gawk -F '\t' '/^NC_045512/{print}' | cut -f 1-5 | gawk -F '\t' '/MantaDEL/{a1=$2;w=$4;n=length(w)-1;printf("Variant NC_045512:%d:Del_%d\nIntMap NC_045512 %d %d \"%d bases %d to %d deleted\"\nVCF %d %s %s\nParent_sequence  NC_045512\nDeletion\nRicha\n\n", a1,n,a1,a1+n+1,n,a1+1,a1+n,a1,$4,$5);}' 

# check is the splicing of all 8long deletions are correlated: one sample = height high deletaion rate ?
tace MetaDB << EOF
  select v from v in ?variant where v->intmap[2] < 90 && v->intmap[3] > 20000
  select -o tatou v,s,r,f,m,w,w1,w2 from v in @, r in v->validation, s in r->sorting_title, f in r[1],m in f[1],w in m[2],w1 in w[2],w2 in w1[2]
EOF
cat tatou | gawk '{v=$1;r=$2"\t"$3;f=$4;vv[v]=1;rr[r]=1;z[v,r]=f;if(f>zm[v])zm[v]=f;}END{printf("\t");for(v in vv)if(zm[v]>5)printf("\t%s",v);for (r in rr){printf("\n%s",r);for(v in vv)if(zm[v]>5)printf("\t%s",z[v,r]);}printf("\n");}' > tatou.txt
cat tatou.txt | head -1 > RESULTS/BigDelFreq.txt
cat tatou.txt | tail -n +2 | sort -V >> RESULTS/BigDelFreq.txt





# 300 base from Dan modified
Dan attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaac aactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgtgtgcagccgatcatcagcacatctaggtttcgtccgggtgtgaccgaaaggtaagatggagagccttgtccctggtttcaacgagaaaac
Cor attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgt tgcagccgatcatcagcacatctaggtttcgtccgggtgtgaccgaaaggtaagatggagagccttgtccctggtttcaacgagaaaac
                                                                           G                                                                                                                                           G


  scripts/submit  tmp/Transloc/$run/t2.snpSub29c  "hello"


set run=COV-20200315-P10-A12-P
echo "bb45_70\tgatctcttgtagatctgttctctaaa" >  tmp/Transloc/t2.extend.motifs 
bin/variant_caller  -extend -wLn 26 -snpFile  tmp/Transloc/t2.extend.motifs -t toto -run $run -o   tmp/Transloc/$run/t2.extend -maxLanes 1  



goto phaseLoop

#####################################################################################################
## analysis of every deletions in the virus

# MAGIC d1
# collate all deletions
zcat tmp/OR/*/d1.*gz | gawk -F '\t' '{a=$5;b=$7;if(a>b){c=a;a=b;b=c;}if(b>a+1)n[a"-"b]+=$11;d[a"-"b]=a"\t"b"\t"b-a-1}END{for(k in n)printf("%s\t%s\t%d\n",k,d[k],n[k]);}' | sort -k 5nr  > _a

# Count all starting in the zone 1 - 100
gawk '{if($2<100)n[int($3/100)]+=$5;}END{for(k in n)printf("%d\t%d\n",100*k,n[k]);}' | sort -k 1nr > RESULTS/DelStaringBefore100.Arrival.txt

# density of all bridges at several dentities
cat _a | gawk '{n[int($2/10),int($3/100)]+=$5;}END{for(i=0;i<=300;i++){printf("%d",100*i);for(j=0;j<=300;j++)printf("\t%d",n[i,j]);printf("\n");}}' > RESULTS/DelDensite.300.txt
cat _a | gawk '{n[int($2/10),int($3/10)]+=$5;}END{for(i=0;i<=3000;i++){printf("%d",10*i);for(j=0;j<=3000;j++)printf("\t%d",n[i,j]);printf("\n");}}' > RESULTS/DelDensite.3000.txt

# confirmations: enter all the deletions in the database

Variant "NC_045512:64:Del_28190::" // seg->a1/a2=61/28265, b1/b2=53/28273 iSeg=15
VCF 64  
Multi_deletion 28190 
Observed_genomic_sequence GTAGATCTGTTCTCTAAACGAACAAACTAAA
Method Magic4
IntMap "NC_045512" 64 28255 "Bases 65 to 28254 (28190 bases) are deleted" 
Found_in_genome
Parent_sequence "NC_045512"
fCounts COV-20200314-P9-A05-P 135 1095 1230 Frequency 11.0
rCounts COV-20200314-P9-A05-P 44 607  651 Frequency 6.8
nsCounts COV-20200314-P9-A05-P 179 1702  1881 Frequency 9.5


cat _a | gawk '{a1=$2+0;a2=$3+0;da=$4+0;nn=$5+0;printf("Variant \"%s:%d:Del_%d::\"\nVCF %d\nMulti_deletion %d\nMethod MagicWG\nIntMap \"%s\" %d %d\nFound_in_genome\nParent_sequence \"%s\"\nnsCounts Any %d\n\n",t,a1,da,a1,da,t,a1,a2,t,nn);}' t=NC_045512 > allDel.ace


#####################################################################################################
## analysis of every deletions in the MALAT NEAT reagion
cat ../tmp/METADATA/av.ns.gene.sponge | grep NEAT
cat ../tmp/METADATA/av.ns.gene.sponge | grep MALAT
NEAT    11 65190245 65213016
MALAT   11 65265233 65273982
set chrom=11
set a1=65190000
set a2=65280000

# MAGIC d1
# collate all deletionsin this 90kb region
set toto=RESULTS/MALAT_NEAT.deletions.txt
echo -n "### $toto : " > $toto
date >> $toto
echo "## Deletion in MM project in chromosome 11, 65.190 M to 65.280 M, 90kb region" >> $toto
echo "## The deletions reported here are found by the aligner inside single molecules, if the 2 reads of the pair support the deletion, they count as 2." >> $toto
zcat tmp/OR/*/d1.*gz | gawk -F '\t' '{k=0;}/^DELETION/{k=1;}/INTRON/{k=1;}{x=$5;y=$7;if(x>y){z=x;x=y;y=z;}if(k==1 && $3==chrom && $6 == chrom && x>a1 && x<a2 && y>a1 && y<a2){d[x"\t"y]=$10;n[x"\t"y]+= $11;}}END{for(k in n)printf("%s\t%s\t%d\n",k,n[k],d[k]);}' chrom=$chrom a1=$a1 a2=$a2 | sort -k 1n >> $toto

# Count all starting in the zone 1 - 100
gawk '{if($2<100)n[int($3/100)]+=$5;}END{for(k in n)printf("%d\t%d\n",100*k,n[k]);}' | sort -k 1nr > RESULTS/DelStaringBefore100.Arrival.txt

# density of all bridges at several dentities
cat _a | gawk '{n[int($2/10),int($3/100)]+=$5;}END{for(i=0;i<=300;i++){printf("%d",100*i);for(j=0;j<=300;j++)printf("\t%d",n[i,j]);printf("\n");}}' > RESULTS/DelDensite.300.txt
cat _a | gawk '{n[int($2/10),int($3/10)]+=$5;}END{for(i=0;i<=3000;i++){printf("%d",10*i);for(j=0;j<=3000;j++)printf("\t%d",n[i,j]);printf("\n");}}' > RESULTS/DelDensite.3000.txt

# confirmations: enter all the deletions in the database

cat _a | gawk '{a1=$2+0;a2=$3+0;da=$4+0;nn=$5+0;printf("Variant \"%s:%d:Del_%d::\"\nVCF %d\nMulti_deletion %d\nMethod MagicWG\nIntMap \"%s\" %d %d\nFound_in_genome\nParent_sequence \"%s\"\nnsCounts Any %d\n\n",t,a1,da,a1,da,t,a1,a2,t,nn);}' t=NC_045512 > allDel.ace





#########################
## parse the official list of snps
## get a list of 2699 substtituion format : T2445G
 cat ../RESULTS/Apr23_UCSC_nextstrainSamples.vcf.txt | cut -f 3 | sed -e 's/,/\n/g' | sort -u > _s
# verify
cat _s | gawk '{a=substr($1,1,1);n=length($1);g=substr($1,n);x=0+substr($1,2,n-2); printf("%s\t%s\t%s\t%d\n",$1,a,g,x);}' > _s1
head _s1
# create an ace file

cat _s1 | gawk -F '\t' '/^#/{next;}{x=$4;a=$2;g=$3;nam=t ":" x ":Sub:" a ":" g ;printf("Variant %s\nIntMap %s %d %d \"Base %d:%s>%s is modified\"\nVCF %d %s %s\nUCSC\nParent_sequence %s\n%s2%s\n\n",nam,t,x,x+1,x,a,g,x,a,g,t,a,g);}' t=NC_045512 > Apr23_UCSC_official_Substitutions.ace



## parse the snps from Chris mason 2020_06_06
cat variants_June5_vcf_jor_vars.clades.table.txt | gawk -F '\t' '/SNV/{x=$4;a=$8;g=$9;printf("Variant NC_045512:%s:Sub:%s:%s\nIntMap %s %d %d \"Base %s:%s>%s is modified\"\nVCF %s %s %s\nRicha\nParent_sequence %s\n%s2%s\n\n",x,a,g,c,x,x+1,x,a,g,x,a,g,c,a,g);}' c=NC_045512 > variants.sub.mason.june5.ace
cat variants_June5_vcf_jor_vars.clades.table.txt | gawk -F '\t' '/DEL/{x=$4;a=$8;g=$9;da=length(a)-length(g);a=substr(a,1,da+1);g=substr(g,1,1);if(da<2)next;printf("Variant NC_045512:%s:Del_%d:%s:%s\nIntMap %s %d %d \"Base %d to %d (%d bases) are deleted\"\nVCF %s %s %s\nRicha\nParent_sequence %s\nMulti_deletion %d\n\n",x,da,a,g,c,x,x+da,x+1,x+da+1, da,x,a,g,c,da);}' c=NC_045512 > variants.del.mason.june5.ace
cat variants_June5_vcf_jor_vars.clades.table.txt | gawk -F '\t' '/DEL/{x=$4;a=$8;g=$9;da=length(a)-length(g);a=substr(a,1,da+1);g=substr(g,1,1);if(da!=1)next;printf("Variant NC_045512:%s:Del:%s:%s\nIntMap %s %d %d \"Base %d:%s is deleted\"\nVCF %s %s %s\nRicha\nParent_sequence %s\nDel%s\n\n",x,a,g,c,x,x+da+1,x+1,substr(a,2),x,a,g,c,substr(a,2));}' c=NC_045512 > variants.del1.mason.june5.ace


















phaseLoop:
 echo done

