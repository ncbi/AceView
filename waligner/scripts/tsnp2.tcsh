#!bin/tcsh -f

set phase=$1
set zone=$2
set run=$3

if ($phase == tsnp2a) goto tsnp2a
if ($phase == tsnp2b) goto tsnp2b
if ($phase == tsnp2c) goto tsnp2c

echo "bad phase in tsnp2.tcsh $1 $2 $3"
date
exit 1

tsnp2a:

if (! -e tmp/TSNP_DB/$zone/database) then
   scripts/tsnp_DB.tcsh $Strategy $zone
endif





  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/TSNP/$run/$zone/tsnp1.2.0.deUno.tsf && ! -e tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf) then
      cat  tmp/TSNP/$run/$zone/tsnp1.2.*.deUno.tsf | bin/tsf --merge > tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp1.MB.0.deUno.tsf && ! -e tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf) then
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.deUno.tsf | bin/tsf --merge > tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf
    endif    
  end

  bin/tace tmp/TSNP_DB/$zone <<EOF
    pparse MetaDB/$MAGIC/runs.ace
    // pparse tmp/SNP_ZONE/$zone.fasta.gz
   save
   quit
EOF

  echo "-R Map NC_045512.2  NC_045512\n" >  tmp/TSNP/map.rename.ace

  set toto=tmp/TSNP_DB/$zone/tsnp2._r
    echo ' ' > $toto.ace
  echo "read-models" > $toto
  echo "query find variant" >> $toto
  echo "kill"  >> $toto
  echo "pparse $toto.ace" >> $toto
  #echo "pparse MetaDB/$MAGIC/runs.ace" >> $toto
  #echo "pparse tmp/SNP_ZONE/$zone.fasta.gz" >> $toto

  set foundIn=Found_in_genome
  set mapIn=IntMap
  if ($Strategy == RNA_seq)  then
    set foundIn=Found_in_mRNA
    set mapIn=mRNA
  endif
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set minSnpFrequency2=$minSnpFrequency
    if (-e tmp/TSNP/$MAGIC.minFrequency.txt) then
      # reset minSnpFrequency to at least twice the percent error rate, the error_profile is given in error per million 
      set z=`cat tmp/TSNP/$MAGIC.minFrequency.txt | gawk '{if ($1 == run) z=($4 + 0)/5000 ; if (mf<z)mf=z;}END{print int(mf);}' mf=$minSnpFrequency run=$run`
      set minSnpFrequency2=$z
    endif

    if (-e tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf) then
      cat tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf | gawk -F '\t' '/^#/{next;}{v=$1;run=$2;a1=$4;a2=$5;m=$7;c=$9;tag=$12;if(substr($12,1,6)=="Multi_")tag=$12" " $13" "$14;if(m>c)c=m;r=c-m;if(c>minC)f=100*m/c;else f=-10;if(f>=minF){if(v!=oldV){split(v,aa,":");seq=aa[1];printf("\nVariant %s\n%s\nParent_sequence %s\n%s %s %d %d\n%s\n",v,foundIn,seq,mapIn,seq,a1,a2,tag);oldV=v;printf("MCounts %s %d %d %d Frequency %.2f\n",run,m,r,c,f);}}}' minF=$minSnpFrequency2 minC=$minSnpCover foundIn=$foundIn mapIn=$mapIn >> $toto.ace
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf) then
      cat tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf | gawk -F '\t' '/^#/{next;}{v=$1;run=$2;a1=$4;a2=$5;m=$7;c=$9;if(m>c)c=m;tag=$12;if(substr($12,1,6)=="Multi_")tag=$12" " $13" "$14;if(m>c)c=m;r=c-m;if(c>minC)f=100*m/c;else f=-10;if(f>=minF){if(v!=oldV){split(v,aa,":");seq=aa[1];printf("\nVariant %s\n%s\nParent_sequence %s\%s %s %d %d\n%s\n",v,foundIn,seq,mapIn,seq,a1,a2,tag);oldV=v;printf("MBCounts %s %d %d %d Frequency %.2f\n",run,m,r,c,f);}}}' minF=$minSnpFrequency2  minC=$minSnpCover foundIn=$foundIn  mapIn=$mapIn >> $toto.ace
    endif    
  end

  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/TSNP/$run/$zone/tsnp1.2.0.bridges.ace) then
      cat  tmp/TSNP/$run/$zone/tsnp1.2.*.bridges.ace | gawk '/^Observed_genomic_sequence/{next;}/^Reference_genomic_sequence/{next;}/^[fr]Counts/{next;}/n[su]Counts/{next;}{print}' >> $toto.ace
      cat  tmp/TSNP/$run/$zone/tsnp1.2.*.bridges.ace | gawk '/^Variant/{v=$2;}/^nsCounts/{m[v]+=$3;r[v]+=$4;c[v]+=$5;}END{for(v in c){f=100*m[v]/c[v];if (f>=-95)printf("Variant %s\n%s %s %d %d %d Frequency %.2f\n\n",v,cc,run,m[v],r[v],c[v],100.0*m[v]/c[v]);}}' run=$run cc=MCounts  >> $toto.ace
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp1.MB.0.bridges.ace) then
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.bridges.ace | gawk '/^Observed_genomic_sequence/{next;}/^Reference_genomic_sequence/{next;}/^[fr]Counts/{next;}/n[su]Counts/{next;}{print}' >> $toto.ace
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.bridges.ace | gawk '/^Variant/{v=$2;}/^nsCounts/{m[v]+=$3;r[v]+=$4;c[v]+=$5;}END{for(v in c){f=100*m[v]/c[v];if (f>=-95)printf("Variant %s\n%s %s %d %d %d Frequency %.2f\n\n",v,cc,run,m[v],r[v],c[v],100.0*m[v]/c[v]);}}' run=$run cc=MBCounts  >> $toto.ace
    endif
  end
  echo "pparse tmp/TSNP/map.rename.ace" >> $toto
  echo 'save' >> $toto
  echo 'quit' >> $toto

  mv $toto.ace $toto.preace
  cat $toto.preace | gawk '/^-/{next}{print}' > $toto.ace

  bin/tace tmp/TSNP_DB/$zone <  $toto

  set remap2g=remap2genes
  if ($Strategy == RNA_seq) set remap2g=remap2genome

  bin/tsnp -db_$remap2g  tmp/METADATA/mrnaRemap.gz  -db tmp/TSNP_DB/$zone 
#  bin/tsnp -db_translate -db tmp/TSNP_DB/$zone 
  bin/tsnp -db_translate -db tmp/TSNP_DB/$zone > tmp/TSNP_DB/$zone/tsnp2.translate.ace
  echo "pparse tmp/TSNP_DB/$zone/tsnp2.translate.ace" |  bin/tacembly tmp/TSNP_DB/$zone -noprompt

touch tmp/TSNP_DB/$zone/tsnp2a.done

goto phaseLoop

########################################################################################
## tsnp2b make words

tsnp2b:
# bin/tsnp -db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone -gzo -o tmp/TSNP_DB/$zone/tsnp2b.$MAGIC --filter "MBcounts && substitution"
bin/tsnp -db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone -gzo -o tmp/TSNP_DB/$zone/tsnp2b.$MAGIC 
touch tmp/TSNP_DB/$zone/tsnp2b.done

goto phaseLoop

########################################################################################
## tsnp2b count words in run

tsnp2c:
  bin/tricoteur -count -wLn 31 -wordFile tmp/TSNP_DB/$zone/tsnp2b.$MAGIC.w31.gz -run $run -gzo -o tmp/TSNP/$run/$zone/tsnp2c.$MAGIC

goto phaseLoop

########################################################################################
## junk
if (0) then



    bin/tace tmp/TSNP_DB/$zone <<EOF
      pparse MetaDB/$MAGIC/runs.ace
      pparse tmp/SNP_ZONE/$zone.fasta.gz
      save
      quit     
EOF
    if (-d tmp/SNP_DB/$zone/database && ! -e tmp/TSNP_DB/$zone/tsnp2.$MAGIC.w31.gz) then
      bin/tsnp -db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone -gzo -o tmp/TSNP_DB/$zone/tsnp2.$MAGIC --filter "MBcounts && substitution"
    endif
 

endif

########################################################################################

phaseLoop:
  echo "phase $phase done"

########################################################################################
########################################################################################
