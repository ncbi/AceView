#!/bin/tcsh

set phase=$1
set run=$2
set lane=$3

#limit stacksize unlimited
#limit datasize unlimited

if (! -d TARGET/Targets) exit 1
echo "# running: a123M.magicblast.tcsh $1 $2 $3"
if (1 && $phase == makeDB) then
  pushd TARGET
  if (! -d BLASTDB) mkdir BLASTDB

  foreach target ($RNAtargets  $DNAtargets)
    if ($target == gdecoy) continue
    set ff=$species.$target
    if (-e Targets/$ff.fasta.gz && ! -e BLASTDB/$ff.nsq) then
      gunzip -c Targets/$ff.fasta.gz | gawk '/^>/{n=split($1,aa,"|");ok=0;if(length(aa[1])<48){ok=1;print aa[1];}next;}{gsub("-","n",$0);if(ok==1)print;}' > BLASTDB/$ff
      pushd BLASTDB
        echo "makeblastdb -in $ff -dbtype nucl -parse_seqids #target=$target"
        makeblastdb -in $ff -dbtype nucl -parse_seqids
        if ($target == genome) then
          dna2dna -decoy -i $ff -O fasta | gawk '/^>/{split($1,aa,"|");ok=0;if(length(aa[1])<48){ok=1;print aa[1];}next;}{gsub("-","n",$0);if(ok==1)print;}' > $species.gdecoy
          echo "makeblastdb -in  $species.gdecoy -dbtype nucl -parse_seqids"
          makeblastdb -in  $species.gdecoy -dbtype nucl -parse_seqids
          \rm $species.gdecoy
        endif
        popd
      \rm  BLASTDB/$ff
    endif
  end 

  popd
  exit 0
endif

if ($phase == postMagic) then
   set ttt=`ls tmp/MAGICBLAST/$lane.*.mbhits.gz | gawk '{n++;if(n>1)t=t ","; t=t $1;}END{print t}'`
   echo $lane
   echo "bin/postMagicBlast -i $ttt -run $lane -introns -pair -tabular -gzo -o tmp/MAGICBLAST/$lane -expression -info  tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt" 
         bin/postMagicBlast -i $ttt -run $lane -introns -pair -tabular -gzo -o tmp/MAGICBLAST/$lane -expression -info  tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt
  exit 0
endif

set inok=0
foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
    if ($run == $run2) then
      set inok=1
      set infile="-query Fastc/$lane.fastc.gz -paired -infmt fastc"
    endif
end
if ($inok == 0) then
    
    set infile="-query Fastc/$lane.fastc.gz -infmt fasta"
    foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
      if ($run == $run2) then
        set infile="-query Fastc/$lane.fastc.gz -infmt fastc"
      endif
    end
endif


foreach target ($RNAtargets $DNAtargets)
  if ($target == gdecoy) continue 
  # default on fasta but prefer using a blast database
  set subject="-subject TARGET/Targets/$species.$target.fasta.gz"
  if (-e TARGET/BLASTDB/$species.$target.nsq) then
    set subject="-db  TARGET/BLASTDB/$species.$target"
  endif
  source scripts/target2target_class.txt

# if the bonus is changed, synchronize in snp.confirmation.tcsh the line
# if($8=="A_mito" || $8 == "B_rrna")s--;
# and synchronize a123M.MagicBlast.tcsh
  set bonus=0
  if ($target == gdecoy) set bonus=2
  if ($target == genome) set bonus=0
  if ($target =~ chrom*) set bonus=0
  if ($target == cloud)  set bonus=-6
  if ($target == snp)    set bonus=-6
  if ($target == introns)    set bonus=-6
  if ($target == rnaGene) set bonus=0
  if ($target == smallRNA) set bonus=0
  if ($target == mito) set bonus=1
  if ($target == chloro) set bonus=1
  if ($target == rrna) set bonus=1
  if ($target == transposon) set bonus=1
  if ($target == av2008) set bonus=-1
  if ($target == RefSeqCurrent) set bonus=-2
  if ($target == Gaj) set bonus=-3
  if ($target == FBK) set bonus=-4
  if ($target == virus) set bonus=-30
  if ($target == bacteria) set bonus=-30

# exemple de sortie de ./bin/time
# 65.28user 7.82system 1:13.40elapsed 99%CPU (0avgtext+0avgdata 1933504maxresident)k

# -outfmt tabular
  if (-e  tmp/MAGICBLAST/$lane.$target.mbhits.gz) continue
   set gt="-reftype transcriptome -splice F -limit_lookup F "
   foreach target2 ($DNAtargets)
     if ($target == $target2) then
       set gt="-reftype genome -splice T -max_intron_length $intronMaxLength"
     endif
   end

  echo "./bin/time bin/magicblast $subject $infile -out tmp/MAGICBLAST/$lane.$target.mbhits -outfmt tabular -gzo -no_unaligned $gt -num_threads 1 -tag $target_class\t$bonus"
       echo -n "$target\tStart\t"
       date
       echo -n "$target\t"  >>  tmp/MAGICBLAST/$lane.err
       (./bin/time bin/magicblast $subject $infile -out tmp/MAGICBLAST/$lane.$target.mbhits.gz -outfmt tabular -gzo -no_unaligned $gt -num_threads 1 -tag $target_class\t$bonus) >>&  tmp/MAGICBLAST/$lane.err
       echo -n "$target\tEnd\t"
       date
       ls -ls tmp/MAGICBLAST/$lane.$target.mbhits.gz       
end

echo  "a3: jobstats "
scripts/jobstats.tcsh $run $lane 1 1

gunzip -c tmp/MAGICBLAST/$lane.*.mbhits.gz | bin/postMagicBlast -run $lane -introns -expression -pair -gzo -o tmp/MAGICBLAST/$lane -info tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt
touch tmp/MAGICBLAST/$lane.align.done
echo -n "done "
date

  if (-e tmp/MAGICBLAST/$MAGIC.job_statistics.txt) \rm tmp/MAGICBLAST/$MAGIC.job_statistics.txt
  if (-d tmp/GENERUNS/$run) \rm -rf  tmp/GENERUNS/$run
  if (-d tmp/GENELANES/$run) \rm tmp/GENELANES/$lane.*
  if (-e tmp/MAGICBLAST/$run/c2.alistats.ace) \rm tmp/MAGICBLAST/$run/c2.alistats.ace

exit 0



foreach lane (`cat MetaDB/TestM/LaneList `)
   scripts/submit tmp/COUNT/$lane.toto55 "bin/postMagicBlast -run $lane -pair -clipali -i tmp/COUNT/$lane.hits.gz -o tmp/COUNT/$lane -expression"
end
qusage 1
cat tmp/COUNT/RNA_AGLR1_A_1.*/f2.*.stats.tsf | bin/tsf -g RNA_AGLR1_A_1 --sumAll -o tmp/COUNT/RNA_AGLR1_A_1.stats


foreach lane (`cat MetaDB/TestMB/LaneList `)
  scripts/submit toto55 "zcat tmp/MAGICBLAST/$lane.*.mbhits.gz | bin/postMagicBlast -run $lane -introns -pair -tabular -gzo -o tmp/MAGICBLAST/$lane -expression -info  tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt"
end


foreach lane (`cat MetaDB/TestMB/LaneList `)
  echo $lane
  scripts/submit tmp/MAGICBLAST/$lane.postMB "scripts/a123M.magicblast.tcsh postMagic xxx $lane"
end
