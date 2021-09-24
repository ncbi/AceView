#!/bin/tcsh

set phase=$1
set run=$2
set lane=$3

if (! -d TARGET/Targets) exit 1
echo "# running: a123M.magicblast.tcsh $1 $2 $3"
if (1 && $phase == makeDB) then
  pushd TARGET
  if (! -d BLASTDB) mkdir BLASTDB

  foreach target ($RNAtargets  $DNAtargets)
    if ($target == gdecoy) continue
    set ff=$species.$target
    if (-e Targets/$ff.fasta.gz && ! -e BLASTDB/$ff.nsq) then
      gunzip -c Targets/$ff.fasta.gz | gawk '/^>/{split($1,aa,"|");ok=0;if(length(aa[1])<8 && length(aa[1])+length(aa[2])<47){ok=1;print aa[1]"|"aa[2];ok=1;next;}if(length(aa[1])<48){ok=1;print aa[1];}next;}{gsub("-","n",$0);if(ok==1)print;}' > BLASTDB/$ff
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


foreach target ($RNAtargets)
  # default on fasta but prefer using a blast database
  set subject="-subject TARGET/Targets/$species.$target.fasta.gz"
  if (-e TARGET/BLASTDB/$species.$target.nsq) then
    set subject="-db  TARGET/BLASTDB/$species.$target"
  endif

# -outfmt tabular
  if (-e  tmp/MAGICBLAST/$lane.$target.sam) continue
  if (-e  tmp/MAGICBLAST/$lane.$target.sam.gz) continue
  echo "./bin/time bin/magicblast $subject $infile -gzo -out tmp/MAGICBLAST/$lane.$target.sam.gz  -no_unaligned -reftype transcriptome -splice F -limit_lookup F -num_threads 1"
       echo -n "$target\tStart\t"
       date
       echo -n "$target\t"  >>  tmp/MAGICBLAST/$lane.err
       (./bin/time bin/magicblast $subject $infile -gzo -out tmp/MAGICBLAST/$lane.$target.sam.gz  -no_unaligned -reftype transcriptome -splice F -limit_lookup F -num_threads 1) >>&  tmp/MAGICBLAST/$lane.err
       echo -n "$target\tEnd\t"
       date
       ls -ls tmp/MAGICBLAST/$lane.$target.sam.gz       
end

foreach target ($DNAtargets)
  echo "PPPPPPP $run $target"
  # default on fasta but prefer using a blast database
    if ($target == gdecoy) continue 
  set subject="-subject TARGET/Targets/$species.$target.fasta.gz"
  if (-e TARGET/BLASTDB/$species.$target.nsq) then
    set subject="-db  TARGET/BLASTDB/$species.$target"
  endif

# exemple de sortie de ./bin/time
# 65.28user 7.82system 1:13.40elapsed 99%CPU (0avgtext+0avgdata 1933504maxresident)k
  
 if (! -e  tmp/MAGICBLAST/$lane.$target.sam && ! -e  tmp/MAGICBLAST/$lane.$target.sam.gz) then
    echo "./bin/time bin/magicblast $subject $infile -gzo -out tmp/MAGICBLAST/$lane.$target.sam.gz -no_unaligned -reftype genome -max_intron_length $intronMaxLength -num_threads 1" 
         echo -n "$target\tStart\t"
         date
         echo -n "$target\t"  >>  tmp/MAGICBLAST/$lane.err
         (./bin/time bin/magicblast $subject $infile -gzo -out tmp/MAGICBLAST/$lane.$target.sam.gz -no_unaligned -reftype genome -max_intron_length $intronMaxLength -num_threads 1) >>&  tmp/MAGICBLAST/$lane.err
         echo -n "$target\tEnd\t"
	 date
         ls -ls tmp/MAGICBLAST/$lane.$target.sam.gz
 endif

# ATTENTION, we assume that the run is stranded plus
# because otherwise we of not know the strand of the intron
 if (-e tmp/MAGICBLAST/$lane.$target.sam.gz && ! -e  tmp/MAGICBLAST/$lane.$target.introns.gz) then 
    # zcat tmp/MAGICBLAST/$lane.$target.sam.gz | gawk -F '\t' -f scripts/sam2introns.awk | gzip > tmp/MAGICBLAST/$lane.$target.introns.gz
 endif
end

touch  tmp/MAGICBLAST/$lane.done

## a kind of bestali
foreach target ($DNAtargets $RNAtargets) 
  if (-e  tmp/MAGICBLAST/$lane.$target.hits.gz && ! -e tmp/MAGICBLAST/$lane.$target.sorted_hits.gz) then 
    if (! -e  tmp/MAGICBLAST/caption) then
      gunzip -c  tmp/MAGICBLAST/$lane.$target.hits.gz  | head -3 | tail -1 >  tmp/MAGICBLAST/caption
    endif
    gunzip -c  tmp/MAGICBLAST/$lane.$target.hits.gz  | gawk -F '\t' '/^#/{next;}{printf("%s\t%09d\t%s\t", $1,$25,target); print;}' | sort | gzip >  tmp/MAGICBLAST/$lane.$target.sorted_hits.gz
  endif
end

zcat  tmp/MAGICBLAST/$lane.*.sorted_hits.gz | sort | gawk -F '\t' '{if($1 == p && $2 < s)next ; p = $1;  s = $2; print}' | gzip >  tmp/MAGICBLAST/$lane.bestali.gz


echo  "a3: jobstats "
scripts/jobstats.tcsh $run $lane 1
touch tmp/COUNT/$lane.align.done
echo -n "done "
date

  if (-e tmp/COUNT/$MAGIC.job_statistics.txt) \rm tmp/COUNT/$MAGIC.job_statistics.txt
  if (-d tmp/GENERUNS/$run) \rm -rf  tmp/GENERUNS/$run
  if (-d tmp/GENELANES/$run) \rm tmp/GENELANES/$lane.*
  if (-e tmp/COUNT/$run/c2.alistats.ace) \rm tmp/COUNT/$run/c2.alistats.ace
exit 0

