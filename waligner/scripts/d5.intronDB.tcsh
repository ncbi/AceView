#!bin/tcsh -f

set phase=$1
set project=$2
echo "phase=$phase  project=$project"
if ($phase == init) then
  if (! -d GeneIndexDB) then
    scripts/geneindex.tcsh g1a
  endif

  if (! -e GeneIndexDB/d5.init.done) then
    set toto=GeneIndexDB/d5.init
    echo 'read-models' > $toto
    foreach target ($Etargets)
      echo "pparse tmp/METADATA/$target.counts.ace" >> $toto
      echo "pparse tmp/METADATA/$target.GENE.info.ace" >> $toto
      echo "pparse tmp/METADATA/$target.GENE.ln.ace" >> $toto
      echo "pparse tmp/METADATA/$target.MRNA.info.ace" >> $toto
      echo "pparse tmp/METADATA/$target.MRNA.ln.ace" >> $toto
      echo "pparse tmp/METADATA/gtf.$target.goodProduct.ace" >> $toto
      echo "pparse tmp/METADATA/gtf.$target.f.intron.ace" >> $toto
      echo "pparse tmp/METADATA/gtf.$target.r.intron.ace" >> $toto
    end
    echo save >> $toto
    echo quit >> $toto

   bin/tacembly GeneIndexDB < $toto 

   touch GeneIndexDB/d5.init.done 
  endif

  goto phaseLoop
endif

if ($phase == cumul) then
  set toto=GeneIndexDB/$project.intron_confirmation.ace
  echo $toto
  if (! -e $toto) then
    if (-e  RESULTS/Expression/unique/introns/$project.introns.INTRON.u.ace.gz2) then
      zcat RESULTS/Expression/unique/introns/$project.introns.INTRON.u.ace.gz | gawk '/^Intron/{z=$0}/_SumOfAllReadsInProject/{if ($4 >0) printf("%s\nRun_U %s %s %s %s %s %s\n\n",z,magic,$3,$4,$5,$6,$7);}' magic=$project > $toto
      echo "pparse $toto" | bin/tacembly GeneIndexDB -no_prompt
    else
      echo "Missing file  RESULTS/Expression/unique/introns/$project.introns.INTRON.u.ace.gz"
    endif
  endif
  goto phaseLoop
endif


phaseLoop:
  echo d5.intronDB phase $phase  done
