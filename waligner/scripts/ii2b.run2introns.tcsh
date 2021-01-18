#!bin/tcsh -f

set run=$1
set uu=$2

# centralize the introm support from lane to run
if (! -e tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) then
  bin/bestali  -intronSupport2ace  -gzo -o tmp/INTRONRUNS/$run/$run.$uu -inFileList tmp/INTRONRUNS/$run/$run.$uu.list
endif


# prepare a short file used by d4 comparison to known RefSeq/av introns
if ($uu == u) then
   gunzip -c tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz | gawk '/^Intron/{split($2,aa,"__");split(aa[2],bb,"_");printf("%s\t%09d\t%09d",aa[1],bb[1],bb[2]);}/^Run_U/{printf("\t%d\n",$6);}END{printf("\n");}' | gzip >  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz

      foreach target ($Etargets)
        gunzip -c tmp/METADATA/$target.fr.introns.gz ZZZZZ.gz  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1,$2,$3]=1;next;}if (ok[$1,$2,$3]==1) {n++;nn+=$4;}}END{printf("Ali %s\nKnown_introns %s %d introns_supported_by %d reads\n\n", run, target, n,nn);}' run=$run target=$target  >  tmp/INTRONRUNS/$run/known_introns.$target.ace
      end

endif

exit 0
