#!bin/tcsh -ef

set group=$1
set chrom=$2
set justMito=$3

# construct the combined wiggle of runs with low genomic contamination

# prepare the list of BF runs

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

set out=tmp/WIGGLEGROUP/$group/$chrom
if ($justMito == 1) set out=tmp/WIGGLEGROUP/$group/$chrom.step1
if (! -d  tmp/WIGGLEGROUP/$group) mkdir tmp/WIGGLEGROUP/$group
if (! -d  tmp/WIGGLEGROUP/$group/$chrom) mkdir tmp/WIGGLEGROUP/$group/$chrom
if ($justMito == 1 && ! -d tmp/WIGGLEGROUP/$group/$chrom.step1) mkdir tmp/WIGGLEGROUP/$group/$chrom.step1

if (! -e $out/wg2b.done) then
# check that all runs are available

  foreach uu (u nu pp)
    foreach fr (f r  ELF ELR ERF ERR)
      if (-e $out/$fr.$uu.chrom.list) \rm  $out/$fr.$uu.chrom.list
      if (-e $out/$fr.$uu.genes.list) \rm $out/$fr.$uu.genes.list
      set ok=0
      foreach run (`cat tmp/WIGGLEGROUP/$group/RunList`)
        set drop=0
        if ($fr ==  ELF || $fr ==  ELR || $fr == ERF || $fr == ERR) then
          foreach run2 (`cat MetaDB/$MAGIC/WiggleDropEndList`)
            if ($run2 == $run) set drop=1
          end
          if ($drop == 1) continue 
        endif
        if (-e tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
          echo -n " tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz " >>   $out/$fr.$uu.chrom.list
          set ok=1
        endif
        if (-e tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
          echo -n " tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz " >>   $out/$fr.$uu.chrom.list
          set ok=1
        endif

        if ($Strategy != Exome && $justMito == 1) then
          if ($fr == f && -e tmp/SNP/$run/SpikeIn.plus.u.BV.gz) then
            echo -n " tmp/SNP/$run/SpikeIn.plus.u.BV.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
          endif
          if ($fr == r && -e tmp/SNP/$run/SpikeIn.minus.u.BV.gz) then
            echo -n " tmp/SNP/$run/SpikeIn.minus.u.BV.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
          endif
        endif
      end
      echo "ok=$ok"
      if ($ok == 1) then
          # contruct the combined wiggles
        if (-e  $out/$fr.$uu.genes.list) then
          if ($Strategy != Exome && $justMito == 0 && ! -e tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.BF.gz) then
            echo  $out/$fr.$uu.genes
            gunzip -c `cat $out/$fr.$uu.genes.list` | bin/wiggle  -I BF -gzo -O BF $out_step  -o tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr -cumul >&   tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.cumul
          endif
        endif
        cat $out/$fr.$uu.chrom.list
        echo
        if ( ! -e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz) then 
          echo  $out/$fr.$uu.chrom
          cat $out/$fr.$uu.chrom.list
          echo aa
          gunzip -c `cat $out/$fr.$uu.chrom.list` | bin/wiggle  -I BV -gzo -O BF $out_step  -o  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr -cumul >&   tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.cumul 
        endif
      endif
    end
   
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.f.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.r.BF.gz) then
      echo "bb "
      ls -ls tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.[fr].BF.gz
      gunzip -c tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.[fr].BF.gz | bin/wiggle -I BF -O BF  $out_step -gzo -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu -cumul
    endif
    
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ELF.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ELR.BF.gz) then
       gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL
       if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF.gz) \rm tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF.gz
       gzip  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF
    endif
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ERF.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ERR.BF.gz) then
       gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER
       if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF.gz) \rm tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF.gz
       gzip  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF
    endif
  
  end

endif

touch $out/wg2b.done
exit 0

          # notice that /R.genes.frns has not been constructed
          # all the files called R.genes.*.*.cumul, i.e. R.genes.u.ERF.cumul, are empty
          if (-e tmp/WIGGLEGROUP/$group/$chrom/R.genes.frns.u.BF.gz && ! -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ns.sponge.total.1.txt) then 

            set fr=ns
            set mask=""
            if(-e tmp/METADATA/RefSeq.$fr.sponge) set mask=$mask",RefSeq_transcripts:tmp/METADATA/RefSeq.$fr.sponge"
            if(-e tmp/METADATA/av.$fr.sponge) set mask=$mask",AceView_transcripts:tmp/METADATA/av.$fr.sponge"
            if(-e tmp/METADATA/seqc.$fr.sponge) set mask=$mask",SEQC_transcripts:tmp/METADATA/seqc.$fr.sponge"
            if(-e tmp/METADATA/EBI.$fr.sponge ) set mask=$mask",EBI_transcripts:tmp/METADATA/EBI.$fr.sponge"
            if(-e tmp/METADATA/RefSeq.$fr.gene.sponge ) set mask=$mask",RefSeq_genes:tmp/METADATA/RefSeq.$fr.gene.sponge"
            if(-e tmp/METADATA/av.$fr.gene.sponge ) set mask=$mask",AceView_genes:tmp/METADATA/av.$fr.gene.sponge"
            if(-e tmp/METADATA/seqc.$fr.gene.sponge ) set mask=$mask",SEQC_genes:tmp/METADATA/seqc.$fr.gene.sponge"
            if(-e tmp/METADATA/EBI.$fr.gene.sponge ) set mask=$mask",EBI_genes:tmp/METADATA/EBI.$fr.gene.sponge"
            if(-e tmp/METADATA/genome.ns.sponge) set mask=$mask",Genome:tmp/METADATA/genome.ns.sponge"

            gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.genes.frns.u.BF.gz | bin/wiggle -I BF -O BF  $out_step -gzo -o tmp/WIGGLEGROUP/$group/$chrom/R.genes.u.ns

            bin/geneelements -sponge 1 -spongeFile $mask  -sxxChromosome $chrom -wiggle  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu.BF.gz | grep Total  |   sed -e "s/^Total/$chrom\t$fr/" | gawk -F '\t' '{r=run;}/^#/{r="#Run\tChrom"}{printf("%s\t",r);print;}' run=$group >  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ns.sponge.total.1.txt

          endif

exit 0

