#!bin/tcsh -ef

set chrom=$1

  if (! -d tmp/YH$chrom) mkdir tmp/YH$chrom
  if (! -d tmp/YH$chrom/database) then
    pushd tmp/YH$chrom
    mkdir database
    ln -s ../../metaData/wspec.aceview_web_site wspec
    popd
    echo y | bin/tacembly  tmp/YH$chrom -no_prompt
  endif


  if (! -e tmp/YH$chrom/ff1.genome.done) then

    bin/tacembly tmp/YH$chrom <<EOF
      pparse TARGET/CHROMS/$species.chrom_$chrom.fasta.gz
      pparse  tmp/f1.strategy.ace
      pparse  tmp/pA/$chrom/pA.$chrom.feature.ace 
      pparse  tmp/TABIX/$chrom.tabix.ace
      pparse MetaDB/$MAGIC/runs.ace
      query find sequence $chrom
      kstore ss
      acem
        make_subseq -dna c t$chrom 600000 10000 // this breaks my 6 contigs into tiles
        quit                    // oct 15 2001, i changed from 400 kb to 600 kb
      kget ss
      Follow DNA
      kill 
      query find sequence YBR_contig // c_NT_sequences 
      kill
      query find sequence genomic
      follow source
      edit  YBR_contig  // c_NT_sequences 
 
      save
      quit
EOF
 
    bin/gene2chrom2 -any -gs -i tmp/YH$chrom  >! tmp/YH$chrom/g2c.gsi.ace

    bin/tacembly  tmp/YH$chrom << EOF
      pparse  tmp/YH$chrom/g2c.gsi.ace
      save
      quit
EOF

  touch  tmp/YH$chrom/ff1.genome.done
endif

  if (-d   tmp/YH$chrom/data) \rm -rf  tmp/YH$chrom/data
  mkdir  tmp/YH$chrom/data
  echo "pparse tmp/XH$chrom/f3.genes.ace" > tmp/YH$chrom/data/_r

  echo "pparse tmp/EHITS.*/$chrom/f3.RefSeq.fasta.gz" >> tmp/YH$chrom/data/_r
  echo "pparse tmp/EHITS.*/$chrom/f3.RefSeq.intmap.ace.gz" >> tmp/YH$chrom/data/_r

  foreach pp (Fatigue NB3 SEQC)
    set ff=/home/mieg/$pp/tmp/XH$chrom/f5.mrna.fasta
    if (-e $ff && ! -e tmp/YH$chrom/data/$pp.mrna.fasta) then
      cat $ff | gawk '/^>/{split($1,aa,"|");s=aa[1];i=index(s,"_");printf("%s_%s_%s\n", substr(s,1,i-1),pp,substr(s,i+1));next;}{print}' pp=$pp > tmp/YH$chrom/data/$pp.mrna.fasta
    endif
    echo "pparse  tmp/YH$chrom/data/$pp.mrna.fasta" >> tmp/YH$chrom/data/_r

    set ff=/home/mieg/$pp/tmp/XH$chrom/f5.mrna2map
    cat $ff.pretxt  | gawk -F '\t' '/\"/{gsub(/\"/,"",$0); printf("S%s_%s\t%s\t%d\t%d\n", chrom, $1,$2,$3,$4);}' chrom=$chrom > $ff.txt
    if (-e $ff.txt) then
      cat $ff.txt | gawk '{gsub(/\"/,"",$0); i = index($1,"_");s=substr($1,i+1);printf("Sequence S%s_%s_%s\nIntMap %s %d %d\n\n", chrom, pp,s,$2,$3,$4);}' chrom=$chrom pp=$pp    > tmp/YH$chrom/ff1.$pp.mrna.intmap.ace
      echo "pparse  tmp/YH$chrom/ff1.$pp.mrna.intmap.ace" >> tmp/YH$chrom/data/_r
    endif

    set ff=/home/mieg/$pp/tmp/EHITS.*/$chrom/f3.pA.ace
    if (-e $ff) then
      echo "pparse $ff" >> tmp/YH$chrom/data/_r
    endif
  end

  foreach ff (`ls  tmp/YH$chrom/data/*.fasta`)
    cat $ff | gawk '/^>/{s=substr($1,2);printf("Sequence %s\nForward\nIs_read\nComposite 10\ncDNA_clone %s\n\n",s,s);}' > $ff.ace
    echo "pparse   $ff.ace" >> tmp/YH$chrom/data/_r
  end
  
  echo save >> tmp/YH$chrom/data/_r
  echo quit >> tmp/YH$chrom/data/_r

  bin/tacembly tmp/YH$chrom <  tmp/YH$chrom/data/_r


  bin/tacembly tmp/YH$chrom << EOF
    query find sequence Is_read && ! IntMap
    edit -D Is_read
    acem
      cdna_1
      quit
    save
    find model
    list -a -f tmp/YH$chrom/ff1.done
    quit
EOF




