#!bin/tcsh -f

set chrom=$1
bin/tacembly tmp/XH$chrom << EOF
  query find read XY* ; ct_ac || other ; ! gt_ag
  edit -D Is_read
  follow from_gene
  acembly
    cdna_73 -clean_killed_mRNA
    quit
  query find tg (other || ct_ac) && ! gt_ag
  follow read
  query ct_ac || other
  edit -D Is_read
  follow from_gene
  acembly
    cdna_73 -clean_killed_mRNA
    quit
  save
  quit
EOF

