#!bin/tcsh -f

set chrom=$1

bin/tacembly tmp/XH$chrom << EOF
  read-models
  query find gene capture == A2
  edit Colour YELLOW
  query find gene capture == R3
  edit Colour RED3
  query capture == A2
  edit Colour ORANGE
  query find est XI_*
  query Composite > 1
  edit COLOUR BLUE2
  query Composite > 2
  edit COLOUR BLUE3
  query Composite > 5
  edit COLOUR BLUE4
  query Composite > 10
  edit COLOUR BLUE5
  query Composite > 20
  edit COLOUR BLUE6
  query Composite > 50
  edit COLOUR BLUE7
  query Composite > 100
  edit COLOUR BLUE8
  save
  quit
EOF

exit 0




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

