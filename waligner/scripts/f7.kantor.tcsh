#!bin/tcsh -f

set chrom=$1

tbly  tmp/XH$chrom << EOF
  query find kantor product
  list -a -f tmp/XH$chrom/f7.kantor.list
  quit
EOF

tbly Kantor << EOF
  key  tmp/XH$chrom/f7.kantor.list
  show -a -f  tmp/XH$chrom/f7.kantor.ace
  quit
EOF

tbly  tmp/XH$chrom << EOF
  pparse   tmp/XH$chrom/f7.kantor.ace
  save
  find mrna
  acembly
    cdna_kantor
    quit
  save
  quit
EOF

touch tmp/XH$chrom/f7.kantor.done

