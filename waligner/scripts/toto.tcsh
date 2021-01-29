#!/bin/tcsh -f
set zone=$1
tbly tmp/SNP_DB/$zone << EOF
  read-models
  // parse tmp/METADATA/RefSeq.MRNA.splicing.ace
  save
  quit
EOF
bin/snp -db_translate -db tmp/SNP_DB/$zone
      bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC -Reference_genome  "$Reference_genome" 
      bin/snp -i tmp/SNPH/$zone/$MAGIC.snp.sorted -db tmp/SNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -differential -o  `pwd`/tmp/SNP_DB/$zone/$MAGIC.differential -Reference_genome  "$Reference_genome "

