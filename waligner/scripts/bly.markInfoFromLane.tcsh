#!/bin/tcsh -f

set run=$1
set lane=$2

pushd $run.$lane
    zcat ../$run/f2.$lane.fastc.gz | gawk '/^>/{gsub("n.",r,$1);}{print ;}' r=$run. | gzip > est.dna.gz
    tacembly . <<EOF
y
      pparse ../acedata/methods.ace
      parse est.dna
      find sequence
      edit Ref_Mrna
      edit Is_read
      acem
        cdna_80
        quit
      save
      list -a -f est.list
      quit
EOF
  popd
end


