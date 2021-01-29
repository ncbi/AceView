#!/bin/tcsh -f

# This is  a very specific script written nov 28 2019 for the project Silver
# This is an UMI project badly annotated in SRA
# each read is UMI + insert
# in SRA it was declared as 2 sublibs of the same SRX

# export a table with the file names

tace MetaDB <<  EOF
  query find run file AND project == $MAGIC
  follow sublibrary_of
  bql -o $MAGIC.umi.files.txt r,f1,f2 from r in @,s1 in  r->sublibraries,s2 in r->sublibraries where s1<s2, f1 in s1->file[2], f2 in s2->file[2]
  quit
EOF

foreach run (`cat  $MAGIC.umi.files.txt | cut -f 1`)
  echo $run
  set f1=`cat  $MAGIC.umi.files.txt | gawk '{gsub(".gz","",$2);if($1==run)print $2}' run=$run` 
  set f2=`cat  $MAGIC.umi.files.txt | gawk '{gsub(".gz","",$3);if($1==run)print $3}' run=$run` 
  gunzip $f1 $f2 
  mkdir Fastc/$run
  cat $f1 | gawk '{nn++;z1=$1;n=getline <f2;z2=$1; if(substr(z1,1,1)==">"){if(0&&z1 != z2){print "error line " nn; exit(1);}next;}if(length(z1)>length(z2)){z0=z1;z1=z2;z2=z0;}printf("%s\t%s\n",z1,z2);}' f2=$f2 | sort > Fastc/$run/all.raw2
end
foreach run (`cat  $MAGIC.umi.files.txt | cut -f 1`)
   cat Fastc/$run/all.raw2 | sort -u | cut -f 2 | sort | gawk '{n++; printf (">n.%d\n%s\n",n,$1);}' | dna2dna -I fasta -O fastc -gzo -splitMb 200 -o Fastc/$run/f2
end
# associate the SRX to its longer SRR
\rm x2r.txt
foreach run (`cat  $MAGIC.umi.files.txt | cut -f 1`)
  set f1=`cat  $MAGIC.umi.files.txt | gawk '{gsub(".gz","",$2);if($1==run)print $2}' run=$run` 
  set f2=`cat  $MAGIC.umi.files.txt | gawk '{gsub(".gz","",$3);if($1==run)print $3}' run=$run` 
  cat $f1 | head -2 | gawk '{nn++;z1=$1;n=getline <f2;z2=$1; if(substr(z1,1,1)==">")next;if(length(z1)<length(z2)){f0=f1;f1=f2;f2=f0;}printf("%s\t%s\n",run,f1);}' f1=$f1 f2=$f2 run=$run | gawk '/DATA/{split($2,aa,"/");split(aa[4],bb,".");printf("%s\t%s\n",$1,bb[1]);}' >> x2r.txt
end
# export the interesting SRR .ace file
cat x2r.txt | gawk '{printf("Run %s\n",$2);}' > x2r.list
tace MetaDB << EOF
  key x2r.list
  show -a -f x2r.preace
EOF
# rename and remove the Sublib tag
cat x2r.txt ZZZZZ x2r.preace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){r2x[$2]=$1;next;}}/^Sublib/{next;}/^Ali/{next;}/^Run/{gsub(/\"/,"",$2);printf("Run %s\n",r2x[$2]);next;}{print;}' > x2r.ace

# save the previous database, reconstruct one with the SRX runs pointing to the UMI rationalized SRX Fastc files
foreach  ff (`ls Fastc/SRX*/f2.*.fasta`)
  set f2=`echo $ff | sed -e 's/\.fasta//'`
  mv $f2.fasta $f2.fastc
end
# \rm Fastc/SRX*/*.raw2.gz