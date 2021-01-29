#!bin/tcsh
set type=$1
set m1=$2
set m2=$3
set byKb=$4

  set totov=tmp/VIRUS/$MAGIC
  set totob=tmp/BACTERIA/$MAGIC
  set totot=tmp/BACTERIA/$MAGIC.tr


  set MG=MetaDB/$MAGIC/Mb_in_genes_with_GeneId_minus_high_genes.txt
  set toto2c=tmp/vir2.$type.$byKb

  echo -n "## $toto2c " > $toto2c.v
  cat $totov.kb_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c.v
  echo >> $toto2c.v
  cat $toto2c.v | gawk '/#Run/{print;last;}' > $totov.999
  cat $totov.kb_aligned_per_gene.txt ZZZZZ $totov.reads_aligned_per_gene.txt | gawk -F '\t' -f scripts/vir2.k1.awk m1=$m1 m2=$m2 byKb=$byKb | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ $MG ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus byKb=$byKb | scripts/tab_sort -k 8,8nr >> $toto2c.v
  echo "\n\n" >> $toto2c.v

  cat $totob.kb_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  > $toto2c.b
  echo >> $toto2c.b
  cat $toto2c.b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totob.kb_aligned_per_gene.txt ZZZZZ $totob.reads_aligned_per_gene.txt | gawk -F '\t' -f scripts/vir2.k1.awk m1=$m1 m2=$m2 byKb=$byKb | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ $MG ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes byKb=$byKb | scripts/tab_sort -k 8,8nr >> $toto2c.b
  echo "\n\n" >> $toto2c.b

  cat $totot.kb_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  > $toto2c.t
  echo >> $toto2c.t
  cat $toto2c.t | gawk '/#Run/{print;last;}' > $totov.999
  cat $totot.kb_aligned_per_gene.txt ZZZZZ $totot.reads_aligned_per_gene.txt | gawk -F '\t' -f scripts/vir2.k1.awk m1=$m1 m2=$m2 byKb=$byKb | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ $MG ZZZZZ VirusDB/transposon.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Transposon byKb=$byKb | scripts/tab_sort -k 8,8nr >> $toto2c.t
  echo "\n\n" >> $toto2c.t

# new code: reconstrruct the Sp(ecies)_group on the fly inside the given quality theshold

# scan the filtered tables just established, then add up the values by remapping the species into species groups
cat $toto2c.v | head -60 | gawk -F '\t' '/^#/{print}' > $toto2c.vg
cat VirusDB/virus_groups.metadata.txt     ZZZZZ $toto2c.v | gawk -F '\t' -f scripts/vir2.k3.awk | scripts/tab_sort -k 8,8nr  >> $toto2c.vg 
cat $toto2c.b | head -60 | gawk -F '\t' '/^#/{print}' > $toto2c.bg
cat VirusDB/bacteria_groups.metadata.txt  ZZZZZ $toto2c.b | gawk -F '\t' -f scripts/vir2.k3.awk | scripts/tab_sort -k 8,8nr  >> $toto2c.bg 

# concatenate in the order we happen to prefer today 2020_06_25
echo "\n\n\n" > $toto2c.sp
cat $toto2c.vg $toto2c.sp $toto2c.v $toto2c.sp  $toto2c.bg $toto2c.sp $toto2c.b $toto2c.sp  $toto2c.t > $toto2c
# \rm $toto2c.vg $toto2c.v  $toto2c.bg $toto2c.b $toto2c.t

exit 0
# below aold code relying on Sp(ecies)_group beeing known to geneindex.c
  echo -n "## $toto2c " > $toto2c
  cat $totov.gene_group.kb_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  echo >> $toto2c
  cat $toto2c | gawk '/#Run/{print;last;}' > $totov.999
  cat $totov.gene_group.kb_aligned_per_gene.txt ZZZZZ $totov.gene_group.reads_aligned_per_gene.txt | gawk -F '\t' -f scripts/vir2.k1.awk m1=$m1 m2=$m2 byKb=$byKb | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ $MG ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus byKb=$byKb | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c


  cat $totob.gene_group.kb_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  echo >> $toto2c
  cat $toto2c | gawk '/#Run/{print;last;}' > $totov.999
  cat $totob.gene_group.kb_aligned_per_gene.txt ZZZZZ $totob.gene_group.reads_aligned_per_gene.txt | gawk -F '\t' -f scripts/vir2.k1.awk m1=$m1 m2=$m2 byKb=$byKb | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ $MG ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes byKb=$byKb | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c


