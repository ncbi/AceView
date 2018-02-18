#!bin/tcsh -f

set ff=$1

# this code can be used to fix an old  MetaDB dumpdir and parse it in the new schema
# gsub("","",\$0) ;

cat << EOF >  modelFix.awk 
{ gsub("Average_fragment_length","Fragment_length_average",\$0) ;
  gsub("Links_gene_genome","Compatible_gene_extension",\$0) ;
  gsub("Covariance_analysis","Autosort",\$0) ;
  print ;
}
EOF


echo "cat $ff | gawk -f modelFix.awk > $ff.modelFixed"
cat $ff | gawk -f modelFix.awk > $ff.modelFixed

