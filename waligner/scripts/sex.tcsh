#!bin/tcsh -f

goto dosage

set toto=chrom_ratio.txt
date >! $toto

cat << EOF > $toto.awk
{ r = \$1 ; c = \$2 ; n = 0 + \$3 ; a = \$4 ; b = \$5 ;
  rr[r] = 1 ; cc[c] = 1 ; 
  if (! c2i[c]){cMax++; c2i[c] = cMax ; i2c[cMax] = c;}
  rca[n,r,c] = a ; rcb[n,r,c] = b ; 
  ca[n,c] += a ; cb[n,c] += b ; 
  if (c != "X" && c != "Y") { ra[n,r] += a ; rb[n,r] += b ; tta[n] += a ; ttb[n] += b ;  }
}
END {
   printf ("\tRun") ;
   for (n = 5 ; n <= 20 ; n = 2*n)
     {
	printf ("\t") ;
        for (i=1 ; i<= cMax ; i++)
	  { c = i2c[i] ; printf("\t%s(%d)",c,n) ; }
     }
   for (r in rr)
     {
       printf("\n%s\t%s",substr(r,4),r) ; 
       for (n = 5 ; n <= 20 ; n = 2*n)
         {
	   printf("\t%d", 0) ;
           for (i=1 ; i<= cMax ; i++)
	     { 
               c = i2c[i] ;
               if (method == 1)
                 printf("\t%.2f", 100 * rca[n,r,c] * tta[n] / (1 + ra[n,r] * ca[n,c])) ; 
               if (method == 2)
                 printf("\t%.2f", 100 * rcb[n,r,c] * ttb[n] / (1 + rb[n,r] * cb[n,c])) ; 
             }
         }
     
     }
   printf ("\n") ; 
}
EOF

  echo "Relative coverage per run of different chromosomes at threshold 5" >> $toto
    cat toto.sex  | gawk -F '\t' -f $toto.awk method=1  > $toto.1
    head -1 $toto.1 >> $toto
    tail -n +2 $toto.1 | sort -k 1,1n >> $toto

  echo >> $toto
  echo >> $toto

  echo "Relative length covered per run of different chromosomes at threshold 5" >> $toto
    cat toto.sex | gawk -F '\t' -f $toto.awk method=2   > $toto.1
    head -1 $toto.1 >> $toto
    tail -n +2 $toto.1 | sort -k 1,1n >> $toto

  echo >> $toto
  echo >> $toto


echo $toto
\cp $toto RESULTS
goto phaseLoop 



#############

dosage:

# plot the dosage per run 
  cat tmp/SPONGE/EXOME/exomeDosage.txt | bin/histo -plot -o tmp/SPONGE/EXOME/geneDosagePlot -columns 4-12 -plain

# transpose the matrix
  cat tmp/SPONGE/EXOME/exomeDosage.txt | gawk -f scripts/transpose.awk > tmp/SPONGE/EXOME/exomeDosageT.txt

  set n=`grep -n Chromosome tmp/SPONGE/EXOME/exomeDosage.txt | gawk '{n=$1}END{i=index(n,":");print substr(n,1,i-1)}'`
tail -n +$n geneDosage.txt | bin/histo -plot -o geneDosagePlot -columns 4-12 -plain

cat  geneDosage.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<2)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > geneDosageR.txt

cat  exomeDosage.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<1)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > exomeDosageR.txt

# PCA and heat map, construction of the file exomeDosageLegend.txt
cat  tmp/SPONGE/EXOME/exomeDosage.index.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<1)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > tmp/SPONGE/EXOME/exomeDosageR.index.txt

R <<EOF
genes = read.table("tmp/SPONGE/EXOME/exomeDosageR.index.txt",sep='\t')
gg=as.matrix(genes)
d2=dim(gg)
gg2=gg[3:d2[1],6:d2[2]]
rownames(gg2)=gg[3:d2[1],1]
colnames(gg2)=gg[1,6:d2[2]]
cc=cor(gg2)

colnames(gg2)[hh$rowInd]
pdf ("exomeDosage.correl.pdf")
hh=heatmap(cc)

write(colnames(gg2)[hh$rowInd],"exomeDosageLegend.txt")
quit()
EOF

\cp exomeDosage.correl.pdf exomeDosageLegend.txt ~/ACEVIEWHELP


######

cat TARGET/GENES/RefSeq.gene_model.txt | gawk -F '\t' '{g=$2;chrom=$6;pos=($7+$8)/2;if(gg[g])next;gg[g]=pos;printf("%s\t%d\t%s\n",chrom,pos,g);}' | sort -k 1,1 -k2,2n |  gzip > geneDosage.intmap.gz

gunzip -c geneDosage.intmap.gz tmp/GENERUNS/RefSeq.u.geneSupport.ace.gz | bin/geneelements -geneDosage -o geneDosage

# plot the dosage per run 
cat geneDosage.txt | bin/histo -plot -o geneDosagePlot -columns 4-12 -plain

dosage2:
# plot the dosage per gene
cat geneDosage.txt | gawk -f scripts/transpose.awk > geneDosageT.txt

set n=`grep -n Chromosome geneDosage.txt | gawk '{n=$1}END{i=index(n,":");print substr(n,1,i-1)}'`
tail -n +$n geneDosage.txt | bin/histo -plot -o geneDosagePlot -columns 4-12 -plain

cat  geneDosage.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<2)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > geneDosageR.txt

cat  exomeDosage.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<1)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > exomeDosageR.txt

cat  exomeDosage.index.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<1)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > exomeDosageR.index.txt


R <<EOF
genes = read.table("tmp/SPONGE/EXOME/exomeDosageR.index.txt",sep='\t')
gg=as.matrix(genes)
d2=dim(gg)
gg2=gg[3:d2[1],6:d2[2]]
rownames(gg2)=gg[3:d2[1],1]
colnames(gg2)=gg[1,6:d2[2]]
cc=cor(gg2)

colnames(gg2)[hh$rowInd]
pdf ("exomeDosage.correl.pdf")
hh=heatmap(cc)

write(colnames(gg2)[hh$rowInd],"exomeDosageLegend.txt")
quit()

EOF


cat variant2strand.txt | gawk '/^#/{next}{printf("%d\t%d\n",$3,$7);}' | bin/histo -smooth -plot -o anisotropy

###############
## select a list 0f 73210 indicative variants with high coverage > 200 and intermediate frequency in the pool: 60 to 140
cat variant2strand.txt | gawk '/^#/{next}{if($7>200 && $6 > 60 && $6 < 140)printf("Variant %s\n",$1);}' > middleVariant.list

## for that list export the mutant/cover counts for the sum of the 2 strands
$tacembly tmp/VariantDB << EOF
  key middleVariant.list
  table -active -o middleVariant.anisotropy.txt -f tables/middleVariant.anisotropy.def
EOF

# mimic the format of the  tmp/SPONGE/EXOME/*/Ghs?.u.$limit.txt files before using 'geneelements -exomeDosage'
# we reformat to just export the 2 strands snp and cover

cat  middleVariant.anisotropy.txt | gawk -F '\t' '/^#/{next}{gsub(/\"/,"",$0);printf("%s\t%s\t%d\t%d\n",$1,$2,$8,$9);}' >  middleVariant.frequency.txt

# export a variantDosage table per Run/Variant, in the Run order given by the exomeDosage dendrogram
cd ..
cat  exomeDosageLegend.txt  tmp/VariantDB/middleVariant.frequency.txt | bin/geneelements -variantDosage -o variantDosage

cat  variantDosage.txt | gawk -F '\t' '/^Chromosome/{n++;}{if(n<1)next;n1=3;}/^Chromosome/{n1=4;}/Number of genes/{next;}{if(NF<3)next;}{if(n1==3){printf("%s",$3);}for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}' > variantDosageR.txt





phaseLoop:
echo done

