BEGIN { gene = "-" ; type = "-" ; gid = "-" ; gf = "-" ; refseq = 0 ; est = 0 ; }
/^#/{next;}
/^ZZZZZ/ {
  zz++ ;
  if (zz == 2)
    {
      if (uu = "u") printf ("# Introns and their unique support per group and run.") ;
      else printf ("# Introns and their support per group and run.") ;
      if (uu = "u") printf ("\n# Only uniquely mapped reads contribute to this table ") ;
      printf ("\n# On the indicated chromosome, a1 and a2 are the coordinates of the first and last base of the intron ") ;
      printf ("in the reference genome %s", RG) ;
      printf ("\n# For example in a gt_ag intron, a1 is the coordinate of the g of gt_, a2 of the g of _ag.") ;
      printf ("\n# If a1 < a2, the intron is on the + strand, if a1 > a2 on the - strand.") ;
      printf ("\n# All coordinates are 1-based, i.e. the first base of the chromosome is numbered 1.") ;
      printf ("\n# For each group and each run, we count the number of reads with best match ") ;
      printf ("across the exon-junction, covering at least 8 nucleotides of each bordering exon.") ;
      printf ("\n# Considering that MAGIC clips the RNA alignments ") ;
      printf ("right before mismatches closer than 8 nucleotides from the edge of the alignment, all reads supporting an intron ") ;
      printf ("are either locally exact or extend further than 8 nucleotides inside the bordering exons.") ;
      printf ("\n# The RefSeq and EST columns show the number of Genbank/dbEST supporting accessions.") ; 
      printf ("\n# The RNA_seq support column shows the number of reads supporting the intron in all RNA-seq analyzed so far in MAGIC.") ; 

      printf ("\n# \n# Chromosome\ta1\ta2\tStrand\tType\tIntron Length\tGene\tRefSeq\tEST\tDeep RNA-seq\tNumber of supports in all runs in project") ;
      printf ("\tNumber of runs with at least 1 support") ;
      printf ("\tNumber of runs with at least 2 supports") ;
      printf ("\tNumber of runs with at least 5 supports") ;
      printf ("\tNumber of runs with at least 10 supports") ;
      printf ("\tNumber of runs with at least 50 supports") ;
      printf ("\tNumber of runs with at least 100 supports") ;
      for (ir = 1 ; ir <= nr ; ir++)
	printf ("\t%s", i2r[ir]) ;
      printf ("\tRefSeq model") ;
      printf ("\n") ;
    }
  next ;
}
{
  gsub (/\"/,"", $0) ;
}
{
  if (zz < 2)
    {
      nr++ ; i2r[nr] = $1 ; r2i[$1]= nr ; nn[ir] = 0 ;
      if (zz == 1) isRun[nr] = 1 ;
      next ;
    }
}
/^Intron /{
  if (a1 > 0 && type != "-")
    {
      nt = 0 ; n1 = 0 ; n2 = 0 ; n5 = 0 ; n10 = 0 ; n50 = 0 ; n100 = 0 ;
      for (ir = 1 ; ir <= nr ; ir++)
	if (isRun[ir] == 1)
	  {
	    nt += nn[ir] ;
	    if (nn[ir] >=   1)   n1++ ;
	    if (nn[ir] >=   2)   n2++ ;
	    if (nn[ir] >=   5)   n5++ ;
	    if (nn[ir] >=  10)  n10++ ;
	    if (nn[ir] >=  50)  n50++ ;
	    if (nn[ir] >= 100) n100++ ;
	  }
      if (nt > 0)
	{
	  printf ("%s\t%d\t%d", chrom, a1, a2) ;
	  da = a2 - a1 ; strand = "+" ; if (da < 0) { strand = "-" ; da = -da ;} da++ ;
	  printf ("\t%s\t%s", strand, type) ;
	  printf ("\t%d\t%s", da, gene) ;
	  printf ("\t%d\t%d\t%d", refseq, est, deep) ; 
	  printf ("\t%d\t%d\t%d\t%d\t%d\t%d\t%d", nt, n1, n2, n5, n10, n50, n100) ;
	  for (ir = 1 ; ir <= nr ; ir++)
	    { printf ("\t%d", nn[ir]) ; }
	  printf ("\t%s", gf) ; 
	  printf ("\n") ;
	}
    }
  a1 = 0 ; a2 = 0 ; type = "-" ; gene = "-" ; gf = "-" ;
  deep = 0 ; refseq = 0 ; est = 0 ; 
  for (ir = 1 ; ir <= nr ; ir++)
    { nn[ir] = 0 ; }
  id = $2 ;
  next ;
}
/^IntMap/ { chrom = $2 ; a1 = $3 ; a2 = $4 ; next ; }
/^Group_U/ {
  r = $2 ; ir = r2i[r] ; 
  if (ir) { nn[ir] += $6 ;  }
  next ;
}
/^Run_U/ {
  r = $2 ; ir = r2i[r] ; 
  if (ir) { dn = $6 - nn[ir] ; if (dn > 0) nn[ir] += dn ; }
  next ;
}
/^Run_nU/ {
  r = $2 ; ir = r2i[r] ; 
  if (ir) { dn = $6 - nn[ir] ; if (dn > 0) nn[ir] += dn ; }
  next ;
}
/^Validated_u/ {
  r = $2 ; ir = r2i[r] ; 
  if (ir) { dn = $3 - nn[ir] ; if (dn > 0) nn[ir] += dn ; }
  next ;
}
/^Other/ { type = $2 ; next ; }
/^g._ag/ { type = $1 ; next ; }
/^.t_ag/ { type = $1 ; next ; }
/^GeneId/ { if (gid != "-") gid = gid "," $2 ; else gid = $2 ; next ; }
/^Gene/ { gene = $2 ; next ; }
/^From_gene/ { if (! gene) { g= $2; i=index(g,"Aug10");if(i>0){j=index(g,".");if(j>i-5)g=substr(g,1,j-1);} gene =g ;} next ; }
/^EST/{ est = $2 ; next ; }
/^RNA_seq/{ deep = $2 ; next ; }
/^RefSeq/{ refseq = $2 ; next ; }
/^From_genefinder/ { if (gf != "-") gf = gf "," $2 ; else gf = $2 ; next ; }
