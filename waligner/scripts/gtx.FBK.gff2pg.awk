
{ if ($3 != "exon") next ;
  chrom = $1 ; a1 = $4 ; a2 = $5 ;

  gsub (/\"/,"",$9) ;
  i=index($9,"gene_id") ;
  split (substr($9,i), aa, ";") ; 
  split (aa[1],bb, " ") ;
  gene = bb[2] ;

  i=index($9,"transcript_id") ;
  split (substr($9,i), aa, ";") ; 
  split (aa[1],bb, " ") ;
  m = bb[2] ;
  if (cm[m] && cm[m] != chrom)
    {
      printf ("\n// ERROR mRNA %s on chrom %s already on %s\n\n", m, chrom, cm[m]) ;
      next ;
    }
  mm[m] = m ; gm[m] = gene ; cm[m] = chrom ;

  if (ag1[gene] < 1 || ag1[gene] > a1) ag1[gene] = a1 ;
  if (ag2[gene] < 1 || ag2[gene] < a2) ag2[gene] = a2 ;

  if (am1[m] < 1 || am1[m] > a1) am1[m] = a1 ;
  if (am2[m] < 1 || am2[m] < a2) am2[m] = a2 ;

  nmx[m]++ ; amx1[m,nmx[m]] = a1 ; amx2[m,nmx[m]] = a2 ;
  mstrand[m] = $7 ; 
}
END {
  for (m in mm)
    {
      m1 = am1[m] ; m2 = am2[m] ;
      s =  mstrand[m] ;
      if (s == "-") { m0 = m1 ; m1 = m2 ; m2 = m0 ; }
      chrom = cm[m] ;
      
      printf("Sequence %s\n", chrom) ;
      printf ("Subsequence %s__%s %d %d\n\n", lab,m,m1,m2) ;
      
      printf("Sequence %s__%s\n", lab,m) ;
      printf("Method %s\n", lab) ;
      printf ("Is_predicted_gene\n -D Source_exons\n") ;
      printf("Source %s\n", chrom) ;
      printf("IntMap %s %d %d\n", chrom, m1, m2) ;
      printf ("Model_of X__%s__%s\n", lab, gm[m]) ;
      
      if (s == "+")
	for (x = 1 ; x<= nmx[m] ; x++)
	  printf("Source_exons %d %d // x=%d\n", amx1[m,x]- m1 + 1, amx2[m,x]- m1 + 1, x) ; 
      if (s == "-")
	for (x = nmx[m] ; x >= 1 ; x--)
	  printf("Source_exons %d %d // x=%d\n", m1 - amx2[m,x] + 1, m1 - amx1[m,x] + 1, x) ;
      printf ("\n") ;
    }
}
