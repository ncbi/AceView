
{
  gsub (/\"/,"",$0) ;
  if (NF >= 9 && index($9,"xon") < 1) next ;
  m = $1 ; gene=$11 ;
  chrom = $2 ; a1 = $3 ; a2 = $4 ; x1 = $5; x2 = $6 ;
  if (substr(chrom,1,3)=="chr")chrom = substr (chrom,4) ;

  if (cm[m] && cm[m] != chrom)
    {
      printf ("\n// ERROR mRNA %s on chrom %s already on %s\n\n", m, chrom, cm[m]) ;
      next ;
    }
  mm[m] = m ; cm[m] = chrom ; if (gene) gm[m] = gene ; 

  if (am1[m] < 1) am1[m] = a1 ;
  am2[m] = a2 ;

  nmx[m]++ ; amx1[m,nmx[m]] = x1 ; amx2[m,nmx[m]] = x2 ;
}
END {
  for (m in mm)
    {
      m1 = am1[m] ; m2 = am2[m] ;
      s =  mstrand[m] ;
      chrom = cm[m] ;
      
      printf("Sequence %s\n", chrom) ;
      printf ("Subsequence %s__%s %d %d\n\n", lab,m,m1,m2) ;
      
      printf("Sequence %s__%s\n", lab,m) ;
      printf("Method %s\n", lab) ;
      printf ("Is_predicted_gene\n -D Source_exons\n") ;
      printf("Source %s\n", chrom) ;
      printf("IntMap %s %d %d\n", chrom, m1, m2) ;
      if (gm[m]) printf("Model_of X__%s__%s\n",lab,gm[m]) ;
      for (x = 1 ; x<= nmx[m] ; x++)
	printf("Source_exons %d %d // x=%d\n", amx1[m,x], amx2[m,x], x) ; 
      printf ("\n") ;
    }
}
