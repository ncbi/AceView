
{ if ($3 != "exon") next ;
  chrom = $1 ; a1 = $4 ; a2 = $5 ; 

  i=index($9,"parent_type=") ;
  split (substr($9,i), aa, ";") ; 
  split (aa[1],bb, "=") ;
  ptype = bb[2] ;

  i=index($9,"Parent=") ;
  split (substr($9,i), aa, ";") ; 
  split (aa[1],bb, "=") ;
  nParents = split (bb[2],parents,",") ;
  for (i = 1 ; i<= nParents ; i++)
    {
      parent = parents[i] ;


  print aa[1] ; 
  next ; 
  i=index($9,"parent") ; 
  split (substr($9,i), aa, ";") ; split (aa[1], bb, ":") ;
  mrna = substr(bb[1],13) ;
  if (old && mrna != old) 
    {
      if (strand == "+")
	{
	  x0 = 1 ; dx = 0 ;
	  for (i = 1 ; i<= nx ; i++)
	    {
	      a1 = aa1[i] ; a2 = aa2[i] ;
	      dx += a2 - a1 + 1 ;
	      printf ("%s\t%d\t%d\t%s\t%d\t%d\n", old, x0, x0 + dx - 1, chrom, a1, a2) ;
	      x0 += dx ;
	    }
	}
      else
	{
	  x0 = 1 ; dx = 0 ;
	  for (i = nx ; i >= 1 ; i--)
	    {
	      a1 = aa1[i] ; a2 = aa2[i] ;
	      dx += a2 - a1 + 1 ;
	      printf ("%s\t%d\t%d\t%s\t%d\t%d\n", old, x0, x0 + dx - 1, chrom, a2, a1) ;
	      x0 += dx ;
	    }
	}
      x0 = $4 - 1 ; 
      nx = 0 ;  
    }
  old = mrna ;
  nx = nx+1 ;
  printf ("#### %s #%d#\t", mrna, nx) ; print ;
  aa1[nx] = a1 ;
  aa2[nx] = a2 ;
  strand = $7 ;
}

END    {
      if (strand == "+")
	{
	  x0 = 1 ; dx = 0 ;
	  for (i = 1 ; i<= nx ; i++)
	    {
	      a1 = aa1[i] ; a2 = aa2[i] ;
	      dx += a2 - a1 + 1 ;
	      printf ("%s\t%d\t%d\t%s\t%d\t%d\n", old, x0, x0 + dx - 1, chrom, a1, a2) ;
	      x0 += dx ;
	    }
	}
      else
	{
	  x0 = 1 ; dx = 0 ;
	  for (i = nx ; i >= 1 ; i--)
	    {
	      a1 = aa1[i] ; a2 = aa2[i] ;
	      dx += a2 - a1 + 1 ;
	      printf ("%s\t%d\t%d\t%s\t%d\t%d\n", old, x0, x0 + dx - 1, chrom, a2, a1) ;
	      x0 += dx ;
	    }
	}
      nx = 0 ;
    }
