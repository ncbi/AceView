{  if ($1 == "Autosome") c = $1 ; else  c ="chr" $1 ;
if (1)
     {
       g =  $2 ;
       ig = g2i[g] ;
       if (0 + ig < 1)
	 {
	   ng++ ;
	   i2g[ng] = g ;
	   g2i[g] = ng ;
	   ig = ng ;
	 }
       x = $3 ;
       a = $4 ;
       b = $5 ;
       ic = c2i[c] ;
       if (0 + ic < 1)
	 {
	   nc++ ;
	   i2c[nc] = c ;
	   c2i[c] = nc ;
	   ic = nc ; 
	 }
       ac[ic,ig,x] = a ;
       bc[ic,ig,x] = b ;
       aa[x] += a ;
     }
}
END {
  for (ii = 0 ;ii < 2 ; ii++)
    {
      if (ii == 0) 
	printf ("\nCorrelation between strands, the peak indicates the length of the insert minus the length of a read\n") ;
      else
	printf ("\nAutocorrelation\n") ;

      for (ic = 1 ; ic <= nc ; ic++)
	{
	  printf("\tDistance") ;
	  for (ig = 1 ; ig <= ng ; ig++)
	    {
	      printf("\t%s:%s", i2c[ic], i2g[ig]) ;
	    }
	  printf("\t\t") ;     
	}
      
      for(x = 0 ; x <= 500 ; x += 10)
	{
	  printf("\n") ;
	  
	  for (ic = 1 ; ic <= nc ; ic++)
	    {
	      printf("\t%d",x) ;
	      for (ig = 1 ; ig <= ng ; ig++)
		{
		  if (ii == 0)
		    z = bc[ic,ig,x] ;
		  else
		    z = ac[ic,ig,x] ;
		  printf("\t%.2f", z) ;
		}
	      printf("\t\t") ;
	    }
	}
      printf ("\n") ;
    }
}
