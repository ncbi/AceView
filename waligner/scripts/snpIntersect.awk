/^#/ { next }
{
  if (irMax == 0)
    {
      run = "Union" ; irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; 
    }
  if (format == "snp")
    {
      chrom = $1 ; pos = $2;
      run = $3 ;
      z = chrom ":" pos ":" a ":" b
    }
  else if (format == "countEdited")
    {
      chrom = $1 ; pos = $2;       chrom2 = $3 ; pos2 = $4; run = $8 ;
      z = chrom ":" pos ":" chrom2 ":" pos2 ":" $5 ":" $6 ":" $7 ;
    }
  else
    {
      print "Bad format parameter in snpSupportStat.awk format =#" format "#" ;
      exit (1) ;
    } 
  ir = r2ir[run] ; if (ir<1) { irMax++ ; ir = irMax ; r2ir[run] = ir ; ir2r[ir] = run ; }

  c[z,ir]++ ; 
  anyC[ir]++ ;

  c[z,1]++ ; anyz[z]++ ;
  anyC[1]++ ;
}
END  {
  nMax = 0 ;
  for (z in anyz)
    {   # histo of how many sample see a given variant */
      n = 0 ; for (ir = 1 ; ir <= irMax ; ir++)
	if (c[z,ir] > 0) n++ ;
      nR[n]++ ; if (n> nMax) nMax = n ;
    }
  printf ("Number of samples") ;
  for (n = 1 ; n <= nMax ; n++) printf ("\t%d", n) ;
  printf ("\nNumber of variants seen in n samples") ;
  for (n = 1 ; n <= nMax ; n++) printf ("\t%d", nR[n]) ;

  for (ir = 1 ; ir <= irMax ; ir++)
    {
      for (jr = 1 ; jr <= irMax ; jr++)
	{
	  for (z in anyz)
	    {
	      if (c[z,ir] > 0 && c[z,jr] > 0) intersect[ir, jr]++ ;
	      if (c[z,ir] > 0 || c[z,jr] > 0) union[ir, jr]++ ;
	    }
	}
    }
   printf ("\n\nThe number of variants seen in both runs (i.e. the cardinal of intersection) is presented over the diagonal\n") ;
   printf ("The number of variants seen in either run (i.e. the cardinal of union) is presented below the diagonal\n") ;
   printf ("The number of variants in a given run is given on the diagonal, it is also the self-union or self-intersectionof the run\n") ;

   printf ("Run") ;
   for (jr = 1 ; jr <= irMax ; jr++) printf ("\t%s", ir2r[jr]) ;
   for (ir = 1 ; ir <= irMax ; ir++)
     {
       printf ("\n%s", ir2r[ir]) ;
       for (jr = 1 ; jr <= irMax ; jr++)
	 {
	   if (ir < jr) printf ("\t%s", intersect[ir,jr]) ;
	   if (ir >= jr) printf ("\t%s", union[ir,jr]) ;
	 }
     }
   printf("\n") ;
}
