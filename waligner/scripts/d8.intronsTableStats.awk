/^#/{next;}
{  n1 = $12 ; 
   if (0+n1 < 1) next ; 
   if (n1 > nMax) 
     nMax = n1 ; 
   k = 0 ;
   for (i = 17 ; i >= 12 ; i--) {
     x = $i - k ; k += $i ; 
     nn[n1,i] += $i/n1 ; # 1 2 5 10 50 100 supports in that many runs
   }
}
END {
  printf ("Number of persons") ;
  for (n1 = 1 ; n1 <= nMax ; n1++)
    {
      printf ("\n%d",n1) ;
      for (i = 12 ; i<= 17 ; i++) 
	printf ("%d", nn[n1,i]) ;
    }
  printf ("\n") ;
}
  
