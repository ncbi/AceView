/^#/ {next ;}
{ if (NN<1)NN=100 ;
  s = $2 ; # sample
  split (s,aa,"_") ;
  died = substr(aa[1], length(aa[1])) ;
  duration = aa[5] ;
  duration = int (duration/NN) ;
  if (died == 1) 
    {
      dd[nn] = duration ;
      if (duration > dMax) dMax = duration ;
    }
  else dd[nn] = -1 ;

  nn++ ; # total number of patients
}
END {
  if (JJ<1) JJ = 3 ;
  n1 = int (nn/JJ) ; n2 = nn - JJ * n1 ;
  # count all samples with duration at least d
  for (d = 0 ; d <= dMax ; d++)
    for (n = 1 ; n<= nn ; n++)
      {
	j = int ((n - n2)/n1) ; # split the patients in 3 classes 
	if (j < 0) j == 0 ;
	if (dd[n] == -1 || dd[n] >= d) 
	  { njd[j,d]++ ; nd[d]++ ; }
      }
  # export
  printf ("# Kaplan Meier survival curve %s :\n", title) ;
  printf ("# Number of surviving patient after d days, spliting the patients in %d quantiles\n", JJ) ;
  printf ("# Survival time in days") ;
  for (j = 1 ; j <= JJ ; j++)
    printf ("\tQuantile %d", j) ;
  printf ("\tAll patients") ;
  for (d = 0 ; d <= dMax ; d++)
    {
      printf ("\n%d", d * NN) ;
      for (j = 0 ; j < JJ ; j++)
	printf ("\t%d", njd[j,d]) ;
      printf ("\t%d", nd[d]) ;
    }
   printf("\n");
}
