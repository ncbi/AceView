/^Insert_length_in_pairs/ { 
    r=$2 ;
  for (i = 3 + 4 ; i <= NF  ; i++)
    {
      # smooth using Pascal coeff 
      # do not do this, we lose the periodicity of the Illumina method
      z = $i ; 
      if (1)
	{
	  delta = 1 ; nn[r, i] = z ;  
	  if(i > imax) imax = i ;
	}
      else
	{
	  delta = 256 ; # (1+1)^8
	  if (z > 0)
	    {
	      nn[r,i-4] += z ; nn[r,i-3] += 8*z ; nn[r,i-2] += 28*z ; nn[r,i-1] += 56*z ; nn[r,i] += 70*z ;
	      nn[r,i+4] += z ; nn[r,i+3] += 8*z ; nn[r,i+2] += 28*z ; nn[r,i+1] += 56*z ;
	      
	      if(i+4 > imax) imax = i + 4 ;
	    }
	}
    }
}
END {
  MX = 1000
  nr2 = 0 ;
  for (r = 1; r <= nr ; r++)  if (ok[r]==1)nr2++ ;
  printf("## Insert_length_in_pairs in %d/%d runs", nr2,nr) ;
  if (nr2 == 0) nr2 = 1 ;

  printf ("\n# Run") ;  for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", i2r[r]);

  totMax = 0 ; printf ("\tAverage") ;
  if (delta < 1) delta = 1 ;

  for (i = 3 ; i <= imax ; i++)
    {
      tot = 0 ;
      for (r = 1 ; r <= nr ; r++)
	{
	  if (ok[r]) 
	    tot += nn[r,i]/delta ;
	}
      tot /= nr2 ;
      if (tot > totMax) totMax=tot ;
      if (i>MX && MX * tot < totMax)
	{ imax = i ; }
    }

  for (r = 1 ; r <= nr ; r++)
    {
      if (ok[r] < 1) continue ;
      mx = 0 ; t = 0 ; a = 0 ;
      for (i = 3 ; i <= imax ; i++)
	{
	  z = nn[r,i]/delta ;
	  t += z ;
	  if (z > mx)
	    { mx = z ; ix = i ; }
	  a += z * ( i - 3) ;
	}
      if (t < 1) t = 1 ;
      av[r] = a/t ;
      mxr[r] = mx ;
      ixr[r] = ix ;
      totr[r] = t ;
      t0 = 0 ; t1 = 0 ; t5 = 0 ; t95 = 0 ; t99 = 0 ; tm = 0 ;
      
      for (i = 3 ; i <= imax ; i++)
	{
	  z = nn[r,i]/delta ;
	  t0 += z ;
	  if (100 * t0 > t && t1 == 0)
	  { t1r[r] = i - 3 ; t1 = i - 3 ; }
	  if (100 * t0 > 5 * t && t5 == 0)
	    { t5r[r] = i - 3 ; t5 = i - 3  }
	  if (2 * t0 > t && tm == 0)
	    { median[r] = i - 3  ; tm = i - 3  ; }
	  if (100 * t0 > 95 * t && t95 == 0)
	  { t95r[r] = i - 3 ; t95 = i - 3 ; }
	  if (100 * t0 > 99 * t && t99 == 0)
	    { t99r[r] = i - 3 ; t99 = i - 3 }
	}
    }
  printf ("\nHigh 99%%") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", t99r[r] + .5 ) ;
  printf ("\nHigh 95%%") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", t95r[r] + .5 ) ;
  printf ("\nAverage") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", av[r] + .5 ) ;
  printf ("\nMedian") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", median[r] + .5 ) ;
  printf ("\nPosition of maximum") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", ixr[r] - 3  + .5 ) ;
  printf ("\nLow 5%%") ;
  for (r = 1 ; r <= nr ; r++)
    if (ok[r]) printf ("\t%d", t5r[r] + .5 ) ;
  printf ("\nLow 1%%") ;
  for (r = 1 ; r <= nr ; r++) 
    if (ok[r]) printf ("\t%d", t1r[r] + .5 ) ;


  printf ("\nRun") ; 
  for (r = 1; r <= nr ; r++) if (ok[r])  printf ("\t%s", i2r[r]);
  for (i = 3 ; i <= imax ; i++)
    {
      tot = 0 ;
      printf ("\n%d",i-3) ; 
      for (r = 1 ; r <= nr ; r++)
	if (ok[r])
	  {
	    printf("\t%d", nn[r,i]/delta) ;
	    tot +=nn[r,i]/delta ;
	  }
      tot /= nr2 ;
      if (tot > totMax) totMax=tot ;
      printf("\t%d", tot) ;
      if (0) printf("\t%d", totMax) ;
      
      if (i>MX + 2 && MX * tot < totMax)
	{ printf("\n"); exit ;}
    }
  printf ("\n") ;
}
