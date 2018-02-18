/^ZZZZZ/{zz = 1; next; }
{ if (zz < 1) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r] = 1 ; next; }}
/^\"/{ 
  if (zz == 1)
    {
      gsub("\"","",$0) ; gsub("NULL","",$0) ;
      split($0,aa,"\t") ;
      z=aa[1];r = r2i[z] ; if (r < 1) next ;
      if (length (aa[2]) > 0) { nMachine++ ; machine[r]=aa[2]; }
      if (length (aa[3]) > 0) { nSample++ ; sample[r]=aa[3];}
      if (length (aa[4]) > 0) { nSystm++ ; systm[r]=aa[4]; }
      if (length (aa[5]) > 0) { nTissue++ ; tissue[r]=aa[5]; }
      if (length (aa[6]) > 0) { nTitre++ ; titre[r]=aa[6]; }
      if (length (aa[7]) > 0) { nSystm++ ; systm[r] = systm[r] " " aa[7]; }
      if (length (aa[8]) > 0) { nRunid++ ; runid[r]=aa[8]; }
      if (length (aa[11]) > 0) { nsTitre++ ; stitre[r]=aa[11]; }
      if (length (aa[12]) > 0) { nsTitre2++ ; stitre2[r]=aa[12]; }
      if (length (aa[13]) > 0) { noTitre++ ; otitre[r]=aa[13]; }
    }
  next ;
}
/^Insert_length_in_pairs/ { 
  z=$2 ; r = r2i[z] ;  if(r<1) next ;
  if(good[r]==1) ok[r] = 1 ;
  for (i = 3 + 4 ; i <= NF  ; i++)
    {
      # smooth using Pascal coeff 
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
    print sample[1], 10+nSample ;
  MX = 1000
  nr2 = 0 ;
  for (r = 1; r <= nr ; r++)  if (ok[r]==1)nr2++ ;
  printf("## Insert_length_in_pairs in %d/%d runs", nr2,nr) ;
  if (nr2 == 0) nr2 = 1 ;
  if (nRunid >= 0) { printf("\n## RunId"); for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", runid[r]); }
  if (nMachine >= 0) { printf("\n## Machine"); for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", machine[r]); }
  if (nSample >= 0) { printf("\n## Sample"); for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", sample[r]); }
  if (nTissue >= 0) { printf("\n## Tissue"); for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", tissue[r]); }
  if (nSystm >= 0) { printf("\n## System"); for (r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", systm[r]); }
  if (nTitre >= 0) { printf("\n## Title");  for( r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", titre[r]); }
  if (nsTitre >= 0) { printf("\n## Sorting_Title");  for( r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", stitre[r]); }
  if (nsTitre2 >= 0) { printf("\n## Sorting_Title_2");  for( r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", stitre2[r]); }
  if (noTitre >= 0) { printf("\n## Other_Title");  for( r = 1; r <= nr ; r++)  if (ok[r]) printf ("\t%s", otitre[r]); }
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
