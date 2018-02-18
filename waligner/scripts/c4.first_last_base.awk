/^ZZZZZ/{zz++ ; next ;}
{
  if (zz < 1)
    { 
      r = $1 ; 
      ir=r2ir[r]; if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}
      good[ir] = 1 ;
      next ;
    }
}
{ 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
}

{
  if (zz == 1) 
    { 
      gsub("\"","",$0) ; gsub("NULL","",$0) ;
      split($0,aa,"\t") ;
      z=aa[1];r = r2ir[z] ; if (r < 1) next ;
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
}
{ 
  if (zz == 2 || zz == 3)
    { 
	r = $1 ; ir=r2ir[r];if(ir<1) next ;
      x = $2 ; 
      if (zz == 2)
	{
	  f1[ir,x] = $3 ; f2[ir,x] = $4 ;
	  f1t[ir] += $3 ; f2t[ir] += $4 ;
	}
      if (zz == 3)
	{
	  l1[ir,x] = $3 ; l2[ir,x] = $4 ;
	  l1t[ir] += $3 ; l2t[ir] += $4 ;
	}
      if (0+x > 0+xmax) xmax = x ; 
      # printf("zz=%s xmax=%s r=#%s# x=#%s#\n", zz, xmax,$1,$2) ; print ;
    }
  next ;
}
END { 
    printf ("# xmax=%d\n", xmax) ;
    
    for (r = 1 ; r <= nr ; r++)
    {
	if (good[r] < 1) continue ;
	printf("\t%s first base", rr[r]) ;
	printf("\t%s last base", rr[r]) ;
	printf("\t%s first base", rr[r]) ;
	printf("\t%s last base", rr[r]) ;
    }
    
    printf ("\t\tRun") ;
    for (r = 1 ; r <= nr ; r++)
    {
	if (good[r] < 1) continue ;
	printf("\t%s first base", rr[r]) ;
	printf("\t%s last base", rr[r]) ;
	printf("\t%s first base", rr[r]) ;
	printf("\t%s last base", rr[r]) ;
    }
    if (nRunid)
    {
	printf ("\nRunId") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	}
	printf ("\t\tRunId") ;
	for (r = 1 ; r <= nr ; r++)
	{ 
	    if (good[r] < 1) continue ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	    printf("\t%s", runid[r]) ;
	}
    }
    
    if (nMachine >= 1)
    {  
	printf ("\nMachine") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	}
	
	printf ("\t\tMachine") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	    printf("\t%s", machine[r]) ;
	}
    }
    
    if (nSample >= 1)
    {
	printf ("\nSample") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	}
	printf ("\t\tSample") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	    printf("\t%s", sample[r]) ;
	}
    }
    
    if (nSystm >= 1)
    {
	printf ("\nSystem") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	}
	printf ("\t\tSystem") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	    printf("\t%s", systm[r]) ;
	}
    }
    
    if (nTissue >= 1)
    {
	printf ("\nTissue") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	}
	printf ("\t\tTissue") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	    printf("\t%s", tissue[r]) ;
	}
    }
    
    if (nTitre >= 1)
    {
	printf ("\nTitle") ;
	for (r = 1 ; r <= nr ; r++)
	{
	    if (good[r] < 1) continue ;
	    printf("\t%s", titre[r]) ;
	    printf("\t%s", titre[r]) ;
	    printf("\t%s", titre[r]) ;
	    printf("\t%s", titre[r]) ;
	}
	printf ("\t\tTitle") ;
      for (r = 1 ; r <= nr ; r++)
      {
	  if (good[r] < 1) continue ;
	  printf("\t%s", titre[r]) ;
	  printf("\t%s", titre[r]) ;
	  printf("\t%s", titre[r]) ;
	  printf("\t%s", titre[r]) ;
      }
    }
    
  if (nsTitre >= 1)
  {
      printf ("\nSorting title") ;
      for (r = 1 ; r <= nr ; r++)
      {
	  if (good[r] < 1) continue ;
	  printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	}
      printf ("\t\tTitle") ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	   printf("\t%s", stitre[r]) ;
	}
    }

  if (nsTitre2 >= 1)
    {
      printf ("\nSorting title 2") ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	}
      printf ("\t\tTitle") ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	   printf("\t%s", stitre2[r]) ;
	}
    }

   if (noTitre >= 1)
    {
      printf ("\nTitle") ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	}
      printf ("\t\tTitle") ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	   printf("\t%s", otitre[r]) ;
	}
    }

   printf ("\nPosition") ;
  for (r = 1 ; r <= nr ; r++)
    {
	   if (good[r] < 1) continue ;
	   printf("\t%s first base 1", rr[r]) ;
	   printf("\t%s last base 1", rr[r]) ;
	   printf("\t%s first base 2", rr[r]) ;
	   printf("\t%s last base 2", rr[r]) ;
    }
  printf ("\t\tPosition") ;
  for (r = 1 ; r <= nr ; r++)
    {
	   if (good[r] < 1) continue ;
	   printf("\t%s first base 1", rr[r]) ;
	   printf("\t%s last base 1", rr[r]) ;
	   printf("\t%s first base 2", rr[r]) ;
	   printf("\t%s last base 2", rr[r]) ;
    }

  for (x = 1 ; x <= xmax ; x++)
    {
      printf ("\n%d", x) 
	for (r = 1 ; r <= nr ; r++)
	  {
	      if (good[r] < 1) continue ;
	      printf("\t%d", 0+f1[r,x]) ;
	      printf("\t%d", 0+l1[r,x]) ;
	      printf("\t%d", 0+f2[r,x]) ;
	      printf("\t%d", 0+l2[r,x]) ;
	  }
      printf ("\t\t%d", x) ;
      for (r = 1 ; r <= nr ; r++)
	{
	   if (good[r] < 1) continue ;
	   if (f1t[r] > 0)
	       printf("\t%.2f", 100.0 * f1[r,x]/f1t[r]) ;
	   else
	    printf ("\t0") ;
	   if (l1t[r] > 0)
	       printf("\t%.2f", 100.0 * l1[r,x]/l1t[r]) ;
	   else
	       printf ("\t0") ;

	   if (f2t[r] > 0)
	       printf("\t%.2f", 100.0 * f2[r,x]/f2t[r]) ;
	   else
	    printf ("\t0") ;
	   if (l2t[r] > 0)
	       printf("\t%.2f", 100.0 * l2[r,x]/l2t[r]) ;
	   else
	       printf ("\t0") ;	}
    }
  printf ("\nTotal", x) 
    for (r = 1 ; r <= nr ; r++)
      {
	   if (good[r] < 1) continue ;
	   printf("\t%d\t%d", 0+f1t[r], 0+l1t[r]) ;
	   printf("\t%d\t%d", 0+f2t[r], 0+l2t[r]) ;
      }
    printf ("\t\t%d", x) ;
    for (r = 1 ; r <= nr ; r++)
    {
	   if (good[r] < 1) continue ;
	   printf ("\t100\t100") ;
	   printf ("\t100\t100") ;
    }
    printf ("\n") ;
}

