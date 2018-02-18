/^#/{ nTitle++ ; 
    if (nTitle >= 3 && nTitle <= 5)
	print ; 
    next ;
}
{
    if (index($2, "Kept ") > 0) next ;
    if (NF > 10) run = $1 ;
    it = t2it[$2] ; 
    if (it+ 0  < 1)
    {
	itMax++ ; it = itMax ;
	t2it[$2] = it ; it2t[it] = $2 ;	
    }
    for (i = 2 ; i <= NF ; i++)
    {
	x = 10 * (i - 3) ;
	zz [it, x] += $i ; zzz [x] += $i ; zzt [it] += $i ; cumul += $i ;
    }
    if (x > xMax)
	xMax = x ;
}
END {
    printf ("# Run\tTranscript\tAverage coverage of the first 8kb") ;
  
    for (x = 0 ; x <= xMax ; x += 50)
	printf ("\t%d bp", x) ;
    kept = 0 ; cumul = 0 ;
    for (it = 1 ; it <= itMax ; it++)
    {
         # zzt[it] > 100000 || zz[it,0] > 1000
	if (1)
	{
	    for (x = 0 ; x <= xMax ; x += 50)
	    {
		u = 0 ;
		for (x1 = x - 20 ; x1 <= x+20 ; x1 += 10)
		{
		    x2 = x1 ; if (x2 < 0) x2 = -x2 ;
		    u += zz [it, x2]/100 ;
		}
		u = int (.4999 + u/5) ;
		uu[x] = u ;
	    }
	    for (x = 0 ; x <= 50 ; x += 50)
		uu[x] = uu[50] ;
	    uu0 = 0 ; k = 0 ;
	    for (x = 0 ; x <= xMax && x <= 8000 ; x += 50)
	    { uu0 += uu[x] ; k++ ; }
	    uu0 = int (.4999 + uu0/k) ;
	    if (uu0 > 1 && uu[0] >= 10)
	    {
		kept++ ;
		printf ("\n%s\t%s", run, it2t[it]) ;
		printf ("\t%d", uu0) ; cumul += uu0 ;
		for (x = 0 ; x <= xMax ; x += 50)
		{
		    printf ("\t%d", uu[x]) ;
		    xxk[x] += uu[x] ; 
		}
	    }
	}
    }
    printf ("\n%s\tKept %d/%d transcripts\t%d", run, kept, nMax, cumul) ;
    for (x = 0 ; x <= xMax ; x += 50)
	printf ("\t%d", xxk[x]) ;
    printf ("\n") ;
}
