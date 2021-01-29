/^#/{ next; }
/^ZZZZZ/ {zz++ ; next ; }
{
    if (zz+0 < 1)
    {
	g = $4 ; s = $5 ; 
	ig = 0 + g2ig[g] ; 
	if (ig == 0)
	{
	    igMax++ ; ig = igMax ;
	    ig2g[ig] = g ; g2ig[g] = ig ;
	}
	s2ig[s] = ig ;
	next ;
    }
}
{
    s = $1 ;
    ig = 0 + s2ig[s] ; if (ig < 1) next ; 
    for (i = 2 ; i <= NF ; i++)
    {
	if (i < 8) zg[ig, i] = $i ;
	if (i >= 8) zg[ig, i] += $i ;
    }
    zg[ig,1]++ ;
    if (NF > iMax) iMax = NF ;
    next ;
}
END {

    for (ig = 1 ; ig <= igMax ; ig++)
	if (zg[ig, 1] > 0)
	{
	    printf ("%s", ig2g[ig]) ;
	    mx = 0 ; 
	    for (i = 2 ; i < 8 ; i++) 
		printf ("\t%s", zg[ig,i])  ;
	    printf ("\t%.2f", zg[ig,8])  ;
	    for (i = 9 ; i <= iMax - 6 ; i++) 
	    {
		printf ("\t%.2f", zg[ig,i])  ;
		w[i] = zg[ig,i] ;
		if (mx < zg[ig,i]+ 0)
		    mx = zg[ig,i] ;
	    }
	    printf ("\t%.2f", mx) ;
	    z = zg[ig,8];if (z==0)z=1;printf ("\t%.2f", 100.0*mx/z) ;
	    

	    # export the number of runs needed to reach 25% 50% 75% of the total 
	    s = 0 ; z25 = z/4.0 ; z50 = 2.0 * z25 ; z75 = 3.0 * z25 ; z100 = .9999 * z ;
	    i25 = 0 ; i50 = 0 ; i75 = 0 ; i100 = 0 ;
	    k = 0 ; s = 0 ;
	    ix = asort (w) ;
	    for (i = ix ; i >=1 ; i--)
	    {
		k++ ;
		s += w[i] ;
		if (s >= z25 && ! i25) i25 = k ;
		if (s >= z50 && ! i50) i50 = k ;
		if (s >= z75 && ! i75) i75 = k ;
		if (s >= z100 && ! i100) i100 = k ; 
	    }
	    printf ("\t%d\t%d\t%d\t%d", i25, i50, i75, i100) ;
	    
	    printf("\n") ;
	}
}
