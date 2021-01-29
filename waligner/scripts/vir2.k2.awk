/^ZZZZZ/ { zz++ ; next ; }
{
    if (zz < 1) 
    {
	isRun[$1] = 1 ;
    }
    if (zz < 2)
    {
	run = $1 ;
	if (isRun[$1] == 1 && $2 > r2Mb[$1]) r2MB[$1] = $2 ;
	next ;
    }
    if (zz < 3) 
    {
	g = $2 ;
	if ($1 == t)
	    ok[g] = 1 ;
	next;
    }
}
/^#Run/ {   # This must be applied to the last#Runnline of the table because some inputs contain several asynchronde definitions, say virus, then bacteria
    for(i = 2 ; i <= NF ; i++)
    {
	colRun[i] = isRun[$i] ;
	if (byKb == 1) col2Mb[i] = r2MB[$i] ; 
    }
    next ;
}
{
    if (ok[$1] == 1)
    {
	printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s",$1,$2,$3,$4,$5,t,$5) ;
	z = 0 ; mx = 0 ; ix = 0 ; w[1] = 0 ;
	for (i = 2 ; i <= NF ; i++)
	{
	    w[i] = 0 ;
	    # if ($i == 6644926) { printf("z=%s i=%d colRun=%s %s %d\n",$i,i,colRun[i],$1,ok[$1]) ; exit(0);}
	    if (colRun[i] == 1 && $i > 0)
	    {
		w[i] = $i ;
		if (col2MB[i] > 0 && $i > 0) w[i] = $i/col2MB[i]  ;
		if (w[i] > 0) z += w[i] ; ;
		if (w[i]> mx) mx = w[i] ; 
	    }
	}
	printf ("\t%.2f",z) ;  # replace column 8
	for (i = 9 ; i <= NF - 6 ; i++)
	    printf ("\t%.2f",w[i]) ;
	printf ("\t%.2f", mx) ;
	if (z==0)z=1;printf ("\t%.2f", 100.0*mx/z) ;

	# export the number of runs needed to reach 25% 50% 75% of the total 
        s = 0 ; z25 = z/4.0 ; z50 = 2.0 * z25 ; z75 = 3.0 * z25 ; z100 = .9999 * z ;
	i25 = 0 ; i50 = 0 ; i75 = 0 ; i100 = 0 ;
	k = 0 ; s = 0 ;
	ix = asort (w) ;
	for (i = ix ; i >=9 ; i--)
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

