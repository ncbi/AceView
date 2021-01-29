/^#/ { next ; }
/^ZZZZZ/ { zz++ ; next ; }
{ 
    if (zz<1)
    {
	for (i = 5 ; i <= NF ; i++)
	    u[$1,i] = $i ;
	next ;
    } 
    if ($4 + 0 == 0) next ;
    printf ("%s",$1) ;
    for (i = 2 ; i <= NF ; i++)
    {
	z = 1000 * u[$1,i] + 0 ;
	if (i<8 || i >= NF-5)
	    printf ("\t%s", $i) ;
	else
	{
	    if (z >= m1 * $i && z < $i*m2)
	    {
		if (byKb == 2)
		{
		    ln = 0 + $4/1000.0 ;
		    if (ln == 0) ln = 1 ;
		    printf("\t%.3f",u[$1,i]/ln) ;
		}
		else
		    printf("\t%s",u[$1,i]) ;
	    }
	    else 
		printf("\t-1");
	}
    }
    printf("\n") ;
}

