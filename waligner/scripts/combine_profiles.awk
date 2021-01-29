/^#/{
    nc++;
    if (nc==2)
    {     # find the column of the total absolute score 
	iabs = 0 ;
	for (i=1;i<=NF;i++)
	{
	    if(index($i,"absolute")>0)
		iabs=i;
	}
	printf ("iabs=%d\n",iabs) ;
	next;
    }
    if (nc == 4)
    {   # locate the score columns
	ii1 = 0 ; ii2 = 0 ;
	for (i=1;i<=NF;i++)
	{
	    if(index($i,"<") + index($i,">") > 0)
	    {
		if (ii1 == 0)
		    ii1 = i ;
		ii2 = i
	    }
	}
	next;
    }
    next ;
}
{ # gene columns
    nc = 0 ;
    gene = $2 ;
    n = 0 ;
    for (i = ii1 ; i <= ii2 ; i++)
	if ($i > 10)
	    n++ ;
    s = $iabs + 0 ;
    score[gene] += s ;
    ng[gene] += n ;
}
END {
    for (g in ng)
    {
	if (ng[g] > 0)
	    printf ("%s\t%.2f\t%d\n", g, score[g], ng[g]) ; 
    }
}
