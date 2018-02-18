/^#/ {
    if (mx < 1)
    {
	for (i = 1 ; i <= NF ;i++)
	    if (tolower(substr($i,1,24)) == "maximal allele frequency")
		mx = i;
	run [mx-1] =  magic ".union" ;
        run[mx] =  magic ".sum" ;
        for (i = mx+1 ; i<=NF ; i++)
	    run[i] = $i ; nRuns=NF;
    }
    next ;
}
/^ZZZZZ/ {rejected++; next; }
{ 
    type = tolower ($11) ;
    k=0 ; kmx = 0 ;
    for (i = mx+1 ; i<=NF; i++)
    {
	if ($i > 5)
	{
	    x = 2 * $i / 100.0 ; 
	    k += x ;
	    if ($i > kmx) kmx = $i ;
	    zz[rejected,type,i] += x ;
	}
    }
    if (k > 0)
    {
	union[type]++ ;
	if (kmx > 20) zz[rejected,type,mx-1]++ ; # snp > 20% at least once in this project
	zz[rejected,type,mx] += k ;
    }
}
END {
    for (t in union)
    {
	v = t ;
	slide = "" ;
	if (substr (v,1,1) == "*") 
	{
	    slide = "*" ;
	    v = substr (v, 2) ;
	}
	if (substr(v,1,3) == "del")
	{
	    p = length(v) - 3 ;
	    w = substr("-----" ,1 , p) substr(v,4) ;
	    v = w;
	}
	if (substr(v,1,3) == "ins")
	{
	    p = length(v) - 3 ;
	    w = substr("+++++",1,p) substr(v,4) ;
	    v = w;
	}
	aceName [t] = slide v;
    }
    for (i = mx -1 ; i <= nRuns ; i++)
    {
	for(t in union)
	{
	    if (zz[1,t,i]>0)
	    {
		printf("%s\t%d\t%s\t%d\n", run[i], zz[1,t,i], aceName[t], zz[2,t,i]) ;
	    }
	}
    }
}
