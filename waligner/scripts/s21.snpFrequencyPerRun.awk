{ 
    line++ ;
    if (line == 1)
    {
	ii0 = 0 ;
	for (i = 1 ; ii0 == 0 && i < NF ; i++)
	{
	    if ($i == "Run") 
		ii0 = i + 2 ;
	}
	for (i = ii0 ;  i <= NF ; i++)
	{
	    i2r[i] = $i ;
	}
	next ;
    }
}
/^#/ {next ; }
{
    for (i = ii0 ;  i <= NF ; i++)
    {
	if ($i >= 0)
	{
	    x = int($i) ;
            zz[i, x



/^
