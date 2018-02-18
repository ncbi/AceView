/^#/{next;}
/^ZZZZZ/{zz++; j = 1 ; next;}
{
  if(zz<1)
    {
      if ($3 > $2) { n++; ca[n] = $1 ; a1[n] = $2 ; a2[n] = $3; }
      next;
    }
  cb = $1 ; b1 = $2 ; b2 = $3 ;
  for (i = j ; i>1 ; i--)
    if (ca[i] != cb ||  a2[i] < b1) 
      break ;
  b11 = b1 ; b12 = b2 ;
  for (i = i ; i <= n + 1 ; i++)
    {
      if (i > n || ca[i] > cb ||  a1[i] > b2) 
	{
	  if (a1[i] < b11 && a2[i] > b11)
	    b11 = a2[i] ;
	  if (b11 + 100 < b2)
	    printf ("%s\t%d\t%d\tb1=%d\tb2=%d\ta1=%d\ta2=%d\ti=%d\tn=%d\n",cb,b11,b2,b1,b2,a1[i],a2[i],i,n) ;
	  break ;
	}
      if (ca[i] < cb || a2[i] < b1)
	continue ;
      if (a1[i] <= b11 && a2[i] > b11)
	b11 = a2[i] ;
      if (a1[i] > b11 + 100)
	printf ("%s\t%d\t%d\n",cb,b11,a1[i]) ;
      if (a2[i] > b11)
	b11 = a2[i] ;
    }
  if (i > n) i = n ;
  j = i ;
}

