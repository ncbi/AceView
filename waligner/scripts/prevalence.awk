/^ZZZZZ/{
  zz++;
  next;
}

/^#Run/{nnr++;if(nnr==1){for (i=1;i<=30 && i<=NF;i++){if($i=="_SumOfAllReadsInProject")imin=i+1;}for(i=imin;i<=NF;i++)if(substr($i,1,3)!="Rhs"){imax=i;next}imax=NF;}}
/^#/{next}

{if(zz<1) {
  selected_gene[$1] = 1 ;
  next ;
}}

{if(zz==1) {
  m = $1 ; g = $5 ; 
  m2g[m] = g ;  
  if (selected_gene[g] + selected_gene[m] < 1) next ;
  if((hasGid=="with" && $6=="")||(hasGid=="no" && $6)) next;
  if(length(g)>=1) 
    {  
      geneok[g] = 1 ;
      if (GM == "GENE" || GM2 == "m2g")
	{
	  ggok[g] = 1 ; gg[g] = -1 ; 
	}
      else
	{ ggok[m] = 1 ; gg[m] = -1 ; }
    }
  next ;
}}
{if(zz==2) {
  g=$1 ;
  if (geneok[g] == 1)
    {
      g2t[g]=$4 ; tt[$4]=1;
    }
  next;
}}


/AceView/{ if (index($2,"AceView")<1) next ;
  g=$1;
  if (GM2 == "m2g") 
    {
      m = g ;
      g = m2g[m] ;
      if (length(g)<1) printf("m=%s no Gene\n", m) ;
    }
  if (0+ggok[g] < 1) 
    {
      # printf("%s\t%s\trejected\n", g, m) ;
      next ;
    }

  n = 0 ;
  for (i = imin ; i <= NF && i < imax ; i++)
    {
      if ($i + 0 > 0)
	n++ ;
    }
  # printf("%s\t%s\taccepted\t%d\t%d\t%d\t%d\n", g, m, zz, n, imin, imax) ;

  if (n > nmax) nmax = n ;
  if (zz == 3 && n > 0+gg[g]) gg[g] = n ; 
  if (zz == 4 && n > 0 + touched[g]) touched[g] = n ; 

  next ;
}
END {
  for(g in gg)
    for(n=-1;n<=nmax;n++)
      {
	if (GM == "GENE" || GM2 == "m2g") 
	  t=g2t[g];
	else 
	  t = g2t[m2g[g]] ;
	if (length(t)<1) 
	  {
	    t="no_type" ; 
	    if(1)printf("ZZZ %s\n",g) ;
	  }

	if (n == -1 || (n == 0 && 0+touched[g]>0) || (n > 0 && 0+gg[g]>=n))
	  { 
	    tt[t]=1;t2n[t,n]++;ttn[n]++;
	  }
      }
  
  printf("A\tType\t%s\tTouched",title);
  for (n=1;n<=nmax;n++)printf("\t%d",n) ;
  printf("\t\tType\t%s\tTouched",title);
  for (n=1;n<=nmax;n++)printf("\t%.1f",100.0*n/(nmax+.000001)) ;
  for(t in tt)
    {
      if(0+t2n[t,0]<1)continue;
      printf("\nC\t%s",t);
      for(n=-1;n<=nmax;n++)
	printf("\t%d",t2n[t,n]);
      printf ("\t\t%s",t) ;
      for(n=-1;n<=nmax;n++)
	{
	  z=100*t2n[t,n]/(.0001+t2n[t,-1]) ;
	  printf("\t%.2f",z);
	}
    }
  printf("\nT\tAny type %s %s", GM2, title);
  for(n=-1;n<=nmax;n++)
    printf("\t%d",ttn[n]);
  printf ("\t\tAny type %s %s", GM2, title) ;
  for(n=-1;n<=nmax;n++)
    {
       z=100*ttn[n]/(.0001+ttn[-1]) ;
        printf("\t%.2f",z);
    }
  printf("\n");
  printf("a\nb\tType\tAnnotated\tTouched");
  for(n=1;n<=nmax;n++)
    printf("\t%d",n);
  for(t in tt)
    {
      printf("\np\t%% %s",t);
      for(n=-1;n<=nmax;n++)
	printf("\t%.1f",100*t2n[t,n]/(.0001+ttn[n]));
    }
  printf("\n");
}
