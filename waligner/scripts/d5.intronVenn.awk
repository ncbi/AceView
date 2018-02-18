/./ { 
  z=$1; 
  for (j = 2 ; j <= NF ; j++)
   {
     split($j,aa,":") ; 
     
     t = aa[2] ; # sample A C D B
     tt[t]++ ;
     
     n = 0 + aa[3];
     if (n >= limit)
       {    
         nn[z] += n ; 
	 nt[z,t,j] = n ;
	 nt[z,"x",j] += n ;
       }
   }
}

function export(t) {
  printf ("\n%s", t) ;
  for (j = 2 ; j <= maxNam ; j++)
    printf ("\t%d", nnn[t,j]);
}

END {
  nam[2] = "Any" ; 
  nam[3] = "New" ;
  nam[4] = "RefSeq-2012" ;
  nam[5] = "AceView-2010" ;
  nam[6] = "Encode-2012" ;
  nam[7] = "WEHI" ;
  maxNam = 7 ;

  for (z in nn)
    {
      for (j = 2 ; j <= maxNam ; j++)
	{
	  if (nt[z,"A",j] + nt[z,"C",j] + nt[z,"D",j] + nt[z,"B",j] > 0)
	      nnn["Union",j]++ ; 

	  za = nt[z,"A",j] ;
	  zb = nt[z,"B",j] ;
	  zc = nt[z,"C",j] ;
	  zd = nt[z,"D",j] ;

	  if (za + zb + zc + zd >= 5)
	    {
	      za = nt[z,"A",j]/(1+1*tt[t]) ;
	      zb = nt[z,"B",j]/(1+1*tt[t]) ;
	      zc = nt[z,"C",j]/(1+1*tt[t]) ;
	      zd = nt[z,"D",j]/(1+1*tt[t]) ;
	      
	      if ((za > 1.05 * zc && zc > 1.05 * zd && zd >= zb) ||  (za >  1.05 * zc && zc >= zd && zd >  1.05 * zb) || (za >= zc && zc >   1.05 * zd && zd >   1.05 * zb))
		nnn["Titrating",j]++ ; 
	      if ((za < .95 * zc && zc < .95 * zd && zd <= zb) ||  (za <  .95 * zc && zc <= zd && zd <  .95 * zb) || (za <= zc && zc <   .95 * zd && zd <  .95 * zb))
		nnn["Titrating",j]++ ; 
	    }
	  
	  for (t1 in tt) 
	    if (nt[z,t1,j] > 0) 
	      {
		nnn[t1,j]++ ;
		for (t2 in tt)
		  if (t2 > t1 && nt[z,t2,j] > 0)
		    {
		      nnn[t1 t2,j]++ ;
		      for (t3 in tt)
			if (t3 > t2 && nt[z,t3,j] > 0)
			  {
			    nnn[t1 t2 t3,j]++ ;
			    for (t4 in tt)
			      if (t4 > t3 && nt[z,t4,j] > 0)
				nnn[t1 t2 t3 t4,j]++ ;
			  }
		    }
	      }
	  
	  # triple Venn a (CD) B
	  if (nt[z,"A",j] > 0  && nt[z,"C",j] + nt[z,"D",j] == 0 && nt[z,"B",j] == 0) nnn["V1A",j]++ ;
	  if (nt[z,"A",j] == 0 && nt[z,"C",j] + nt[z,"D",j]  > 0 && nt[z,"B",j] == 0) nnn["V1CD",j]++ ;
	  if (nt[z,"A",j] == 0 && nt[z,"C",j] + nt[z,"D",j] == 0 && nt[z,"B",j]  > 0) nnn["V1B",j]++ ;
	  if (nt[z,"A",j]  > 0 && nt[z,"C",j] + nt[z,"D",j] == 0 && nt[z,"B",j]  > 0) nnn["V2AB",j]++ ;
	  if (nt[z,"A",j]  > 0 && nt[z,"C",j] + nt[z,"D",j]  > 0 && nt[z,"B",j] == 0) nnn["V2ACD",j]++ ;
	  if (nt[z,"A",j] == 0 && nt[z,"C",j] + nt[z,"D",j]  > 0 && nt[z,"B",j]  > 0) nnn["V2BCD",j]++ ;
	  if (nt[z,"A",j]  > 0 && nt[z,"C",j] + nt[z,"D",j]  > 0 && nt[z,"B",j]  > 0) nnn["V3ABCD",j]++ ;
	}
    }

 
 printf ("Venn diagram considering %d or more supports\n",  limit) ;
 for (j = 2 ; j <= maxNam ; j++)
   printf("\t%s", nam[j]) ;
 
 export("Union") ;
 export("Titrating") ;
 for (t1 in tt) 
   export(t1) ;

 for (t1 in tt) 
   for (t2 in tt)
     if (t2 > t1)
       export(t1 t2) ;

 for (t1 in tt) 
   for (t2 in tt)
     for (t3 in tt)
       if (t3 > t2 && t2 > t1)
	 export(t1 t2 t3) ;

 for (t1 in tt) 
   for (t2 in tt)
     for (t3 in tt)
       for (t4 in tt)
	 if (t4 > t3 && t3 > t2 && t2 > t1)
	   export(t1 t2 t3 t4) ;

 printf("\n") ;
 printf("\n") ;
 export("V1A") ;
 export("V1B") ;
 export("V1CD") ;
 printf("\n") ;
 export("V2AB") ;
 export("V2ACD") ;
 export("V2BCD") ;
 printf("\n") ;
 export("V3ABCD") ;
 printf("\n") ;
 printf("\n") ;
 printf("\n") ;
}

