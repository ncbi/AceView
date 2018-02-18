BEGIN { 
  nBad = 0 ;
   itMax = split("Sexy_Brain Sexy_UHR Sexy_any pA pAnorm any adps adrn brain breast colon heart kidney liver lung lymph ovary prostate skelet testes thyrd wbc", i2t, / /) ;

   # itMax = split("adps adrn", i2t, / /) ;

  #itMax = split(nams, i2t, /,/) ;
  #printf ("itMax= %d\n", itMax) ;
  for (it = 1 ; it <= itMax ; it++)
    {
      t2i[i2t[it]] = it ; # so zero is rejected
      # printf ("%d $%s$ \n", it, i2t[it]) ;
    }
}
{  
  gsub(/\"/,"",$0) ; gsub(/G_/,"",$0) ; 
  g = $1 ; t = $2 ; z = $3 ;
  if (z < limit) {  next ;}
  
  it1 = t2i[t] ; 
  if (it1 < 1)
    {
      bad[t] = 1 ; nBad++ ; next ;  # we cannot handle new names
    }
  if ( gg[g] != 1) { gg[g] = 1 ; ngT++ ; }
  gt[g,it1] = 1 ;
}
END {
  ngI = ngT ;
  if (nBad > 0)
    {
      for (t in bad) printf ("UNKNOWN name \"%s\",\n", t) ;
      printf("I only know:") ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
    }
  else
    {
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{
	  for (g in gg)
	    {
	      # printf ("it1=%d g=%s gt=%d\n", it1, g, gt[g, it1]) ;
	      if (gt[g, it1] == 1)
		{
		  nn[it1,it1]++ ;
		  for (it2 = 1 + it1 ; it2 <= itMax ; it2++) # intersection
		    {
		      if (gt[g, it2] == 1)
			nn[it1,it2]++ ;
		    }
		}
	      else
		{
		  if (missing[g] == 0) { ngI-- ; missing[g] = 1 ; }
		}
	      for (it2 = 1 ; it2 < it1 ; it2++) # union
		{
		  if (gt[g, it1] + gt[g, it2] >= 1)
		    nn[it1,it2]++ ;
		}
	    }
	}

      if (ngT == 0) ngT = -1 ;

      printf("%d genes > %d Union/Intersection", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    printf("\t%d", nn[it1,it2]) ;
	}
      printf("\nDiagonale") ;
        for (it2 = 1 ; it2 <= itMax ; it2++)
	  printf("\t%d", nn[it2,it2]) ;

      printf ("\n\nUnion: %d genes in the union of all columns\n", ngT) ;
      printf("%d genes > %d Union", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    {
	      if (it1 <= it2) { jt1 = it2 ; jt2 = it1 ; }
	      else  { jt1 = it1 ; jt2 = it2 ; }
	      printf("\t%d", nn[jt1,jt2]) ;
	    }
	}
      printf("\nDiagonale") ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%d", nn[it2,it2]) ;
      
      printf ("\n\nIntersection: %d genes in the intersection of all columns\n", ngI) ;
      printf("%d genes > %d Intersection", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    {
	      if (it1 >= it2) { jt1 = it2 ; jt2 = it1 ; }
	      else  { jt1 = it1 ; jt2 = it2 ; }
	      printf("\t%d", nn[jt1,jt2]) ;
	    }
	}
      printf("\nDiagonale") ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%d", nn[it2,it2]) ;
      
      printf ("\n\nUnion/Intersection\n") ;
      printf("%d genes > %d Union/Intersection", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    {
	      z = ngT ;
	      printf("\t%.1f", 100 * nn[it1,it2]/z) ;
	    }
	}


      printf ("\n\nUnion: %d genes in the union of all columns\n", ngT) ;
      printf("%d genes > %d Union", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    {
	      if (it1 <= it2) { jt1 = it2 ; jt2 = it1 ; }
	      else  { jt1 = it1 ; jt2 = it2 ; }
	      z = ngT ;
	      printf("\t%.1f", 100 * nn[jt1,jt2]/z) ;
	    }
	}

      
      printf ("\n\nIntersection: %d genes in the intersection of all columns\n", ngI) ;
      printf("%d genes > %d Intersection", ngT, limit) ;
      for (it2 = 1 ; it2 <= itMax ; it2++)
	printf("\t%s", i2t[it2]); 
      for (it1 = 1 ; it1 <= itMax ; it1++)
	{ 
	  printf("\n%s", i2t[it1]) ;
	  for (it2 = 1 ; it2 <= itMax ; it2++)
	    {
	      if (it1 >= it2) { jt1 = it2 ; jt2 = it1 ; }
	      else  { jt1 = it1 ; jt2 = it2 ; }
	      z = ngT ;
	      printf("\t%.1f", 100 * nn[jt1,jt2]/z) ;
	    }
	}
    }
  printf ("\n") ;
}



