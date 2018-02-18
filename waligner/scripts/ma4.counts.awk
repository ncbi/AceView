/^ZZZZZ/{zz++;next;}
/^#/{next;}
{ 
  line++ ; ok = 0 ; 

  g = $7 ; if (g == "") next ;

  gids = $8 ;

  if (avOnly == 1 && gids == "")
    next ;
  if (avOnly == 2 && gids != "")
    next ;

  hasProbe = 0 + $9 ;   
  map = $10 ; if (map == "") next ;
    
  RNAp[4] = 0 ; 
  RNAm[4] = 0 ; 
  AGLp[4] = 0 ;
  AGLm[4] = 0 ;

  for (type = 0 ; type < 4 ; type++)
    {
      RNAp[type] = 0 + $(40 + 6 * type) ;
      AGLp[type] = 0 + $(42 + 6 * type) ;

      RNAp[4] += RNAp[type] ;
      AGLp[4] += AGLp[type] ;


      RNAm[type] = 0 - $(41 + 6 * type) ;
      AGLm[type] = 0 - $(43 + 6 * type) ;

      RNAm[4] += RNAm[type] ;
      AGLm[4] += AGLm[type] ;
    }

  for (type = 0 ; type < 5 ; type++)
    {
      # printf("YYYY %s %d %d\n", g, type,  RNAp[type]) ;
      # never seen

      cc[21] = "Probes" ; 
      nn[21, type] += hasProbe ;
      cc[22] = "Differential probes" ;
      nn[22, type] += AGLp[type] + AGLm[type] ;

      if (0 + AGLp[type] + AGLm[type] + RNAp[type] + RNAm[type] == 0)
	continue ;

      # only RNA
      else if (0 + AGLp[type] + AGLm[type] == 0)
	{
	  cc[2] = "DEG with probe seen only in RNA-seq" ;
	  cc[3] = "DEG seen in RNA-seq, no probe" ;

	  cc[6] = "Complex DEG with probe seen only in RNA-seq" ;
	  cc[7] = "Complex DEG seen in RNA-seq, no probe" ;

	  k = 0 ; if (hasProbe == 0) k = 1 ;
	  if (RNAp[type] > 0 && RNAm[type] > 0)
	    nn [6+k,type]++ ;
	  else
	    nn [2+k,type]++ ;
	}

      # only AGL
      else if (0 + RNAp[type] + RNAm[type] == 0)
	{  
	  cc[4] = "DEG seen only in Agilent" ;
	  cc[8] = "Complex DEG seen only in Agilent" ;

	  # complex
	  if (AGLp[type] > 0 && AGLm[type] > 0)
	    nn [8,type]++ ;
	  else
	    nn [4,type]++ ;
	}
      
      # Both
      else 
	{
	  cc[1] = "DEG seen consistently as up or down in RNA-seq and Agilent" ;
	  cc[5] = "Complex DEG seen consistently in RNA-seq and Agilent" ;
          cc[9] = "Contradiction between RNA-seq and Agilent, one sees the gene up, the other down"

	  # coherence
	  if (AGLp[type] +  RNAp[type] == 0 ||  AGLm[type] + RNAm[type] == 0)
	    nn [1, type]++ ;
	  # complex
	  else if ((AGLp[type] > 0  && AGLm[type] > 0) || (RNAp[type] > 0 && RNAm[type] > 0))
	    nn [5, type]++ ;
	  # contradiction
	  else
	    nn [9, type]++ ;
	}
    }
}
END {
  # create redundant partial sums

  for (type = 0 ; type < 5 ; type++)
    {
      cc[10] = "Any" ; nn [10, type] = 0 ;
      for (i = 1 ; i <= 9 ; i++)
	nn [10, type] += nn[i, type] ;

      cc[11] = "RNA-seq and Agilent" ;
      nn [11, type] = nn [1, type] +  nn[5, type] ;
      
      cc[12] = "RNA-seq only" ; 
      nn [12, type] = nn [2, type] +  nn[3, type] +  nn [6, type] +  nn[7, type] ;
      
      cc[13] = "Agilent only" ; 
      nn [13, type] = nn [4, type] +  nn[8, type] ;
      
      cc[14]= "Contradictions" ;
      nn [14, type] = nn [9, type] ;
      
      cc[15] = cc[3] ;  nn [15, type] = nn [3, type] ;
      cc[16] = cc[2] ;  nn [16, type] = nn [2, type] ;
      cc[17] = cc[1] ;  nn [17, type] = nn [1, type] ;
      cc[18] = cc[4] ;  nn [18, type] = nn [4, type] ;
      cc[19] = cc[4] ;  nn [19, type] = nn [9, type] ;
      
      cc[20] = "Complex genes" ; 
      nn [20, type] = nn [5, type] +  nn[6, type] +  nn [7, type] +  nn[8, type] ;
    }

  split ("Stage1/Stage4/Stage4s/MNA/Any stage",types,"/") ;

  printf ("# Type") ;
  for (n = 1 ; n <= 22 ; n++)
    printf ("\t%s", cc[n]) ;
  
  for (type = 0 ; type < 5 ; type++)
    {
      printf ("\n%s", types[type+1]) ;
      for (n = 1 ; n <= 22 ; n++)
	printf ("\t%d", nn[n,type]) ;
    }
  printf ("\n") ;

}
