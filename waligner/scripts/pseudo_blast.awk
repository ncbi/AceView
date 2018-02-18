/^gF /{ gF = $2 ; next ; }
/^gR /{ gR = $2 ; next ; }

/^>/{ zz = 1 ;
 a = substr($1,2) ;
 gsub(/_F3/,"",a) ;
 next ;
}
{
  # printf("zz=%d ",zz); print;
  if (zz < 1)
    {
      s[$1] = $4 ;  /* clayton version of the a sequence */
      next ;
    }
  nseq++ ;
  # find the best ali relative to the forward genome
  bestScore = 0 ; ggg = gF ;
  for (strand = -1 ; strand < 2 ; strand += 2)
    {
      if (strand == -1) gg = gR ;
      if (strand ==  1) gg = gF ;
      ggL = length(gg) ;
      for (dx = 0 ; dx < ggL ; dx++)
	{
	  n = 0 ;
	  for (i = 1 ; i <= length($1) &&  i + dx <= ggL ; i++)
	    if (substr(gg,i+dx,1) == substr($1,i,1))
	      n++ ;
	  if (n > bestScore)
	    {
	      bestScore = n ; bestDx = dx ; bestStrand = strand ;
              ggg = substr(gg, dx+1, 51) ;
	    }
	}
    }
  if (2 * bestScore < length($1)) next ; 
  qq = 15000000 + 10*nseq + 1000*bestDx*bestStrand +  1000000*bestStrand ;

  printf("%d\t%s\tReceived\t%s",qq,a,s[a]) ;
  gsub(/t/,"T",$1) ;
  # count the errors relative to the original sequence
  n = 0 ;
  for (i = 3 ; i <= length($1) ; i++)
    if (substr(s[a],i,1) != substr($1,i,1))
      n++ ;
  if (n>0)
    {
      printf("\n%d\t%s\tSq2SW  %d\t  ",qq+1,a,n) ;
      for (i = 3 ; i <= length($1) ; i++)
	if (substr(s[a],i,1) != substr($1,i,1))
	  printf("|") ;
	else 
	  printf(" ") ;
      n++ ;
    }
  printf ("\n%d\t%s\tOriginal\t%s", qq+2,a, $1) ;


  # count the errors relative to the genome aligned sequence
  n = 0 ;
  for (i = 3 ;i <= length($1) ; i++)
    if (substr(ggg,i,1) != substr($1,i,1))
      n++ ;
  if (n>0)
    {
      printf("\n%d\t%s\tSq2Ge  %d\t  ",qq+3,a,n) ;
      for (i = 3 ; i <= length($1) ; i++)
	if (substr(ggg,i,1) != substr($1,i,1))
	  printf("|") ;
	else 
	  printf(" ") ;
      n++ ;
    }
  printf ("\n%d\t%s\tGenome  \t%s %d %d %d\n%d\n%d\n", qq+4,a, ggg, bestStrand, bestScore, bestDx,qq+5,qq+6) ;

  # printf ("bestscore = %d bestDx=%d bestStrand=%d \n%s\n\n", bestScore, bestDx, bestStrand, ggg) ;
  # exit ;
}
