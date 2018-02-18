{line++;}
{ gsub (/\"/,"",$0); }
{ 
  if (line == 1)
    {
      print ;
      for (col = 1 ; col <= NF ; col++)
	{
	  n = split ($col, aa, "::") ;
	  fnam [col] = "ZERO" ;	  fparam [col] = "ZERO" ;
	  if (n > 1)
	    {
	      split (aa[2], bb, " ")  ;
	      fnam [col] = bb[1] ;
	      if (bb[2])
		fparam [col] = bb[2] ; 
	       # print col, aa[2], "::", bb[1], ":", bb[2] ;
	    }
	}
      next ;
    }  
}
{ if (line <  -2) next ; }
{
  for (col = 1 ; col <= NF ; col++)
    if (fnam [col] == "set")     # expect  :: set ali
      param[fparam[col]] = $col ;
  for (col = 1 ; col <= NF ; col++)
    {
      if (col > 1) printf ("\t") ;
      
      ff = fnam[col] ;
      # printf ("%s:",ff) ;
      if (ff == "ZERO" || ff == "set")
      {
	  if ($col == "NULL")
	      printf("") ;
	  else
	      printf ("%s", $col) ;
      }
      if (ff == "percent")
	{
          if ($col == "NULL")
	      printf("") ;
	  else
	  {
	      z = 0 + fparam [col] ;
	      if (z + 0 > 100 || (z+0 > 0 && 100 * $col > 98 * z))
		  printf ("%.3f%%", $col/z) ;
	      else if (z + 0 > 10)
		  printf ("%.2f%%", $col/z) ;
	      else if (z + 0 > 1)
		  printf ("%.1f%%", $col/z) ;
	      else
	      {
		  if (fparam [col] == "ZERO")
		      printf ("%.2f%%", $col) ;
		  else
		  {
		      z = param[fparam[col]] ;
		      # printf ("%s*%s:",fparam[col],z);
		      if (z + 0 > 0)
			  printf ("%.2f%%", 100*$col/z) ;
		  }
	      }
	  }
	}
      if (ff == "div")
	{
          if ($col == "NULL")
	      printf("") ;
	  else
	  {
	      z = 0 + fparam [col] ;
	      if (z + 0 == 10)
		  printf ("%.1f", $col/z) ;
	      if (z + 0 == 100)
		  printf ("%.2f", $col/z) ;
	      if (z + 0 == 1000)
		  printf ("%.3f", $col/z) ;
	      if (fparam [col] == "")
	      {
		  z = param[fparam[col]] ;
		  if (z + 0 > 0)
		      printf ("%.2f%%", 100*$col/z) ;
		  else
		      printf("") ; 
	      }
	  }
	}
    }
  printf ("\n") ;

}
