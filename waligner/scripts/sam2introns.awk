{ seq= $1;
  chrom = $3 ;
  a1 = $4 ;
  flag = $2 ;
  isUp = int(flag/16) %2 ;
  isRead2 = int(flag/128) %2 ;
  isIntronUp = (isUp + isRead2) * (1 - isUp * isRead2) ;
  if (index($0, "XS:A:+") > 0) isIntronUp = 0
  if (index($0, "XS:A:-") > 0) isIntronUp = 1
  cigar = $6 ;
  ln = length (cigar) ;
  dx = 0 ;
  for (i = 1 ; i <= ln ; i++)
  {
      ddx = 0 ;
      for (j = 0 ; i + j < ln ; j++)
      {
	  c = substr (cigar, i+j, 1) ;
	  k = index ("0123456789", c) ;
          if (k >= 1)
	      ddx = 10 * ddx + k - 1 ;
	  else
	      break ;
      }
      i = i + j ;
      op = substr (cigar, i, 1) ;
      if (op == "N")
      {
	  r12 = ">" ;
	  if (isRead2 == 1) r12 = "<" ;
	  if (isIntronUp)
	      printf("%s%s\t%d\t%s\t%d\t%d\t%d\n", seq, r12, flag, chrom, ddx, a1+dx+ddx-1, a1+dx) ; 
	  else
	      printf("%s%s\t%d\t%s\t%d\t%d\t%d\n", seq, r12, flag, chrom, ddx, a1+dx-1, a1+dx+ddx) ;
      }
      else if (op == "S" || op == "I")
	  continue ;
      dx += ddx ; # D M X, move on target
  }
}
   
   
