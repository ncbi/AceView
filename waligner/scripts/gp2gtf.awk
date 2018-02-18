{ s = $1; chrom = $2 ; strand = $3 ; 
  a1 = $4 + 1 ; a2 = $5 ; p1 = $6 + 1 ; p2 = $7 ;
  nn = $8 ;
  split ($9,xx1,",") ; split ($10,xx2,",") ;
  
  ln = 0 ; frame = 0 ;
  for (i = 1 ; i <= nn ; i++)
    {
      if (strand == "+") iX = i ;
      else iX = nn - i + 1 ;

      u1 = xx1[i] + 1 ; u2 = xx2[i] ; 
      printf ("%s\t%s", chrom, method) ;
      printf ("\texon\t%d\t%d\t.\t%s\t.", u1, u2, strand) ;
      printf ("\tgene_id=\"%s\";", s) ;
      printf (" transcript_id=\"%s\";", s) ;
      printf (" exon_number %d;", iX) ;
      printf ("\n") ;

      v1 = p1 ; if (u1 > v1) v1 = u1 ;
      v2 = p2 ; if (u2 < p2) v2 = u2 ;
      
      if (v1 <= v2)
	{
          frame = (99999999 + frame - ln) % 3 ;
          printf ("%s\t%s", chrom, method) ;
          printf ("\tCDS\t%d\t%d\t.\t%s\t%d", v1, v2, strand, frame) ;
          printf ("\tgene_id=\"%s\";", s) ;
          printf (" transcript_id=\"%s\";", s) ;
          printf (" exon_number %d;", iX) ;
          printf ("\n") ;
          ln += v2 - v1 + 1 ; 
	}
    }
} 
