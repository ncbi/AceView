BEGIN { printf ("\t# Reads with best alignment to a single site, in transcripts or genome		Reads aligned in 2 to 9 sites, in transcripts or genome	Reads aligned in mitochondria, rRNA, ERCC at a single site	Reads aligned in mitochondria, rRNA, ERCC in 2 to 9 sites (rare)	Reads aligned to other transcripts, at a single site	Reads aligned to other transcripts, at 2 to 9 sites	Genomic reads aligned at a single site	Genomic reads aligned at 2 to 9 sites") ;
 printf ("\tPartition in single gene\tin 2 to 9 genes\tin single genomic site\tin 2 to 9 genomic sites");
 printf ("\tPartition") ;
}
/^Cumulated_mismatches/{for(i = 2 ; i < 9 ;i++)cumulatedMM[i] += $i;} 
{
  if (tag != $1) next ;

  ln=$2 ; if (BP == 0 && ln > 100) ln = 100 ;
  if(ln>lnMax) lnMax = ln ;
  for(i = 3 ;i<10 ;i++)
     { vCumul[i] += $i ; n[ln,i] += $i ; hCumul[ln] += $i ; cumul += $i ;cumulBp += ln*$i ;}
}
END {
  for (ln = 0 ; ln <= lnMax ; ln++)
    {
      if (ln >= 0 && hCumul[ln] == 0) continue ;
      printf("\n%d\t%d\t%d",ln,hCumul[ln],BP*ln*hCumul[ln]) ;
      anyp = 0 ; if (cumul > 0) anyp = 100 * hCumul[ln]/cumul ;
      
      for (i = 3 ; i<7 ; i++)
	printf ("\t%d",n[ln,i]) ;
      if (hCumul[ln] == 0) hCumul[ln] = 1 ;
      for(i = 3 ;i<7 ;i++)
	printf ("\t%.2f", 100 * n[ln,i]/hCumul[ln]) ;
      for( i = 3 ; i<7 ; i++)
	{
	  u = vCumul[i] ; if (u == 0) u = 1 ;
	  printf ("\t%.2f", 100 * n[ln,i]/u) ;
	}
      printf ("\t%.2f", anyp) ;
    }
  printf ("\nTotal\t%d\t%d", cumul, BP*cumulBp) ;
  for (i = 3 ; i<7 ; i++)
    printf("\t%d",vCumul[i]) ;
  if (cumul == 0) cumul = 1 ;
  for(i = 3 ; i<7 ;i++)
    printf("\t%.2f",100*vCumul[i]/cumul) ;
  for (i = 3 ; i<8 ; i++)
    printf("\t100")  ;
  printf("\n") ;
  if (0 && tag == "Per_cent_mismatch")  # cumulatedMM is utterly wrong
    {
      n1 = cumulatedMM[2] ;
      printf("Total number of missmatches\t") ;
      for (i = 2 ; i<3 ; i++)
	printf ("\t%d", cumulatedMM[i]) ; 
      printf("\n") ;
    }
}
