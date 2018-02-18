/^ZZZZZ/{zz++ ; next ;}
{ if (zz < 1) { r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} next ; }}
/^\"/ { 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
  r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}

  if (zz == 2) 
    {
      gsub (/NULL/,"0",$0) ;
      x = $2 ; 
      ss[ir,1,x] = $3 ;sst[ir,1] += $3 ; sstt[1] += 0+$3 ;
      ss[ir,2,x] = $4 ;sst[ir,2] += $4 ; sstt[2] += 0+$4 ;
      ss[ir,3,x] = $5 ;sst[ir,2] += $5 ; sstt[3] += 0+$5 ;
      if (x > xmax) xmax = x ;
    }
  if (zz == 3) 
    {
      hasQual = 1 ;
      gsub (/NULL/,"0",$0) ;
      ii = 0 + substr($2,2) ; if (ii == 0) ii = 1 ; 
      x = $3 ;
      qual[ir,ii,x] = $4 ;
      if (x > xmax) xmax = x ;
    }
  if (zz == 1)
    {
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") nmachine = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ; if ($2 != "NA") nsample = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ; if ($7 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ; if ($7 != "NA") ntissue = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; title[ir]=$6; if ($6 != "NA") ntitle = 1 ;
      gsub(/NULL/,"NA",$8) ; gsub(/ /,"_",$8) ; runid[ir]=$8; if ($8 != "NA") nrunid = 1 ;
      next ;
    }
}
END { 

  printf ("%s", mydate) ;    for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mydate) ;   for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mydate) ;  for (r = 1 ; r <= nr ; r++)printf("\t") ;


  mytitle = "Mismatches observed in uniquely mapped reads, per position in the read, i.e. per sequencing cycle" ;
  printf ("\n%s", mytitle) ;  for (r = 1 ; r <= nr ; r++)printf("\t") ;

  mytitle = "Partition of mismatches, per position in read, i.e. per cycle. Please look for eventual mismatches spikes or defects across runs." ;
  printf ("\t\t%s", mytitle) ; for (r = 1 ; r <= nr ; r++)printf("\t") ;

  mytitle = "Base quality: a posteriori quality per base position, -10 Log10(number of mismatches with the reference per base aligned at each position), 10 is poor, 30 is excellent. This curve is a direct measure of the quality of the bases in the run, and is used in Magic in replacement of the FastQ quality factor. This single small quality file per run and per base position is a more accurate representation of the base quality, it ameliorates the SNV calls and as a side bonus minimizes CPU, memory and disk usage. Note that the third table, which gives the quality, is not the log of the second table, which is a profile normalized at 100%%." ;
  printf ("\t\t%s", mytitle) ; for (r = 1 ; r <= nr ; r++)printf("\t") ;

  printf ("\nRun") ;
  for (r = 1 ; r <= nr ; r++)
    {
      printf("\t%s", rr[r]) ;
    }
  printf ("\t\tRun") ;
  for (r = 1 ; r <= nr ; r++)
    {
      printf("\t%s", rr[r]) ;
    }

  if (hasQual)
    {
      printf ("\t\tRun") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", rr[r]) ;
	}
    }

  if (nrunid)
    {
      printf ("\nRunId") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", runid[r]) ;
	}
      for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tRunId") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", runid[r]) ;
	    }
	}
    }

  if (nmachine)
    {
      printf ("\nMachine") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", machine[r]) ;
	}
      for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tMachine") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", machine[r]) ;
	    }
	}
    }

  if (nsample)
    {
      printf ("\nSample") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", sample[r]) ;
	}
     for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tSample") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", sample[r]) ;
	    }
	}
    }
  
  if (nsystm)
    {
      printf ("\nSystem") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", systm[r]) ;
	}
     for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tSystem") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", systm[r]) ;
	    }
	}
    }

  if (ntissue)
    {
      printf ("\nTissue") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", tissue[r]) ;
	}
     for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tTissue") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", tissue[r]) ;
	    }
	}
    }

  if (ntitle)
    {
      printf ("\nTitle") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\t%s", title[r]) ;
	}
      for (ii = 0 ; ii <= 1 + hasQual ; ii++)
	{
	  printf ("\t\tTitle") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      printf("\t%s", title[r]) ;
	    }
	}
    }

  printf ("\nPosition") ;
  for (r = 1 ; r <= nr ; r++)
    {
      printf("\t%s", rr[r]) ;
    }
  printf ("\t\tPosition") ;
  for (r = 1 ; r <= nr ; r++)
    {
      printf("\t%s", rr[r]) ;
    }
  if (hasQual)
    {
      printf ("\t\tPosition") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  printf("\tA posteriori quality") ;
	}
    }

  for (ii = 1 ; ii < 3 ; ii++)
    {
      if (sstt[ii] < 1) continue ;
      
      if (ii > 1) printf ("\n\n\n\n") ;
      for (x = 1 ; x <= 1000 && x < xmax + 8 ; x++)
	{
	  printf ("\n%d", x) 
	    for (r = 1 ; r <= nr ; r++)
	      {
		printf("\t%d", 0+ss[r,ii,x]) ;
	      }
	  printf ("\t\t%d", x) ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      if (sst[r,ii] > 0)
		printf("\t%.2f", 100.0 * ss[r,ii,x]/sst[r,ii]) ;
	      else
		printf ("\t0") ;
	    }
	  if (hasQual)
	    {
	      printf ("\t\t%d", x) ;
	      for (r = 1 ; r <= nr ; r++)
		{
		  if (qual[r,ii,x] > 0)
		    printf("\t%d", qual[r,ii,x]) ;
		  else
		    printf ("\t0") ;
		}
	    }
	}
    }
  printf ("\nTotal", x) 
    for (r = 1 ; r <= nr ; r++)
      {
	printf("\t%d", 0+sst[r,1]+sst[r,2]+sst[r,3]) ;
      }
  printf ("\t\t%d", x) ;
  for (r = 1 ; r <= nr ; r++)
    {
      if (0 + sstt[1]+sstt[2]+sstt[3] > 0)
	printf("\t%.2f", 100.0 * (0+sst[r,1]+sst[r,2]+sst[r,3])/ (0+sstt[1]+sstt[2]+sstt[3])) ;
      else
	printf ("\t0") ;
    }
  printf ("\n") ;
}

