/^ZZZZZ/{zz++ ; next ; }
{ if (zz < 1) { r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} next ; }}
/^\"/ { 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
  r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} 
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
      gsub(/NULL/,"NA",$11) ; gsub(/ /,"_",$11) ; stitle[ir]=$11; if ($11 != "NA") nstitle = 1 ;
      gsub(/NULL/,"NA",$12) ; gsub(/ /,"_",$12) ; stitle2[ir]=$12; if ($12 != "NA") nstitle2 = 1 ;
      gsub(/NULL/,"NA",$13) ; gsub(/ /,"_",$13) ; stitle[ir]=$13; if ($13 != "NA") notitle = 1 ;
      next ;
    }
  if (zz == 2) { bp[ir] = $3 ; next ; }
  if (zz == 3) { rejected[ir] = $4 ; next ; }
  for (i=3 ; i<=13 ; i++) { ss[ir,i] = $i ; sst[ir] += $i ; }
}
END { 
  printf ("Run") ;
  if (nrunid == 1)  printf ("\tRunId") ;
  if (nmachine == 1)  printf ("\tMachine") ;
  if (nsystm == 1)  printf ("\tSystem") ;
  if (ntissue == 1)  printf ("\tTissue") ;
  if (nsample == 1)  printf ("\tSample") ;
  if (ntitle == 1)  printf ("\tTitle") ;
  if (nstitle == 1)  printf ("\tSorting Title") ;
  if (nstitle2 == 1)  printf ("\tSorting Title 2") ;
  if (notitle == 1)  printf ("\tOther Title") ;

  printf ("\t%d %s", 1, type) ;
  printf ("\t2 antisense genes at 1 genomic site", type) ;

  for (i=4 ; i<=12 ; i++)
    {
      printf ("\t%d %s", i - 2, type) ;
      if (i > 3) printf("s") ;
    }
  printf ("\t10 or more %ss, rejected", type) ;
  printf ("\t\tRun") ; 

  printf ("\t%% %d %s", 1, type) ;
  printf ("\t%% 2 antisense genes at 1 genomic site", type) ;
  for (i=4 ; i<=12 ; i++)
    printf ("\t%% %d %ss", i - 2, type) ;

  printf ("\t%% 10 or more %ss, rejected", type) ;
  for (r = 1 ; r <= nr ; r++)
    {
      printf ("\n%s",  rr[r]) ;
      if (nrunid == 1)  printf ("\t%s", runid[r]) ;
      if (nmachine == 1)  printf ("\t%s", machine[r]) ;
      if (nsystm == 1)    printf ("\t%s", systm[r]) ;
      if (ntissue == 1)   printf ("\t%s", tissue[r]) ;
      if (nsample == 1)   printf ("\t%s", sample[r]) ;
      if (ntitle == 1)    printf ("\t%s",title[r]) ;
      if (nstitle == 1)    printf ("\t%s",stitle[r]) ;
      if (nstitle2 == 1)    printf ("\t%s",stitle2[r]) ;
      if (notitle == 1)    printf ("\t%s",otitle[r]) ;

      printf ("\t%d\t%d", ss[r,3],ss[r,13]) ;
      for (i=4 ; i<=12 ; i++)
	printf ("\t%d", ss[r,i]) ;
      printf("\t%d", rejected[r]) ;
      
      printf ("\t\t%s",  rr[r]) ;

	  if (sst[r] > 0)
	    printf("\t%.2f", 100.0 * ss[r,3]/sst[r]) ;
	  else
	    printf ("\t0") ;

	  if (sst[r] > 0)
	    printf("\t%.2f", 100.0 * ss[r,13]/sst[r]) ;
	  else
	    printf ("\t0") ;

      for (i=4 ; i<=12 ; i++)
      	{
	  if (sst[r] > 0)
	    printf("\t%.2f", 100.0 * ss[r,i]/sst[r]) ;
	  else
	    printf ("\t0") ;
	} 
      if (sst[r] > 0)
	printf("\t%.2f", 100.0 * rejected[r]/sst[r]) ;
      else
	printf ("\t0") ;
    }
  printf ("\n") ;
}

