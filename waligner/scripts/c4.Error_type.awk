/^ZZZZZ/{zz++ ; next ;}
{ if (zz < 1) { r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} next ; }}
/^\"/ { 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
  r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}

  if (zz == 3) 
    {
      nn[ir,1,"Ambiguous"] = $2 ; 
      next ;
    }
  if (zz == 2) 
    { 
      ii = 0 + substr($2,2) ; if (ii == 0) ii = 1 ; 
      ii = 1 ;       # we are not interested in splitting the fragments 
      if (ii > iimax) iimax = ii ;

      r = 0 ; t = $3 ;
      # if (t == "Any" &&  ssc[ir,ii] < 1)  ssc[ir,ii] = $5 ;
      if (t == "Any")  ssc[ir,ii] += $5 ;
      if (index (t, "n") > 0) next ;
      t1 = t ; gsub(/[at*g>c+-]/,"",t1); if (length(t1)>0) next;
      if (substr(t,1,1) == "*")
	{
	  t = substr(t,2) ; r = 1;
	}
      
      isI = 10000 * index (t,"+++")  + 100 * index (t,"++")  + index (t,"+") ;
      isD = 10000 * index (t,"---")  + 100 * index (t,"--")  + index (t,"-") ;
      
      gsub (/a/,"A",t) ;  gsub (/t/,"T",t) ;  gsub (/g/,"G",t) ;  gsub (/c/,"C",t) ;
      gsub (/+++/, "Ins ", t) ; gsub (/++/, "Ins ", t) ; gsub (/+/, "Ins ", t) ;
      gsub (/---/, "Del ", t) ; gsub (/--/, "Del ", t) ; gsub (/-/, "Del ", t) ;
      
      if (r == 1)
	{
	  xxr[ir,ii,t] += $4;
	  xxr[ir,ii, "Any"] += $4 ;  
	  
	  if (isD >= 10000) xxr[ir,ii,"Triple deletion"] += $4 ;
	  else if (isD >= 100) xxr[ir,ii,"Double deletion"] += $4 ;
	  else if (isD >= 1) xxr[ir,ii,"Deletion"] += $4 ;
	  
	  if (isI >= 10000) xxr[ir,ii,"Triple insertion"] += $4 ;
	  else if (isI >= 100) xxr[ir,ii,"Double insertion"] += $4 ;
	  else if (isI >= 1) xxr[ir,ii,"Insertion"] += $4 ;
	  
	  if (isI + isD == 0) 
	    {
	      xxr[ir,ii,"Substitution"] += $4 ;
	      if (t == "A>G" || t == "G>A" || t == "T>C" || t == "C>T")
		xxr[ir,ii,"Transition"] += $4 ;
	      else
		xxr[ir,ii,"Transversion"] += $4 ;
	    }
	}
      
      nn[ir,ii, t] += $4 ;
      nn[ir,ii, "Any"] += $4 ;  
      
      if (isD >= 10000) nn[ir,ii,"Triple deletion"] += $4 ;
      else if (isD >= 100) nn[ir,ii,"Double deletion"] += $4 ;
      else if (isD >= 1) nn[ir,ii,"Deletion"] += $4 ;
      
      if (isI >= 10000) nn[ir,ii,"Triple insertion"] += $4 ;
      else if (isI >= 100) nn[ir,ii,"Double insertion"] += $4 ;
      else if (isI >= 1) nn[ir,ii,"Insertion"] += $4 ;
      
      if (isI + isD == 0)
	{
	  nn[ir,ii,"Substitution"] += $4 ;
	  if (t == "A>G" || t == "G>A" || t == "T>C" || t == "C>T")
	    nn[ir,ii,"Transition"] += $4 ;
	  else
	    nn[ir,ii,"Transversion"] += $4 ;
	}
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

  printf ("%s", mydate) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mydate) ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mydate) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mydate) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;

  printf ("\n%s", mytitle) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mytitle) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mytitle) ;    for (ii = 1 ; ii <= iimax ; ii++)for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\t%s", mytitle) ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;

  printf ("\nAbsolute observed counts") ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\tCounts per kb aligned") ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\tPercent of each type") ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;
  printf ("\t\tPercent of sliding insertion deletion (e.g. in dimer)") ;   for (ii = 1 ; ii <= iimax ; ii++) for (r = 1 ; r <= nr ; r++)printf("\t") ;


  if (1)
    {

      printf ("\nRun") ; if (iimax > 1) printf("/%d",ii)
      for (r = 1 ; r <= nr ; r++)
	{  
	  for (ii = 1 ; ii <= iimax ; ii++)
	    { printf("\t%s", rr[r]) ; if (iimax > 1) printf("/%d",ii) ; }
	}
      printf ("\t\tRun") ; if (iimax > 1) printf("/%d",ii)
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    { printf("\t%s", rr[r]) ; if (iimax > 1) printf("/%d",ii) ; }
	}
      printf ("\t\tRun") ; if (iimax > 1) printf("/%d",ii)
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    { printf("\t%s", rr[r]) ; if (iimax > 1) printf("/%d",ii) ; }
	}
      printf ("\t\tRun") ; if (iimax > 1) printf("/%d",ii)
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    { printf("\t%s", rr[r]) ; if (iimax > 1) printf("/%d",ii) ; }
	}
      
      if (nrunid)
	{
	  printf ("\nRunId") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", runid[r]) ;
	    }
	  for (jj = 0 ; jj <= 1 + hasQual ; jj++)
	    {
	      printf ("\t\tRunId") ;
	      for (r = 1 ; r <= nr ; r++)
		{
		  for (ii = 1 ; ii <= iimax ; ii++)
		    printf("\t%s", runid[r]) ;
		}
	    }
	}

      if (nmachine)
	{
	  printf ("\nMachine") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", machine[r]) ;
	    }
	  printf ("\t\tMachine") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", machine[r]) ;
	    }
	  printf ("\t\tMachine") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", machine[r]) ;
	    }
	  printf ("\t\tMachine") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", machine[r]) ;
	    }
	}
      
      if (nsample)
	{
	  printf ("\nSample") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", sample[r]) ;
	    }
	  printf ("\t\tSample") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", sample[r]) ;
	    }
	  printf ("\t\tSample") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", sample[r]) ;
	    }
	  printf ("\t\tSample") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", sample[r]) ;
	    }
	}
      
      if (nsystm)
	{
	  printf ("\nSystem") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", systm[r]) ;
	    }
	  printf ("\t\tSystem") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", systm[r]) ;
	    }
	  printf ("\t\tSystem") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", systm[r]) ;
	    }
	  printf ("\t\tSystem") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", systm[r]) ;
	    }
	}
      
      if (ntissue)
	{
	  printf ("\nTissue") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", tissue[r]) ;
	    }
	  printf ("\t\tTissue") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", tissue[r]) ;
	    }
	  printf ("\t\tTissue") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", tissue[r]) ;
	    }
	  printf ("\t\tTissue") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", tissue[r]) ;
	    }
	}
      
      
      printf ("\nMb uniquely aligned") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    printf("\t%s", ssc[r,ii]) ;
	}
      printf ("\t\tMb uniquely aligned") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    printf("\t%s", ssc[r,ii]) ;
	}
      printf ("\t\tMb uniquely aligned") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    printf("\t%s", ssc[r,ii]) ;
	}
      printf ("\t\tMb uniquely aligned") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
	    printf("\t%s", ssc[r,ii]) ;
	}
      
      printf ("\nTotal number of mismatches") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
            printf("\t%d", nn[r,ii,"Any"]) ;
	}
      
      printf ("\t\tTotal number of mismatches") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
            printf("\t%d", nn[r,ii,"Any"]) ;
	}
      
      printf ("\t\tTotal number of mismatches") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
            printf("\t%d", nn[r,ii,"Any"]) ;
	}
      
      printf ("\t\tTotal number of mismatches") ;
      for (r = 1 ; r <= nr ; r++)
	{
	  for (ii = 1 ; ii <= iimax ; ii++)
            printf("\t%d", nn[r,ii,"Any"]) ;
	}
      
      if (ntitle)
	{
	  printf ("\nMismatch type") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", title[r]) ;
	    }
	  printf ("\t\tMismatch type") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", title[r]) ;
	    }
	  printf ("\t\tMismatch type") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", title[r]) ;
	    }
	  printf ("\t\tMismatch type") ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%s", title[r]) ;
	    }
	}
      
      
      Types = "Any,Substitution,Transition,Transversion,Insertion,Deletion,Double insertion,Double deletion,Triple insertion,Triple deletion,A>G,T>C,G>A,C>T,A>T,T>A,G>C,C>G,A>C,T>G,G>T,C>A,Ins A,Ins T,Ins G,Ins C,Del A,Del T,Del G,Del C,Ins AA,Ins TT,Ins GG,Ins CC,Ins AG,Ins CT,Ins AC,Ins GT,Ins TG,Ins CA,Ins TC,Ins GA,Ins AT,Ins TA,Ins GC,Ins CG,Del AA,Del TT,Del GG,Del CC,Del AG,Del CT,Del AC,Del GT,Del TG,Del CA,Del TC,Del GA,Del AT,Del TA,Del GC,Del CG,Ins AAA,Ins TTT,Ins GGG,Ins CCC,Ins AAT,Ins ATT,Ins AAG,Ins CTT,Ins AAC,Ins GTT,Ins TTA,Ins TAA,Ins TTG,Ins CAA,Ins TTC,Ins GAA,Ins GGA,Ins TCC,Ins GGT,Ins ACC,Ins GGC,Ins GCC,Ins CCA,Ins TGG,Ins CCT,Ins AGG,Ins CCG,Ins CGG,Ins ATA,Ins TAT,Ins ATG,Ins CAT,Ins ATC,Ins GAT,Ins AGA,Ins TCT,Ins AGT,Ins ACT,Ins AGC,Ins GCT,Ins ACA,Ins TGT,Ins ACG,Ins CGT,Ins TAG,Ins CTA,Ins TAC,Ins GTA,Ins TGA,Ins TCA,Ins TGC,Ins GCA,Ins TCG,Ins CGA,Ins GAG,Ins CTC,Ins GAC,Ins GTC,Ins GTG,Ins CAC,Ins GCG,Ins CGC,Ins CAG,Ins CTG,Del AAA,Del TTT,Del GGG,Del CCC,Del AAT,Del ATT,Del AAG,Del CTT,Del AAC,Del GTT,Del TTA,Del TAA,Del TTG,Del CAA,Del TTC,Del GAA,Del GGA,Del TCC,Del GGT,Del ACC,Del GGC,Del GCC,Del CCA,Del TGG,Del CCT,Del AGG,Del CCG,Del CGG,Del ATA,Del TAT,Del ATG,Del CAT,Del ATC,Del GAT,Del AGA,Del TCT,Del AGT,Del ACT,Del AGC,Del GCT,Del ACA,Del TGT,Del ACG,Del CGT,Del TAG,Del CTA,Del TAC,Del GTA,Del TGA,Del TCA,Del TGC,Del GCA,Del TCG,Del CGA,Del GAG,Del CTC,Del GAC,Del GTC,Del GTG,Del CAC,Del GCG,Del CGC,Del CAG,Del CTG,Ambiguous" ;

	nTypes = split (Types, aaa, ",") ;
	for (i=1;i<=nTypes;i++)
	  {
	    t = aaa[i] ;
	  


# observed numbers
	  printf ("\n%s",t) ;
	  for (r = 1 ; r <= nr ; r++)
	    {
	      for (ii = 1 ; ii <= iimax ; ii++)
		printf("\t%d", nn[r,ii,t]) ;
	    }
	  
# number per megabase/1000 == per kb
	  printf ("\t\t%s",t) ;
	  for (r = 1 ; r <= nr ; r++)
	    for (ii = 1 ; ii <= iimax ; ii++)
	      {
		if (ssc[r,ii] > 0)
		  printf("\t%.5f", .001 * nn[r,ii,t]/ssc[r,ii]) ;
		else
		  printf ("\t0") ;
	      }
	  
# percent among all types
	  printf ("\t\t%s",t) ;
	  for (r = 1 ; r <= nr ; r++)
	    for (ii = 1 ; ii <= iimax ; ii++)
	      {
		if (nn[r,ii,"Any"] > 0)
		  printf("\t%.2f", 100.0 * nn[r,ii,t]/nn[r,ii,"Any"]) ;
		else
		  printf ("\t0") ;
	    }
	  
# percent same as neighbour
	  printf ("\t\t%s",t) ;
	  for (r = 1 ; r <= nr ; r++)
	    for (ii = 1 ; ii <= iimax ; ii++)
	      {
		if (xxr[r,ii,t] > 0)
		  printf("\t%.2f", 100.0 * xxr[r,ii,t]/nn[r,ii,t]) ;
		else
		  printf ("\t0") ;
	      }
	}
      printf ("\n") ;
    }
}

