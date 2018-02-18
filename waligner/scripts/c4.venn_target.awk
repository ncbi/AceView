/^ZZZZZ/{zz++ ; next ;}

{
  if (zz < 1)
    { 
      r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}
      next ;
    }
}
{ 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
}

{
  if (zz == 1) 
    { 
      r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") nmachine = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ;if ($3 != "NA") nsample = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ;if ($7 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ;if ($5 != "NA") ntissue = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; title[ir]=$6;if ($6 != "NA") ntitle = 1 ;
      gsub(/NULL/,"NA",$8) ; gsub(/ /,"_",$8) ; runid[ir]=$8; if ($8 != "NA") nrunid = 1 ;
      next ;
    }
}

/^RUN/{ r = $2 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;}next;}

{ 
  if (zz == 2)
    { 
      t1 = $2 ;
      j = split($3,aa,";") ;
      z="" ; z1 = "" ;
      for(i=1;i<=j;i++)
	{
	  split(aa[i],bb,"_") ;
	  z = z z1 bb[2] ;
	  z1 = "_" ;
	}

      t2 = $2 "\t" z ; 

      f1[ir,t2] = $1 ;

      f1[-1,t2] += $1 ;
      if (t1 == "Details") 
	{ 
	  ffr[-1] += $1 ;
	  ffr[ir] += $1 ;
	  n0t[z] += $1 ;
	}
      if (t1 == "Singlet") 
	n1t[z] += $1 ;
      if (t1 == "Doublet") 
	n2t[z] += $1 ;
      if (t1 == "Triplet") 
	n3t[z] += $1 ;

    }
  next ;
}
END { 
  printf ("Run\tTitle\tSample\tSinglet\tRun\tAny") ;
  for (t in n1t)
    printf("\tSinglet_%s", t) ;
  printf ("\tDoublet\tRun") ;
  for (t in n2t)
    printf("\tDoublet_%s", t) ;
  printf ("\tTriplet\tRun") ;
  for (t in n3t)
    printf("\tTriplet_%s", t) ;
  printf ("\tDetails\tRun") ;
  for (t in n0t)
    printf("\tJust_%s", t) ;

  rr[-1] = "Any" ; title[-1] = "Any" ; sample[-1] = "-" ;
  for (r = -1 ; r <= nr ; r++)
    {
      if (r == 0) continue ;

      printf("\n%s\t%s\t%s", rr[r], title[r], sample[r]) ;

      printf("\t\t%s", rr[r]) ;
      printf("\t%d", ffr[r]) ;
      for (t in n1t)
	printf("\t%d", f1[r,"Singlet\t" t]) ;

      printf("\t\t%s", rr[r]) ;
      for (t in n2t)
	printf("\t%d", f1[r,"Doublet\t" t]) ;

      printf("\t\t%s", rr[r]) ;
      for (t in n3t)
	printf("\t%d", f1[r,"Triplet\t" t]) ;

      printf("\t\t%s", rr[r]) ;
      for (t in n0t)
	printf("\t%d", f1[r,"Details\t" t]) ;
    }
  printf ("\n") ;
}

