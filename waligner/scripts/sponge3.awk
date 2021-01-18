/^#/{next;}
/^ZZZZZ/{zz++ ; next ;}
{
  if (zz <5)
    { 
      if(length($1)>0)
	{
	  r = $1 ;
	  if (zz == 0) { isSublib[r] = 1 ; next ; }
          if (isSublib[r]) next ;
	  ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;ir2r[nr]=r;ir=nr;}
	  isgroup[ir] = zz ;	  
	}
      next ;
    }
}
{ 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
}

{
  if (zz == 5) 
    { 
      r = $1 ; ir=0+r2ir[r];if(ir<1)next;
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") nmachine = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ;if ($3 != "NA") nsample = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ;if ($7 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ;if ($5 != "NA") ntissue = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; titre[ir]=$6;if ($6 != "NA") ntitle = 1 ;
      gsub(/NULL/,"NA",$11) ; gsub(/ /,"_",$11) ; stitre[ir]=$11;if ($11 != "NA") nsorting_title = 1 ;
      gsub(/NULL/,"NA",$12) ; gsub(/ /,"_",$12) ; stitre2[ir]=$12;if ($12 != "NA") nsorting_title2 = 1 ;
      gsub(/NULL/,"NA",$13) ; gsub(/ /,"_",$13) ; otitre[ir]=$13;if ($13 != "NA") nother_title = 1 ;
      
      gsub(/NULL/,"NA",$8) ; gsub(/ /,"_",$8) ; runid[ir]=$8; if ($8 != "NA") nrunid = 1 ;

      next ;
    }
}
{ 
    r=$1; ir = r2ir[r] ; if (ir < 1) next ; 
 
   if (level+0 > 0 && $5 != level) next ;
   t = $6 ; it = t2it[t] ; if(it<1){itmax++;t2it[t]=itmax;it=itmax; it2t[it]=t;}
   chrom = $2 ; chroms[chrom] = 1 ;
   levels[$5] = 1 ;
   xx[1,ir,it,$5] += $7 ; xx[2,ir,it,$5] += $8 ; xx[3,ir,it,$5] += $9 ;
   if(0 && $6=="Genome")print " === Gege:: ", xx[3,ir,it,$5], "===",ir,ir2r[ir],it,$5,"=== "
}
END {
  jjmax = 0 ;
  jjmax++ ; nam[jjmax] = "Extent on genome" ;
  jjmax++ ; nam[jjmax] = "Above threshold (bp)" ;
  jjmax++ ; nam[jjmax] = "RNA in extent (bp)" ;
  jjmax++ ; nam[jjmax] = "% genome" ;
  jjmax++ ; nam[jjmax] = "% genome above threshold" ;
  jjmax++ ; nam[jjmax] = "% extent above threshold" ;
  jjmax++ ; nam[jjmax] = "% bases absorbed in extent above threshold " ;
  jjmax++ ; nam[jjmax] = "fold coverage of extent" ;
  jjmax++ ; nam[jjmax] = "fold coverage of supported region" ;

  itgenome = t2it["Genome"] ; 
  if (0)
    {
      for (g in t2it)
	printf ("%s %s\n", g, t2it[g]) ;
    }
  printf("#Run\t1 Threshold") ;
  if (nrunid>0) printf("\t2 RunId") ;
  if (nsample>0) printf("\t3 Sample") ;
  if (ntissue>0) printf("\t4 Tissue") ;
  if (ntitle>0) printf("\t5 Title") ;
  if (nsorting_title>0) printf("\t6 Sorting Title") ;
  if (nsorting_title2>0) printf("\t7 Sorting Title2") ;
  if (nother_title>0) printf("\t8 Other Title") ; 
  if (nmachine > 0)  printf ("\t9 Machine") ;
  if (nsystm > 0)printf ("\t9 System") ;
  printf ("\t10 Genome\t11 Mapped bp on genome") ;

  for (it = 1 ; it <= itmax ; it++)
    {
      for (jj = 1 ; jj <= jjmax ; jj++)
	printf ("\t%06d %s %s", 1000*jj + it, it2t[it], nam[jj]) ;
    } 

  for (level in levels)
    {
      tot = 1 ;
      for (ir = 1 ; ir <= nr ; ir++)
	if (tot < xx[1,ir,itgenome,level])
	  tot = xx[1,ir,itgenome,level] ;

      for (ir = 1 ; ir <= nr ; ir++)
	{
	    if (isgroup[ir] < 3) continue ;
	  printf ("\n%s\t%d", ir2r[ir],level) ;
	  
	  if (nrunid>0) printf("\t%s", runid[ir]) ;
	  if (nsample>0) printf("\t%s",sample[ir]) ;
	  if (ntissue>0) printf("\t%s",tissue[ir]) ;
	  if (ntitle>0) printf("\t%s", titre[ir]) ;
	  if (nsorting_title>0) printf("\t%s", stitre[ir]) ;
	  if (nsorting_title2>0) printf("\t%s", stitre2[ir]) ;
	  if (nother_title>0) printf("\t%s", otitre[ir]) ;
	  if (nmachine > 0)  printf ("\t%s", machine[ir]); 
	  if (nsystm > 0)  printf ("\t%s", systm[ir]); 

	  tot2 = xx[3,ir,itgenome,level] ;
	  printf ("\t%d",tot) ;    
	  printf ("\t%d",tot2) ;
	  if (0) print "===",ir,ir2r[ir],itgenome,level,"=== "
	  for (it = 1 ; it <= itmax ; it++)
	    {
	      for (jj = 1 ; jj <= 3 ; jj++)
		printf ("\t%d",  xx[jj,ir,it,level]) ;
	      
	      printf ("\t%.3f",  100*xx[1,ir,it,level]/(1+tot)) ;
	      printf ("\t%.3f",  100*xx[2,ir,it,level]/(1+tot)) ;
	      printf ("\t%.3f",  100*xx[2,ir,it,level]/(1+xx[1,ir,it,level])) ;
	      printf ("\t%.3f",  100*xx[3,ir,it,level]/(1+tot2)) ;
	      printf ("\t%.3f",  xx[3,ir,it,level]/(1+xx[1,ir,it,level])) ;
	      printf ("\t%.3f",  xx[3,ir,it,level]/(1+xx[2,ir,it,level])) ;
	    }
	}
    }
  printf ("\n") ;      
}
