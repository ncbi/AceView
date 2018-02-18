
{ 
  gsub (/c_NC_001807.4\|human_mithochondria/,"mito",$0) ;
  mrna = $1; x = $2 ; n = $3 ; t = $4; m = $5 ; ts=$6;
  gsub(/HELdge/,"HEL",m) ; gsub(/S_brain/,"ILM_nS",m) ;
  gsub(/ILM_R/,"ILM_",m) ; if(m~"R454") m="R454" ;
  tissue = $6 ;
  
  if (0 && t != "a>g") next ;
  if (x < x0 - 6 || x > x0 + 6)next;
  if (x>xmax[mrna]) xmax[mrna] = x ;
  if (!xmin[mrna]) xmin[mrna] = x ;
  if (x<xmin[mrna]) xmin[mrna] = x ;
  zz [mrna,x,m,ts,t] += n;  
  mm[m] += n ; 
  tt[t] = t ;
  tissues[ts] += n ;
  mrnas[mrna]+=n;
}
END {
  printf ("Target\tPosition\tManip\tTissue\ttype\tCount\n") ; 
  
  for(mrna in mrnas)
    for (x=xmin[mrna]; x<=xmax[mrna]; x++)
      for (m in mm)
	for (ts in tissues)
	  for (t in tt)
	    {
	      z = zz[mrna,x,m,ts,t] ;
	      if (z > zmin)
		printf ("%s\t%d\t%s\t%s\t%s\t%d\n", mrna,x,m,ts,t,z);
	    }
}
