
BEGIN { mm["LIF_S"]=0;}
/S_brain/{next}
{ gsub (/c_NC_001807.4\|human_mithochondria/,"mito",$0) ;
  mrna = $1; x = $2 ; n = $3 ; t = $4; m = $5 ;
  gsub(/HELdge/,"HEL",m) ; gsub(/S_brain/,"ILM_nS",m) ;
  gsub(/ILM_R/,"ILM_",m) ; if(m~"R454") m="R454" ; 
  gsub(/LIF_R/,"LIF_",m) ; 
  tissue = $6 ;
  gsub(/\*/,"",t);
  gsub(/>/,"2",t);
  gsub(/+/,"p2",t);
  gsub(/-/,"m2",t);
  if (t != type) next ;
  if (x>xmax[mrna]) xmax[mrna] = x ;
  nt[mrna,x,tissue] += n ;
  mmx[mrna,x] += n; 
  nmt[mrna,x,m,tissue] += n;  
  mm[m] += n ;
  mt[m,tissue] += n ;
  tt[tissue] += n ;
  mrnas[mrna]+=n;
}
END {
  j = 1 ;  for (i in mm) { mindex[j]=i ; j++; } nms = asort(mindex) ;

  printf ("Target\tType\tPosition") ;
  for(t in tt) printf("\tAny %s",t) ; 
  for (i = 1 ; i <= nms ; i++) 
    {
      m = mm [mindex[i]] ;
    }
  for (m in mm)
    {
      for(t in tt) if(mt[m,t]>0)printf("\t%s %s",m,t) ; 
    }
  printf ("\tNumber of platforms detecting the change\tTotal number of occurences\n") ;
  
  for(mrna in mrnas)
    {
      for (x=1; x<=xmax[mrna]; x++)
	{
	  z=0;ztt=0; 
	  for(m in mm)
	    {
	      zt=0;
	      for(t in tt)
		{
		  ztt += nmt[mrna,x,m,t] ;
		  mm1[substr(m,1,3)] += nmt[mrna,x,m,t] ;
		}
	    } 
	  for(m in mm1) if (mm1[m]>0) {z++ ; mm1[m] = 0 ;}
	  if (ztt>=minZ && z>=minProtocol)
	    {
	      printf("%s\t%s\t%d",mrna,type2,x) ;
	      for(t in tt) printf("\t%d", nt[mrna,x,t]);
	      for (i = 1 ; i <= nms ; i++) 
		{
		  m = mm [mindex[i]] ;
		}
	      for (m in mm)
		{
		  for(t in tt)
		    if (mt[m,t]>0)
		      printf("\t%d", nmt[mrna,x,m,t]);
		}
	      printf("\t%d",z) ;
	      printf("\t%d\n", mmx[mrna,x]) ;
	    }
	}
    }
}
