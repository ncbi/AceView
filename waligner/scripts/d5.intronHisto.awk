BEGIN { n=0; 
 for(i = 1 ; i<= 7 ; i++) 
   {
     if(i == 1) 
       {
	 x = 1 ;
	 n++ ; tt[n] = 1 ; 
	 n++ ; tt[n] = 2 ; 
	 n++ ; tt[n] = 3 ; 
	 n++ ; tt[n] = 4 ; 
	 n++ ; tt[n] = 5 ; 
	 n++ ; tt[n] = 6 ; 
	 n++ ; tt[n] = 7 ; 
	 n++ ; tt[n] = 8 ; 
	 n++ ; tt[n] = 9 ; 
        }
     else
       {
	 x = 10*x ;
	 n++ ; tt[n] = x ; 
	 n++ ; tt[n] = 2*x ; 
	 n++ ; tt[n] = 5*x ; 
       }
   }
 tMax = n ;

 nam[0] = "Any" ; 
 nam[1] = "New" ;
 nam[2] = "RefSeq-2012" ;
 nam[3] = "AceView-2010" ;
 nam[4] = "Encode.37.70-2013" ;
 nam[5] = "Magic-NB-2013" ;
 maxNam = 5 ;

 printf ("Run\tRun\tAligned reads\tRead Cumul\tAligned kb\tkb cumul\tspanning reads") ;
 for (j = 0 ; j <= maxNam ; j++)
   {
     printf ("\t%s spanning reads", nam[j]) ;
     for (t = 1 ; t <= tMax ; t++)
       printf ("\t%s exon junctions dicovered in at least %d reads", nam[j], tt[t]) ;
   }
 printf("\tIntrons") ;
 for (j = 0 ; j <= maxNam ; j++)
   printf("\tTitrating_A>B in %s\tTitrating_A<B in %s", nam[j], nam[j]) ;

 printf ("\n") ;
 n1 = 0 ;
}

/^ZZZZZ/ {zz++; next ; }
{ if (zz < 4)
  {
# expect RefSeq ZZZZZ AceView ZZZZZ EnsemblIntrons
    k = 1 ; if (zz == 1) k = 2 ; if (zz == 2) k = 4 ; if (zz == 3) k = 8 ;
    z = $1 ":" 0+$2 ":" 0+$3 ;
    if (donor == 1) z = $1 ":" 0+$2 ;
    if (donor == 2) z = $1 ":" 0+$3 ;
    if (fuzzy > 0)  { _x3 = (0 + $2) % fuzzy ; _x4 = (0 + $3) % fuzzy ; z = $2 ":" _x3 ":" _x4 ; }
    isKnown[z] += k ;
    # printf("Adding %s\n",z);
  }
}

function exportTitration(out,p) {
  ntrp = 0 ; ntrm = 0 ; nii = 0 ;
  for (i = 0 ; i<= maxNam ; i++)
    { ttrp[i] = 0 ;  ttrm[i] = 0 ; }

  for (z in isKnown) 
    {
      if (nnj[0,z]>0)
	{
	  tr=0 ;nii++ ;
	  za = nnj["A",z] ;
	  zb = nnj["B",z] ;
	  zc = nnj["C",z] ;
	  zd = nnj["D",z] ;

	   if (za + zb + zc + zd >= 5) 
	     {
	       za = nnj["A",z]/(1+acdbReads["A"])  ;
	       zb = nnj["B",z]/(1+acdbReads["B"])  ;
	       zc = nnj["C",z]/(1+acdbReads["C"])  ;
	       zd = nnj["D",z]/(1+acdbReads["D"])  ;
	       
	      if ((za > 1.05 * zc && zc > 1.05 * zd && zd >= zb) ||  (za >  1.05 * zc && zc >= zd && zd >  1.05 * zb) || (za >= zc && zc >   1.05 * zd && zd >   1.05 * zb))
		 { tr = 1 ; ntrp++ ; }
	      if ((za < .95 * zc && zc < .95 * zd && zd <= zb) ||  (za <  .95 * zc && zc <= zd && zd <  .95 * zb) || (za <= zc && zc <   .95 * zd && zd <  .95 * zb))
		 { tr = -1 ; ntrm++ ; }
	     }
	  if (p > 0) 
	    {
	      printf ("%s", z) > out ;
	      for (j = 0 ; j <= maxNam ; j++)
		printf ("\t%s:%s:%d", nam[j],tissue,nnj[j,z]) > out ;
	      
	      printf ("\tA;A:%d", nnj["A",z]) > out ;
	      printf ("\tC;C:%d", nnj["C",z]) > out ;
	      printf ("\tD;D:%d", nnj["D",z]) > out ;
	      printf ("\tB;B:%d", nnj["B",z]) > out ;
	      if (tr > 0) printf("\tTritrating_plus") > out ;
	      if (tr < 0) printf("\tTritrating_minus") > out ;
	      if (tr < 0) printf("\tFluctuating") > out ;
	    }			    
	  k = isKnown[z] ; if (tr > 0) ttrp[0]++ ;  if (tr < 0) ttrm[0]++ ; 
	  if (k % 8 == 0)
	    { if (p > 0) printf("\tNew")  > out ; if (tr > 0) ttrp[1]++ ;  if (tr < 0) ttrm[1]++ ; } 
	  else
	    if (p > 0) printf("\t") > out ;
	  
	  if ((int(k) % 2) == 1)
	    { if (p > 0) printf("\t%s", nam[2])  > out; if (tr > 0) ttrp[2]++ ;  if (tr < 0) ttrm[2]++ ; } 
	  else
	    if (p > 0) printf("\t") > out ;
	  
	  if ((int(k/2) % 2) == 1)
	    { if (p > 0) printf("\t%s", nam[3])  > out; if (tr > 0) ttrp[3]++ ;  if (tr < 0) ttrm[3]++ ; } 
	  else
	    if (p > 0) printf("\t")  > out;
	  
	  if ((int(k/4) % 2) == 1)
	    { if (p > 0) printf("\t%s", nam[4]) > out ; if (tr > 0) ttrp[4]++ ;  if (tr < 0) ttrm[4]++ ; } 
	  else
	    if (p > 0) printf("\t")  > out;
	  
	  if ((int(k/8) % 2) == 1)
	    { if (p > 0) printf("\t%s", nam[5]) > out ; if (tr > 0) ttrp[5]++ ;  if (tr < 0) ttrm[5]++ ; } 
	  else
	    if (p > 0) printf("\t")  > out;
	  
	  if (p > 0) printf ("\n") > out ;
	}
    }
  if (p == 0)
    {
      printf("\t%d",nii) ;
      for (i = 0 ; i<= maxNam ; i++)
	printf ("\t%d\t%d", ttrp[i], ttrm[i]) ;
    }
}

function addUp(j,z,n) {
  nnj[j,z] +=  n ;
  sjReads[j] += n ;
  x = nnj[j,z] ;
  for (t = 1 ; t <= tMax ; t++)
    if (x - n < tt[t] &&  x >= tt[t]) 
      { 
	nnnj[j,t]++ ;
	# print "### ",j,z,n,t,nnnj[j,t] ;
      }
}
/^RUN/ { 
  n1 = 0 ;
  run = $2 ; sample = $3 ; ali = $4 ; ali2 += ali ; kb = $5 ; kb2 += kb ;

  # print "\n\nreading from " inF ;
  for (pp=0 ; pp<1 ; pp++)
    {
      if (substr(sample,1,2) == "NB")
	inF = "NB/tmp/OR/" run "d1." run ".de_uno.txt" ;
      else
	{
	  inF = "tmp/OR/" run "d1." run ".de_uno.txt" ;
	  split (sample, aa, "_") ;
	  acdb = aa[3] ;
	}
      while ( (getline < inF) > 0)
	{
# print $0 ;
	  if (NF >= 7 &&  (gtag < 1 || (gtag==10 && $5 == "gt_ag") || (gtag==1 && $5 != "gt_ag") || (gtag==2 && $5 != "gt_ag" && $5 != "gc_ag")))
	    { 
	      n1++ ;
	      if (n1 >= 0)
		{
		  z = $2 ":" $3 ":" $4 ;
		  if (donor == 1) z = $2 ":" $3 ;
		  if (donor == 2) z = $2 ":" $4 ;
		  if (fuzzy > 0)  { _x3 = (0 + $3) % fuzzy ; _x4 = (0 + $4) % fuzzy ; z = $2 ":" _x3 ":" _x4 ; }
		  n = $7 ;
		  if (countOnce == 1) n = 1
		  addUp(0,z,n) ;
		  if (acdb && length(acdb) == 1) addUp(acdb,z,n) ;
		  sReads += n ;
		  acdbReads[acdb] += n ;
		  k = isKnown[z] ; 
		  isKnown[z] = k ; # ensure that all introns are listed in isKnown
# printf("Searching %s : %s\n",z,k); 
		  if (k % 8 == 0)
		    addUp(1,z,n) ;
		  if ((int(k) % 2) == 1)
		    addUp(2,z,n) ;
		  if ((int(k/2) % 2) == 1)
		    addUp(3,z,n) ;
		  if ((int(k/4) % 2) == 1)
		    addUp(4,z,n) ;
		  if ((int(k/8) % 2) == 1)
		    addUp(5,z,n) ;
		}
	    }
	}
    }

  # printf ("Found %d records\n", n1) ;
  if (n1 > 0)
    {
      printf ("%s\t%s\t%.0f\t%.0f\t%d\t%d\t%d",run,sample,ali,ali2,kb,kb2,sReads) ;
      for (j = 0 ; j <= maxNam ; j++)
	{
	  printf ("\t%d", sjReads[j]) ;
	  for (t = 1 ; t <= tMax ; t++)
	    printf ("\t%d", nnnj[j,t]) ;
	}
      if (titration > 0) exportTitration(out, 0) ;
      printf ("\n") ;
    }
  
  close (inF) ;
}
END {  if (titration > 0) exportTitration(out, 1) ; }

