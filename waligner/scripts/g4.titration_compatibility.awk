BEGIN { minA = 10000; minB = 10000 ;}
/^ZZZZZ/ {zz = 1 ; next ; }
# gather the index in the any-any group
/^Genes/{next;}
/^Gene/ { if (zz < 1)
  { split ($1,aa,"\""); g = aa[2]; gg[g] = 1 ; print g; next ;}
}
/^Group_nU/ {
  i=index ($1,"A1_any_sampleA_any_591runs") ; if (i>0) {split($1,bb," "); g2A[g] = bb[3]; if(bb[3]>1 && bb[3]<minA) minA=bb[3];}
  i=index ($1,"A1_any_sampleB_any_591runs") ; if (i>0) {split($1,bb," "); g2B[g] = bb[3]; if(bb[3]>1 && bb[3]<minB) minB=bb[3];}
  i=index ($1,"B_I_any_A_any_395runs") ; if (i>0) {split($1,bb," "); if (bb[3]>1) g2IA[g] = bb[3]; if (bb[10]>=0) vg2IA[g] = bb[10]; }
  i=index ($1,"B_I_any_B_any_395runs") ; if (i>0) {split($1,bb," "); if (bb[3]>1) g2IB[g] = bb[3]; if (bb[10]>=0) vg2IB[g] = bb[10]; }
  i=index ($1,"B_L_any_A_any_190runs") ; if (i>0) {split($1,bb," "); if (bb[3]>1) g2LA[g] = bb[3]; if (bb[10]>=0) vg2LA[g] = bb[10]; }
  i=index ($1,"B_L_any_B_any_190runs") ; if (i>0) {split($1,bb," "); if (bb[3]>1) g2LB[g] = bb[3]; if (bb[10]>=0) vg2LB[g] = bb[10]; }
  next;
}

{ if (zz < 1) next ; }
# gather the columns corresponding to ILM or LIF, A or B
/^ERCC-/{next}
/^ERCC_/{next}
{ line++; }
{ if (line == 2)
  {
    for(i = 1 ; i <= NF ; i++) 
      {
	n = split ($i, aaa, "_") ;
	p = aaa[1] ; lab = aaa[2] ; a = aaa[3] ; 
	if (p == "B" || p == "D")
	  { p = aaa[2] ; lab = aaa[3] ; a = aaa[4] ; }
	if (p == "A1" && a == "sampleA")
	  { p = "any" ; a = "A" ; }
	if (p == "A1" && a == "sampleB")
	  { p = "any" ; a = "B" ; }
	if ((a == "A" || a == "B"))
	  {
	    i2p[i] = p ; i2a[i] = a ; i2lab[lab] = lab ; pp[p] = 1 ; aa[a] = 1 ; labs[lab] = 1 ;  nok++;
	  }
	i2r[i] = $i ;
	if (i > imax+0) imax = i ;
      }
  #  print "Recognized ",nok, "groups" ;
    next ;
  }
}
# jump the other title lines
{ if (line < 6) next ; }
# memorize the genes
{
 # print $0 ;
  for (i = 2 ; i <= NF ; i++)
    {    
      if (! i2r[i]) continue ;
# print "i=",i,"line=",line,"g=",$i ;
      g = $i ; 
      if (g && length(g)>1)
	{
	  ig[i,g] = 1 ; nig[i]++ ;
	  if (0 && g == "zofyby"){ printf("%s %s -> %s %s\n",g,i2r[i],i2p[i],i2a[i]);}
	  if (i2a[i])
	    { 
	      if (0) print "ZZZ",i,"__",g,"__",i2r[i],i2p[i],i2a[i] ;
	      gg[g] = 1 ;
	      p = i2p[i] ; a = i2a[i] ; g2ap [g,a,p]++ ;
	    }      
	}
    }
# print "line=",line,NF,$4,$12,$NF ;
#  exit (0) ;
}
# analyse the compatibility table
# A1, B1 unmistakable: seens 14 ILM and 6 LIF among the groups 
# A2, B2 harder 4 ILM 2 LIF
# A3, B3 dubious > 4 total
END {
# print "endline=",line;
  minA = int (minA + 1.5) ; minB = int (minB + 1.5) ;
  if (minA < minB) minA = minB ;
  if (minA > minB) minB = minA ;
  printf ("#Run\tCondition\tTitrating genes\tA\tB\tBad\tA_Perfect\tA_Medium\tA_Dubious\tB_Perfect\tB_Medium\tB_Dubious") ;
  for (i = 1 ; i <= imax ; i++)
    {
      if (! i2r[i]) continue ;
      if (index(i2r[i],"/")>0) continue ;
      for (a in aa)
	for (k = 1 ; k <= 3 ; k++)
	  nk[k,a] = 0 ;
      for (g in gg)
	if (ig[i,g])
	  for (a in aa)
	    {
              nI = 0 ; if (g2ap[g,a,"I"]) nI = 0+g2ap[g,a,"I"] ;
              nL = 0 ; if (g2ap[g,a,"L"]) nL = 0+g2ap[g,a,"L"] ;
	      if (nI + nL > 0)
		{
		  if (a == "A") b = "B" ; else b = "A" ;
		  nIb = 0 ; if (g2ap[g,b,"I"]) nIb = 0+g2ap[g,b,"I"] ;
		  nLb = 0 ; if (g2ap[g,b,"L"]) nLb = 0+g2ap[g,b,"L"] ;

		  if (nIb + nLb == 0)
		    {
		      if (nI >= 14 && nL >= 6) { nk[1,a]++ ;  }
		      else if (nI >= 4 && nL >= 2) { nk[2,a]++ ;  }
		      else if (nI + nL >= 5) { nk[3,a]++ ;  }
		    }
		  else if (0 && i2a[i])
		    printf ("ERROR run %s gene %s  nI=%d nL=%d   nIb=%d nLb=%d\n", i2r[i], g, nI, nL, nIb, nLb) ;
		}
	    }


      z = i2r[i] "\t" title "\t" nig[i] ; zu = "" ;
      for (a in aa)
	{
	  n1[a] = 0 ;
	  for (k = 1 ; k <= 3 ; k++)
	    {
	      n = 0 + nk[k,a] ;
	      zu = zu "\t" n ;
	      n1[a] += n ; 
	    }
	  z = z "\t" n1[a] ;
	}
      bad = "" ;
      if (n1["A"] > n1["B"])
	{
	  if (n1["B"] > 0) bad = n1["B"] ; 
	  printf ("\nA\t%s\t%s\t%s", z, bad,zu) ;
	}
      else
	{
	  if (n1["A"] > 0) bad = n1["A"] ; 
	  printf ("\nB\t%s\t%s\t%s", z, bad,zu) ;
	}
    }
  printf ("\n") ;

  printf ("#Gene\tA index in 591 runs\tB index in 591 runs\tMax\tDelta\tA ILM\tB ILM\tA LIF\tB LIF\tsigma A ILM\tsigma B ILM\tsigma A LIF\tsigma B LIF\tA>B\tA>B @ ILM\tA>B @ LIF\tB>A\tB>A @ ILM\tB>A @ LIF\tTitrating A\tTitrating B")  > OGG ;
  for (i = 1 ; i < imax ; i++)
    {
      if (! i2p[i]) continue ;
      printf ("\t%s", i2r[i])  > OGG ;
    }
  for (g in gg)
    {
      a = g2A[g] ; if (a<minA) a = minA ;
      b=g2B[g] ; if (b<minB) b = minB ;
      if (a == minA && b == minB) continue ;
      x = a ; y = b - a ; if(b>x){ x = b  ;} 
      if (a == minA) a = "" ; if (b == minB) b = "" ; 
      printf ("\n%s\t%s\t%s\t%.2f\t%.2f",g,a,b,x,y)   > OGG ;
      printf ("\t") > OGG ; if (g2IA[g]) printf ("%.2f", g2IA[g]) > OGG ;
      printf ("\t") > OGG ; if (g2IB[g]) printf ("%.2f", g2IB[g]) > OGG ;
      printf ("\t") > OGG ; if (g2LA[g]) printf ("%.2f", g2LA[g]) > OGG ;
      printf ("\t") > OGG ; if (g2LB[g]) printf ("%.2f", g2LB[g]) > OGG ;
      printf ("\t") > OGG ; if (vg2IA[g]) printf ("%.2f", vg2IA[g]) > OGG ;
      printf ("\t") > OGG ; if (vg2IB[g]) printf ("%.2f", vg2IB[g]) > OGG ;
      printf ("\t") > OGG ; if (vg2LA[g]) printf ("%.2f", vg2LA[g]) > OGG ;
      printf ("\t") > OGG ; if (vg2LB[g]) printf ("%.2f", vg2LB[g]) > OGG ;

      z = "" ;
      for (a in aa)
	{
	  if (a != "A" && a != "B")
	    continue ;
	  cg = "" ; cgI = "" ; cgL = "" ;
	  nI = 0 ; if (g2ap[g,a,"I"]) nI = 0+g2ap[g,a,"I"] ;
	  nL = 0 ; if (g2ap[g,a,"L"]) nL = 0+g2ap[g,a,"L"] ;
	  
	  if (a == "A") b = "B" ; else b = "A" ;
	  nIb = 0 ; if (g2ap[g,b,"I"]) nIb = 0+g2ap[g,b,"I"] ;
	  nLb = 0 ; if (g2ap[g,b,"L"]) nLb = 0+g2ap[g,b,"L"] ;

	  if (nIb + nLb + nI + nL < 6)
	    cg = "G" ;
	  else
	    {
	      if (nI + nL > 0)
		{
		  
		  if (nIb + nLb == 0)
		    {
		      if (nI >= 14 && nL >= 6)     cg = "A" ;
		      else if (nI >= 4 && nL >= 2) cg = "B" ;
		      else if (nI + nL >= 5)       cg = "C" ;
		    }
		  else
		    {
		      if (nIb + nLb >= 10 * (nI + nL)) ;
		      else
			{
			  if (nI + nL >= 10 * (nIb + nLb)) cg = "D" ;
			  else 
			    {
			      if (nI + nL > 5 &&  nIb + nLb > 5) cg = "E" ;
			      else 
				cg = "F" ;
			    }
			}
		    }
		}
	    }
	  if (nI + nIb < 3)
	    cgI = "GI" ;
	  else
	    {
	      if (nIb == 0)
		{
		  if (nI >= 14)     cgI = "AI" ;
		  else if (nI >= 5) cgI = "BI" ;
		  else if (nI >= 3) cgI = "CI" ;
		}
	      else
		{
		  if (nIb >= 10 * (nI)) ; 
		  else
		    {
		      if (nI >= 10 * (nIb)) cgI = "DI" ;
		      else 
			{
			  if (nI >= 3 &&  nIb >= 3) cgI = "EI" ;
			  else 
			    cgI = "FI" ;
			}
		    }
		}
	    }
	
	  if (nLb + nL < 3)
	    cgL = "GL" ;
	  else
	    {
	      if (nLb == 0)
		{
		  if (nL >= 6)      cgL = "AL" ;
		  else if (nL >= 4) cgL = "BL" ;
		  else if (nL >= 3) cgL = "CL" ;
		}
	      else
		{
		  if (nLb >= 10 * (nL)) ;
		  else
		    {
		      if (nL >= 10 * (nLb)) cgL = "DL" ;
		      else 
			{
			  if (nL >= 3 &&  nLb >= 3) cgL = "EL" ;
			  else 
			    cgL = "FL" ;
			}
		    }
		}
	    }

	  z = z "\t" cg "\t" cgI "\t" cgL ;
	}
      printf ("%s",z)  > OGG ;  
      z = "" ;
      for (a in aa) uu[a] = 0 ;
      for (i = 1 ; i < imax ; i++)
	{
	  if (! i2p[i]) continue ;
	  x = "" ;
	  if (ig[i,g]) { x = 1 ; uu[i2a[i]]++ ;}
	  z = z "\t" x ;
	}

	printf ("\t%d", uu["A"])  > OGG ;  
	printf ("\t%d", uu["B"])  > OGG ;  
      printf ("%s", z)  > OGG ;  
    }
}
  
 


