/^MiMxFc/ { 
    g=$2 ; gg[g]=1 ; gMin[g]=$4 ; gMax[g]=$3 ; gFc[g]=$5 ; next ; 
}
{
    if ($5+$6+1 == 0) next ; 
    if ($2+0>0) next ; 
    f = $1 ; ff[f] = 1 ; g = $2 ; gg[g] = 1 ; ln[g] = $3 ; cap[g] = $4 ; 
    g1[g] += $5 ; g2[g] += $6 ; z1[f,g] = $5 ; z2[f,g] = $6 ; 
}
END {
    nf = split("RNA_Total,RNA_PolyA,ILMR3,AGLR2,ROCR2,Nanopore.titr_AGLR2,PacBio2.titr.ccs3_AGLR2,Nanopore.titr_ROCR3,PacBio2.titr.ccs3_ROCR3,ILMR2,ILMR2_lowQ,AGLR1,ROCR1,ILMR1,BSPR1",ff,",")  ; 
    split("Total,PolyA,I3,A2,R2,A2,A2,R3,R3,I2,I2,A1,R1,I1,B1",fCap,",")  ; 
    printf ("#Gene\tLength\tMax Index in Total\tMin Index in Total\tFold Change\tCapture\tTruth\tSum B>A\tSum A>B\tInconsistent\tSum B>A no capture\tSum A>B no capture") ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s B>A",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s A>B",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s inconsistent",ff[i]) ; 

    for(g in gg)
    {
	bad = "" ; trueg = "non-DEG" ; 
	gnc1 = z1[ff[1],g] + z1[ff[2],g] ; 
	gnc2 = z2[ff[1],g] + z2[ff[2],g] ; 
	gc1 = g1[g] - gnc1 ; 
	gc2 = g2[g] - gnc2 ; 

	if (g1[g]*g2[g]>0)
	    bad = "Inconsistent" ; 
	
	if (gnc1*gnc2>0)
	{
	    if (gnc1 > gnc2)
	    {
		if ((gnc1 > 400 && gnc1 > 2 * gnc2) || (gc1 > 500 && gc1 > 3 * gc2))
		    trueg = "TrueB_" ;
		else
		    trueg = "TrueBA_" ;
		trueg = trueg  int((g1[g]+50)/100) "_a" int((g2[g]+50)/100) ; 
	    }
	    else
	    {
		if ((gnc2 > 400 && gnc2 > 2 * gnc1) || (gc2 > 500 && gc2 > 3 * gc1))
		    trueg = "TrueA_" ;
		else
		    trueg = "TrueAB_a" ;
		trueg = trueg int((g2[g]+50)/100) "_b" int((g1[g]+50)/100) ; 
	    }
	}
	else if (gnc1>0)
	{
	    if (5*gc2<g1[g])
	    {
		trueg = "TrueB_" ;
		trueg  = trueg int((g1[g]+50)/100) ; 
		if (gc2>0) trueg = trueg "_a" int((gc2+50)/100) ; 
	    }
	    else
	    {
		trueg  =  "InconsistentAB_a"int((g2[g]+50)/100)"_b"int((g1[g]+50)/100) ; 
	    }
	}
	else if (gnc2>0)
	{
	    if (5*gc1<g2[g])
	    {
		trueg = "TrueA_" ;
		trueg = trueg  int((g2[g]+50)/100) ; 
		if (gc1>0) trueg = trueg "_b" int((gcc+50)/100) ; 
	    }
	    else 
	    {
		trueg  =  "InconsistentAB_a"int((g2[g]+50)/100)"_b"int((g1[g]+50)/100) ; 
	    }
	    
	}
	else
	{
	    if (gc1 > 0)
	    {
		if (gc1 > 2*gc2)
		{
		    if (gc1 >= 500)
			trueg  =  "NewTrueB_" ;
		    else if (gc1 >= 500)
			trueg  =  "NewBprobable_" ;
		    else if (gc1 >= 0)
			trueg  =  "NewDubiousB_" ;
		    else
			trueg  =  "NewBignore_" ;

		    trueg  = trueg int((gc1+50)/100) ; 
		    if (gc2 > 0) trueg = trueg "_a" int((gc2+50)/100) ;
		}
		else if (gc2 > 2*gc1)
		{
		    if (gc2 >= 500)
			trueg  =  "NewTrueA_" ;
		    else if (gc2 >= 500)
			trueg  =  "NewAprobable_" ;
		    else if (gc2 >= 0)
			trueg  =  "NewDubiousA_" ;
		    else
			trueg  =  "NewAignore_" ;
		    trueg  = trueg int((gc2+50)/100) ; 
		    if (gc1 > 0) trueg = trueg "_b" int((gc1+50)/100) ;
		}
		else if (gc1 >= gc2)
		{		
		    trueg  =  "InconsistentAB_a"int((gc2+50)/100) "_b" int((gc1+50)/100) ; 
		}
		else 
		{		
		    trueg  =  "InconsistentAB_a"int((gc2+50)/100) "_b" int((gc1+50)/100) ; 
		}
	    }
	    else if (gc2 > 0)
	    {
		if (gc2 >= 500)
		    trueg  =  "NewTrueA_" ;
		else if (gc2 >= 500)
		    trueg  =  "NewAprobable_" ;
		else if (gc2 >= 0)
		    trueg  =  "NewDubiousA_" ;
		else
		    trueg  =  "NewAignore_" ;
		trueg  = trueg int((gc2+50)/100) ; 
	    }

	}

	split (trueg, aa, "_") ;
	ttt[aa[1]] += 1 ;
	gTrue[g] = aa[1] ;

	printf ("\n%s\t%d\t%s\t%s\t%s\t%s\t%s\t%.0f\t%.0f\t%s",g,ln[g],gMax[g],gMin[g],gFc[g],cap[g],trueg,gc1,gc2,bad) ; 
	printf ("\t%d\t%d", gnc1, gnc2) ; 
	for (i = 1 ; i <= nf ; i++)
	    printf ("\t%s",z1[ff[i],g]) ; 
	for (i = 1 ; i <= nf ; i++)
	    printf ("\t%s",z2[ff[i],g]) ; 
	for (i = 1 ; i <= nf ; i++)
	{
	    bad = "" ; 
	    if (z1[ff[i],g]*z2[ff[i],g]>0)
		bad = "Inconsistent" ; 
	    printf ("\t%s",bad) ; 
	}
    }
    printf ("\n") ; 

######  stats
    for(g in gg)
    {
	t = gTrue[g] ;
	gnc1 = z1[ff[1],g] + z1[ff[2],g] ; 
	gnc2 = z2[ff[1],g] + z2[ff[2],g] ; 
	gc1 = g1[g] - gnc1 ; 
	gc2 = g2[g] - gnc2 ; 
        N1 = g1[g] ; N2 = g2[g] ;

	for (i = 1 ; i <= nf ; i++)
	{
	    isCap = 0 ;
	    if (index(cap[g], fCap[i]) > 0)
		isCap = 1 ;
	    n1 = z1[ff[i],g] ; 
	    n2 = z2[ff[i],g] ; 
	    if (n1+n2 == 0)
		ss[t,i,3,isCap]++ ;          # missed
	    else if (n1*n2 >0)
		ss[t,i,2,isCap]++ ;          # false
	    else if (n1>0 && N1 > N2)
		ss[t,i,1,isCap]++ ;          # true
	    else if (n2 >0 && N2 > N1)
		ss[t,i,1,isCap]++ ;          # true
	    else if (n1>0 && N1 > N2)
		ss[t,i,2,isCap]++ ;          # false
	    else if (n2 >0 && N1 > N2)
		ss[t,i,2,isCap]++ ;          # false
	}
    }
    
    i2tMax = split ("TrueA,TrueB,NewTrueA,NewTrueB,TrueAB,TrueBA,NewAB,NewBA,NewDubiousA,NewDubiousB,InconsistentAB,InconsistentBA,non-DEG", i2t, ",") ;
    i2tMax = split ("True,True,NewTrue,NewTrue,TrueAB,TrueAB,NewAB,NewAB,NewDubious,NewDubious,Inconsistent,Inconsistent,non-DEG", i2t2, ",") ;

    out = outf".deg_truth.txt" ;
    printf ("### Sensitivity and specificity of the captured platforms, in each case the numbers are seen:contradicted:not seen\n") > out ;    
    printf ("# Type\tCumul") > out ;
    for (i = 1 ; i <= nf ; i++)
	printf ("\tTargeted %s", ff[i]) > out ;
    for (i = 1 ; i <= nf ; i++)
	printf ("\tNot targeted %s", ff[i])  > out ;
    split ("true,false,missed", nam,",") ;
    for (j = 1 ; j <= i2tMax ; j++)
    {
	t = i2t [j] ;
	done[t] = 1 ;
	if (ttt[t] < 1) continue ;
	for (k = 1 ; k <= 3 ; k++)
	{
	    printf ("\n%s.%s\t%d", t,nam[k],ttt[t])  > out ;
	    for (i = 1 ; i <= 2 ; i++)
		printf ("\t%d", ss[t,i,k,0])  > out ;
	    for (i = 3 ; i <= nf ; i++)
		printf ("\t%d", ss[t,i,k,1])  > out ;
	    printf ("\t") ;
	    for (i = 1 ; i <= nf ; i++)
		printf ("\t%d", ss[t,i,k,0])  > out ;
	}
    }
    printf ("\n")  > out ;
    for (j = 1 ; j <= i2tMax ; j+=2)
    {
	t = i2t [j] ;
	t2 = i2t[j+1] ;
	t3 = i2t2[j] ;
	done[t] = 1 ;
	if (ttt[t] < 1) continue ;
	for (k = 1 ; k <= 3 ; k++)
	{
	    printf ("\n%s.%s\t%d", t3,nam[k],ttt[t]+ttt[t2])  > out ;
	    for (i = 1 ; i <= 2 ; i++)
		printf ("\t%d", ss[t,i,k,0]+ss[t2,i,k,0])  > out ;
	    for (i = 3 ; i <= nf ; i++)
		printf ("\t%d", ss[t,i,k,1]+ss[t2,i,k,1])  > out ; 
	    printf ("\t") ;
	    for (i = 1 ; i <= nf ; i++)
		printf ("\t%d",  ss[t,i,k,0]+ss[t2,i,k,0])  > out ;
	}
    }
    printf ("\n")  > out ;

    for (t in ttt)
	if (ttt[t] + 0 > 0 && done[t] + 0 < 1) 
	    printf ("\t%s %d \n", t, ttt[t])   > out ;


################################
    out = outf ".deg_truth.histos.txt" ; 
    printf ("### Histograms per type\n# type") > out ; 

    for (g in gg)
    {
	t = gTrue [g] ;
	i = int(gMax[g] + .5) ; wMax[t,i] += 1 ;
	i = int(gMin[g]+ .5) ; wMin[t,i] += 1 ;
	i = int((gMax[g]-gMin[g])*5 + .5) ;if (i>30)i=30; wFc[t,i] += 1 ;
	# printf ("%s %s %s %d\n", g, gMax[g], gMin[g], i) > out
    }
    

    printf("\tx") > out ;
    for (j = 1 ; j <= i2tMax ; j++)
	printf ("\t%s Max ", i2t[j])  > out;
    printf("\t\t")  > out ;
    for (j = 1 ; j <= i2tMax ; j++)
	printf ("\t%s Min ", i2t[j])  > out ;
    printf("\t\t")  > out ;
    for (j = 1 ; j <= i2tMax ; j++)
	printf ("\t%s FC ", i2t[j])  > out ;


    for (i = 0 ; i <= 30 ; i++)
    {
	printf ("\n\t%d", i)  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    printf ("\t%d", wMax[i2t[j], i])  > out ;
	printf ("\t\t%d", i)  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    printf ("\t%d", wMin[i2t[j], i])  > out ;
	printf ("\t\t%.1f", 2**(i/5.0))  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    printf ("\t%d", wFc[i2t[j], i])  > out ;
    }
    printf ("\n\t%d", i) > out ;

}
