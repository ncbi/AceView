/^MiMxFc/ { 
    g=$2 ; gg[g]=1 ; gMin[g]=$3 ; gMax[g]=$4 ; gFc[g]=$5 ; next ; 
}
{
    if ($5+$6+1 == 0) next ; 
    if ($2+0>0) next ; 
    f = $1 ; ff[f] = 1 ; g = $2 ; gg[g] = 1 ; ln[g] = $3 ; cap[g] = $4 ; 
    g1[g] += $5 ; g2[g] += $6 ; z1[f,g] = $5 ; z2[f,g] = $6 ; 
}
END {
    nf = split("RNA_Total,RNA_PolyA,ILMR3,AGLR2,ROCR2,Nanopore.titr_AGLR2,PacBio2.titr.ccs3_AGLR2,Nanopore.titr_ROCR3,PacBio2.titr.ccs3_ROCR3,ILMR2,ILMR2_lowQ,AGLR1,ROCR1,ILMR1,BSPR1",ff,",")  ; 
    printf ("#Gene\tLength\tMax Index in Total\tMin Index in Total\tFold Change\tCapture\tTruth\tSum B>A\tSum A>B\tInconsistent\tSum B>A no capture\tSum A>B no capture") ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s B>A",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s A>B",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s inconsistent",ff[i]) ; 

    for(g in gg)
    {
	bad = "" ; trueg = "" ; 
	gnc1 = z1[ff[1],g] + z1[ff[2],g] ; 
	gnc2 = z2[ff[1],g] + z2[ff[2],g] ; 
	gc1 = g1[g] - gnc1 ; 
	gc2 = g2[g] - gnc2 ; 

	if (g1[g]*g2[g]>0)
	    bad = "Inconsistent" ; 
	
	if (gnc1*gnc2>0)
	{
	    if (gnc1 > gnc2)
		trueg = "TrueBA_b"int((g1[g]+50)/100)"_a"int((g2[g]+50)/100) ; 
	    else
		trueg = "TrueAB_a"int((g2[g]+50)/100)"_b"int((g1[g]+50)/100) ; 
	}
	else if (gnc1>0)
	{
	    if (5*gc2<g1[g])
	    {
		trueg = "TrueB_" int((g1[g]+50)/100) ; if (gc2>0) trueg = trueg "_a" int((gc2+50)/100) ; 
	    }
	    else
	    {
		trueg  =  "IncB_"int((g1[g]+50)/100)"_a"int((g2[g]+50)/100) ; 
	    }
	}
	else if (gnc2>0)
	{
	    if (5*gc1<g2[g])
	    {
		trueg = "TrueA_" int((g2[g]+50)/100) ; if (gc1>0) trueg = trueg "_b" int((gcc+50)/100) ; 
	    }
	    else 
	    {
		trueg  =  "IncA_"int((g2[g]+50)/100)"_b"int((g1[g]+50)/100) ; 
	    }
	    
	}
	else
	{
	    if (gc1 > 0)
	    {
		if (gc1 > 5*gc2)
		{
		    trueg  =  "NewB_"int((gc1+50)/100) ; if (gc2 > 0) trueg = trueg "_a" int((gc2+50)/100) ;
		}
		else if (gc2 > 5*gc1)
		{
		    trueg  =  "NewA_"int((gc2+50)/100) ; if (gc2 > 0) trueg = trueg "_b" int((gc1+50)/100) ;
		}
		else if (gc1 >= gc2)
		{		
		    trueg  =  "NewBA_b"int((gc1+50)/100)" "_a" "int((gc2+50)/100) ; 
		}
		else 
		{		
		    trueg  =  "NewAB_a"int((gc2+50)/100)" "_b" "int((gc1+50)/100) ; 
		}
	    }
	    else if (gc2 > 0)
	    {
		trueg  =  "NewA_"int((gc2+50)/100) ; 
	    }

	}
	
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
}

