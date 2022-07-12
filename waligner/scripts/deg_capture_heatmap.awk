/^LnMiMxFc/ { 
    g=$2 ; gg[g]=1 ; g1[g] = 0 ; g2[g] = 0 ; ln[g] = $3 ; gMin[g]=$5 ; gMax[g]=$4 ; gFc[g]=$6 ; next ; 
}

/^Truth/ {
    g = $2 ;
    nTruth++ ; truth[g] = $3 ;
    next ;
}

{
    if ($5+$6+0 == -1) next ; 
    if ($2+0>0) next ; 
    f = $1 ; ff[f] = 1 ; g = $2 ; gg[g] = 1 ; cap[g] = $4 ; 
    g1[g] += $5 ; g2[g] += $6 ; z1[f,g] = $5 ; z2[f,g] = $6 ; 
}
END {
    nf = split("RNA_Total,RNA_PolyA,ILMR3,AGLR2,ROCR2,Nanopore.titr_AGLR2,PacBio2.titr.ccs3_AGLR2,Nanopore.titr_ROCR3,PacBio2.titr.ccs3_ROCR3,ILMR2,ILMR2_lowQ,AGLR1,ROCR1,ILMR1,BSPR1",ff,",")  ; 
    nf = split("RNA_Total,RNA_PolyA,AGLR1,AGLR2,ROCR1,ROCR2,ILMR1,ILMR2_lowQ,ILMR2,Nanopore.titr_AGLR2,PacBio2.titr.ccs3_AGLR2,Nanopore.titr_ROCR3,PacBio2.titr.ccs3_ROCR3,BSPR1,ILMR3",ff,",")  ; 
    split("Total,PolyA,A1,A2,R1,R2,I1,I2,I2,A2,A2,R3,R3,B1,I3",fCap,",")  ; 
    printf ("#Gene\tLength\tMax Index in Total\tMin Index in Total\tFold Change\tCapture\tTruth\tInconcistency\tSum B>A capture\tSum A>B capture\tSum B>A no capture\tSum A>B no capture") ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s B>A",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s A>B",ff[i]) ; 
    for (i = 1 ; i <= nf ; i++)
	printf ("\t%s problem",ff[i]) ; 

    for (g in gg)
    {
	if (length(g) == 0 || g == 0)
	    continue ;

	bad = "" ; 
	ok1 = 0 ;
	ok2 = 0 ;
	okI = 0 ;
	gnc1 = z1[ff[1],g] + z1[ff[2],g] ; 
	gnc2 = z2[ff[1],g] + z2[ff[2],g] ; 
	gc1 = g1[g] - gnc1 ; 
	gc2 = g2[g] - gnc2 ; 
	
	if (0) printf ("\nzzz %s", g) ;
	trueg = "non-DEG" ; 

	if (nTruth < 1)
	{
	    if (g1[g] + g2[g] > 0)
	    {
		if (g2[g] >= 3 *g1[g] && g2[g] >= 0)       
		{
		    if (gnc2 > 0 && gnc1 == 0)
		    {
			trueg = sprintf ("TrueA_%03d",  int((g2[g]+50)/100)) ;
			if (g1[g] > 50)
			    trueg = trueg "_b" int((g1[g]+50)/100) ; 
			ok2 = 1 ;
		    }
		    else if (g2[g] > 350)
		    {   # at least 1 capture platforms or 3 platforms  excluding BSPR1 and ILMR3
			k = 0 ; k2 = 0 ;
			for (i = 1 ; i <= nf - 2 ; i++)
			{
			    if (z2[ff[i],g] > 1 * z1[ff[i],g] && z2[ff[i],g]>0)
			    {
				if (index(cap[g], fCap[i]) > 0)
				    k++ ;
				else
				    k2++ ;
			    }
			}
			if (k >= 1)
			{
			    trueg = sprintf ("NewTrueA_%03d",  int((g2[g]+50)/100)) ;
			    if (g1[g] > 50)
				trueg = trueg "_b" int((g1[g]+50)/100) ; 
			    ok2 = 1 ;
			}
			if (k2 >= 3)
			{
			    trueg = sprintf ("TrueA_%03d",  int((g2[g]+50)/100)) ;
			    if (g1[g] > 50)
				trueg = trueg "_b" int((g1[g]+50)/100) ; 
			    ok2 = 1 ;
			}
		    }
		}
		else if (g1[g] >= 3 *g2[g] && g1[g] >= 0)       
		{
		    if (gnc1 > 0 && gnc2 == 0)
		    {
			trueg = sprintf ("TrueB_%03d",  int((g1[g]+50)/100)) ;
			if (g2[g] > 50)
			    trueg = trueg "_a" int((g2[g]+50)/100) ; 
			ok1 = 1 ;
		    }
		    else if (g1[g] > 350)
		    {   # at least 1 capture platforms or 3 platforms  excluding BSPR1 and ILMR3
			k = 0 ; k2 = 0 ;
			for (i = 1 ; i <= nf - 2 ; i++)
			{
			    if (z1[ff[i],g] > 1 * z2[ff[i],g] && z1[ff[i],g]>0)	
			    {
				if (index(cap[g], fCap[i]) > 0)
				    k = k + 1 ;
				else
				    k2 = k2 + 1 ;
			    }
			    if (0 && 1 + z1[ff[i],g] +  z2[ff[i],g] > 0)
				printf ("######## %s %s %f %f k=%d k2=%d\n" , g, ff[i], z1[ff[i],g], z2[ff[i],g], k, k2) ;
			}
			if (k >= 1)
			{
			    trueg = sprintf ("NewTrueB_%03d", int((g1[g]+50)/100)) ;
			    if (g2[g] > 50)
				trueg = trueg "_a" int((g2[g]+50)/100) ; 
			    ok1 = 1 ;
			}
			if (k2 >= 3)
			{
			    trueg = sprintf ("TrueB_%03d",  int((g1[g]+50)/100)) ;
			    if (g2[g] > 50)
				trueg = trueg "_a" int((g2[g]+50)/100) ; 
			    ok1 = 1 ;
			}
			if (0) printf ("######### %s k=%d k2=%d\n",g,k,k2);
		    }
		}
	    }
	    if (ok1 + ok2 == 0 && g1[g] * g2[g] > 0)
	    { 
		trueg = "Inconsistent-undecidable" ; 	    
		okI = 1 ;
	    }
	}
	else
	{
	    trueg = truth[g] ;
	    split (trueg, aa, "_") ;
	    gTrue[g] = aa[1] ;
	    if (aa[1] == "TrueB" || aa[1] == "NewTrueB") ok1 = 1 ;
	    if (aa[1] == "TrueA" || aa[1] == "NewTrueA") ok2 = 1 ;
	    if (aa[1] == "Inconsistent-undecidable") okI = 1 ;
	}
	    
	split (trueg, aa, "_") ;
	gTrue[g] = aa[1] ;


######################################################
######  stats

	t = gTrue[g] ;
	for (i = 1 ; i <= nf ; i++)
	{
	    bads[i] = "" ;
	    isCap = 0 ;
	    nt = "nt " ;
	    if (index(cap[g], fCap[i]) > 0)
	    { isCap = 1 ; nt = "" ; }
	    isCap2 = isCap ;
	    if (i < 3) { isCap2 = 1 ; nt = "" ; };
	    n1 = z1[ff[i],g] ; 
	    n2 = z2[ff[i],g] ; 
	    if (n1 + n2 == 0)
	    {
		# print "g=" g " t=" t " ok1=" ok1 " ok2=" ok2 " ff=" ff[i] " cap2=" isCap2 ;
		if (ok1 > 0 && isCap2)
		{  t= "B-Missed" ; ss[t,i,isCap]++ ; BMissed[g] = 1 ; bads[i] = t ; dd[g,i,"Missed"] = 1 ; }          # missed
		else if (ok2 > 0 && isCap2)
		{  t= "A-Missed" ; ss[t,i,isCap]++ ; AMissed[g] = 1 ; bads[i] = t ; dd[g,i,"Missed"] = 1 ; }          # missed
	    }
	    else if (ok1 + ok2 && n1 > 3*n2)
	    {
		if (ok1 > 0)
		{  t= gTrue[g] ; ss[t,i,isCap]++ ; bads[i] = nt t " seen" ; dd[g,i, nt "Seen"] = 1 ; }          # true
		else if (ok2 > 0) 
		{ t = "Opposite" ; ss[t,i,isCap]++ ; Opposite[g] = 1 ; bads[i] = t ; dd[g,i, nt "Problem"] = 1 ; }         # false
	    }
	    else if (ok1 + ok2 && n2 > 3*n1)
	    {
		if (ok2 > 0)
		{  t= gTrue[g] ; ss[t,i,isCap]++ ;  bads[i] = nt t " seen" ; dd[g,i, nt "Seen"] = 1 ; }          # true
		else if (ok1 > 0) 
		{ t = "Opposite" ; ss[t,i,isCap]++ ; Opposite[g] = 1 ; bads[i] = t ;  dd[g,i, nt "Problem"] = 1 ; }         # false
	    }
	    else if (ok1 > 0 && isCap2)
	    {  t= "B-Missed" ; ss[t,i,isCap]++ ; BMissed[g] = 1 ; bads[i] = t ; dd[g,i,"Missed"] = 1 ; }          # missed
	    else if (ok2 > 0 && isCap2)
	    {  t= "A-Missed" ; ss[t,i,isCap]++ ; AMissed[g] = 1 ; bads[i] = t ; dd[g,i,"Missed"] = 1 ; }          # missed
	    else if (ok1 + ok2 + okI >= 1 && n1 * n2 > 0)
	    { t = "Internal-inconsistency" ; ss[t,i,isCap]++ ; IntInc[g] = 1 ; bads[i] = t ;  dd[g,i, nt "Problem"] = 1 ; }     
	    else if (okI == 1 && n1 * n2 == 0)
	    { t = "Inconsistent-undecidable" ; ss[t,i,isCap]++ ; bads[i] = t ; }     
	    else
	    { t = "Weak" ; ss[t,i,isCap]++ ; bads[i] = t ; }
	    
	}


	if (Opposite[g] == 1)
	    ttt["Opposite"]++ ;
	if (AMissed[g] == 1)
	    ttt["A-Missed"]++ ;
	if (BMissed[g] == 1)
	    ttt["B-Missed"]++ ;
	if (IntInc[g] == 1)
	    ttt["Internal-inconsistency"]++ ;
	if (ExtInc[g] == 1)
	    ttt["External-inconsistency"]++ ; 
	
	newT = 0 ; bad= "-" ;
	if (IntInc[g] == 1)
	{
	    bad = "Internal-inconsistency" ;
	}
	else if (Opposite[g] == 1)
	{
	    bad = "External-inconsistency" ; 
	}
	else if (okI == 1)
	{ 
	    trueg = "Inconsistent-undecidable" ; 	    
	    bad = trueg ;
	    okI = 1 ;
	    newT = 1 ; 
	}
	else if (ok1 + ok2 == 0 && g1[g] + g2[g] > 0)
	{ 
	    trueg = "Weak" ; 
	    newT = 1 ; 
	}

	if (nTruth < 100 && newT + okI >= 1)
	{
	    gTrue[g] = trueg ;
	    if (g2[g] > 50)
		trueg = trueg "_a" int((g2[g]+50)/100) ; 
	    if (g1[g] > 50)
		trueg = trueg "_b" int((g1[g]+50)/100) ; 
	}
    
	split (trueg, aa, "_") ;
	t = aa[1] ;
	gTrue[g] = aa[1] ;
	ttt[t] ++ ;
    
##################################################

	if (cap[g]=="") cap[g]="-";

	printf ("\n%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s",g,ln[g],gMax[g],gMin[g],gFc[g],cap[g],trueg, bad) ; 

	printf ("\t%d\t%d", gc1, gc2) ; 
	printf ("\t%d\t%d", gnc1, gnc2) ; 
	for (i = 1 ; i <= nf ; i++)
	    printf ("\t%.0f",z1[ff[i],g]) ; 
	for (i = 1 ; i <= nf ; i++)
	    printf ("\t%.0f",z2[ff[i],g]) ; 
	for (i = 1 ; i <= nf ; i++)
	    printf ("\t%s", bads[i]) ; 
    }
    printf ("\n") ;

################################
    out = outf".deg_truth.txt" ;
    i2tMax = split ("TrueA,TrueB,NewTrueA,NewTrueB,A-Missed,B-Missed,Opposite,Internal-inconsistency,Inconsistent-undecidable,Weak,non-DEG", i2t, ",") ;

    printf ("### Sensitivity and specificity of the captured platforms, in each case the numbers are seen:contradicted:not seen\n") > out ;    
    printf ("### File %s : %s\n", outf".deg_truth.txt", strftime())  > out ; 
    printf ("# Type\tCumul") > out ;
    for (i = 1 ; i <= nf ; i++)
	printf ("\tTargeted %s", ff[i]) > out ;
    for (i = 1 ; i <= nf ; i++)
	printf ("\tNot targeted %s", ff[i])  > out ;
    split ("seen,seen opposite,missed", nam,",") ;
    for (j = 1 ; j <= i2tMax ; j++)
    {
	t = i2t [j] ;
	done[t] = 1 ;
	if (ttt[t]+0 < 0) continue ;

	    printf ("\n%s\t%d", t, ttt[t])  > out ;
	    for (i = 1 ; i <= 2 ; i++)   # for the truth columns we export the non-captured values
		printf ("\t%d", ss[t,i,0])  > out ;
	    for (i = 3 ; i <= nf ; i++)
		printf ("\t%d", ss[t,i,1])  > out ;
	    printf ("\t") ;
	    for (i = 1 ; i <= nf ; i++)
		printf ("\t%d", ss[t,i,0])  > out ;

    }
    printf ("\n")  > out ;

    for (t in ttt)
	if (ttt[t] + 0 > 0 && done[t] + 0 < 1) 
	    printf ("\t%s %d \n", t, ttt[t])   > out ;


################################
    out = outf ".deg_truth.histos.txt" ; 
    printf ("### File %s : %s\n", outf".deg_truth.histo.txt", strftime())  > out ; 
    printf ("### Histograms per type\n# type") > out ; 

    for (g in gg)
    {
	t = gTrue [g] ;
	i = int(gMax[g] + .5) ; if (i > 22) i = 22 ;  wMax[t,i] += 1 ;
	i = int(gMin[g]+ .5) ; wMin[t,i] += 1 ;
	i = int(ln[g]/200 + .5) ; if (i>30)i=30;wLn[t,i] += 1 ;
	dz = 2**(gMax[g]-gMin[g]) ; i = int(5*(dz-1) + .5) ; if (i>22)i=22; wFc[t,i] += 1 ;
	# printf ("%s %s %s %d\n", g, gMax[g], gMin[g], i) > out
    }

    printf("\tx") > out ;
    for (j = 1 ; j <= i2tMax ; j++)
	if (ttt[i2t[j]] > 0)
	    printf ("\t%s Max ", i2t[j])  > out;
    printf("\t\t")  > out ;
    if (0)
    {
	for (j = 1 ; j <= i2tMax ; j++)
	    if (ttt[i2t[j]] > 0)
		printf ("\t%s Min ", i2t[j])  > out ;
	printf("\t\t")  > out ;
    }
    for (j = 1 ; j <= i2tMax ; j++)
	if (ttt[i2t[j]] > 0)
	    printf ("\t%s FC ", i2t[j])  > out ;
    printf("\t\t")  > out ;
    for (j = 1 ; j <= i2tMax ; j++)
	if (ttt[i2t[j]] > 0)
	    printf ("\t%s Ln ", i2t[j])  > out ;
    
    for (i = 0 ; i <= 22 ; i++)
    {
	printf ("\n\t%d", i)  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    if (ttt[i2t[j]] > 0)
		printf ("\t%d", wMax[i2t[j], i])  > out ;
	if (0)
	{
	    printf ("\t\t%d", i)  > out ;
	    for (j = 1 ; j <= i2tMax ; j++)
		if (ttt[i2t[j]] > 0)
		    printf ("\t%d", wMin[i2t[j], i])  > out ;
	}
	printf ("\t\t%.2f", 1+i/5.0)  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    if (ttt[i2t[j]] > 0)
		printf ("\t%d", wFc[i2t[j], i])  > out ;
	printf ("\t\t%d", 200 * i)  > out ;
	for (j = 1 ; j <= i2tMax ; j++)
	    if (ttt[i2t[j]] > 0)
		printf ("\t%d", wLn[i2t[j], i])  > out ;
    }
    printf ("\n\t%d", i) > out ;


################################

    d2tMax = split ("Seen,Missed,Problem,nt Seen,nt Problem", d2t,",") ;

    out = outf ".deg_truth.histos_per_platform.txt" ; 
    printf ("### File %s : %s\n", outf".deg_truth.histo.txt", strftime())  > out ; 
    printf ("### Histograms per type per platform\n# type") > out ; 


    for (d = 1 ; d <= d2tMax ; d++)
    {
	t = d2t[d] ;
	for (f = 1 ; f <= nf ; f++)
	    for (g in gg)
	    {
		if (dd[g,f,t] == 1)
		{
		    i = int(gMax[g] + .5) ; if (i > 22) i = 22 ;  wdMax[t,f,i] += 1 ;
		    i = int(ln[g]/200 + .5) ; if (i>30)i=30;wdLn[t,f,i] += 1 ;
		    dz = 2**(gMax[g]-gMin[g]) ; i = int(5*(dz-1) + .5) ; if (i>22)i=22; wdFc[t,f,i] += 1 ;
		}
	    }

    
	printf ("\n### %s", t) > out

	printf("\tx") > out ;
	for (f = 1 ; f <= nf ; f++)
	    printf ("\t%s Max %s", t, ff[f])  > out;
	printf("\t\t")  > out ;
	for (f = 1 ; f <= nf ; f++)
	    printf ("\t%s FC %s", t, ff[f])  > out;
	printf("\t\t")  > out ;
	for (f = 1 ; f <= nf ; f++)
	    printf ("\t%s LN %s", t, ff[f])  > out;
	printf("\t\t")  > out ;
	
	
	for (i = 0 ; i <= 22 ; i++)
	{
	    printf ("\n\t%d", i)  > out ;
	    for (f = 1 ; f <= nf ; f++)
		printf ("\t%d", wdMax[t,f,i]) > out ;


	    printf ("\t\t%.2f", 1+i/5.0)  > out ;
	    for (f = 1 ; f <= nf ; f++)
		printf ("\t%d", wdFc[t,f,i]) > out ;

	    printf ("\t\t%d", 200*i)  > out ;
	    for (f = 1 ; f <= nf ; f++)
		printf ("\t%d", wdLn[t,f,i]) > out ;
	}

	printf ("\n\n\n") > out ;
    }
    
    
    
}

################################
