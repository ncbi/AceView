{ 
    m = $1; mm[m] = 1 ;
    for (j = 4 ; j <= 6 ; j++)
	zz[m,$2,$3, j] = $j ; 
}
END {
    for (m in mm)
    {
	m1 = m ; 
	gsub (".intron_support", "", m1);
	gsub ("GREG_", "", m1);
	gsub ("HG19", "Human", m1);
	gsub ("PFAL", "Malaria", m1);
    
	split (m1, aa, "/") ;
	method = aa[1] ;
    	split (aa[2], bb, ".") ;
	species = bb[1] ; 
	run = bb[2] bb[3] 
        printf ("%s\t%s\t%s", species, run, method) ;
	printf ("\t1") ;
	for (j = 4 ; j <= 6 ; j++)
	{
	    if (j > 4) printf ("\t\t%s", method) ;
	    gold = zz[m, "TRUE", "GOLD", j] ;
	    tp = zz[m, "TRUE", "RUN", j] ;
	    fp = zz[m, "FALSE", "RUN", j] ;		
	    fn = gold - tp ;
	    p = 0 ; r = 0 ; f = 0 ;
	    if (tp > 0)
	    {
		p = tp/(tp+fp) ; r = tp/(tp+fn) ; f = 2 * p * r / (p+r) ;
	    }
	    printf ("\t%d\t%d\t%d\t%d\t%d", gold, tp+fp, fp,  tp, fn) ;
	    printf ("\t%.4f\t%.4f\t%.4f", p, r, f) ;
	}
	printf ("\n") ;
    }
}


