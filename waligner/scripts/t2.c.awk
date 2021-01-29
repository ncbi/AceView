{
    v = $1 ;  k=split(v,aa,":") ; ok = 0 ;
    if(length($4)<1) next ;
    if (k == 4 && index(aa[3], "Del_") > 0)  { ok = 1 ; vv[v] = 1 ; split (aa[3],bb,"_") ; x1[v] = aa[2]+1 ; if (forceX1 > 0) x1[v] = forceX1 ; dx[v] = bb[2] ; x2[v] = x1[v] + dx[v] - 1 ;nam[v]=x1[v]-1"Del_"dx[v];delins="dele";}
    if (k == 4 && length(aa[3]) > length(aa[4]) && length(aa[4])==1) { ok=1 ;  vv[v] = 1 ; x1[v] = aa[2] + 1 ; dx[v] = length(aa[3]) - length(aa[4]) ; x2[v] = x1[v] + dx[v] - 1 ; nam[v]=x1[v]-1"Del_"dx[v]":"aa[3]":"aa[4];delins="dele";}
    if (k == 4 && length(aa[3]) < length(aa[4]) && length(aa[3])==1) { ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; dx[v] = length(aa[4]) - length(aa[3]) ; x2[v] = x1[v] + 1 ; nam[v]=x1[v]"Ins_"dx[v]":"aa[3]":"aa[4]; delins="inser";}
    if (k == 4 && type>=100 && type < 110 && length(aa[3]) == 1 && length(aa[4])==1) { ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; dx[v] = 0 ; x2[v] = x1[v] ; nam[v]=aa[3] x1[v] aa[4]; delins="substitu";}
    if (ok == 0) next ;
    st = $2 ; t = $3 ; r = $4 ; i = r2i[r] ; if (i<1) { iMax++; i = iMax ; i2r[i] = r; r2i[r] = i ;} i2st[i] = st ;i2t[i] = t ;
    if (iStart < 1 && substr(t,1,1) != "A") iStart = i ;  # groups are supposed to start with A so they come first at this stage of the hack 
    m[v,i] = $5 ; w1[v,i] = $6 ; w2[v,i] = $7 ; nw[i] ++ ; ww1[i] += $6 ; ww2[i] += $7 ; mm[i] += $5 ; if (i == 1)mmv[v] = $5 ;
}
END {
    for (i = 1 ; i <= iMax ; i++) 
    { # take the average at the common donor site  
	if (forceX1>0) { if (nw[i] < 1) nw[i] = 1; ww1[i] /= nw[i] ;  }
	denom[i] = ww1[i] + mm[i] ; if (denom[i] == 0) denom[i] = 1 ;
    }

    tt = "" ; if (type == 9) tt = "\t\t" ;if (type == 109) tt = "\t" ;
    if (type%100 == 0) printf("# Relative frequency of the %stions", delins) ;
    if (type%100 == 1) printf("# Number of %sted reads", delins) ;
    if (type == 2) printf("# Number of reads continuous accross the proximal site") ;
    if (type == 102) printf("# Number of reference reads") ;
    if (type == 3) printf("# Number of reads continuous accross the distal site") ;
    if (type == 9) printf("# %stion, proximal, distal counts", delins) ;

    printf("\t\t\t\t\tSorting_title") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2st[i], tt) ;
    printf("\n# \t\t\t\t\tTitle") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2t[i], tt) ;
    printf("\n# \t\t\t\t\tRun") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2r[i]) ;
    if ((type % 100) < 9)
    {
	if (forceX1) printf("\n# \t\t\tReads supporting any %stions initiated at position %d\t\t", delins, forceX1) ;
	else  printf("\n# \t\t\tReads supporting any %stions in this table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d", mm[i]) ;
	printf("\n# \t\t\tReads supporting the normal structure at proximal position\t\t") ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d", int(ww1[i])) ; 
	printf("\n# \t\t\t% Reads supporting the normal structure at proximal position\t\t") ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%.3f", 100*ww1[i]/denom[i]) ;
	printf("\n# \t\t\t% Reads supporting any %stion in the table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%.3f", 100*mm[i]/denom[i]) ;
    }
    if (type%100 == 9)
    {
	printf("\n# \t\t\tReads supporting any %stions in the table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d\t%d\t%d", mm[i], ww1[i], ww2[i]) ;
    }
    if (delins == "substitu") printf("\n# Substitution\tModified base\tModified base\tNumber of modified reads") ;
    if (delins == "dele") printf("\n# Deletion\tFirst deleted base\tLast Deleted base\tNumber of deleted reads") ;
    if (delins == "inser") printf("\n# Insertion\tBase before insertion\tBase after insertion\tNumber of inserted reads") ;
    printf ("\tMaximum\tSeen at least in in") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2t[i]) ;

    for (v in vv)
    {
        printf ("\nNC_045512.2:%s",nam[v]) ;
        printf("\t%d\t%d\t%d", x1[v], x2[v], mmv[v]) ;

	zMax = 0 ; iBest = 0 ;
	if (forceX1 > 0 && type%100 == 0) { for (i = iStart ; i <= iMax ; i++) {{k = 0 + denom[i]; if (k < 20) z = -20 ; else z = 100.0 * m[v,i]/denom[i] ; if (z > zMax){ zMax = z ; iBest = i ;}}}}
	if (forceX1 == 0  && type%100 == 0) { for (i = iStart ; i <= iMax ; i++) {k = 0 + m[v,i] + w1[v,i] ; z = 0 ; if (k > 20) z = 100.0 * m[v,i]/k ; if (z > zMax){ zMax = z ; iBest = i ;}}}
	if (type%100 == 1) { for (i = iStart ; i <= iMax ; i++) { z = m[v,i] ; if (z > zMax) { zMax = z ; iBest = i ;} }}
	if (type%100 == 2) { for (i = iStart ; i <= iMax ; i++) { z = w1[v,i] ;  if (z > zMax) { zMax = z ; iBest = i ;} }}
	if (type == 3) { for (i = iStart ; i <= iMax ; i++) { z = w2[v,i] ;  if (z > zMax) { zMax = z ; iBest = i ;} }}
	printf ("\t%.2f\t%s", zMax, i2r[iBest]) ;

	if (forceX1 > 0 && type%100 == 0) { for (i = 1 ; i <= iMax ; i++) {k = 0 + denom[i]; if (k < 20) printf( "\t%d", -20) ; else printf( "\t%.3f", 100.0 * m[v,i]/denom[i]) ; }}
	if (forceX1 == 0  && type%100 == 0) { for (i = 1 ; i <= iMax ; i++) {k = 0 + m[v,i] + w1[v,i] ; if (k < 20) printf( "\t%d", -20) ; else printf( "\t%.3f", 100.0 * m[v,i]/k) ; }}
        if (type%100 == 1) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", m[v,i]) ; }
        if (type%100 == 2) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", w1[v,i]) ; }
        if (type == 3) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", w2[v,i]) ; }

        if (type == 9) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d\t%d\t%d", m[v,i], w1[v,i], w2[v,i]) ; }
        if (type == 109) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d\t%d\t%d", m[v,i], w1[v,i], w2[v,i]) ; }
    }
    printf("\n") ;
}

