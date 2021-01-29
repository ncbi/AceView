{
    v = $1 ;  k=split(v,aa,":") ; ok = 0 ;
    target[v] = aa[1] ;
    if (k<4) next ;
    if (aa[3] == "Del")  { ok = 1 ; vv[v] = 1 ; x1[v] = aa[2]+1 ; dx[v] = 1 ; x2[v] = x1[v] + dx[v] - 1 ; nam[v]=x1[v]-1"Del"substr(aa[4],1);delins="dele";s1[v]=substr(aa[4],2);s2[v]="";}
    if (substr(aa[3],1,4) == "Del_")  { ok = 1 ; vv[v] = 1 ; x1[v] = aa[2]+1 ;if (forceX1 > 0) x1[v] = forceX1 ; split(aa[3],bb,"_");dx[v] = 0+bb[2] ; x2[v] = x1[v] + dx[v] - 1 ; nam[v]=x1[v]-1 aa[3] ":" substr(aa[4],2);delins="dele";s1[v]=dx[v]"bp";if(length(aa[4])>0)s1[v]=substr(aa[4],2);s2[v]="";}
    if (aa[3] == "Ins") { ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; dx[v] = 1 ; x2[v] = x1[v] + 1 ; nam[v]=x1[v]"Ins"substr(aa[5],1); delins="inser";s1[v]="";s2[v]=substr(aa[5],2);}
    if (substr(aa[3],1,4) == "Ins_") { ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; split(aa[3],bb,"_");dx[v] = 0+bb[2] ; x2[v] = x1[v] + 1 ; nam[v]=x1[v] aa[3] ":" substr(aa[5],2); delins="inser";s1[v]="";s2[v]=dx[v]"bp";if(length(aa[5])>0)s2[v]=substr(aa[5],2);}
    if (aa[3] == "Sub")  {ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; dx[v] = 1 ; x2[v] = x1[v] ; nam[v]=aa[4] x1[v] aa[5]; delins="substitu";s1[v]=aa[4];s2[v]=aa[5];printf ("@@@@ %s %d %d\n",v, x1[v],x2[v]);}
    if (substr(aa[3],1,4) == "Sub_") {ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; split(aa[3],bb,"_");dx[v] = 0+bb[2] ;  x2[v] = x1[v]+dx[v]-1 ; nam[v]=aa[4] x1[v] aa[5]; delins="substitu";s1[v]=aa[4];s2[v]=aa[5];printf ("XXXX %s %d %d\n",v, x1[v],x2[v]);}
    if (substr(aa[3],1,7)"DelIns_") {ok=2 ; vv[v] = 1 ; x1[v] = aa[2] ; split(aa[3],bb,"_");dx[v] = 0+bb[2] ; dx2[v] = 0+bb[3] ;  x2[v] = x1[v]+dx[v]-1 ; nam[v]=aa[4] x1[v] aa[5]; delins="substitu";s1[v]=aa[4];s2[v]=aa[5];}
    if (ok == 0) next ;
    printf ("XXYY %s %s %d %d\n",v,aa[3],x1[v],x2[v]);
    st = $2 ; t = $3 ; r = $4 ; i = r2i[r] ; if (i<1) { iMax++; i = iMax ; i2r[i] = r; r2i[r] = i ;} i2st[i] = st ;i2t[i] = t ;
    if (iStart < 1 && substr(t,1,1) != "A") iStart = i ;  # groups are supposed to start with A so they come first at this stage of the hack 
    m[v,i] = $5 ; w1[v,i] = $6 ; w2[v,i] = $7 ; nw[i] ++ ; ww1[i] += $6 ; ww2[i] += $7 ; mm[i] += $5 ; if (i == 1)mmv[v] = $5 ;
}
END {
    minF = 4 ;
    
    for (i = 1 ; i <= iMax ; i++) 
    { # take the average at the common donor site  
	if (forceX1>0) { if (nw[i] < 1) nw[i] = 1; ww1[i] /= nw[i] ;  }
	denom[i] = ww1[i] + mm[i] ; if (denom[i] == 0) denom[i] = 1 ;
    }

    tt = "" ; if (type == 9) tt = "\t\t" ;if (type == 109) tt = "\t" ;
    if (type%100 == 0) printf("# Frequency of the %stions", delins) ;
    if (type%100 == 1) printf("# Number of %sted reads", delins) ;
    if (type == 2) printf("# Number of reads continuous accross the proximal site") ;
    if (type == 102) printf("# Number of reference reads") ;
    if (type == 3) printf("# Number of reads continuous accross the distal site") ;
    if (type == 9) printf("# %stion, proximal, distal counts", delins) ;

    printf(" \t\t\t\t\t\t\t\tSorting_title") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2st[i], tt) ;
    printf("\n# \t\t\t\t\t\t\t\tTitle") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2t[i], tt) ;
    printf("\n# \t\t\t\t\t\t\t\tRun") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2r[i]) ;
    if ((type % 100) < 9)
    {
	if (forceX1) printf("\n# \t\t\t\t\t\tReads supporting any %stions initiated at position %d\t\t", delins, forceX1) ;
	else  printf("\n# \t\t\t\t\t\tReads supporting any %stions in this table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d", mm[i]) ;
	printf("\n# \t\t\t\t\t\tReads supporting the normal structure at proximal position\t\t") ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d", int(ww1[i])) ; 
	printf("\n# \t\t\t\t\t\t% Reads supporting the normal structure at proximal position\t\t") ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%.3f", 100*ww1[i]/denom[i]) ;
	printf("\n# \t\t\t\t\t\t% Reads supporting any %stion in the table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%.3f", 100*mm[i]/denom[i]) ;
    }
    if (type%100 == 9)
    {
	printf("\n# \t\t\t\t\t\tReads supporting any %stions in the table\t\t", delins) ;
	for (i = 1 ; i <= iMax ; i++) printf( "\t%d\t%d\t%d", mm[i], ww1[i], ww2[i]) ;
    }
    if (delins == "substitu") printf("\n# Substitution\tReference\tVariant\tModified base\tModified base\tLength\tNumber of modified reads") ;
    if (delins == "dele") printf("\n# Deletion\tReference\tVariant\tFirst deleted base\tLast Deleted base\tLength\tNumber of deleted reads") ;
    if (delins == "inser") printf("\n# Insertion\tReference\tVariant\tBase before insertion\tBase after insertion\tLength\tNumber of inserted reads") ;
    printf ("\tMaximum\tSeen at least in in") ;
    for (i = 1 ; i <= iMax ; i++) printf( "\t%s"tt, i2t[i]) ;

    for (v in vv)
    {
	zMax = 0 ; iBest = 0 ;
	if (forceX1 > 0 && (type%100)%9 == 0) { for (i = iStart ; i <= iMax ; i++) {{k = 0 + denom[i]; if (k < 20) z = -20 ; else z = 100.0 * m[v,i]/denom[i] ; if (z > zMax){ zMax = z ; iBest = i ;}}}}
	if (forceX1 == 0  && (type%100)%9 == 0) { for (i = iStart ; i <= iMax ; i++) {k = 0 + m[v,i] + w1[v,i] ; z = 0 ; if (k > 20) z = 100.0 * m[v,i]/k ; if (z > zMax){ zMax = z ; iBest = i ;}}}
	if (type%100 == 1) { for (i = iStart ; i <= iMax ; i++) { z = m[v,i] ; if (z > zMax) { zMax = z ; iBest = i ;} }}
	if (type%100 == 2) { for (i = iStart ; i <= iMax ; i++) { z = w1[v,i] ;  if (z > zMax) { zMax = z ; iBest = i ;} }}
	if (type == 3) { for (i = iStart ; i <= iMax ; i++) { z = w2[v,i] ;  if (z > zMax) { zMax = z ; iBest = i ;} }}

	if (zMax < minF) continue ;
        printf ("\n%s\t%d\t%s\t%s\t%s",target[v],x1[v],nam[v],s1[v],s2[v]) ;
        printf("\t%d\tx2=%d\t%d\t%d", x1[v], x2[v], dx[v],mmv[v]) ;

	printf ("\t%.2f\t%s", zMax, i2r[iBest]) ;

	if (forceX1 > 0 && type%100 == 0) { for (i = 1 ; i <= iMax ; i++) {k = 0 + denom[i]; if (k < 20) printf( "\t%d", -20) ; else printf( "\t%.3f", 100.0 * m[v,i]/denom[i]) ; }}
	if (forceX1 == 0  && type%100 == 0) { for (i = 1 ; i <= iMax ; i++) {k = 0 + m[v,i] + w1[v,i] ; if (k < 20) printf( "\t%d", -20) ; else printf( "\t%.3f", 100.0 * m[v,i]/k) ; }}
        if (type%100 == 1) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", m[v,i]) ; }
        if (type%100 == 2) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", w1[v,i]) ; }
        if (type == 3) { for (i = 1 ; i <= iMax ; i++) printf( "\t%d", w2[v,i]) ; }

        if (type%100 == 9)   { for (i = 1 ; i <= iMax ; i++) printf( "\t%d\t%d\t%d", m[v,i], w1[v,i], w2[v,i]) ; }
    }
    printf("\n") ;
}

