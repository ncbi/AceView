/^ZZZZZ/{zz = 1; next; }
{ if (zz < 1) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r]=1; next; }}

/^\"/{ gsub("\"","",$0);split($0,aa,"\t");rg=aa[1];machine[rg]=aa[2];sample[rg]=aa[3];systm[rg]=aa[4];tissue[rg]=aa[5];titre[rg]=aa[6];stitre[rg]=aa[11];otitre[rg]=aa[12];systm[rg] = systm[rg] " " aa[7];runid[rg]=$8;next;}
/^RUN/ {z=$2; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} next;}
{  z=$1;  t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;} nnrt[r,t] += $2 ; nnr[r]+= $2 ; nnt[t] += $2; next; }
END {
    printf ("Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (titre[i2r[r]] != "NULL") z = titre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSorting_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (stitre[i2r[r]] != "NULL") z = stitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nOther_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (otitre[i2r[r]] != "NULL") z = otitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSample") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (sample[i2r[r]] != "NULL") z = sample[i2r[r]] ; printf ("\t%s", z) ; } printf ("\t-") ;
    printf ("\nPlatform") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (platform[i2r[r]] != "NULL") z = platform[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf("\nRun") ; for (r = 1 ; r <= nr ; r++) printf ("\t%s",i2r[r]) ; printf ("\tAny run") ;
    printf("\nRunId") ; for (r = 1 ; r <= nr ; r++) printf ("\t%s",runid[i2r[r]]) ; printf ("\tAny run") ;

    for (t = 1 ; t <= nt ; t++)
      {
	printf ("\n%s", i2t[t]) ;
	for (r = 1 ; r <= nr ; r++)
	  printf ("\t%d",nnrt[r,t]) ;
	printf ("\t%d",nnt[t]) ;
      }
   printf ("\n") ;

    printf ("Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (titre[i2r[r]] != "NULL") z = titre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSorting_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (stitre[i2r[r]] != "NULL") z = stitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nOther_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (otitre[i2r[r]] != "NULL") z = otitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSample") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (sample[i2r[r]] != "NULL") z = sample[i2r[r]] ; printf ("\t%s", z) ; } printf ("\t-") ;
    printf ("\nPlatform") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (platform[i2r[r]] != "NULL") z = platform[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf("\nRun") ; for (r = 1 ; r <= nr ; r++) printf ("\t%s",i2r[r]) ; printf ("\tAny run") ;

    for (t = 1 ; t <= nt ; t++)
      {
	printf ("\n%s", i2t[t]) ;
	for (r = 1 ; r <= nr ; r++)
	  {
	    u = nnt[t]; if (u == 0) u = 1 ; 
	    printf ("\t%.2f",100.0 * nnrt[r,t]/u) ;
	  }
	printf ("\t100") ;
      }
   printf ("\n") ;

    printf ("\nTitle") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (titre[i2r[r]] != "NULL") z = titre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSorting_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (stitre[i2r[r]] != "NULL") z = stitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nOther_Title") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (otitre[i2r[r]] != "NULL") z = otitre[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf ("\nSample") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (sample[i2r[r]] != "NULL") z = sample[i2r[r]] ; printf ("\t%s", z) ; } printf ("\t-") ;
    printf ("\nPlatform") ;  for (r = 1 ; r <= nr ; r++) { z = "-"; if (platform[i2r[r]] != "NULL") z = platform[i2r[r]] ; printf ("\t%s",z) ;}  printf ("\t-") ;
    printf("\nRun") ; for (r = 1 ; r <= nr ; r++) printf ("\t%s",i2r[r]) ; printf ("\tAny run") ;

    for (t = 1 ; t <= nt ; t++)
      {
	printf ("\n%s", i2t[t]) ; u2 = 0 ; v2 = 0 ;
	for (r = 1 ; r <= nr ; r++)
	  {
	    u = nnr[r]; u2 += u ; if (u == 0) u = 1 ; 
	    printf ("\t%.2f",100.0 * nnrt[r,t]/u) ; v2 += nnrt[r,t] ;
	  }
	if (u2 == 0) u2 = 1 ; 
	printf ("\t%.2f",100.0 * v2/u2) ;
      }
   printf ("\n") ;

}

           
  
