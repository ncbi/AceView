BEGIN {} 
/^ZZZZZ/{zz++ ; next ;}
{ if (zz < 1) { r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} next ; }}
{gsub (/\"/,"",$0);}
{  if (zz == 1)
    {
      r = $1 ; ir=r2ir[r]; if(ir < 1) next ; 
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") nmachine = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ; if ($2 != "NA") nsample = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ; if ($7 != "NA") nsystm = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ; if ($7 != "NA") ntissue = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; title[ir]=$6; if ($6 != "NA") ntitle = 1 ;
      next ;
    }
}
/^# Chromosome/ { 
  for (i = 1 ; i <= NF ; i++)
    { ir = r2ir[$i] ; i2ir[ir] = ir ; nam[i] = $i ;}
  imax=i ;
  next;
}
/^#/{ next ; }
{ 
  for (i = 9 ; i <= NF ; i++)
    {
      if ($i < 1) continue ; 
      if ($8 + 0 < 1) { nnNoRefSeq[i]++ ; }
      if ($9 + 0 < 1) { nnNoEst[i]++ ; }
      z = 1 ; j = 0 ;
      for ( k = 1 ; k <= 8 ; k++)
        {
          if (k>1) z = 10 * z ;
          for (a = 1 ; a <= 5 ; a++)
            if (a == 1 || a == 2 || a == 5)
              {
                j++ ; 
                if (0 + $i > a * z) 
                  { 
                    nn[j]++ ; n[i,j]++;
                    if (0)  printf("j=%d a=%d z = %d i=%d x=%d\n",j,a,z,i,$i); 
                    if ($8 + 0 < 1) { nNoRefSeq[i,j]++; nnNoRefSeq[i]++ ; }
                    if ($9 + 0 < 1) { nNoEst[i,j]++; nnNoEst[i]++ ; }
                  }
              }
         }
     }
}
END {
      printf ("# Support") ;
      for (i = 9 ; i <= imax ; i++)
        printf ("\t%s", nam[i]) ;
      if (nsample)
       { 
          printf ("\nSample") ;
          for (i = 9 ; i <= imax ; i++)
           { ir = i2ir[i] ; printf ("\t%s", sample[ir]) ; }
       }
      if (nsystm)
       { 
          printf ("\nSystem") ;
          for (i = 9 ; i <= imax ; i++)
           { ir = i2ir[i] ; printf ("\t%s", systm[ir]) ; }
       }
      if (ntissue)
       { 
          printf ("\nTissue") ;
          for (i = 9 ; i <= imax ; i++)
           { ir = i2ir[i] ; printf ("\t%s", tissue[ir]) ; }
       }
      if (ntitle)
       { 
          printf ("\nTitle") ;
          for (i = 9 ; i <= imax ; i++)
           { ir = i2ir[i] ; printf ("\t%s", title[ir]) ; }
       }
      if (1)
       { 
          printf ("\nNo RefSeq") ;
          for (i = 9 ; i <= imax ; i++)
           { printf ("\t%d", nnNoRefSeq[i]) ; }
       }
      if (1)
       { 
          printf ("\nNo AceView") ;
          for (i = 9 ; i <= imax ; i++)
           { printf ("\t%d", nnNoEst[i]) ; }
       }
      z = 1 ; j = 0 ;
      for ( k = 1 ; k <= 8 ; k++)
        {
          if (k>1) z = 10 * z ;
          for (a = 1 ; a <= 5 ; a++)
            if (a==1 || a==2 || a==5)
              {
                j++ ; 
                if (nn[j] > 0)
                  {
                    printf ("\n%d", a*z) ;
                    for (i = 9 ; i <= imax ; i++)
                      printf ("\t%d", n[i,j]);
                   }
               }
        }
  printf("\n");
}
