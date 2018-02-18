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

/^Gene\ / {split ($1,aa," ") ; gene = aa[2];}
{ gsub(/Group_U/,"Run_U",$0);}
/^Run_U/{
  split ($2,aa," ") ;
  r=aa[1] ; z=aa[2] ; v=aa[9] ;
  ir=r2ir[r];if(ir < 1)next ;
  if (ir == 1 && z > 4 && v > 0.01) gvz[gene] = v * z ;
  gz[gene,ir] = z ;
  #  printf("g=%s ir=%s r=%s z=%s gvz=%s\n",gene,ir,r,z,gvz[gene]);
  next;
}
END {
      printf ("# Gene\tVariability") ;
      for (ir = 1 ; ir <= nr ; ir++) 
        printf ("\t%s", rr[ir]) ;
      if (nsample)
       { 
          printf ("\nSample") ;
          for (ir = 1 ; ir <= nr ; ir++)
           { printf ("\t%s", sample[ir]) ; }
       }
      if (nsystm)
       { 
          printf ("\nSystem") ;
          for (ir = 1 ; ir <= nr ; ir++)
           { printf ("\t%s", systm[ir]) ; }
       }
      if (ntissue)
       { 
          printf ("\nTissue") ;
          for (ir = 1 ; ir <= nr ; ir++)
           { printf ("\t%s", tissue[ir]) ; }
       }
      if (ntitle)
       { 
          printf ("\nTitle") ;
          for (ir = 1 ; ir <= nr ; ir++)
           { printf ("\t%s", title[ir]) ; }
       }
      printf ("\nZZZZZ") ;
      for (gene in gvz)
        {
          printf ("\n%s\t%d", gene, 1000 * (gvz[gene]+0)) ;
          for (ir = 1 ; ir <= nr ; ir++)
            printf ("\t%.2f", gz[gene,ir]);
        }
  printf("\n");
}
