/^nAli/{next}
/^nBp/{next}
/^ERR/{next}
{ if($1 != "~")ok = 0 ; if($1==chrom)ok=1 ;
 x = $2 ; if (x<xmin || x>xmax) next ; 
 if(length($3) != 1) next ;
 bp[x] = $3 ;	
 if (format == 0)
   {
     if($4 == "+") { okp[x] += $6; oop[x] += $9 ; }
     else if($4 == "-")  { okm[x] += $6; oom[x] += $9 ; }
   }
 else  # merge previous results 
   {
     okp[x] += $4 ; okm[x] += $5 ; oop[x] += $6 ; oom[x] += $7 ;
   }
}
END {
  for (x = xmin ; x <= xmax ; x++)
    if (bp[x]) printf ("%s\t%s\t%s\t%d\t%d\t%d\t%d\n", chrom, x, bp[x], okp[x], okm[x], oop[x], oom[x]) ;
}
