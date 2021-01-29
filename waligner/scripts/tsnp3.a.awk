 /^ZZZZZ/{zz++;next;}
{ if(zz==1){r2t[$2]=$1;next;}}
{ if(zz<1){r2st[$2]=$1;next;}}
{
   run = $3 ;
   if (any == 1) run = "any" ;
   if (any == 2) {split(r2st[run],aa,"_");run=aa[1];if(run=="Q")next;}
   v = $1"\t" run ; vv[v]=1 ;
 }
 /Variant/  { vm[v]  += $4 ; next ; }
 /RefDonor/ { vw1[v] += $4 ; next ; }
 /RefAccep/ { vw2[v] += $4 ; next ; }
 /VarDonor/ { vm1[v] += $4 ; next ; }
 /VarAccep/ { vm2[v] += $4 ; next ; }
 /Ref\t/    { vw[v]  += $4 ; next ; }
END {
  for (v in vv)
     printf ("%s\t6\t%d\t%d\t%d\t%d\t%d\t%d\n",v,vm[v],vw[v],vm1[v],vm2[v],vw1[v],vw2[v]);
}
