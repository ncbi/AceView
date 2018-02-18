/^#/{next;}
{ chrom=$1;a1=$2;a2=$3 ;
 type=$4;

 if (type == "Gene")
   {
     ng++ ; 
     gene[ng] = $5 ; geneStrand[ng] = $6 ; geneTitle[ng] = $7 ;
     g1[ng] = a1 ; g2[ng] = a2 ;
     if (a2 > g2Max) { g2Max2 = g2Max ; g2Max = a2 ; }
     next ;
   }
 if (type == "mRNA")
   {
     if (a2 < m2Max) next ;
     nm++ ; 
     mrna[nm] = $5 ; 
     ma1[nm] = a1 ; ma2[nm] = a2 ;      mx1[nm] = $6 ; mx2[nm] = $7 ;
     if (a2 > m2Max) { m2Max = a2 ; }
     next ;
   }
 if (type == "Variant")
   {
     nok = 0 ;
     if (a2 < g2Max)
       for (g = ng ; g > 0   ; g--)
	 {
	   if (g2[g] >= a1 && g1[g] <= a1)
	     {
	       nok++ ;
	       printf ("Variant\t%s\tGene\t%s\t%s\t%s\n", $5, gene[g], geneStrand[g], geneTitle[g]) ;
	     }
	   if (nok> 0 && a2 > g2Max2) break ;
	 }
     nok = 0 ;
     if (a2 < m2Max)
       for (m = nm ; m > 0   ; m--)
	 {
	   if (ma2[m] >= a1 && ma1[m] <= a1)
	     {
	       nok++ ;
	       if (mx1[m] < mx2[m])
		 { u1 = mx1[m] + a1 - ma1[m] ; u2 =  1 ; }
	       else
		 { u1 = mx1[m] - a1 + ma1[m] ; u2 = - 1 ; }
	       printf ("Variant\t%s\tmRNA\t%s\t%d\t%d\n", $5, mrna[m], u1, u2) ;
	     }
	   if (nok> 0) break ;
	 }
   }
}

 
