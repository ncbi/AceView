/^#/{ line++ ; if(aceOut < 1) print ; if (iMin > 1) next ; for (i = 1 ; i <= NF ; i++) { if ($i == "Run") iMin = i + 2 ; if ($i == "SNV type") iType = i; run[i] = $i; if(i > iMin) print iMin,i,run[i];}iMax = NF ; }
	{ last;
    if (protein)
    {
	ok = 0 ;
	split ($iType, bb, " ") ; w = bb[1] ; if(substr(w,1,3)=="Ter") w="Ter";
	k = index ("AA_substitution AA_to_Stop Frame_preserving_indel Extension Frameshift Met1del Met1_lost Met1_gained Stop_to_AA Ter", w) ;
	if (k < 1)
	    next ;
    }
    k = split ($10, aa, ":") ; 
    gene = aa[1] ; if (gene == "") next; 
    genes[gene] = 1; 
    af = $(iMin-3) ; # allele frequency in cohort
    for (i = iMin ; i <= iMax ; i++) 
    {
	z = $i ; 
	if (z > -10)
	{
	    if ( zgene[gene, i] < 1000) zgene[gene,i] = 1000 ;
	    if (af > 50) z = 100 - z ; 
	    zgene[gene, i] += z ;
	    genes[gene] += z ;
	}
    }
    next;
}
END {

    for (gene in genes)
    {
	if (genes[gene] < 1*(iMax - iMin)) continue ;
	if (aceOut < 1)
	{
	    printf("%s", gene) ;
	    for (i = 2 ; i < iMin - 1 ; i++)
		printf ("\t") ;
	    for (i = iMin ; i <= iMax ; i++)
	    {
		z =  zgene[gene, i] + 0 ;
		if (z < 1000) z = -10 ;
		else z -= 1000 ;
		printf ("\t%d", z) ;
	    }
	    printf ("\n") ;
	}
	else
	{
	    printf("Gene \"%s\"\n", gene) ;
	    for (i = iMin ; i <= iMax ; i++)
	    {
		z =  zgene[gene, i] + 0 ;
		if (z < 1000) continue ;
		printf ("Run_U %s 0 %d seqs %d tags %d kb\n", run[i], z, z, z) ;
	    }
	    printf ("\n") ;
	}
    }

}
