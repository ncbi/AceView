/^#/{ next ; }
{if (chrom && $1 != chrom) next ;}
{
    # search for monomodals in columns
    chi2 = $23 ;
   # ref=26 low 27 mid 28 high = 29 pure = 30
    n = $26 + $27 + $28 + $29 + $30 ;
    bad = $27 ;
   # print n, bad ; next ;

    if (10 * $27 > n) next ;

    print ;
}
