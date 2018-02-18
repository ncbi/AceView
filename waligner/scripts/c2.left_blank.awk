/^First_base_aligned/{
    x = 0 + $2 ; 
    k1 = 0 + $3 ; if (x > 0 && k1 > 0) n1[x] = k1 ; if (k1 > k1Max) {x1Max = x ; k1Max = k1 ; } nn1 += k1 ;
    k2 = 0 + $4 ; if (x > 0 && k2 > 0) n2[x] = k2 ; if (k2 > k2Max) {x2Max = x ; k2Max = k2 ; } nn2 += k2 ;
    }
END {
    z1 = 0 ; if(100 * k1Max > 40 * nn1) z1 = x1Max ; 
    z2 = 0 ; if(100 * k2Max > 40 * nn2) z2 = x2Max ; 
    printf ("Run %s\n", run) ;  # do NOT -D left_blank, because after it is used the reads are reported aligned from base 1
    if (z1>1 || z2 > 1) printf("Jump5 %d", z1) ; 
    if (z2 > 1) printf("  %d",z2) ; 
    printf("\n\n") ;
}


