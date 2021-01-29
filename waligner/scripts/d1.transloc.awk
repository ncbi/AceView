BEGIN {printf("#Run\tRead\tChrom A\tfrom\tto\tChrom B\tfrom\tto\tscore A\tscore B\tAli A\tAli B\tx1 A\tx2 A\tx1 B\tx2 B\tc1 A\tc2 A\tc1 B\tc2 B\n");}
/Z_genome/ {
    gsub (/>/,"",$1) ;
    gsub (/</,"",$1) ;
    seq = $1 ;    score = $2 ;    ali = $4 ;
    x1 = $6 ; x2 = $7 ;    c1 = $27 ; c2 = $28 ;
    chr = $11 ; a1 = $12 ; a2 = $13 ;
    cc1 = c1 ; if (oldC1 > cc1) cc1 = oldC1 ;
    cc2 = c2 ; if (oldC2 < cc2) cc2 = oldC2 ;
    dc = c2 - c1 ; dcc = 3 * (cc2 - cc1) ;
    
    if (seq == oldSeq && (chr != oldChr || a1 > oldA1 + 1000000 || a1 < oldA1 - 1000000) && ali > minAli && oldAli > minAli && dcc < dc && dcc < oldDc) 
    {
	printf("%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t", run, seq, oldChr, oldA1, oldA2, chr, a1, a2) ;
	printf("%d\t%d\t%d\t%d\t", oldScore, score, oldAli, ali) ;
	printf("%d\t%d\t%d\t%d\t", oldX1, oldX2, x1, x2) ;
	printf("%d\t%d\t%d\t%d\n", oldC1, oldC2, c1, c2) ;

	printf("%s\t%s\t%s\t%d]\t%d\t%s\t%d\t%d\t", run, seq, chr, a1, a2, oldChr, oldA1, oldA2) ;
	printf("%d\t%d\t%d\t%d\t", score, oldScore, ali, oldAli) ;
	printf("%d\t%d\t%d\t%d\t", x1, x2, oldX1, oldX2) ;
	printf("%d\t%d\t%d\t%d\n",c1, c2,  oldC1, oldC2) ;
    }
    oldSeq = seq ;    oldScore = score ; oldAli = ali ;
    oldX1= x1 ; oldX2 = x2 ; oldC1 = c1 ; oldC2 = c2 ; oldDc = dc ;
    oldChr = chr ; oldA1 = a1 ; oldA2 = a2 ;
}

