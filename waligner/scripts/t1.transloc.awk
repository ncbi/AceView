BEGIN {printf("#Run\tRead\tGene_A\tGene_B\tmRNA_S\tfrom\tto\tmRNA_B\tfrom\tto\tscore A\tscore B\tAli A\tAli B\tx1 A\tx2 A\tx1 B\tx2 B\tc1 A\tc2 A\tc1 B\tc2 B\n");}
/^#/{next ;}
{
    seq = $1 ;    score = $2 ;    ali = $4 ; tc = $8 ; # target_class
    x1 = $6 ; x2 = $7 ;    # exon
    c1 = $27 ; c2 = $28 ;  # chain in x orientation
    gene = $9 ; mrna = $11; a1 = $12 ; a2 = $13 ;
    gsub (/^X__/,"",gene) ;
    cc1 = c1 ; if (oldC1 > cc1) cc1 = oldC1 ;   # overlap start
    cc2 = c2 ; if (oldC2 < cc2) cc2 = oldC2 ;   # overlap end
    dc = c2 - c1 ;           # chain length
    dcc = 3 * (cc2 - cc1) ;  # 3 * length of overlap, must be smaller than old and new chain length

    sens = (oldA2 - oldA1) * (a2 - a1) ;
    
    if (sens > 0)
    { 
	wo =  "++" ;
	if ((a2 - a1) * (x2 - x1) * (x1 + x2 - oldX1 - oldX2) > 0 )
	    wab = 0 ; # oldGene ---> ---> gene 
	else
	    wab = 1 ; # gene  ---> ---> oldGene
    }
    else
    { 
	if ((oldX2 > oldX1) * (x1 + x2 - oldX1 - oldX2) > 0)
	{ 
	    wo =  "+-" ;
	    if (oldGene < gene)
		wab = 0 ;  # oldGene ---> <--- gene
	    else
		wab = 1 ;  # gene  ---> <--> oldGene
	}
	else
	{
	    wo =  "-+" ;
	    if (oldGene < gene)
		wab = 0 ;  # oldGene <--- ---> gene
	    else
		wab = 1 ;  # gene  <--- ---> oldGene
	}
    }

    #print "......", wab, wo, oldGene, gene, sens, oldA1,oldA2,a1,a2
    if (seq == oldSeq && gene != oldGene && ali > minAli && oldAli > minAli && dcc < dc && dcc < oldDc) 
    {
	if (wab == 0)
	{
	    printf("%s\t%s\t%s__%s%s\t%s\t%d\t%d\t%s\t%d\t%d\t", run, seq, oldGene, gene, wo, oldMrna, oldA1, oldA2, mrna, a1, a2) ;
	    printf("%d\t%d\t%d\t%d\t", oldScore, score, oldAli, ali) ;
	    printf("%d\t%d\t%d\t%d\t", oldX1, oldX2, x1, x2) ;
	    printf("%d\t%d\t%d\t%d\n", oldC1, oldC2, c1, c2) ;
	}
	else
	{
	    printf("%s\t%s\t%s__%s%s\t%s\t%d\t%d\t%s\t%d\t%d\t", run, seq, gene, oldGene, wo, mrna, a1, a2, oldMrna, oldA1, oldA2) ;
	    printf("%d\t%d\t%d\t%d\t", score, oldScore, ali, oldAli) ;
	    printf("%d\t%d\t%d\t%d\t", x1, x2, oldX1, oldX2) ;
	    printf("%d\t%d\t%d\t%d\n",c1, c2,  oldC1, oldC2) ;
	}
    }
    oldSeq = seq ;    oldScore = score ; oldAli = ali ;
    oldX1= x1 ; oldX2 = x2 ; oldC1 = c1 ; oldC2 = c2 ; oldDc = dc ;
    oldGene = gene ; oldMrna = mrna ; oldA1 = a1 ; oldA2 = a2 ;

    link = $22 ;
    i = index(link, "links_to:") ;
    if (i > 0)
    {
	c = substr (seq, length(seq)) ;
	frag = substr (seq, 1, length(seq) - 1) ;
	if (c == "<")   # comes first in alphabetic order
	{
	    zOldFrag = frag ;
	    zOldGene = gene ;
	    zOldA1 = a1 ;  
	    zOldA2 = a2 ;  
	    next ;
	}
	lGene = substr (link, i+9) ;
	gsub (/^X__/,"",lGene) ;
	
	if (c != ">" || frag != zOldFrag || lGene != zOldGene)   # not a pair
	    next ;  
	if (a1 > a2)
	{ x0 = x1 ; x1 = x2 ; x2 = x0 ; }
	sens = - (zOldA2 - zOldA1) * (a2 - a1) ; # flip since this is a pair
	# print zOldX1, zOldX2, x1, x2, sens, a1, a2, strand ;
	if (sens > 0)
	{
	    wo =  "++" ;
	    if (strand * (a2 - a1) > 0) # assume the > read is strand
		wab = 0 ; # zOldGene ---> ---> gene 
	    else
		wab = 1 ; # gene  ---> ---> zOldGene
	}
	else
	{
	    if (strand * (a2 - a1) > 0) # assume the > read is strand
	    { 
		if (strand == -1)
		{ wo =  "-+" ; wab = 0 ; } # zOldGene <--- ---> gene
		else
		{ wo =  "+-" ; wab = 1 ; }  # gene  <--- ---> zOldGene
	    }
	    else
	    {
		if (strand == -1)
		{ wo =  "-+" ; wab = 0 ; } # zOldGene <--- ---> gene
		else
		{ wo =  "+-" ; wab = 1 ; }  # gene  <--- ---> zOldGene
	    }
	}
	
	if (wab == 0)
	{
	    z = frag "#" zOldGene "#" gene wo ;
	    if (z != oldZ)
		printf ("%s\t%s\t%s__%s%s\tPAIR\n", run, seq, zOldGene, gene, wo) ;
	}
	else
	{
	    z = frag "#" gene "#" zOldGene wo ;
	    if (z != oldZ)
		printf ("%s\t%s\t%s__%s%s\tPAIR\n", run, seq, gene, zOldGene, wo) ;
	}
	oldZ = z ;
    }
}
