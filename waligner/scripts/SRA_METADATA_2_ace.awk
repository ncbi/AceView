# params: prefix=Gce nn=21 will create as first run Gce22
/^Run/{next;}

{   nn++ ; 
    printf("\nRun %s%d\nSRR %s\nRunId %s\n",prefix,nn, $1,$1) ; 
    n=split($2,aa,"/") ; if (n==3 && aa[3]+0>0 && aa[3]<1900) aa[3]+=2000 ; if (n >1) printf("Submission_date %s-%s-%s\n",aa[3],aa[1],aa[2]) ;   
    else { n=split($2,aa,"-") ; if (n==3 && aa[3]+0>0 && aa[3]<1900) aa[3]+=2000 ; printf("Submission_date %s-%s-%s\n",aa[3],aa[2],aa[1]) ;    }
    printf("Spots %d  Mb_in_SRA %d Average_length %d\n", $4, int($5/1000000), $7) ; 
    if ($6>1)
	printf("Paired_end\nFile fasta/1 DATA/%s_1.fasta.gz\nFile fasta/2 DATA/%s_2.fasta.gz\n",$1,$1) ; 
    else 
	printf("File fasta DATA/%s.fasta.gz\n",$1) ; 
    printf("!Details %s\n",$10) ;  
    nsrx[$11]++ ;  srx2srr[$11] = srx2srr[$11] "," prefix nn ; 
    if($12) 
    {
	gsub(/ /,"_",$12) ; printf("Sample \"%s\"\n",$12) ;
    }
    if (substr($13,1,5) == "miRNA" || $14 == "size fractionation")
	printf("Small_RNA\n") ; 
    else if ($13 == "RIP-Seq")
	printf("RIP\n") ; 
    else if (index($14,"cDNA")>0 || index($14,"PolyA")>0  || index($13,"cDNA")>0  || index($13,"EST")>0 )
	printf("PolyA\n") ; 
    else 
	printf("RNA\n") ;  
    if ($13 || $14 || $17 || $18)
    {
	printf ("Library \"") ;
	if ($13 || $14)
	    printf ("SRA_Library_strategy_and_selection %s %s ", $13, $14) ; 
	if ($17 > 0)
	    printf (" SRA_Insert_size %s ", $17) ;
	if ($18 > 0) 
	    printf (" SRA_Insert_sigma %s ", $18) ;
	printf ("\"\n") ;
    }

    gsub("Illumina ","",$20);
    if ($19 == "ILLUMINA") printf ("Illumina \"%s\"\n", $20) ;
    if ($19 == "LS454") printf ("Roche_454 \"%s\"\n", $20) 
    if ($19 == "ION_TORRENT") printf (" Ion_Torrent \"%s\"\n", $20) ;
    if ($19 == "ABI_SOLID") printf (" SOLiD \"%s\"\n", $20) ;

    if ($21 || $42)
    {
	printf("Author \"") ;
	if ($21) printf ("%s", $21) ;
	if ($42) printf (" %s", $42) ;
	printf("\"\n") ;
    }
    if ($23)
	printf ("Reference pm%s\n", $23) ;
    if ($30)
    {
	gsub ("Caenorhabditis elegans ", "", $30) ;
	printf ("Title \"%s\"\n", $30) ;
    }
}

END { 
    printf("\n") ; 
    for(x in nsrx)
    {
	if (nsrx[x] > 1)
	{
	    n = split(srx2srr[x],aa,",") ;
	    printf("Run %s\nAdd_counts\n",x) ; 
	    for(i = 2 ; i <= n ; i++)
		printf("Runs %s\n",aa[i]) ; 
	    printf("\n") ; 
	}
    }
}
