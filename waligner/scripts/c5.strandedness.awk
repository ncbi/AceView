{
 gsub(/\\\//,"/",$0);
}

/^ZZZZZ/{zz++; next; }
{ if (zz < 1) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r]=1; isGroup[r]=0; next; }}
{ if (zz < 2) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r]=1; isGroup[r]=1; ngroups++; next; }}
{ if (zz < 3) { z=$1 ; r = r2i[z] ; if(r>=1) private[r] = 1 ; next ;}}
    
{
  gsub("\"","",$0);
}

{ if (zz == 3)
    {
	nr0 = nr ;
	gsub (/System:/,"",$0);
	gsub (/NULL/,"NA",$0);
	gsub (/NA,/,"",$0);
	gsub (/ /,"_",$0);
	
	split($0,aa,"\t");
	z = aa[1] ; ir = r2i[z]; 
	# print z, ir ;
	if (ir<1) next ;
	
	for (i = 1 ; i <= NF ; i++) $i = aa[i] ; 
	# print z,ir,$11 ;
	machine[ir]=$2; if ($2 != "NA") nMachine = 1 ;
	if ($2 != "NA") { machine[ir] = $2;  nmachine = 1 ; }
	if ($9 != "NA") { machine[ir] = machine[ir] ", " $9 ;  nMachine = 1 ; }
	
	if ($3 != "NA") { if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ;  nSample = 1 ; }
	if ($4 != "NA") { if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; nSystm = 1 ;}
	if ($5 != "NA") { if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ; nTissue = 1 ; }
	if ($7 != "NA") { if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ; nSystm = 1 ; }
	if ($8 != "NA") {runid[ir] = $8 ;  nrunid = 1 ; }
	
	if ($6 != "NA") { titre[ir]=$6; nTitre = 1 ;}
	if ($11 != "NA") { stitre[ir]=$11; nsTitre = 1 ;}
	if ($12 != "NA") { stitre2[ir]=$12; nsTitre2 = 1 ;}
	if ($13 != "NA") { otitre[ir]=$13; noTitre = 1 ;}
     next ;
    }
}
{
  split ($1,aa, "~") ;
  r = aa[1] ; ir=r2i[r];if(ir<1){nr++;r2i[r]=nr;i2r[nr]=r;ir=nr;}
  if (zz == 4) 
    {
	if (aa[2] != "Exons_juntions")
	{
	    f = aa[2] ; 
	    if (! f) f = "f" ;
	    if (f == "f") f = "single" ;
	    if (f == "f2") f = "pair" ;
	    ff[f] = 1 ;
	}
      t = $2 ; ss[ir,f,t] = $3 ;  nn[ir,f] += $4 + $5 + $6 ; 
      # print "###",f,ir,i2r[ir],nn[ir,f] ;
      
      rrp[ir,f,t] = $4 ; rrm[ir,f,t] = $5 ; rra[ir,f,t] = $6 ; 
      nnt[t] += $4 + $5 + $6 ;
      if (0 && f != "ns")
	{
	  f = "ns" ; ff[f] = 1 ; nn[ir,f] += $4 + $5 + $6 ;  
	  t = $2 ; ss[ir,f,t] += $3 ; 
	  rrp[ir,f,t] += $4 ; rrm[ir,f,t] += $5 ; rra[ir,f,t] += $6 ; 
	} 
      # printf( "### r=%s ir=%s f=%s t=%s :: %d %d\n", r, ir, f, t,rrp[ir,f,t],rrm[ir,f,t]) ;
      # gsub(/Forward/,"Stranded",$7) ; 
      proto[ir,f] = $7 ; 
    }
}
END { 
     nr = nr0 ;
    want_all = 1 ;

    for (pass = 1 ; pass <= 1 ; pass++)
    {
	if (pass == 1) printf ("\n# Run") ;  
	else printf ("\t\t\t# Run") ;  
	
	if (nSample + want_all >= 0) printf("\tSample"); 
	if (nrunid + want_all >= 0) printf("\tRunId"); 

	if (nTitre + want_all >= 0) printf("\tTitle"); 
	if (nsTitre + want_all >= 0) printf("\tSorting_Title");  
	if (nsTitre2 + want_all >= 0) printf("\tSorting_Title_2");  
	if (noTitre + want_all >= 0) printf("\tOther_Title");  

	if (nTissue + want_all >= 0) printf("\tTissue"); 
	if (nSystm + want_all >= 0) printf("\tSystem"); 
	# if (nMachine + want_all >= 0) printf("\tMachine");
	
	t = "B_rrna"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    
	    printf ("\trRNA on plus strand") ;
	    printf ("\trRNA on minus strand") ;
	    printf ("\trRNA ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of rRNA") ;
	    printf ("\t%% on minus strand on rRNA") ;
	    printf ("\t%% ambiguous strand on rRNA") ;
	}
	t = "0_SpikeIn"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    
	    printf ("\tSpike-In on plus strand") ;
	    printf ("\tSpike-In on minus strand") ;
	    printf ("\tSpike-In ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of Spike-In") ;
	    printf ("\t%% on minus strand of Spike-In") ;
	    printf ("\t%% ambiguous strand on rRNA") ;
	}
	t = "KT_RefSeq"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    if (ntitle > 0) printf("\tTitle") ;
	    if (nstitle == 1)  printf("\tSorting title") ;
	    if (nstitle2 == 1)  printf("\tSorting title 2") ;
	    if (notitle == 1)  printf("\tOther title") ;
	    
	    printf ("\tRefSeq on plus strand") ;
	    printf ("\tRefSeq on minus strand") ;
	    printf ("\tRefSeq ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of RefSeq") ;
	    printf ("\t%% on minus strand of RefSeq") ;
	    printf ("\t%% ambiguous strand on RefSeq") ;
	}
	t = "ET_av" ;
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    printf ("\tAceView on plus strand") ;
	    printf ("\tAceView on minus strand") ;
	    printf ("\tAceView ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of AceView") ;
	    printf ("\t%% on minus strand of AceView") ;
	    printf ("\t%% ambiguous strand on AceView") ;
	}
	t = "LT_UCSC" ;
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    printf ("\tucscGenes on plus strand") ;
	    printf ("\tucscGenes on minus strand") ;
	    printf ("\tucscGenes ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of ucscGenes") ;
	    printf ("\t%% on minus strand of ucscGenes") ;
	    printf ("\t%% ambiguous strand on ucscGenes") ;
	}
	t = "NT_MiT"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    if (ntitle > 0) printf("\tTitle") ;
	    if (nstitle == 1)  printf("\tSorting title") ;
	    if (nstitle2 == 1)  printf("\tSorting title 2") ;
	    if (notitle == 1)  printf("\tOther title") ;
	    
	    printf ("\tMiT on plus strand") ;
	    printf ("\tMiT on minus strand") ;
	    printf ("\tMiT ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of MiT") ;
	    printf ("\t%% on minus strand of MiT") ;
	    printf ("\t%% ambiguous strand on MiT") ;
	}
	t = "MT_EBI"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    if (ntitle > 0) printf("\tTitle") ;
	    if (nstitle == 1)  printf("\tSorting title") ;
	    if (nstitle2 == 1)  printf("\tSorting title 2") ;
	    if (notitle == 1)  printf("\tOther title") ;
	    
	    printf ("\tEBI on plus strand") ;
	    printf ("\tEBI on minus strand") ;
	    printf ("\tEBI ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of EBI") ;
	    printf ("\t%% on minus strand of EBI") ;
	    printf ("\t%% ambiguous strand on EBI") ;
	}
	t = "Exons_juntions" ;
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    printf ("\tJunctions on plus strand") ;
	    printf ("\tJunctions on minus strand") ;
	    printf ("\tJunctions ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of Junctions") ;
	    printf ("\t%% on minus strand of Junctions") ;
	    printf ("\t%% ambiguous strand on Junctions") ;
	}
	t = "Z_genome"  ; 
	if (nnt[t] + 0 > 0)
	{
	    printf ("\t\tRun") ;
	    printf ("\tGenome on plus strand") ;
	    printf ("\tGenome on minus strand") ;
	    printf ("\tGenome  ambiguous strand") ;
	    printf ("\t%% aligning on plus strand of Genome") ;
	    printf ("\t%% aligning on minus strand of Genome") ;
	    printf ("\t%% ambiguous strand on Genome") ;
	}
    }
  printf ("\n") ;
  for (r = 1 ; r <= nr ; r++)
  {
      if (private[r] == 1) continue ;
      if (nr >= 24)
      {
	  split (title[r], gnam, "_") ;
	  if (ng == 0) 
	    printf("\n") ;
	  else
	    {
	      if ((gnam[1] != oldgnam1 || (ngroups > 24 && gnam[2] != oldgnam2)) && ! (gnam[3] == "any" && oldgnam3 == "any"))
		printf("\n") ;
	    }
	  oldgnam1 = gnam[1] ; oldgnam2 = gnam[2] ; oldgnam3 = gnam[3] ;
	  ng++ ;
	}
      printf("%s", i2r[r]) ;
      if (nSample + want_all >= 0) printf("\t%s", sample[r]); 
      if (nrunid + want_all >= 0) printf("\t%s", runid[r]); 
      
      if (nTitre + want_all >= 0) printf("\t%s", titre[r]);
      if (nsTitre + want_all >= 0) printf("\t%s", stitre[r]);
      if (nsTitre2 + want_all >= 0) printf("\t%s", stitre2[r]);
      if (noTitre + want_all >= 0) printf("\t%s", otitre[r]);
      
      if (nTissue + want_all >= 0) printf("\t%s", tissue[r]);
      if (nSystm + want_all >= 0) printf("\t%s", systm[r]); 
      # if (nMachine + want_all >= 0) printf("\t%s",  machine[r]( ;

      for (f in ff)
	{ 
	  #  print f,r,i2r[r],nn[r,f] ;
	  if (nn[r,f] + 0 < 1) continue ;

	  kkn = split ("B_rrna:0_SpikeIn:KT_RefSeq:ET_av:LT_UCSC:NT_MiT:MT_EBI:Exons_juntions:Z_genome",kkk,":") ;
	  for(kk = 1 ; kk <= kkn ; kk++)
	  {
	      t = kkk[kk] ;
	      if (nnt[t] + 0 < 1) continue ;

	      printf("\t\t%s", i2r[r]) ;
	      printf("\t%d\t%d\t%d", rrp[r,f,t], rrm[r,f,t], rra[r,f,t]) ;
	      
	      w1 = 0 + rrp[r,f,t] + rrm[r,f,t] ; if (w1 == 0) w1 = 1 ;
	      w = 0 + 100.0 * rrp[r,f,t] / w1 ;
	      if (rrp[r,f,t] + rrm[r,f,t] > 1000)
		  printf("\t%.3f", w) ;
	      else
		  printf("\t", w) ;
	      
	      w = 0 + 100.0 * rrm[r,f,t] / w1 ;
	      if (rrp[r,f,t] + rrm[r,f,t] > 1000)
		  printf("\t%.3f", w) ;
	      else
		  printf("\t", w) ;
	      
	      w2  = 0 + w1 + rra[r,f,t] ;
	      w = 0 +100.0*rra[r,f,t] / w2
	      if (rrp[r,f,t] + rrm[r,f,t] > 1000)
		  printf("\t%.3f", w) ;
	      else
		  printf("\t", w) ;
	  }
	}	  
      printf ("\n") ;
    }
  printf ("\n") ;
}

