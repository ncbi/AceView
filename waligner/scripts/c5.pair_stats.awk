BEGIN{
  z = "Aligned_fragments" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Non_compatible_pairs" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Orphans_Any" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs_inside_gene" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_but_links_2_transcripts_of_a_gene" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_gene_extension" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs_in_genome" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs_in_mito" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs_in_other_target" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Compatible_pairs_in_spikeIn" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Incompatible_topology" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}

  z = "Links_2_genes" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Links_gene_to_distant_genome" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Too_distant_on_genome" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Links_gene_to_rRNA" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Links_gene_to_mito" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  z = "Orphans_1" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  if (0)
  {
      z = "Orphans_2" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_3" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_4" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_5" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_6" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_7" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_8" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
      z = "Orphans_9" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
  } 
  z = "Orphans_multi" ; t = t2i[z] ; if(t<1){nt++;t=nt;i2t[nt]=z;t2i[z]=nt;}
}

{
 gsub(/\\\//,"/",$0);
}

/^ZZZZZ/{zz++; next; }
{ if (zz < 1) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r]=1; isGroup[r]=1; ngroups++; next; }}
{ if (zz < 2) { z=$1 ; r = r2i[z] ; if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; good[r]=1; isGroup[r]=0; next; }}
{ if (zz < 3) { z=$1 ; r = r2i[z] ; if(r>=1) private[r] = 1 ; next ;}}

{
  gsub("\"","",$0);
}

{ if (zz == 3)
    {
	nr0 = nr ;
	gsub (/System:/,"",$0);
	gsub (/NULL/,"NA",$0);
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

{r = 0 ; }

/^Ali /{ z=$2; ir = r2i[z] ; next; }
{ if (ir < 1) next ;}
/^Orphans/{ t= "Orphans_" $2;  if($2+0 > 1) t = "Orphans_multi" ; it = t2i[t] ; nnrt[ir,it] += $3 ; if (0) print "KKKKKK ",t,it,$3; next ; }
{ t = $1 ; it = t2i[t] ; if (it >0) {nnrt[ir,it] = $2 ; next ; }}


END {
    nr = nr0 ;
    want_all = 1 ;

    for (pass = 1 ; pass <= 2 ; pass++)
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
	
       	for (t = 1 ; t <= nt ; t++)
	    printf ("\t%s", i2t[t]) ;
    }

    for (r = 1 ; r <= nr ; r++)
    {
	if (private[r] == 1) continue ;
	for (pass = 1 ; pass <= 2 ; pass++)
	{
	    if (pass == 1)  printf ("\n%s",  i2r[r]) ;
	    else  printf ("\t\t\t%s",  i2r[r]) ;
	    
	    if (nSample  + want_all >= 0) printf("\t%s", sample[r]) ;
	    if (nRunid   + want_all >= 0)printf ("\t%s", runid[r]); 
	    
	    if (nTitre   + want_all >= 0)  printf ("\t%s", titre[r]); 
	    if (nsTitre  + want_all >= 0) printf ("\t%s", stitre[r]); 
	    if (nsTitre2 + want_all >= 0)  printf ("\t%s", stitre2[r]); 
	    if (noTitre  + want_all >= 0) printf ("\t%s", otitre[r]); 
	    
	    if (nTissue  + want_all >= 0) printf ("\t%s", tissue[r]); 
	    if (nSystm   + want_all >= 0) printf ("\t%s", systm[r]); 
	    # if (nMachine + want_all >= 0) printf("\t%s", machine[r]) ;
	    
	    for (t = 1 ; t <= nt ; t++)
	    {
		if (pass == 1)
		    printf ("\t%.0f",1.0*nnrt[r,t]) ;
		else
		{
		    z = nnrt[r,1] ; if (z == 0) z = 1 ;
		    printf ("\t%.2f", 100*nnrt[r,t]/z) ;
		}   
	    }
	}  
    }
}


           
  
