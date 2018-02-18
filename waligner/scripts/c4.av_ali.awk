BEGIN { 
  nt = 0 ; nr=0;
  m = "any" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Average_length_aligned" ;
  m = "any1" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Average_length_aligned_on_forward_strand_of_first_best_target"
  m = "any2" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Average_length_aligned_on_reverse_strand_of_first_best_target"
  m = "B_rrna" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "rRNA" ;
  m = "A_mito" ; 
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "mitochondria" ;
  m = "C_chloro" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "chloroplast" ;
  m = "ET_av" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "AceView_2010" ;
  m = "KT_RefSeq" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "RefSeq_37.104_models" ;
  m = "LT_seqc" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "SEQC_Transcriptome_2013" ;
  m = "LT_UCSC" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "ucscGene_Nov_2011" ;
  m = "MT_EBI" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Encode.37.70_Jan2013" ;
  m = "NT_MiT" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "MichiganFeb2015" ;
  m = "OT_rnaGene" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "RNA_gene" ;
  m = "PT_tRNA" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "tRNA" ;
  m = "QT_smallRNA" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "smallRNA" ;
  m = "1_DNASpikeIn" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "DNA_SpikeIn"
  m = "0_SpikeIn" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "RNA_SpikeIn"
  m = "Z_genome" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Continuous genome alignment" ;
  m = "z_gdecoy" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "Decoy_genome_specificity_control" ;
}
/^ZZZZZ/{zz++ ; next ; }
{ if (zz < 1) {  r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} next; }}
/^\"/ { 
  gsub(/\"/,"",$0) ;
  gsub(/\\\//,"/",$0) ;
  r = $1 ; ir=r2ir[r];if(ir<1){nr++;r2ir[r]=nr;rr[nr]=r;ir=nr;} 
  if (zz == 2)
    {
      if (0+$5 < 1000) next ;
      t = $2 ; 
      it=t2i[t];if(it<1){nt++;t2i[t]=nt;tt[nt]=t;it=nt;} 
      nrt[ir,it] = $3 ; ntt[it] += $3;
      ln[ir] = $4 ;
    }
  else
    {
      gsub (/System:/,"",$0);
      gsub(/NULL/,"NA",$2) ; gsub(/ /,"_",$2) ; machine[ir]=$2; if ($2 != "NA") ok["machine"] = 1 ;
      gsub(/NULL/,"NA",$3) ; gsub(/ /,"_",$3) ; if(sample[ir]) { if (index(sample[ir],$3) < 1) sample[ir] = sample[ir] ", " $3 ; } else  sample[ir]= $3 ; if ($3 != "NA") ok["sample"] = 1 ;
      gsub(/NULL/,"NA",$4) ; gsub(/ /,"_",$4) ; if(systm[ir])  { if (index(systm[ir],$4) < 1) systm[ir]  = systm[ir]  ", " $4 ; } else  systm[ir] = $4 ; if ($4 != "NA") ok["systm"] = 1 ;
      gsub(/NULL/,"NA",$7) ; gsub(/ /,"_",$7) ; if(systm[ir])  { if (index(systm[ir],$7) < 1) systm[ir]  = systm[ir]  ", " $7 ; } else  systm[ir] = $7 ; if ($7 != "NA") ok["systm"] = 1 ;
      gsub(/NULL/,"NA",$5) ; gsub(/ /,"_",$5) ; if(tissue[ir]) { if (index(tissue[ir],$5) < 1) tissue[ir]  = tissue[ir]  ", " $5 ; } else  tissue[ir] = $5 ; if ($5 != "NA") ok["tissue"] = 1 ;
      gsub(/NULL/,"NA",$6) ; gsub(/ /,"_",$6) ; title[ir]=$6; if ($6 != "NA") ntitle = 1 ;
      gsub(/NULL/,"NA",$8) ; gsub(/ /,"_",$8) ; runid[ir]=$8; if ($8 != "NA") ok["runid"] = 1 ;
      gsub(/NULL/,"NA",$11) ; gsub(/ /,"_",$11) ; stitle[ir]=$11; if ($11 != "NA") nstitle = 1 ;
      gsub(/NULL/,"NA",$12) ; gsub(/ /,"_",$12) ; stitle2[ir]=$12; if ($12 != "NA") nstitle2 = 1 ;
      gsub(/NULL/,"NA",$13) ; gsub(/ /,"_",$13) ; stitle[ir]=$13; if ($13 != "NA") notitle = 1 ;
      next ;
    }
}
END { 
  printf ("\nRun") ;
  if (ok["runid"] == 1) printf("\tRunId") ;
  if (ok["machine"] == 1) printf("\tMachine") ;
  if (ok["systm"] == 1) printf("\tSystem") ;
  if (ok["tissue"] == 1) printf("\tTissue") ;
  if (ok["sample"] == 1) printf("\tSample") ;
  if (ntitle == 1) printf("\tTitle") ;
  if (nstitle == 1)  printf ("\tSorting Title") ;
  if (nstitle2 == 1)  printf ("\tSorting Title 2") ;
  if (notitle == 1)  printf ("\tOther Title") ;

  printf ("\tAverage adaptor clipped length") ;
  for (t = 1 ; t <= nt ; t++)
    {
      if (ntt[t] > 0)
	printf ("\t%s", i2nam[t]) ;
    }
  printf ("\t\tRun") ;
  if (ntitle == 1) printf("\tTitle") ;
  if (nstitle == 1)  printf ("\tSorting Title") ;
  if (nstitle2 == 1)  printf ("\tSorting Title 2") ;
  if (notitle == 1)  printf ("\tOther Title") ;

  for (t = 1 ; t <= nt ; t++)
    {
      if (ntt[t] > 0)
	printf ("\t%s %%", i2nam[t]) ;
    }
  for (r = 1 ; r <= nr ; r++) 
    {
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
      printf("\n%s", rr[r]) ;
      if (ok["runid"] == 1) printf("\t%s", runid[r]) ;
      if (ok["machine"] == 1) printf("\t%s", machine[r]) ;
      if (ok["systm"] == 1) printf("\t%s", systm[r]) ;
      if (ok["tissue"] == 1) printf("\t%s", tissue[r]) ;
      if (ok["sample"] == 1) printf("\t%s", sample[r]) ;
      if (ntitle == 1)    printf ("\t%s",title[r]) ;
      if (nstitle == 1)    printf ("\t%s",stitle[r]) ;
      if (nstitle2 == 1)    printf ("\t%s",stitle2[r]) ;
      if (notitle == 1)    printf ("\t%s",otitle[r]) ;

      printf("\t%s", ln[r]) ;
      for (t = 1 ; t <= nt ; t++)
	{
	  if (ntt[t] > 0)
	   {
	     if (nrt[r,t] > ln[r])
	       nrt[r,t]= ln[r] ;
	     if (nrt[r,t] > 0)
	       printf("\t%s", nrt[r,t]) ;  
	     else
	       printf ("\t") ;
	   }
	}
      printf("\t\t%s", rr[r]) ;
      if (ntitle == 1)    printf ("\t%s",title[r]) ;
      if (nstitle == 1)    printf ("\t%s",stitle[r]) ;
      if (nstitle2 == 1)    printf ("\t%s",stitle2[r]) ;
      if (notitle == 1)    printf ("\t%s",otitle[r]) ;

      for (t = 1 ; t <= nt ; t++)
	{
	  if (ntt[t] > 0)
	   {
	     if (nrt[r,t] > 0)
	       {
		 z=ln[r]; if (0+z<1)z=1;printf("\t%.2f", 100*nrt[r,t]/z) ;
	       }
	     else
	       printf ("\t") ;
	   }
	}
    }
  printf ("\n") ;
}
