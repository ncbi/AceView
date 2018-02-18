/^ZZZZZ/ { 
  zz++ ; next ; 
}
{ if (zz<1 && $2 + 0 > 5)
  {
    gm = $1 ;
    scoregm[gm] = $2 ;
    next ;
  }
}
{ if (zz == 1)  # register all gene and mrna coordinates 
  {
    m = $1 ; g = $5 ; 

    if (GM == "GENE" || GM2 == "m2g") 
      gm = g ;
    else 
      { 
	gm = m ;
	m2g[m] = g ;
      }
    
    s = scoregm[g] + scoregm[m] ;
    if (s > score[gm])
      score[gm] = s ;

    split ($2, aa, ":") ; 
    chrom = aa[1] ; 
    gg[gm] = 1 ;
    g2chrom[gm] = chrom ;
    if (length($6)>0) gid2g[g] = $6 ;

    split (aa[2], bb, "-") ; a1 = bb[1] ; a2 = bb[2] ;
    if (a1 < a2) a1 = a2 ;

    if (a1 > 0 + g2a[gm])
      g2a[gm] = a1 ;

    next ;
  }
}
{ if (zz == 2)
  {
    g = $1 ; gid = $2 ;
    gCancer[g] = 1 ;
    gidCancer[gid] = 1 ;
    next ;
  }
}
END {
  # count the genes per chromosomes
  for (g in gg)
    {
      nnn++ ;
      chroms[g2chrom[g]]++ ;
      if (0 + score[g] > 0)
	{
	  nnk++ ;                # total number of genes with score
	  dTscore += score[g] ;  # cumulated score 
	}
    }
# the desired number of genes per section is simply
  nng = nnn / 1000 ;              # 300 sections is a good guess
  nng = 50 * int(nng/50) ;     # rounding
  if (nng < 50) nng = 50 ;     # cannot compute stats on less than 100 objects
  if (nng > 300) nng = 100 * int(nng/100) ; 

# in each chromosome, order the genes per coordinates
# allocate each gene to its block

# modify the number to string conversion, otherwise asort does NOT sort numerically
  old=CONVFMT ;  CONVFMT="%09d" ;

  for (chrom in chroms)
    {
      delete(kk) ; delete(nn) ;
      
      delete(sa1) ;  delete(sa2) ; delete(sn) ; # sections coordinates and counts
      delete(sk) ;                 # section number of DEGs 
      delete(sScore) ;             # section cumulated score
      delete(sGenes) ;             # list of genes::scores in section
      delete(sCancer) ;
      i = 1 ;
      for (g in gg)
	{
	  if (g2chrom[g] == chrom)
	    {
	      a = g2a[g] ; atxt = ""  (a + .1) ;
	      kk[i++] = atxt ;
	      # if (chrom==22) printf( "UUUUU\t%s\t%d\t%d\t%d\n", g,g2a[g],i,kk[i-1]);
	    }
	}
      nGenes = asort(kk) ;  chromStart = kk[1] ; chromLength = kk[nGenes] ;

      # establish the sections, count up to nng genes, break on holes > 2Mb
      jMax = 0 ; 
      k = 1 ; a = chromStart ;  a1old = a ; a2old = a
      for (i = 2 ; i <= nGenes ; i++)
	{
	  a = kk[i] ; k++ ;
	  if (i == nGenes || k >= nng || a - a2old > 2000000)
	    { 
	      if (2*k >= nng || k >= 50)  # create a new section
		{
		  jMax++ ; sa1[jMax] = a1old ; sa2[jMax] = a2old ; sn[jMax] = k ;
		}
	      else if (a1old - sa2[jMax] <  2000000)   # increment the previous local section
		{
		  sa2[jMax] = a2old ; sn[jMax] += k ;
		}
	      else   # drop these few genes
		{

		}
	      k = 1 ; a1old = a ; 
	    }
	  a2old = a ;
	}
      
      # count the DEGs in each section
      for (g in gg)
	if (g2chrom[g] == chrom)
	  {
#		      if (g == "MYCN") print "FFFFFFF g=" g "  score=" score[g] "    a=" 0 + g2a[g] ;
	    if (score[g] + 0 > 0)
	      {
		a = 0 + g2a[g] ;
		for (j = 1 ; j <= jMax ; j++)
		  if (a >= 0 + sa1[j] && a <= 0 + sa2[j])
		    {
		      sk[j]++ ;           # DEG in block
		      s = score[g] ;
		      sScore[j] += s ;
		      sGenes[j] = sGenes[j] "::" g ":" 1000 * s ;
		      
# if (g == "MYCNOS" ) printf("XXXX%s\t%d\tblock\t%d\t%d\t%d\t%s\tscore=%s\tnn\t%d\tnk\t%d\t%s\n",chrom,g2a[g],j,sa1[j],sa2[j],g,score[g],sn[j],sk[j],sGenes[j]);
		      break ;
		    }
	      }
	    if (GM == "GENE" || GM2 == "m2g") 
	      gm = g ;
	    else 
	      gm = m2g[g] ; 
	    gid = gid2g[gm] ;
	    # if (g == "MYCN.aAug10") print "PPPPPPP g=" g "   gm="gm  ;
	    if (gCancer[gm] == 1 || gidCancer[gid] == 1)
	      {
		a = 0 + g2a[g] ;
		for (j = 1 ; j <= jMax ; j++)
		  if (a >= 0 + sa1[j] && a <= 0 + sa2[j])
		    {
		      tt = ":" gm ":" ;
		      if (index(sCancer[j],tt) == 0)
			sCancer[j] = sCancer[j] tt ;
		      break ;
		    }
	      }
	  }
    
# report the distribution for this chromosome
      nnb = nnn - nnk ;
      for (jj = 1 ; jj <= jMax ; jj++)  # for each block in the chromosome 
	{
	  a = sk[jj] + 0 ;
	  b = sn[jj] - a ;
	  c = nnk - a ;
	  d = nnb - b ;
	  at = (nnk*sn[jj]) /nnn ;
	  bt = (nnb*sn[jj]) /nnn ;
	  ct = (nnk* (nnn-sn[jj]) ) /nnn ;
	  dt = (nnb* (nnn-sn[jj]) ) /nnn ;
	  chi2= (a-at) * (a-at) / (1+at) + (b-bt) * (b-bt) / (1+bt) + (c-ct) * (c-ct) / (1+ct) + (d-dt) * (d-dt) / (1+dt)  ;

	  a1 = sa1[jj] ;
	  a2 = sa2[jj] ;
	  printf ("%s\t%s\t%d\t%d\t%s\t%d",target,chrom,a1,a2,title,sign)  ;
	  printf ("\t%d\t%d\t%.2f\t%.2f\t%d",a,b,at,bt,a+b) ;
	  printf ("\t%.2f\t%.2f\t%.2f\t%d\t%s", 100*a/ (a+b+.1),100*nnk/nnn,chi2,sScore[jj],sCancer[jj])  ;

	  if (((GM == "GENE" || GM2 == "m2g") && chi2 > 19) || chi2 > 24) 
	    {
	      n = split(sGenes[jj], aa, "::") ; 
	      delete(cc) ; delete(dd) ;
	      for (i = 1 ; i <= n ; i++)
		{
		  split(aa[i], bb, ":") ;
		  g = bb[1] ; s = bb[2] ; 
                  stxt = "a"  (s + .1) g ; cc[stxt] = g ;
		 # print "DDDDDDD " stxt ;
		}
	      n = asorti(cc,dd) ;
	      j = dd[n] ; g = cc[j] ; s0 = score[g] ;
	      printf("\t%s(%.1f)", g, s0) ;
	      for (i = n-1  ; i >= 2 ; i--)
		{
		  j = dd[i] ; g = cc[j] ; s = score[g] ;
		  if (2 * s >= s0)
		    printf(", %s(%.1f)", g, s) ;
		  else
		    break ;
		}
	    }
	  else
	    printf("\t") ;

	  printf("\n") ; 
	      # printf("UUUU %.1f // %s\n", chi2, bGenes[jj]) ;
	}
    }
}
