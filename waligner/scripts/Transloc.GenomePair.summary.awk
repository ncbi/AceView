# construct a plot of the local coverage of the exons and introns around the break point

BEGIN { delta = 150 ; }
/^ZZZZZ/{zz++;next;}
# zz == 0 all.genome.introns, grab orientation and coordinates of the last line
{ if (zz<1) { nnn++ ; if(nnn>0){ fr = $2 ; u1 = $7 ; u2 = $9 ; if(u1>u2){u0=u1;u1=u2;u2=u0;}support = $13 ;} }
}

# zz == 1 tmp/METADATA/mrnaRemap.av.txt, grab the gene structure
{ 
  if (zz == 1) 
    {
      if ($8 == gene1) 
	{
	  ng[1]++ ; chroms[1] = $5 ; 
	  ga1[1,ng[1]] = $6 ; ga2[1,ng[1]] = $7 ;  
	}     
     if ($8 == gene2) 
	{
	  ng[2]++ ; chroms[2] = $5 ; 
	  ga1[2,ng[2]] = $6 ; ga2[2,ng[2]] = $7 ;  
	}
    }
 
}

# zz == 2 wiggle in BV format
{
  if (zz == 2)
    {
      x = int ($1) ;
      if (x > 0 && ((x >= u1 - delta && x <= u1 + delta)  || (x >= u2 - delta && x <= u2 + delta) ))
	{
	  ww[x] = $2 ;
	  if (xmin == 0) xmin = x ;
	  xmax = x ;
	}
    }
}

# zz == 3 fasta file of the recombined region
{
  if (zz == 3)
    {
      if(substr($1,1,1)== ">") next ;
      dna = dna $1 ;
    }
}

# export
END {
  split (map1,mm,":") ; split(mm[4],bb,"_") ; gg1a1 = bb[1] ; gg1a2 = bb[2] ;
  split (map2,mm,":") ; split(mm[4],bb,"_") ; gg2a1 = bb[1] ; gg2a2 = bb[2] ;
  s1 = s2 = 1 ;
  if (gg1a1 > gg1a2) s1 = -1 ; 
  if (gg2a1 > gg2a2) s2 = -1 ; 

  if (fr == "F1.F2")
    {  ig = 1 ; x01 = gg1a1 ; x02 = gg2a1 ; s1 = s1 ; s2 = s2 ; NN1 = N1 ; }
  if (fr == "F2.F1")
    {  ig = 2 ; x01 = gg2a1 ; x02 = gg1a1 ; s1 = s1 ; s2 = s2 ; NN1 = N2 ; }
  if (fr == "F1.R2")
    {  ig = 1 ; x01 = gg1a1 ; x02 = gg2a2 ; s1 = s1 ; s2 = -s2 ; NN1 = N1 ; }
  if (fr == "R1.F2")
    {  ig = 1 ; x01 = gg1a2 ; x02 = gg2a1 ; s1 = -s1 ; s2 = s2 ; NN1 = N1 ; }

  genes[1] = gene1 ; genes[2] = gene2 ;
  strands[1] = s1 ; strands[2] = s2 ;
  stype[1] = ">" ; stype[-1] = "<" ;
  printf ("gene1=%s gene2=%s FR=%s u1=%d u2=%d NN1=%d x01=%d x02=%d gg1a1=%d/%d gg2a1=%d/%d s1=%d s2=%d xmin=%d xmax=%d\n",gene1,gene2,fr,u1,u2,NN1,x01,x02,gg1a1,gg1a2,gg2a1,gg2a2,s1,s2,xmin,xmax);
  
  ig1 = ig ; ig2 = 3 - ig ;
  for (x = xmin ; x <= xmax ; x += step)
    {
      if (x > u1 + delta && x < u2 - delta)
	{ kk = 1 ; continue ; }
      if (kk == 1) { kk = 0 ; for (i = 1 ; i <= 5 ; i++) print ("\n") ;}
      exon1 = exon2 = intron = 0 ;
      if (x < ND)
	{
	  chrom = chroms[ig1] ;
	  gene = genes[ig1] ;
	  strand = strands[ig1] ;
	  y = x01 + (x-300) * s1 ;
	  if (s1 == -1) y++ ;
	  for (i = 1 ; i <= ng[ig1] ; i++)
	    if ((y >= ga1[ig1,i]  && y <= ga2[ig1,i]) || (y >= ga2[ig1,i]  && y <= ga1[ig1,i]))
	      {
		exon1 += 10 ;
	      }
	}
      else
	{
	  chrom = chroms[ig2] ;
	  gene = genes[ig2] ;
	  strand = strands[ig2] ;
	  y = x02 + (x - 300 - NN1 - Nn) * s2 ;
          if (s2 == -1) y ++ ;
	  for (i = 1 ; i <= ng[ig2] ; i++)
	    {
	      if ((y >= ga1[ig2,i]  && y <= ga2[ig2,i]) || (y >= ga2[ig2,i]  && y <= ga1[ig2,i]))
		{
		  exon2 += 10 ;
		  if(0)printf("x=%d y=%d exon2=%d i=%d ng=%d a1=%d a2=%d\n",x,y,exon2,i,ng[ig2], ga1[ig2,i],ga2[ig2,i]) ;
		}
	    }
	}
      if (x > u1 && x < u2)
	intron = support ;
      printf ("%s%s\t%s\t%d\tc=%s\t%d\t%d\t%d\t%d\t%d\n", gene, stype[strand], substr(dna,x,1),x, chrom,y, ww[x],exon1,exon2, intron) ;
    }
}
