{gsub (/\"/,"",$0);}
{gsub(/gi\|251831106\|ref\|NC_012920.1\|/,"chrM",$0);}
/^Variant/{
  v=$2; 

  split (v, aa, ":") ;
  c = "chr" ; ch = aa[1] ;
  if (substr (ch, 1, 3) == "chr") c = "" ; 

  z = aa[2] ; imax = length(z) ;
  for (i = 1 ; substr(z, i, 1) >= "0" && substr(z, i, 1) <= "9" && i < imax ; i++) ;
  pos = substr (z, 1, i - 1) ;
  chrom[v]= c ch "\t" pos ; 

  z = substr (z, i) ;
  if (substr(z,2,1) == "2")
    {
      a[v] = substr(z,1,1) ; b[v] = substr(z,3,1) ;
    }
  if (substr(z,1,3) == "Ins")
    {
      a[v] = "." ; b[v] = substr(z,4) ; 
    }
  if (substr(z,1,3)== "Del")
    {
      a[v] = substr(z,4) ; b[v] = "." ;
    }
  next;
}
/^fCountsNU/{mutpNU[v]=$3;refpNU[v]=$4; ambpNU[v]=$5;next;} # ATTENTION there may be multiple fcounts lines
/^rCountsNU/{mutmNU[v]=$3;refmNU[v]=$4; ambmNU[v]=$5;next;} # because of the constructed type behind the numbers
/^fCounts/{mutp[v]=$3;refp[v]=$4; ambp[v]=$5;next;} # ATTENTION there may be multiple fcounts lines
/^rCounts/{mutm[v]=$3;refm[v]=$4; ambm[v]=$5;next;} # because of the constructed type behind the numbers
END {
  printf("# Chromosome\tPosition\tSample\tMethod\tReference\tVariant\tGenotype\tFrequency\tCoverage\tVariant=\tRef=\tV+\tR+\tV-\tR-\tAmbiguous\n");  
  hello = snpListFileName ;
  for(v in chrom)
    { 
      if (length(b[v]) == 0)
	  continue ;

      mutNU = 0 + mutpNU[v] + mutmNU[v] ;
      wildNU = 0 + refpNU[v] + refmNU[v]
      totNU = mutNU + wildNU ;

      if (totNU >=  mincover && mutNU >= 5 && 100 * mutNU >= 40 * totNU)
	{
	  mutp[v] += mutpNU[v] ; refp[v] += refpNU[v] ; ambp[v] += ambpNU[v] ; 
	  mutm[v] += mutmNU[v] ; refm[v] += refmNU[v] ; ambm[v] += ambmNU[v] ; 
	}

      mut = 0 + mutp[v] + mutm[v] ;
      wild = 0 + refp[v] + refm[v]
      tot = mut + wild ;
      if (tot >= mincover && mut >= minmutant && 100 * mut >= limit * tot)
      {
        printf ("Variant \"%s\"\n", v) > snpListFileName ;


        printf ("%s\t%s\tAceView\t%s\t%s\t%s", chrom[v],sample,a[v],b[v],"-") ;
        printf ("\t%d\t%d", int(100*(mut/tot)), tot) ;
        printf ("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", mut, wild, 0 + mutp[v], 0 + refp[v], 0 + mutm[v], 0 + refm[v], 0 + ambp[v] + ambm[v], 0 + ambp[v], 0 + ambm[v]) ;
        if (3 * (0+ambp[v]+ambm[v]) > tot)
          printf ("\tO_O") ;
	else 
          printf("\t-") ;
        if (tot > 50)
          printf ("\tBig") ;
	else 
          printf("\t-") ;
        printf ("\n") ;
      }
    }
}
