/^Average allele frequency/{for(i=2;i<=NF;i++)runs[i]=$i;next;}
{ t = 0 ; k = 0 ;}

/^All sites/{t="All_sites";k=0;}
/^Sites not measurable/{t="Not_measurable_sites";k=1; t0 = t ;}
/^Rejected/{t="Rejected_sites";k=1;}
/^Well measured/{t="Genomic";k=1;}
/^Reference/{t="Genomic";k=2;}
/^Low/{t="Genomic";k=3;}
/^Mid/{t="Genomic";k=4;}
/^High/{t="Genomic";k=5;}
/^Pure/{t="Genomic";k=6;}

/^Protein changing Sites/{t="Protein_changing";k=1;}
/^Protein changing Reference/{t="Protein_changing";k=2;}
/^Protein changing Low/{t="Protein_changing";k=3;}
/^Protein changing Mid/{t="Protein_changing";k=4;}
/^Protein changing High/{t="Protein_changing";k=5;}
/^Protein changing Pure/{t="Protein_changing";k=6;}

/^exonic Sites/{t="Exonic";k=1;}
/^exonic Reference/{t="Exonic";k=2;}
/^exonic Low/{t="Exonic";k=3;}
/^exonic Mid/{t="Exonic";k=4;}
/^exonic High/{t="Exonic";k=5;}
/^exonic Pure/{t="Exonic";k=6;}

{if (t != 0){ if(k > tags[t]) tags[t]=k;for(i=2;i<=NF;i++)z[t,k,i]=$i;} next;}
END {
  kk[1] = "any" ; kk[2] = "reference" ; kk[3] = "low" ; kk[4] = "mid" ; kk[5] = "high" ; kk[6] = "pure" ;
  for (i in runs) 
    {
      printf("Ali \"%s\"\n", runs[i]); 
      for (t in tags)
        {
	  if (t == "All_sites") 
            {
              printf ("Tested_sites") ;
                printf (" %d ", z[t0,1,i] + z[t,0,i] + 0) ;
            }
          else
            {
	     printf ("%s ", t) ;
             for (k = 1 ; k <= tags[t] ; k++)
               printf (" %d %s ", z[t,k,i] + 0, kk[k]) ;
             } 
           printf ("\n") ;
        }
      
      printf ("\n") ;
    }
}
