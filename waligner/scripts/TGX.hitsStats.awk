# parse a magic hits file, to obtain basic statistics
/^#/ { next}

{ p = $1 ;              #column 1 is the read identifier 
  if (p == oldp) next ;   # scan each read only once (the hits file is sorted 
  oldp = p ;

  score = $2 ;      # score, drop it
  mult = $3 ;       # multiplicity = number of time this exact sequence was seen in the fastq file
  
  nAli += mult ;
  lnAli += $5 ;     # length of the alignment
  nAliBp += mult * lnAli ;    # length of the alignment
  uu = $10 ;              # number of locus the read is aligning to
  if (uu == -2) uu = 1 ;  # special case of one genomic locus 2 genes in antisense

  nerr = 0+$15 ;      # number of mismatches in this alignment
  nM[nerr] += mult ;
  if (uu == 1) {
    nAliBpU += mult * lnAli ; 
    nMU[nerr] += mult ;

  # scan the error types
    if (0 + $15 > 0)
      {
	n = split ($16,aa,",") ;
	for (i = 1 ; i <= n ; i++)
	  {
	    split (aa[i], bb, ":") ;
	    errType[bb[2]] += mult ;
	  }
      }
  }
}
END {       # export all data repeating the run name, to simplify collating
  # Number of aligned reads
  printf ("%s\tAligned_reads\t%d\n", run, nAli) ;
  printf ("%s\tAligned_Bases\t%d\n", run, nAliBp) ;
  printf ("%s\tUniquely_aligned_reads\t%d\n", run, nAliBpU) ;
  # Number of reads with up to 10 mismatches
  printf ("%s\tReads_with_0_to_10_missmatches", run) ;
  for (i=0 ; i<= 10 ; i++) printf ("\t%d", 0+nM[i]) ;  
  printf ("\n%s\tUniquely_aligned_reads_with_0_to_10_missmatches", run) ;
  for (i=0 ; i<= 10 ; i++) printf ("\t%d", 0+nMU[i]) ;  
  printf ("\n") ;

  for (t in errType)
    printf("%s\tErrType\t%s\t%d\n", run, t, errType[t]) ;
}
