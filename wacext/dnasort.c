/*************************************************************************************/

static int parseFastaFile (JP *jp, long unsigned int *nbp, BOOL isRead)
{
  AC_HANDLE h = ac_new_handle () ;
  long unsigned int nTarget, nb = 0 ;
  int newNam ;
  SEQ *seq = 0 ;
  BOOL last = FALSE ;
  char *cp, *cr ;
  Array seqs = 0 ;
  BigArray dna =  isRead ? jp->pDnaArray : jp->tDnaArray ;
  BigArray dna2chrom = 0 ;
  ACEIN ai ;
  BOOL justStrand = jp->justStrand ;
  DICT *dict = isRead ? jp->readDict : jp->targetDict ;
  long unsigned int NN = isRead ? jp->pNN : jp->NN ;
  int nClipped = 0, nbClipped = 0 ;
  
  if (isRead)
    {
      dict = jp->readDict = dictHandleCreate (256, jp->h) ;
      seqs = jp->reads  = arrayHandleCreate (100000, SEQ, jp->h) ;
      dna = jp->pDnaArray = bigArrayHandleCreate (10000000, char, jp->h) ;
      ai = aceInCreate (jp->pFileName, jp->gzi , h) ;
      if (!ai)
	usage ("Sorry i cannot open the read fasta file") ;
    }
  else
    {
      dna2chrom = jp->dna2chrom = bigArrayHandleCreate (1000, int, jp->h) ;
      dict = jp->targetDict = dictHandleCreate (256, jp->h) ;
      seqs = jp->chroms = arrayHandleCreate (100, SEQ, jp->h) ;
      dna = jp->tDnaArray = bigArrayHandleCreate (1000, char, jp->h) ;
      if (NN == 0)
	{
	  jp->NN0 = NN = 1 ;  /* leave a 1b gap where we can stick the individual reads */
	  bigArray (dna, NN << 4, char) = 0 ;
	}
      ai = aceInCreate (jp->tFileName, 0, h) ;
      if (!ai)
	usage ("Sorry i cannot open the target fasta file") ;
    }

  nTarget = arrayMax (seqs) ;
  
  while (TRUE)
    {
      if (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (!cp || !*cp)
	    continue ;
	  {  /* get rid of fastc format */
	    char *cq = strstr (cp, "><") ;
	    if (cq) *cq = 0 ;
	  }
	}
      else
	last = TRUE ;
      if (last || *cp == '>')
	{
	  /* register the read and
	   * add a null base to propagate correctly bacwards
	   * the clause w[i] cannot be repeated if (w-1)[i+1] is not
	   */
	  if (!last)
	    dictAdd (dict, cp+1, &newNam) ;
	  if (seq)
	    {
	      if (isRead)
		{
		  bigArray (dna, NN++, char) = 0 ; /* 00 terminate each sequence */
		  bigArray (dna, NN++, char) = 0 ;
		}
	      else
		{
		  bigArray (dna, NN++, char) = 0 ; /* 00 terminate each sequence */
		  bigArray (dna, NN++, char) = 0 ;
		  if (NN % 1024) /* start on a new page and associate this page to the current newName */
		    NN = NN - (NN % 1024) + 1024 ;
		  bigArray (dna2chrom, NN >> 10, int) = nTarget ;
		}
	      seq->end = NN - 1 ;
	      nb += seq->ln ;
	    }

	  if (seq && seq->ln && isRead)
	    {
	      if (seq->ln <  jp->readLengthMin)
		jp->readLengthMin = seq->ln ;
	      if (seq->ln > jp->readLengthMax)
		jp->readLengthMax = seq->ln ;
	    }
	  
	  if (seq && isRead && jp->clipPolyA)
	    nClipped += jumpAlignSoftClipPolyA (jp, seq, &nbClipped) ;
	  if (seq && (isRead || jp->exportRepeats)  && !justStrand)
	    { /* add the strand inverted read */
	      char *cr ;
	      long unsigned int j, n0 = NN ;
	      

	      bigArray (dna, NN + seq->ln + 2, char) = 0 ; /* make room */
	      for (cr = arrp (dna, seq->start + seq->ln - 1, char), j = seq->ln ; j-- ; cr--)
		bigArray (dna, NN++, char) = complementBase[(int)*cr] ;
	      bigArray (dna, NN++, char) = 0 ;
	      bigArray (dna, NN++, char) = 0 ;
	      if (NN % 1024) /* start on a new page and associate this page to the current newName */
		NN = NN - (NN % 1024) + 1024 ;
	      seq = arrayp (seqs, nTarget++, SEQ) ;
	      seq->nam = (seq-1)->nam ;
	      seq->isDown = FALSE ;
	      seq->ln = (seq-1)->ln ;
	      seq->start = n0 ;
	      seq->end = NN - 1 ;  /* last letter of the read (inclusive ? ) */
	      nb += seq->ln ;
	      bigArray (dna2chrom, NN >> 10, int) = nTarget ;
	    }
	  if (last)
	    break ;
	  seq = arrayp (seqs, nTarget++, SEQ) ;
	  seq->nam = newNam ;
	  seq->isDown = TRUE ;
	  seq->start = NN ;
	}
      else
	{
	  seq->ln += strlen(cp) ;
	  for (cr = cp ; *cr ; cr++)
	    bigArray (dna, NN++, char) = dnaEncodeChar [(int)*cr] ;
	}
    }
  *nbp = nb ;

  if (isRead && jp->reads) 
    arraySort (jp->reads, seqOrder) ; /* sort by read length */

  if (isRead) 
    jp->pNN = NN ;
  else
    {
      long int k, kMax = NN/1024 ;
      int chrom, *ip ;

      jp->NN = NN ;
      
      /* complete the dna2chrom table */
      chrom = bigArray (dna2chrom, kMax - 1, int) ;
      for (k = 0, chrom = 0, ip = bigArrayp (dna2chrom, 0, int) ; k < kMax ; k++, ip++)
	if (*ip) 
	  chrom = *ip ;
	else
	  *ip = chrom ;
    } 
  ac_free (h) ;
  
  fprintf (stderr, "// %s: Parsed %lu %s, %lu bases\n", timeShowNow (), nTarget, isRead ? "fragments" : "targets", nb) ;
  if (jp->clipPolyA)
    fprintf (stderr, "//   clipped %d bases in %d reads\n", nbClipped, nClipped) ;
  return nTarget ;
} /* parseFastaFile  */

