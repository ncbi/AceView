/*
  Version 3 mieg 990419
  - Added capacity to swap a cosmid to the other strand
  Version 2 rd 980524
  - various cleanups
  - add -G option to insert 200 bp between St Louis juxtapositions
  - fix bug: sometimes overlap coords were removed incorrectly
  - remove search loop on overlap: spec says no search
  - introduced local olda2, oldb2 because otherwise the juxtaposition
	signal did not get interpreted if there was a length error
*/

/* note that in this code we use biologists coordinates
   that hh->a1,a2 = 1 length
   */

#include "../wac/ac.h"
#include "mytime.h"

typedef struct{ AceKey seq ; 
  int a1, a2, olda2, xx, dx, missMatchPos1, missMatchPos2, doSwap; 
  BOOL isUp, wrongDna, newMatch, missMatch, gap200 ; } SEQ ;

static AceTag _Subsequence = 0 ;

static BOOL fiddleGaps = FALSE, isMakeChrom = FALSE ;
static BOOL doSave = FALSE, isSTL = FALSE ;
static BOOL exportOverlaps = FALSE ;
static int minOverlap = 12 ;

/***************************************************************/

static void swapA1 (Array seqs)
{
  int j1 = arrayMax(seqs) ;
  SEQ *hh ;

  while (j1--)
    {
      hh = arrayp(seqs, j1, SEQ) ;
      if (hh->a1 > hh->a2)
	{ 
	  int tmp = hh->a1 ; 
	  hh->a1 = hh->a2 ; hh->a2 = tmp ; 
	  hh->isUp = ! hh->isUp ; 
	}
    }
}

/***************************************************************/

static void backSwapA1 (Array seqs)
{
  int j1 = arrayMax(seqs) ;
  SEQ *hh ;

  while (j1--)
    {
      hh = arrayp(seqs, j1, SEQ) ;
      if (hh->isUp)
	{ 
	  int tmp = hh->a1 ; 
	  hh->a1 = hh->a2 ; hh->a2 = tmp ; 
	  hh->isUp = ! hh->isUp ; 
	}
    }
}

/***************************************************************/

static int seqOrder (void *a, void *b) 
{
  SEQ *ga = (SEQ *)a, *gb = (SEQ *)b ;

  return 
    ga->a1 != gb->a1 ?
    ga->a1 - gb->a1 :
    ga->a2 - gb->a2 ;
}

/***************************************************************/

static void checkSeqs (Array seqs)
{
  SEQ *hh ;
  char *cp, *cq ;
  int dx, i, a1, a2, b1, b2, a3, b3, z3, z4, iSwap ;
  Array dna1 = 0, dna2 = 0 ;
  int olda2 = 0, oldb2 = 0 ;

  acePushContext () ;
  a1 = a2 = b1 = b2 = 0 ;
  for (i=0 ; i < arrayMax(seqs) ; i++)
    {
      hh = arrp (seqs, i, SEQ) ;
      if (b2 && b2 > a2)	/* because sometimes get internal subsequence */
	{ aceFree (dna1) ; dna1 = dna2 ; dna2 = 0 ;
	  a1 = b1 ; a2 = b2 ;
	  olda2 = oldb2 ;
	}
      else
	aceFree (dna2) ;
      b1 = hh->a1 ; oldb2 = b2 = hh->a2 ; 
      hh->olda2 = hh->a2 ;
      if (!aceKeyHasTag (hh->seq, aceTag ("DNA", 0)))
	continue ;
      dna2 = aceDnaGet (hh->seq, 0) ;
      if (!dna2)
	continue ;
      if (hh->isUp)
	reverseComplement(dna2) ;
      if (b2 != b1 + arrayMax(dna2) -1) /* safe precaution  !*/
	{ b2 = hh->a2 = b1 + arrayMax(dna2) -1 ;  
	  hh->wrongDna = TRUE ;
	}
      if (!a1 || !a2 || !b1 || !b2)
	continue ; 
      if (!dna1)
	continue ;
      if (b1 == olda2 + 1) /* juxtaposed in original */
	if (fiddleGaps)		/* NB only if DNA in both */
	  { hh->dx = 200 + a2 - olda2 ;
	    hh->gap200 = TRUE ;
	    continue ;
	  }
	else /* try to overlap correctly */
	  for (iSwap = 0 ; iSwap < 2 && !hh->newMatch ; iSwap++)
	    {
	      if (iSwap) reverseComplement (dna2) ;
	      for (dx = -5000 ; dx < -minOverlap ; dx++)
		{
		  if (b1 -a1 + dx < 0 ||
		      b1 -a1 + dx > arrayMax (dna1))
		    continue ;
		  b3 = arrayMax (dna2) ;
		  a3 = arrayMax (dna1) - b1 + a1 - dx ;
		  z3 = a3 < b3 ? a3 : b3 ;
		  cp = arrp (dna1, b1 - a1 + dx, char) ;
		  cq = arrp (dna2, 0, char) ;
		  while (z3--)
		    if (! (*cp++ & *cq++))
		      break ;
		  if (z3 == -1) 
		    { hh->newMatch = TRUE ;
		    hh->dx = dx ;
		    if (iSwap)
		      {hh->doSwap = TRUE ; 
		      }
		    break ;
		    }
		} 
	       if (iSwap && !hh->newMatch) 
		 reverseComplement (dna2) ;
	    }
      else if (b1 <= a2)	/* check if correct */
	{ z3 = (a1 < b1) ? b1 : a1 ;
	  cp = arrp (dna1, z3-a1, char) ;
	  cq = arrp (dna2, z3-b1, char) ;
	  while (z3 <= a2 && z3 <= b2)
	    { 
	      if (!(*cp++ & *cq++))
		{ /* try to swap */
		  reverseComplement (dna2) ;
		  z4 = (a1 < b1) ? b1 : a1 ;
		  cp = arrp (dna1, z4-a1, char) ;
		  cq = arrp (dna2, z4-b1, char) ;
		  while (z4 <= a2 && z4 <= b2)
		    {
		      if (!(*cp++ & *cq++))
			{ 
			  reverseComplement (dna2) ;
			  hh->missMatch = TRUE ;
			  hh->missMatchPos1 = cp - arrp(dna1,0,char) ;
			  hh->missMatchPos2 = cq - arrp(dna2,0,char) ;
			  hh->dx = a2 - b1 + 1 ; /* create a juxtaposition */
			  break ;
			}
		      z4++ ;
		    }
		  if (!hh->missMatch)
		    hh->doSwap = TRUE ;
		  break ;
		}
	      ++z3 ;
	    }
	} 
    }
  
  acePopContext (TRUE) ;
}

/***************************************************************/

static void shiftSeqs (Array seqs)
{
  SEQ *hh ;
  int xx = 0 , i ;

  for (i=0 ; i < arrayMax(seqs) ; i++)
    {
      hh = arrp (seqs, i, SEQ) ;
      xx += hh->dx ;
      hh->xx = xx ;
    }
}

/***************************************************************/

static void printOverlaps (AceKey link, Array seqs)
{
  SEQ *hh, *hh1, *hh2 ;
  int i = arrayMax (seqs), a1, a2, b1, b2, delta, tmp ;
     
  for (i=0 ; i < arrayMax(seqs) - 1 ; i++)
    {
      hh = arrp (seqs, i, SEQ) ;
      hh1 = arrp (seqs, i + 1, SEQ) ;
      a1 = hh->a1; a2 = hh->a2 ;
      b1 = hh1->a1; b2 = hh1->a2 ;
      if (a1 > a2) { tmp = a1 ; a1 = a2 ; a2 = tmp ; }
      if (b1 > b2) { tmp = b1 ; b1 = b2 ; b2 = tmp ; }
      delta =   b1 - a1 + 1 ;
      printf ("Sequence %s\nOverlap_right %s %d",
	       aceName(hh->seq), aceName(hh1->seq), 
	       b1 - a1 + 1) ;
      if (b1 > a2 + 1)
	{
	  hh2 = 0 ;
	  if (i > 0)
	   {
	     hh2 = arrp (seqs, i - 1, SEQ) ;
	     if (hh2->a2 >= b1 - 1)
	       printf (" // Included in %s",  aceName(hh2->seq)) ;
	     else 
	       hh2 = 0 ;
	   }
	  if (!hh2)
	    printf (" // Gap %d bp  (a1=%d a2=%d  b1=%d b2=%d)",
		  b1 - a2 - 1,a1,a2,b1,b2) ;
	}
      else if (b1 <= a2)
	printf (" // DoubleCover %d bp", -b1 + a2 + 1) ;
      else if (b1 == a2 + 1)
	printf (" // Juxtaposition") ;
      printf("\n\n") ;
    }
}

/***************************************************************/
/* if link = 0, as frommakeChrom, the printout adds to a table
   otherwise the printout is an acefile
   */

static void printSeqs (AceKey link, Array seqs)
{
  SEQ *hh ;
  int ii = arrayMax (seqs) ;
  static AceKey oldLink = 0 ;
     
  for (ii=0 ; ii < arrayMax(seqs) ; ii++)
    {
      hh = arrp (seqs, ii, SEQ) ;
      if (!hh->xx && !hh->wrongDna && !hh->doSwap)
	{
	  if (ii > 0 && !link) 
	    {
	      int nn ;
	      SEQ *hh0 = hh - 1 ;

	      
	      if (hh0->a1 < hh0->a2) nn = hh0->a2; else nn = hh0->a1 ;
	      if (hh->a1 < hh->a2) nn -= hh->a1; else nn -= hh->a2 ;
              nn++ ;
	      if (nn > 0)
		printf ("overlap_verified %d",nn) ;
	      else if (nn==0)
		printf ("juxtapositin") ;
	      else if (nn<0)
		printf ("gap %d",-nn) ;
	    }
	  continue ;
	}
      if (link && link != oldLink)
	printf ("\nSequence %s\n", aceName(link)) ;
      if (!hh->a1 || !hh->a2)
	{
	 if (link)
	   printf ("-D Subsequence %s  // was %d %d  shift: %d\n", 
		   aceName(hh->seq), 
		   hh->a1, hh->a2, hh->xx) ;
	 else
	   printf("shift %d",hh->xx) ;
	}
      else
	{ 
	  int a1 = hh->a1, a2 = hh->a2 ;
	  if (hh->doSwap)
	    { int tmp = a1 ; a1 = a2 ; a2 = tmp ; }
	  if (link)
	    printf ("Subsequence %s %d %d // was %d %d  shift: %d", 
		  aceName(hh->seq), a1 + hh->xx, a2 + hh->xx,
		  hh->a1, hh->olda2, hh->xx) ;
	  else
	    printf("shift %d",hh->xx) ;
	  if (hh->wrongDna)
	    printf (" wrong_length (%s%d)", 
		    (hh->a2 > hh->olda2) ? "+" : "", 
		    hh->a2 - hh->olda2) ;
	  if (hh->newMatch)
	    printf (" new_match") ;
	  if (hh->doSwap)
	    printf (" swapping") ;
	  if (hh->gap200)
	    printf (" 200bp_gap") ;
	  if (hh->missMatch)
	    printf (" missmatch_at_pos %d / %d of_previous/this_dna",
		    hh->missMatchPos1, hh->missMatchPos2) ;
	  if (link) 
	    printf ("\n") ;
	}
      oldLink = link ;
    }
}

/***************************************************************/

static void saveSeqs (AceKey link, Array seqs)
{
  SEQ *hh ;
  int i = arrayMax (seqs) ;
  AceInstance Link ;

  acePushContext () ;
  
  if ((Link = aceUpdateKey (link, 0)))
    {
      if (aceGotoTag (Link, _Subsequence))
	aceRemove(Link) ;
      
      for (i=0 ; i < arrayMax(seqs) ; i++)
	{
	  hh = arrp (seqs, i, SEQ) ;
	  if ((hh->a1 || hh->a2) &&
	      aceAddTag (Link, _Subsequence) &&
	      aceAddKey (Link, hh->seq) &&
	      aceAddInteger (Link, hh->a1 + hh->xx))
	    aceAddInteger (Link, hh->a2 + hh->xx) ;
	}
    }
  else 
    fprintf(stderr, "// %s \n", aceErrorMessage (0)) ;
  acePopContext (TRUE) ;
}

/***************************************************************/

static void fixLink (KEY link, BOOL exportOverlap)
{
  AceInstance Link ;
  static Array units = 0 ;
  SEQ *hh ; AceUnit *up ;
  int ii, j1 ;
  Array seqs = 0 ;
  AceKey key ;

  if (!_Subsequence)
    _Subsequence = aceTag ("Subsequence", 0) ;

  acePushContext () ;

  if (!(Link = aceOpenKey(link, 0)))
    goto abort ;

  if (!units) units = arrayCreate (300, AceUnit) ;

  aceGetArray (Link, _Subsequence, units, 3) ;
  seqs = aceArrayCreate (50, SEQ) ;
  for (j1 = ii = 0 ; ii < arrayMax(units) ; ii+= 3)
    {
      up = arrp(units, ii, AceUnit) ;
      key =  up[0].k ;
      if (!isSTL && !aceKeyHasTag (key, aceTag ("Genomic_canonical", 0)))
	continue  ;
      hh = arrayp(seqs, j1++, SEQ) ;
      hh->seq = up[0].k ;
      hh->a1 = up[1].i ;
      hh->a2 = up[2].i ;
      hh->isUp = FALSE ;
    }

  swapA1 (seqs) ; 
  arraySort (seqs, seqOrder) ;

  checkSeqs (seqs) ;
  shiftSeqs (seqs) ;

  backSwapA1 (seqs);

  if (doSave)
    saveSeqs (link, seqs) ;
  else if (exportOverlap)
    printOverlaps (link, seqs) ;
  else
    printSeqs (link, seqs) ;
    
 abort:
  arrayDestroy (seqs) ;
  acePopContext (TRUE) ; /* should save and messfree */
}

/***************************************************************/

static void fixJunction (KEY junction)
{
  AceInstance Link = 0 ;
  static Array units = 0 ;
  AceKey link, seq ;
  AceUnit *up ;
  int ii, j1, a1, a2, b1, b2, u1, u2 ;
  AceKeySet parts = 0, sups ;
  int new = 0 ;

  if (!_Subsequence)
    _Subsequence = aceTag ("Subsequence", 0) ;

  acePushContext () ;

  parts =  aceQueryKey(junction, "Follow parts",0) ;
  sups =  aceQueryKey(junction, "Follow source",0) ;
  if (!keySetMax(parts) || keySetMax(sups) != 1)
    goto abort ;
  link = keySet (sups, 0) ;
  
  if (!(Link = aceOpenKey(link, 0)))
    goto abort ;

  if (!units) units = arrayCreate (300, AceUnit) ;

  b1 = b2 = u1 = u2 = 0 ;
  aceGetArray (Link, _Subsequence, units, 3) ;
  new = 1 ;
  for (ii = 0 ; ii < arrayMax(units) ; ii+= 3)
    {
      up = arrp(units, ii, AceUnit) ;
      seq = up[0].k ;
      a1 = up[1].i ;
      a2 = up[2].i ;
      if (seq == junction)
	{ b1 = a1 ; b2 = a2 ; } /* old values */
      for (j1 = 0 ; j1 < keySetMax(parts) ; j1++)
	if (seq == keySet(parts, j1))
	  {
	    if (new) { new = 0 ; u1 = u2 = a1 ; }
	    if (a1 < u1) u1 = a1 ;
	    if (a2 < u1) u1 = a2 ;
	    if (a1 > u2) u2 = a1 ;
	    if (a2 > u2) u2 = a2 ;
	  }
    }

  if (doSave)
    /* saveSeqs (link, seqs) */ ;
  else
    {
      if (u1 != b1 || u2 != b2)
	{
	  printf ("\nSequence %s\n", aceName(link)) ;
	  printf ("Subsequence %s %d %d // was %d %d  shift: %d\n", 
		  aceName(junction), u1, u2, b1, b2, u1 - b1) ;
	 
	}
    }

    
 abort:
  acePopContext (TRUE) ; /* should save and messfree */
}

/***************************************************************/

static void checkOneOverlap (AceKey k1, int a1, int a2, AceKey k2, int b1, int b2)
{
  Array seqs ; 
  SEQ *hh ;
  acePushContext () ;

  seqs = aceArrayCreate (8, SEQ) ;

  hh = arrayp(seqs, 0, SEQ) ;
  hh->seq = k1 ;
  hh->a1 = a1 ;
  hh->a2 = a2 ;
  hh->isUp = FALSE ;

  hh = arrayp(seqs, 1, SEQ) ;
  hh->seq = k2 ;
  hh->a1 = b1 ;
  hh->a2 = b2 ;
  hh->isUp = FALSE ;

  swapA1 (seqs) ; 
  arraySort (seqs, seqOrder) ;
  checkSeqs(seqs) ;  
  shiftSeqs (seqs) ;
  backSwapA1 (seqs);
  printSeqs (0, seqs) ;
    
  arrayDestroy (seqs) ;
  acePopContext (TRUE) ; /* should save and messfree */
}

/***************************************************************/

static void fixAll (void)
{
  int n ;
  AceKeySet ks = 0 ;

  ks = aceQuery(0, "FIND Sequence  *link* AND Subsequence",0) ;

  if (!keySetMax(ks))
    { printf ("// No Sequence called *link*  in this database, sorry\n%s\n",
	      aceErrorMessage(0)) ;
      return ;
    }
	
  for (n = 0 ; n < keySetMax(ks) ; n++)
    fixLink (keySet(ks,n), FALSE) ;

  aceFree (ks) ;
  ks = aceQuery(0, "FIND Sequence  Junction",0) ;
  for (n = 0 ; n < keySetMax(ks) ; n++)
    fixJunction (keySet(ks,n)) ;
}

/***************************************************************/

static void doExportOverlaps (void)
{
  int n ;
  AceKeySet ks = 0 ;


  ks = aceQuery(0, "FIND Sequence  *link* AND Subsequence",0) ;

  if (!keySetMax(ks))
    { printf ("// No Sequence called *link*  in this database, sorry\n%s\n",
	      aceErrorMessage(0)) ;
      return ;
    }
	
  for (n = 0 ; n < keySetMax(ks) ; n++)
    {
      printf ("\n// Link %s\n\n",aceName(keySet(ks,n))) ;
      fixLink (keySet(ks,n), TRUE) ;
    }

  aceFree (ks) ;
}

/***************************************************************/

static void makeChromosomes (void)
{
  int n, xx, cumul, dnaLength, guessDnaLength, pmap1, pmap2, estimated_length ;
  int a1, a2, b1, b2 ;
  AceKeySet ks = 0 ;
  AceKey k0, k1, k2, so, map, pmap,clone, 
    _Overlap_right, _Estimated_length,
    _Source, _DNA, _Clone, _Map, _Pmap, _Flipped ;
  AceInstance Seq, Clone ;
  BOOL flipped0, flipped1 ;

  ks = aceQuery(0, "FIND Sequence (Genomic_canonical || genomic ) &&  overlap_right && !overlap_left",0) ;


  _Overlap_right = aceTag ("Overlap_right", FALSE) ;
  _Source = aceTag ("Source", FALSE) ;
  _DNA = aceTag ("DNA", FALSE) ;
  _Map = aceTag ("Map", FALSE) ;
  _Clone = aceTag ("Clone", FALSE) ;
  _Pmap = aceTag ("Pmap", FALSE) ;
  _Flipped = aceTag ("Flipped", FALSE) ;
  _Estimated_length = aceTag ("Estimated_length", FALSE) ;

  if (!keySetMax(ks))
    { printf ("// No Sequence with overlap_right and no overlap_left\n%s\n",
	      aceErrorMessage(0)) ;
      return ;
    }

  printf("// contig cosmid a1 a2 length source"
	 "ov_r_cosmid ov_r_value map pamp pmap1 pmap2 overlap_verif\n\n") ;

  a1 = a2 = 0 ;
  for (n = 0 ; n < keySetMax(ks) ; n++)
    {
      k1 = keySet(ks, n) ; cumul = 0 ; k0 = 0 ; flipped0 = flipped1 = FALSE ;
      while (k1)
	{
	  Seq = aceOpenKey (k1,0) ;
	  k2 = so = 0 ; xx = 0 ; dnaLength = 0 ; map = 0 ; pmap = 0 ; clone = 0 ;
	  pmap1 = pmap2 = 0 ; flipped1 = FALSE ;
	  if (Seq)
	    {
	      if (!aceGotoTag (Seq, _Source) || !aceGotoChild (Seq) ||
		  !aceGetKey (Seq, &so))
		so = 0 ;
	      if (!aceGotoTag (Seq, _Clone) || !aceGotoChild (Seq) ||
		  !aceGetKey (Seq, &clone))
		clone = 0 ;
	      if (!aceGotoTag (Seq, _Map) || !aceGotoChild (Seq) ||
		  !aceGetKey (Seq, &map))
		map = 0 ;
	      if (aceGotoTag (Seq, _Flipped))
		flipped1 = TRUE ;
	      if (!aceGotoTag (Seq, _DNA) || !aceGotoChild (Seq) ||!aceGotoChild (Seq) ||
		  !aceGetInteger (Seq, &dnaLength))
		dnaLength = 0 ;
	      if (!aceGotoTag (Seq, _Estimated_length) || !aceGotoChild (Seq) ||
		  !aceGetInteger (Seq, &estimated_length))
		estimated_length = 0 ;
	      if (!aceGotoTag (Seq, _Overlap_right) || !aceGotoChild (Seq) ||
		  !aceGetKey (Seq, &k2))
		k2 = 0 ;
	      else  if (!aceGotoChild (Seq) ||
			!aceGetInteger (Seq, &xx))
		xx = 0 ;
	      aceFree (Seq) ;
	    }
	  if (clone && (Clone = aceOpenKey (clone, 0)))
	    {
	      if (!aceGotoTag (Clone, _Pmap) || !aceGotoChild (Clone) ||
		  !aceGetKey (Clone, &pmap))
		pmap = 0 ;
	      else if (!aceGotoChild (Clone) || !aceGetInteger (Clone, &pmap1) ||
		       !aceGotoChild (Clone) || !aceGetInteger (Clone, &pmap2))
		pmap1 = pmap2 = 0 ;
	      aceFree (Clone) ;
	    }
	  if (dnaLength)
	    guessDnaLength = dnaLength ;
	  else if (estimated_length)
	    guessDnaLength = estimated_length ;
	  else if(xx)
	    guessDnaLength = xx - 1 ;
	  else
	    guessDnaLength = 5000 ;
	  if (flipped1) 
	    { b2 = cumul + 1 ; b1 = cumul+guessDnaLength ; }
	  else
	    { b1 = cumul + 1 ; b2 = cumul+guessDnaLength ; }
	      
	  printf ("%d\t%s\t%d\t%d\t%d\t",n+1, aceName(k1), b1, b2, dnaLength) ;
	  if (so) 
	    printf("%s\t",aceName(so));
	  else
	    printf("No_source\t") ;
	  if (k2) 
	    printf("%s\t",aceName(k2));
	  else
	    printf("End\t") ;
	  printf("%d\t",xx) ;
	  if (xx)
	    cumul += xx - 1 ;
	  else if (dnaLength)
	    cumul += dnaLength ;
	  else
	    cumul += 3000 ;
	  if (map) 
	    printf("%s\t",aceName(map));
	  else
	    printf("No_map\t") ;
	  if (pmap && (pmap1 || pmap2)) 
	    printf("%s\t%d\t%d\t",aceName(pmap),pmap1,pmap2);
	  else
	    printf("No_pmap\t0\t0\t") ;
	  if (k0)
	    checkOneOverlap (k0, a1, a2, k1, b1, b2) ;
	  a1 = b1 ; a2 = b2 ; k0 = k1 ; flipped0 = flipped1 ;
	  printf("\n") ;
	  k1 = k2 ;
	}
      printf("\n\n") ;
    }

  aceFree (ks) ;
}

/***************************************************************/

static void die (void)
{ fprintf (stderr, "// %s\n", aceErrorMessage(0)) ; exit (1) ; }

static void usage (void)
{ fprintf (stderr, "Usage: fixsuperlink [options] $ACEDB\n"
	   "  version 3, 19 apil 1999\n"
           "  fixes the *link* subsequence coordinates\n"
	   "\t-m <min overlap>: min overlap when looking for new match (default 1)\n"
	   "\t-a:  export fixes as an ace file to stdio (default)\n"
	   "\t-w:  write directly in the database\n"
	   "\t-G:  fix St Louis gaps to 200bp\n"
           "\t-O:  create overlap.ace from link data\n"
	   "\t-stl: inside -O use all subsequence not just the Genomic_canonical\n"
	   "\t-makeChrom: verify overlap_right and create a table of cosmids, source, pmap\n") ;
}

int main (int argc, char **argv)
{ 
  /*
extern int utMainPart (float x) ;
extern int utMainRoundPart (float x) ;
  for ( x = -21; x < 21  ; x++)
    printf ("%g %4d %4d\n", x, utMainPart (x), utMainRoundPart(x)) ;
  exit (0) ;
  */

  for (--argc, ++argv ; argc > 1 ; --argc, ++argv)
    if (!strcmp (*argv, "-w"))
      doSave = TRUE ;
    else if (!strcmp (*argv, "-a"))
      doSave = FALSE ;
    else if (!strcmp (*argv, "-G"))
      fiddleGaps = TRUE ;
    else if (!strcmp (*argv, "-O"))
      exportOverlaps = TRUE ;
    else if (!strcmp (*argv, "-stl"))
      isSTL = TRUE ;
    else if (!strcmp (*argv, "-makeChrom")) 
      isMakeChrom = TRUE ;
    else if (!strcmp (*argv, "-m") && argc > 2)
      { --argc ; ++argv ; 
        minOverlap = atoi(*argv) ; 
        --minOverlap ; if (minOverlap < 0) minOverlap = 0 ;
      }
    else
      { fprintf (stderr, "Unrecognized option %s\n", *argv) ;
        usage() ;
	exit (-1) ;
      }
  
  if (!argc)
      { fprintf (stderr, "Database location MUST be specified\n") ;
        usage() ;
	exit (-1) ;
      }

  if (!aceInit(*argv))
    die () ;
  aceOpenContext("fixsuperlink") ;

  if (doSave && !(aceGetWriteAccess(TRUE)))
    die () ;

  if (doSave)
    printf("// Directly editing the database\n") ;
  else
    printf("// This acefile was produced by fixsuperlink version 3 (Avril 99) \n") ;
  
  printf ("// database: $ACEDB: %s\n", *argv) ;
  printf ("// date:\t %s\n", timeShowNow()) ;
  printf ("// please send editions and comments to mieg@kaa.crbm.cnrs-mop.fr\n") ;

  if (exportOverlaps)
    doExportOverlaps() ;
  else if (isMakeChrom)
    makeChromosomes() ;
  else
    fixAll() ;

  /* aceStatus() ; */

  aceQuit(TRUE) ;   /* savesession */
  printf("\n// done (if the file is empty, it means that all test succeeded)\n") ;
  exit (0) ;
}

/***************************************************************/
/***************************************************************/
