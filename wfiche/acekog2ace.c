#include "regular.h"
#include "aceio.h"
#include "mytime.h"

#define PREFIX(b,s) ( strncmp(b,s,sizeof(s)-1) == 0 )

static char *getPrefix (ACEIN fi, char *prefix)
{
  char *cp ;
  
  while (aceInCard (fi))
    {
      cp = aceInPos (fi) ;
      if (!cp || strncmp (cp, prefix, strlen(prefix)))
	continue ;
      return cp + strlen(prefix) ; /* the rest of the line */
    }
  return 0 ;
}

static int aceKog2aceParse (ACEIN fi, BOOL isN)
{
  int nKantor = 0, subjectTitle, gene, link, product ;
  int a1, a2, x1, x2,  b1, y1, eValue, nHit ;
  float bitScore = 0, bestScore = 0 ;
  char *cp, *cq, species ;
  Stack s = stackCreate (1000) ;
  int spec[256], tag[256] ;
  double expect ;

  pushText (s, "-") ;
  for (a1 = 0 ; a1 < 256 ; a1++)
    spec[a1] = 0 ;
  spec['w'] = stackMark (s) ;
  pushText (s, "Caenorhabditis elegans") ;

  spec['a'] = stackMark (s) ;
  pushText (s, "Arabidopsis thaliana") ;

  spec['h'] = stackMark (s) ;
  pushText (s, "Homo sapiens") ;

  spec['m'] = stackMark (s) ;
  pushText (s, "Mus musculus") ;

  spec['r'] = stackMark (s) ;
  pushText (s, "Rattus norvegicus") ;

  spec['d'] = stackMark (s) ;
  pushText (s, "Drosophila melanogaster") ;

  tag['w'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNWorm" : "AceKogWorm") ;

  tag['a'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNAra" : "AceKogAra") ;

  tag['h'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNHuman" : "AceKogHuman") ;

  tag['m'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNMouse" : "AceKogMouse") ;

  tag['r'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNRat" : "AceKogRat") ;

  tag['d'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNDroso" : "AceKogDroso") ;
 lao: 

  bestScore = 0 ;
  cp = getPrefix (fi, "<b>Query=</b>") ;
  if (!cp)
    {
      stackDestroy (s) ;
      return nKantor ;
    }

  if (!cp || !*cp++ || !*cp)
    goto lao ;
  nKantor++ ;
  if (isN)
    {
      printf ("\nmRNA \"%s\"\n", cp+5) ;
      printf("AceKogN_Date %s\n", timeShowNow ());
      if (0) printf ("-D AKGN\n") ;/* no because we merge several species */
    }
  else
    {
      printf ("\nKantor \"%s\"\n", cp) ;
      printf("AceKog_Date %s\n", timeShowNow ());
      printf("Kantor_Date %s\n", timeShowNow ());
      if (0) printf ("-D AKG\n") ;/* no because we merge several species */
    }
  nHit = 0 ;

 nextHit:
  cp = getPrefix (fi, "><a name =") ;
  if (!cp) goto lao ;

 /* one letter species code */
  cp = strstr (cp, "</a>") ;
  if (!cp) goto lao ;
  cp += 4 ;
  species = *cp++ ; 

 /* first word of title is gene name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  gene = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

 /* next word of title is locuslink/newname name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  link = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

  /* next word of title is product name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  product = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

  /* rest if full title on several lines */
  subjectTitle = stackMark (s) ; 
  pushText (s, aceInUnprotect (fi, cp)) ;
  if (aceInCard (fi) &&
      (cp = aceInPos (fi)) &&
      (cp = strstr (cp, "SET\\/")))
    {
      cp += 4 ; *cp = ' ' ;
      cq = strstr (cp, "query") ;
      if (cq ) *cq = 0 ;
      catText (s, cp) ;
    }

  /* Score =  214 bits (545), Expect = 1e-55 */
  cp = getPrefix (fi, "Score") ;
  if (! aceInWord (fi)  || /* Score */
      ! aceInWord (fi)  || /* = */
      ! aceInFloat (fi, &bitScore)  ||  /* 214 */
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! (cp = aceInWord (fi)))
    {
      fprintf (stderr, "// bad score line %d: %s\n", aceInStreamLine (fi), cp) ;
      goto lao ;
    }
  eValue = stackMark (s) ;
  pushText (s, cp) ; /* very often out of range of a float */
  expect = 1 ;
  if (strstr (cp, "e-")) 
    expect = 0 ;
  else
    sscanf (cp, "%lf", &expect) ;
  aceInCard (fi) ;  aceInCard (fi) ;  aceInCard (fi) ;
  b1 = y1 = -1 ;
  while (1) 
    {
      /* Query: 135 IESSRDLLHRIKDEVGAPGIVVGVSVDGKEVWSEGLGYADVENRVPCKPETVMRIASISK 194 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp)
	continue ;
      if (strcmp (cp, "Query:"))
	break ;
      if (!aceInInt (fi, &a1) || !aceInWord (fi) || !aceInInt (fi, &a2))
	break ;
      if (b1 == -1) b1 = a1 ;
      /* I+ +++L+       G PG+ + VS+DGK VW  G GYA++E+   C  ++VMRIASISK */
      aceInCard (fi) ;  /* jump that line */
      
      /* Sbjct: 115 IKKAKELVETTMAIQGIPGLSIAVSLDGKMVWKSGFGYANLESFAKCTGDSVMRIASISK 294 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp || strcmp (cp, "Sbjct:"))
	break ;
      if (!aceInInt (fi, &x1) || !aceInWord (fi) || !aceInInt (fi, &x2))
	break ;
      if (y1 == -1) y1 = x1 ; 
      if (!aceInCard (fi)) break ;
    }

  if (y1 != -1 && 100 * bitScore > 88 * bestScore
      && (!isN || (isN && expect < .001))      
      )      /* success export one hit */
    {
      printf ("AceKog%s %s \"%s\" %g %d %d %d %d %s\n"
	      , isN ? "N" : ""
	      , stackText (s, product)
	      , stackText (s, isN ? tag[(int)species] : spec[(int)species])
	      , bitScore, b1, a2, y1, x2
	      , aceInProtect (fi, messprintf ("%s [%s] eValue = %s",
					      stackText (s, subjectTitle) 
					      , stackText (s, spec[(int)species])
					      , stackText (s, eValue)
					      )
			      )
	      ) ;
      if (1 || !nHit) printf ("%s %s %s\n"
			 , stackText (s, tag[(int)species])
			 , stackText (s,gene)
			 , stackText (s,link)
			 ) ;
      if (! nHit || bitScore > bestScore)
	bestScore = bitScore ;
      nHit++ ;
    }
  if (nHit >= 36)
    goto lao ;
  goto nextHit ;
} /* aceKog2aceParse */

/******************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage acekog2ace < blast.out >! blast.ace \n") ;
  exit (1) ;
} /* usage */

/******************************************************************/
/* split the stdin file on Query= */
static int aceKog2aceSplit (ACEIN fi, BOOL isN)
{
  Stack s = stackCreate (10000) ;
  ACEIN splitFi = 0 ;
  int nn = 0, line = 0 ;
  char *cp, *prefix = "<b>Query=" ;

  while (TRUE) 
    {
      if (aceInCard (fi))
	cp = aceInPos (fi) ;
      else
	cp = 0 ;
      line++ ;
      if (!cp || !strncmp (cp, prefix, strlen(prefix)))
	{
	  if (stackMark (s))
	    {
	      splitFi = aceInCreateFromText (stackText (s, 0), 0, 0) ;
	      nn += aceKog2aceParse (splitFi, isN) ;
	      aceInDestroy (splitFi) ;
	    }
	  stackClear (s) ;
	}
      if (cp)
	{ catText (s, cp) ; catText (s, "\n") ; }
      else
	break ;
    }

  stackDestroy (s) ;
  return nn ;
}

/******************************************************************/

int main (int argc, char *argv[])
{
  if (argc == 1 || (argc == 2 && !strcmp(argv[1],"-n")))
    {
      ACEIN fi = aceInCreateFromStdin (0, 0, 0) ;
      int nn = aceKog2aceSplit (fi,argc==2 ? TRUE : FALSE) ;
      fprintf (stderr, "// Parse %d blast outputs\n", nn) ;
      aceInDestroy (fi) ;
      printf ("\n") ;
    }
  else
    usage () ;
    
  return 0 ;
}

/******************************************************************/
/******************************************************************/
