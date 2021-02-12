/* histo.c
 * Authors: Jean Thierry-Mieg, NCBI, 
 * 1972, C. R. Acad Sci Hebd Seances Acad Aci D 1972 Dec 6; 275(23):2751-4.
 * rewritten in april 2010
 * Histograms and smoothed histograms
*/

#include <ac.h>
#include "bigarray.h"
#include "bitset.h"
#include <math.h>

#define NNh 100
 
/* #define ARRAY_CHECK */

typedef struct sxStruct {  
  AC_HANDLE h ; 
  BOOL gzi, gzo, snp, smooth, plain, plot ;
  int KL ; 
  float min, max ;
  BOOL hasMin, hasMax ;
  const char *inFileName, *outFileName, *title, *columns, *runFileName ;
  KEYSET cols ;
  DICT *runDict ;
  Associator run2title ;
  float mutMedian, mutAverage, mutSigma ;
  float wildMedian, wildAverage, wildSigma ;
  float coverMedian, coverAverage, coverSigma ;

  ACEIN ai ; ACEOUT ao ; int n, NN ; double aa[NNh+1] ; } SX ;

/*************************************************************************************/
/*************************************************************************************/
/* register the samll values up to LL by multiplicity
 * register the larger values as a list
 */
static void histoSmoothRegister (BigArray a, int n, int LL)
{
  if (n >= 0 && n < LL)
    bigArray (a, n, float)++ ;
  else
    {
      long int m = bigArrayMax (a) ;
      bigArray (a, m, float) = n ;
    }
  return ;
} /* histoSmoothRegister */

/*************************************************************************************/

static int histoSmoothGetData (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = sx->ai ;
  int i, k, n, NNz = 100 ;
  int n1 = 0, n2 = 0, n3 = 0 ;
  BigArray bMutant = bigArrayHandleCreate (100000, float, h) ;
  BigArray bWild = bigArrayHandleCreate (100000, float, h) ;
  BigArray bCover = bigArrayHandleCreate (100000, float, h) ;
  Array hh = arrayHandleCreate (101, double, h) ;
  int NN = 200 ;
  int LL = 10000 ;
  Array dd = arrayHandleCreate (NN * NN, int, h) ;

  array (dd, NN * NN - 1, int) = 0 ; /* make space */
  while (aceInCard (ai))
    {
      if (! aceInInt (ai, &n))
	continue ;
      aceInNext (ai) ;
      if (aceInInt (ai, &k))
	{
	  
	  histoSmoothRegister (bMutant, n, LL) ;
	  histoSmoothRegister (bWild, k - n, LL) ;
	  histoSmoothRegister (bCover, k, LL) ;

	  if (k < NN && n < NN)
	    {
	      n2++ ;
	      arr (dd, k * NN + n, int) += 1 ;
	    }
	  else
	    {
	      n1++ ;
	      utSmoothHisto (hh, k, n, 1) ;
	    }
	}
    } 
  fprintf (stderr, "parsed %d lines computed %d values :: %s\n", sx->n, n1, timeShowNow ()) ;
  for (k = 0 ; k < NN ; k++)
    for (n = 1 ; n < NN ; n++)
      {
	i = arr (dd, k * NN + n, int) ;
	if (i > 0)
	  {
	    n3++ ;
	    utSmoothHisto (hh, k, n, i) ;
	  }
      }
  fprintf (stderr, "computed %d/%d additional values :: %s\n", n2, n3, timeShowNow ()) ;
   
  for (i = 0 ; i <= NNz ; i++)
    sx->aa[i] = arr (hh, i, double) ;

  bigFloatVariance (bMutant, LL, &sx->mutMedian, &sx->mutAverage, &sx->mutSigma) ;
  bigFloatVariance (bWild, LL, &sx->wildMedian, &sx->wildAverage, &sx->wildSigma) ;
  bigFloatVariance (bCover, LL, &sx->coverMedian, &sx->coverAverage, &sx->coverSigma) ;

  ac_free (h) ;
  return sx->n ;
} /* histoSmoothGetData */

/*************************************************************************************/
typedef struct hhStruct {Array aa, bb ; BitSet bs ; int n ; const char *title ; } HH ;
static int histoSnpGetData (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ; 
  Array hhh = arrayHandleCreate (8, HH, sx->h) ;
  HH *hh ;
  ACEIN ai = sx->ai ;
  DICT *dict = sx->runDict ;
  float cover, mutant, wild ;
  const char *ccp ;
  int i, run ;
  int LL = 10000 ;
  BigArray bMutant = bigArrayHandleCreate (100000, float, h) ;
  BigArray bWild = bigArrayHandleCreate (100000, float, h) ;
  BigArray bCover = bigArrayHandleCreate (100000, float, h) ;

  if (! dict)
    dict = dictHandleCreate (100, h) ;
  aceInSpecial (ai, "\n") ;

  if (0)   utSmoothHistoTest () ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ; if (! ccp || *ccp == '#') continue ;
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ;
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ;
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ;
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ;
      aceInStep (ai, '\t') ; if (! (ccp = aceInWord (ai))) continue ; /* run name */
      dictAdd (dict, ccp, &run) ;
      hh = arrayp (hhh, run, HH) ;
      if (! hh->title)
	{
	  hh->title = hprintf(sx->h, "%s", ccp) ;
	  hh->aa = arrayHandleCreate (101, double, sx->h) ;
	  hh->bs = bitSetCreate (80000, sx->h) ;
	}
  
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ; /* genotype */
      aceInStep (ai, '\t') ; if (! aceInWord (ai)) continue ;  /* frequency */
      aceInStep (ai, '\t') ; if (! aceInFloat (ai, &cover)) continue ;
      aceInStep (ai, '\t') ; if (! aceInFloat (ai, &mutant)) continue ;
      aceInStep (ai, '\t') ; if (! aceInFloat (ai, &wild)) continue ;

      histoSmoothRegister (bMutant,  mutant, LL) ;
      histoSmoothRegister (bWild, wild, LL) ;
      histoSmoothRegister (bCover, cover, LL) ;

      utSmoothHisto (hh->aa, mutant, cover, 1) ;
    }  
  
  bigFloatVariance (bMutant, LL, &sx->mutMedian, &sx->mutAverage, &sx->mutSigma) ;
  bigFloatVariance (bWild, LL, &sx->wildMedian, &sx->wildAverage, &sx->wildSigma) ;
  bigFloatVariance (bCover, LL, &sx->coverMedian, &sx->coverAverage, &sx->coverSigma) ;

  if (1)
    {
      AC_HANDLE h0 = ac_new_handle () ;
      int hMax = dictMax (dict) + 1 ;
      int jj ;
      ACEOUT ao = 0 ;

      if (0 && sx->outFileName)
	system (messprintf("\\rm %s %s.tcsh %s.ps", sx->outFileName, sx->outFileName, sx->outFileName)) ;
      ao = aceOutCreate (sx->outFileName, ".txt", FALSE, h0) ;
      if (ao)
	{
	  aceOutDate (ao, "##", sx->title) ;
	  aceOut (ao, "Value") ;
	  for (jj = 1 ; jj < hMax ; jj++)
	    {
	      hh = arrayp (hhh, jj, HH) ;
	      if (!hh->title)
		continue ;
	      aceOutf (ao, "\t%s", hh->title) ;
	    } 
	  aceOut (ao, "\n") ;
	  if (sx->run2title)
	    {
	      aceOut (ao, "Title") ;
	      for (jj = 1 ; jj < hMax ; jj++)
		{
		  char *cp ;
		  hh = arrayp (hhh, jj, HH) ;
		  if (!hh->title)
		    continue ;
		  cp = 0 ;
		  assFind (sx->run2title, assVoid(jj), &cp) ;
		  aceOutf (ao, "\t%s", cp ? cp : "-") ;
		} 
	    }
	  aceOut (ao, "\n") ;
	  for (i = 0 ; i <= 100 ; i++)
	    {
	      aceOutf (ao, "%d", i) ;
	      for (jj = 1 ; jj < hMax ; jj++)
		{
		  hh = arrayp (hhh, jj, HH) ;
		  if (!hh->title)
		    continue ;
		  if (i < sx->min || (sx->max && i > sx->max))
		    aceOutf (ao, "\t0") ;
		  else
		    aceOutf (ao, "\t%f", array(hh->aa, i, double)) ;
		}
	      aceOut (ao, "\n") ;
	    }
	}
      if (sx->plot)
	{
	  ao = aceOutCreateToFile (messprintf("%s.tcsh", sx->outFileName), "wx", h0) ;
	  aceOutf(ao, "#!/bin/tcsh -f\n") ;
	  aceOutf(ao, "  gnuplot -bg white << EOF\n") ;
	  aceOutf(ao, "    set style data linespoints\n") ;
	  aceOutf(ao, "    set terminal postscript color \n") ;
	  aceOutf(ao, "    set terminal postscript solid \n") ;
	  aceOutf(ao, "    set terminal postscript landscape\n") ;
	  aceOutf(ao, "    set output '%s.ps'\n", sx->outFileName) ;
	  if (sx->title) 
	    aceOutf (ao, "set title '%s'\n", sx->title) ;

	  for (i = 0, jj = 1 ; jj < hMax ; jj++)
	    {
	      hh = arrayp (hhh, jj, HH) ;
	      if (!hh->title)
		continue ;
	      if (arrayMax (hh->aa))
		aceOutf(ao, "%s '%s.txt' using %d with line title '%s'"
			, i++  ? ",\\\n\t" : "plot "
			, sx->outFileName
			, jj + 1
			, hh->title 
			) ;
	    }
	  aceOutf(ao, "\nEOF\n") ;
	  aceOutDestroy (ao) ; /* will close */
	  system (messprintf("(sleep 1 ; tcsh %s.tcsh ; gv %s.ps ) &", sx->outFileName, sx->outFileName)) ;
	}
      ac_free (h0) ;
    }
  

  ac_free (h) ;
  return sx->n ;
} /* histoSnpGetData */

/*************************************************************************************/

static void histoSmoothExport (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  double *aa = sx->aa ;
  ACEOUT ao = sx->ao ;
  int i, NNz = NNh ;
  if (sx->NN < NNh) NNz = sx->NN ;

  ao = sx->ao = aceOutCreate (sx->outFileName, ".txt", sx->gzo, sx->h) ;
  aceOutf (ao, "# %s\n", timeShowNow()) ;
  if (sx->title) aceOutf (ao, "# %s\n#", sx->title) ;
  aceOutf (ao, "# Analyse %d points\tMedian\tAverage\tSigma\n", sx->n) ;
  aceOutf (ao, "# Type 1\t%f\t%f\t%f\n", sx->mutMedian, sx->mutAverage, sx->mutSigma) ;
  aceOutf (ao, "# Type 2\t%f\t%f\t%f\n", sx->wildMedian, sx->wildAverage, sx->wildSigma) ;
  aceOutf (ao, "# Total\t%f\t%f\t%f\n", sx->coverMedian, sx->coverAverage, sx->coverSigma) ;

  aceOutf (ao, "\n\n\n#\t%s\n", sx->title ? sx->title : "Smoothed histogram") ;
  for (i = 0 ; i <= NNz ; i++)
    aceOutf (ao, "%.2f\t%.5f\n", 100.0 * i/(double)NNz, aa[i]) ;
  
  if (sx->plot)
    {
      system (messprintf("\\rm %s.tcsh %s.ps", sx->outFileName, sx->outFileName)) ;
      ao = aceOutCreateToFile (messprintf("%s.txt", sx->outFileName), "w", h) ;
      if (ao)
	{
	  aceOutf (ao, "// %s\n", timeShowNow()) ;
	  for (i = 0 ; i <= NNz ; i++)
	    aceOutf (ao, "%.2f\t%.3f\n", 1.0 * i/(double)NNz, aa[i]) ;

	  ao = aceOutCreateToFile (messprintf("%s.tcsh", sx->outFileName), "wx", h) ;
	  aceOutf(ao, "#!/bin/tcsh -f\n") ;
	  aceOutf(ao, "  gnuplot -bg white << EOF\n") ;
	  aceOutf(ao, "    set style data linespoints\n") ;
	  aceOutf(ao, "    set terminal postscript color \n") ;
	  aceOutf(ao, "    set terminal postscript solid \n") ;
	  aceOutf(ao, "    set terminal postscript landscape\n") ;
	  aceOutf(ao, "    set output '%s.ps'\n", sx->outFileName) ;
	  if (sx->title) 
	    aceOutf (ao, "set title '%s'\n", sx->title) ;

	  aceOutf(ao, "    plot [0:100] ") ;
	  /* order is important since i do not control the colors */
	  if (0)
	    aceOutf(ao, "  '%s.txt'  using (log(\\$2)) with line title 'distrib'"
		    , sx->outFileName, sx->title ? sx->title : "smoothed SX hito") ;
	  else
	    aceOutf(ao, "  '%s.txt'  using 2 with line title '%s'"
		    , sx->outFileName, sx->title ? sx->title : "smoothed SX hito") ;

	  aceOutf(ao, "\nEOF") ;
	  aceOutDestroy (ao) ; /* will close */

	  system (messprintf("(sleep 1 ; tcsh %s.tcsh ; gv %s.ps ) &", sx->outFileName, sx->outFileName)) ;
	}
    }
  ac_free (h) ;
} /* histoSmoothExport */

/*************************************************************************************/
/*************************************************************************************/
/* plot a histogram of integers given on stdin */

static int histoPlainGetData (SX *sx)
{
  ACEIN ai = sx->ai ;
  ACEOUT ao = sx->ao ;
  KEYSET cols = sx->cols ;
  int h, i, n, col, nn = 0, line = 0 ;
  float f, fmin = 0, fmax = 0, df, ddf ; 
  const char *ccp ;
  int hMax = keySetMax (cols), hMax2 = 0 ;
  Array hhh = arrayHandleCreate (8, HH, sx->h) ;
  Array total = 0, correl = 0, correl2 = 0 ;
  HH *hh = 0 ;
  BOOL isTitle ;
  
  for (h = 0 ; h < hMax ; h++)
    {
      hh = arrayp (hhh, h, HH) ;
      hh->title = hprintf(sx->h, "c%d", keySet (cols, h)) ;
      hh->bb = arrayHandleCreate (sx->NN+1, int, sx->h) ;
      array (hh->bb, sx->NN, int) = 0 ;
      hh->aa = arrayHandleCreate (1000, float, sx->h) ;
      hh->bs = bitSetCreate (80000, sx->h) ;
      if (h == 0 && sx->title)
	hh->title = hprintf(sx->h, "%s", sx->title) ; 
    }
  
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char cutter = 0 ;
      col = 1 ; line++ ;
      if (line == 1 && *aceInPos(ai) == '#')
	isTitle = TRUE ;
      else
	isTitle = FALSE ;
      for (h = 0 ; h < hMax ; h++)
	{
	  hh = arrayp (hhh, h, HH) ;
	  /* locate on the corect column */
	  while (col < keySet (cols, h))
	    {
	      if (col > 1 && cutter != '\t')
		h = hMax ;
	      ccp = aceInWordCut (ai, "\t", &cutter) ;
	      col++ ;
	    }
	  if (h == hMax)
	    break ;

	  if (col > 1 && cutter != '\t')
	    break ;
	  col++ ;
	  if (isTitle)
	    {
	      ccp = aceInWordCut (ai, "\t", &cutter) ;
	      if (ccp && *ccp)
		hh->title = strnew (ccp, sx->h) ;
	      continue ;
	    }
	  if (!aceInFloat (ai, &f))
	    {
	      cutter = 0 ;
	      if (aceInStep (ai, '\t'))
		cutter = '\t' ;
	      continue ; /* missing value */
	    }
	  cutter = 0 ;
	  if (aceInStep (ai, '\t'))
	    cutter = '\t' ;

	  bitSet (hh->bs, line) ;
	  array (hh->aa, line, float) = f ;

	  if (hMax2 <= h) hMax2 = h + 1 ;
	  if (f < sx->min) continue ;
	  if (f > sx->max) continue ;
	  if (nn++ == 0) fmin = fmax = f ;
	  if (f < fmin) fmin = f ;
	  if (f > fmax) fmax = f ;
	}
    }
  if (!nn)
    messcrash ("There was no int number on stdin, sorry") ;
  fprintf (stderr, "// found %d data points\n", nn);

  hMax = hMax2 ;
  if (hMax2)
    {
      total = arrayHandleCreate (hMax, double, sx->h) ;
      correl = arrayHandleCreate (hMax * hMax, double, sx->h) ;
    }
  
  if (sx->NN && sx->hasMin && sx->min < fmin) fmin = sx->min ;
  if (sx->NN && sx->hasMax && sx->max > fmax) fmax = sx->max ;
  
  n = fmin >= 0 ? 1 : -1 ; fmin *= 1000 * n ;
  if (n < 0) 
    {
      for (i = 0, df = utMainPart (fmin), ddf = df/10 ; fmin > df + i * ddf ; i++) ;
      fmin = df + i * ddf ;
    }
  else 
    {
      for (i = 0, df = utMainPart (2*fmin), ddf = df/10 ; fmin < df - i * ddf ; i++) ;
      fmin = df - i * ddf ;
    }
  fmin *= n/1000.0 ;
  
  n = fmax >= 0 ? 1 : -1 ; fmax *= 1000 * n ;
  if (n > 0) 
    {
      for (i = 0, df = utMainPart (fmax), ddf = df/10 ; fmax > df + i * ddf ; i++) ;
      fmax = df + i * ddf ;
    }
  else 
    {
      for (i = 0, df = utMainPart (2*fmax), ddf = df/10 ; fmax < df - i * ddf ; i++) ;
      fmax = df - i * ddf ;
    }
  fmax *= n/1000.0 ;

  if (fmin > 0 && fmin < fmax/4) fmin = 0 ;
  df = fmax - fmin ;
  ddf = 1000 * df/sx->NN ;
  for (df = 0 ; ddf > utMainPart (ddf + df) ; df += ddf/10) ;
  ddf = utMainPart (ddf + df)/1000.0 ;
  if (ddf == 0) ddf = 1 ;

  /* compute the histo and the total */
  for (h = 0 ; h < hMax ; h++)
    {
      hh = arrayp (hhh, h, HH) ;
      if (arrayExists(hh->aa))
	{
	  for (i = 0 ; i < arrayMax (hh->aa) ; i++)
	    if (bit(hh->bs, i))
	      {
		f = arr (hh->aa, i, float) ;
		array (total, h, double) += f ;
		if (f < sx->min) continue ;
		if (f > sx->max) continue ;
		n = ((f + ddf/2) - fmin)/ddf ;
		if (n >= 0 && n <= sx->NN)
		  array (hh->bb, n, int)++ ;
	      }
	}
    }
  if (1) /* correl calculation and export */
    {
      /* compute the average */
      for (h = 0 ; h < hMax ; h++)
	{
	  hh = arrayp (hhh, h, HH) ;
	  n = bitSetCount (hh->bs) ;
	  if (n > 0)
	    array (total, h, double) /= n ;
	}
      /* compute the correl */
      {
	AC_HANDLE handle2 = ac_new_handle () ;
	int h1, h2 ;
	double z ;
	HH *hh1, *hh2 ;
	ACEOUT ao = aceOutCreate (sx->outFileName, ".correl.txt", sx->gzo, handle2) ;
	
	for (h1 = 0 ; h1 < hMax ; h1++)
	  for (h2 = h1 ; h2 < hMax ; h2++)
	    {
	      hh1 = arrayp (hhh, h1, HH) ;
	      hh2 = arrayp (hhh, h2, HH) ;
	      if (arrayExists(hh1->aa) && arrayExists(hh2->aa) )
		for (i = 0 ; i < arrayMax (hh1->aa) && i < arrayMax (hh2->aa) ; i++)
		  if (bit(hh1->bs, i) && bit(hh2->bs, i))
		    {
		      array (correl, h1 * hMax + h2, double) += 
			(array (hh1->aa, i, float) - array (total, h1, double)) *
			(array (hh2->aa, i, float) - array (total, h2, double)) 
			;
		    }
	    }
	/* normalize the correl */
	correl2 = arrayHandleCopy (correl, sx->h) ;
	for (h1 = 0 ; h1 < hMax ; h1++)
	  for (h2 = h1 ; h2 < hMax ; h2++)
	    {
              z = array (correl2, h1 * hMax + h1, double) * array (correl2, h2 * hMax + h2, double) ;
	      if (z == 0) z = 1 ;
	      z = sqrt (z) ;
	      z =  array (correl2, h1 * hMax + h2, double) / z ;
	      array (correl, h1 * hMax + h2, double) = z ;
	    }
	/* print out */
	if (sx->title)
	  aceOutf (ao, "%s\n%s\n", timeShowNow(), sx->title) ;
	else
	  aceOutf (ao, "// %s\n", timeShowNow()) ;
	for (h2 = 0 ; h2 < hMax ; h2++)
	  {
	    hh = arrayp (hhh, h2, HH) ;
	    if (arrayMax (hh->aa))
	      aceOutf (ao, "\t%s", hh->title) ;
	  } 
	z = 0 ; n = 0 ;
	for (h1 = 0 ; h1 < hMax ; h1++)
	  {
	    hh = arrayp (hhh, h1, HH) ;
	    if (arrayMax (hh->aa))
	      {
		aceOutf (ao, "\n%s", hh->title) ;
		for (h2 = 0 ; h2 < hMax ; h2++)
		  {
		    hh = arrayp (hhh, h2, HH) ;
		    if (arrayMax (hh->aa))
		      {
			n++ ; z += array (correl, h1 < h2 ? h1 * hMax + h2 : h2 * hMax + h1, double) ;
			aceOutf (ao, "\t%.2g", array (correl, h1 < h2 ? h1 * hMax + h2 : h2 * hMax + h1, double)) ;
		      }
		  }
	      }
	  }
	aceOutf (ao, "\nAverage correlation %.2g\n", n ? z/n : 0) ;
	ac_free (handle2) ;
      }
    }

  fprintf (stderr, "Found %d values in %d lines %d columns fmin=%.1g fmax=%.1g ddf=%.g\n", nn, line, hMax, fmin, fmax,ddf) ;

  if (1)
    {
      AC_HANDLE h0 = ac_new_handle () ;

      if (sx->outFileName)
	system (messprintf("\\rm %s %s.tcsh %s.ps", sx->outFileName, sx->outFileName, sx->outFileName)) ;
      ao = aceOutCreate (sx->outFileName, ".txt", FALSE, h0) ;
      if (ao)
	{
	  if (sx->title)
	    aceOutf (ao, "%s\n%s\n", timeShowNow(), sx->title) ;
	  aceOut (ao, "Value") ;
	  for (h = 0 ; h < hMax ; h++)
	    {
	      hh = arrayp (hhh, h, HH) ;
	      if (arrayMax (hh->aa))
		aceOutf (ao, "\t%s", hh->title) ;
	    } 
	  aceOut (ao, "\n") ;
	  for (i = 0 ;  fmin + i * ddf <= fmax ; i++)
	    {
	      aceOutf (ao, "%g", fmin + i * ddf) ;
	      for (h = 0 ; h < hMax ; h++)
		{
		  hh = arrayp (hhh, h, HH) ;
		  if (arrayMax (hh->aa) && i < arrayMax (hh->bb))
		    aceOutf (ao, "\t%d", arr(hh->bb, i, int)) ;
		}
	      aceOut (ao, "\n") ;
	    }
	}
      if (sx->plot)
	{
	  ao = aceOutCreateToFile (messprintf("%s.tcsh", sx->outFileName), "wx", h0) ;
	  aceOutf(ao, "#!/bin/tcsh -f\n") ;
	  aceOutf(ao, "  gnuplot -bg white << EOF\n") ;
	  aceOutf(ao, "    set style data linespoints\n") ;
	  aceOutf(ao, "    set terminal postscript color \n") ;
	  aceOutf(ao, "    set terminal postscript solid \n") ;
	  aceOutf(ao, "    set terminal postscript landscape\n") ;
	  aceOutf(ao, "    set output '%s.ps'\n", sx->outFileName) ;
	  if (sx->title) 
	    aceOutf (ao, "set title '%s'\n", sx->title) ;

	  for (h = 0 ; h < hMax ; h++)
	    {
	      hh = arrayp (hhh, h, HH) ;
	      if (arrayMax (hh->aa))
		aceOutf(ao, "%s '%s.txt' using %d with line title '%s'"
			, h ? ",\\\n\t" : "plot "
			, sx->outFileName
			, h + 2
			, hh->title 
			) ;
	    }
	  aceOutf(ao, "\nEOF\n") ;
	  aceOutDestroy (ao) ; /* will close */
	  system (messprintf("(sleep 1 ; tcsh %s.tcsh ; gv %s.ps ) &", sx->outFileName, sx->outFileName)) ;
	}
      ac_free (h0) ;
    }
  
  return sx->n ;
} /* histoPlainGetData */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*  Karhunen Loeve (a form of principal component analysis)
 * Given a 4 column table
 *  Id run count cover : data(id,run)
 * We compute for each id the smooth histo across all runs:  profile(id)[101]
 * Then we construct the K-L base of functions: KL(0...)[101]
 * Then for each id, we decompose it on the KL basis 
 * the we export the first few profiles and the histo of projection on these components
 */
#define KLMX 101
typedef struct klvStruct {int id ; double z[KLMX];} KLV ;

/* given two vectors, a[] and b[], return the scalar product */
static double klScalProd (KLV *a, KLV *b)
{
  int i = KLMX ;
  double s = 0 ;
  while (i--) s += a->z[i] *  b->z[i] ;

  return s ;
} /* klScalProd */

/*************************************************************************************/
/* given a vector a[], and a base vector b[], return the scalar product and sets 
 * the orthogonal remainder vector
 */
static double klProject (KLV *a, KLV *b, KLV *remainder)
{
  double s = klScalProd(a,b) ;
  int i = KLMX ;
  while (i--)  remainder->z[i] = a->z[i] - s * b->z[i] ;
  return s ;
} /* klProject */

/*************************************************************************************/
/* select as new pivot the vector with the highest component in any direction
 */

static int klSelectPivot (Array aa, KLV *pivot)
{
  double z, *zp, zmax = 0 ;
  int i, ii, bestii = 0 ;
  KLV *up ;

  for (ii = 0, up = arrp (aa, 0, KLV) ; ii < arrayMax (aa) ; ii++, up++)
    {
      up = arrp (aa, ii, KLV) ; 
      for (i = 0, zp = up->z ; i < KLMX ; i++, zp++)
	{
	  z = *zp ; if (z < 0) z = -z ;
	  if (z > zmax) { zmax = z ; bestii = ii ; }
	}
    }
  up = arrp (aa, bestii, KLV) ; 
  *pivot = *up ;

  return bestii ;
} /* klSelectPivot */

 /*************************************************************************************/
/* compute the normalized average profile of a set of an Array of vectors
 * anti-symmetrizing with respect to the yperplan orthogonal to the pivot 
 */

static int klAverage (KLV *av, Array aa, KLV *pivot)
{
  double *xp, *zp, zz ;
  int ii = arrayMax (aa) ;
  int i, sign ;
  KLV *up ;

  memset (av, 0, sizeof (KLV)) ;
  while (ii--)
    {
      up = arrp (aa, ii, KLV) ; 
      zp = av->z ; xp = up->z ;
      sign = (pivot && klScalProd (up, pivot) < 0 ? -1 : 1) ;
      if (sign < 0) continue ;
      i = KLMX ;
       while (i--)
	*zp++ += sign * *xp++ ;
    }
  
  zz = klScalProd (av, av) ;

  for (zz = sqrt(zz), i = 0, zp = av->z ; i < KLMX ; i++, zp++)
    *zp /= zz ;
  *pivot = *av ;
  return arrayMax (aa) ;
} /* klAverage */

/*************************************************************************************/
/* given and array of profiles, compute the KL base
 */
static void klBase (Array aa, Array base, Array projections, KLV *energy)
{
  int kk ; /* successive base vectors */
  int ii = arrayMax (aa), nPivot ;
  KLV *proj, pivot ;
  double z, zz = 1 ;
  Array remainder = arrayCopy (aa) ;

  /* on the first pass we compute the tru average around a null pivot
   * on the following pass we select as pivot the position with the
   * highest coverage and pivot around it nPivot times
   */
  memset (&pivot, 0, sizeof(KLV)) ;
  for (kk = 0 ; kk < KLMX && kk < arrayMax(aa) && zz > 1.0E-6 ; kk++)
    {
      if (kk) klSelectPivot (remainder, &pivot) ;
      for (nPivot = 0 ; nPivot < (kk ? 1 : 5) ; nPivot++)
	klAverage (arrayp (base, kk, KLV), remainder, &pivot) ;
      zz = 0 ; ii = arrayMax (aa) ;
      while (ii--)
	{
	  proj = arrayp (projections, ii, KLV) ;
	  z = proj->z[kk]  = klProject (arrayp (remainder, ii, KLV), arrayp (base, kk, KLV), arrayp (remainder, ii, KLV)) ;
	  zz += z * z ;
	}
      energy->z[kk] = zz ;
    }
  ac_free (remainder) ;
  
  return ;
} /* klBase */

/*************************************************************************************/

static void klTransform (SX *sx, Array aa, Array *basep, Array *projectionsp, KLV *energy)
{
  int ii = arrayMax (aa) ;

  *basep = arrayHandleCreate (ii, KLV, sx->h) ;
  *projectionsp = arrayHandleCreate (ii, KLV, sx->h) ;

  klBase (aa, *basep, *projectionsp, energy) ;

  return ;
} /* klTransform */

/*************************************************************************************/
/* locate the 3 main maxima */
static void klExportSmoothedTypology (ACEOUT ao, KLV *up)
{
  int i, xMax[5], oldx, nMin = 0, nMax = 0, slope ;
  double yMax[5], z = 0, *zp, oldz ;
  BOOL isLog = 1 ; 
  double l10 = 1, delta = .01 ;

  memset (xMax, 0, sizeof (xMax)) ;
  memset (yMax, 0, sizeof (yMax)) ;
  if (isLog)
    { 
      l10 = log((double)10.0) ;
      delta = 0.06 ;  /* 1/10  in log10 scale */
    }
  oldx = -1 ; oldz = -10 ; slope = -1 ;
  for (i = 0, zp = up->z ; i <= KLMX ; i++, zp++)
    {
      if (i < KLMX)
	{
	  if (isLog)
	    { 
	      z = *zp ; if (z < 0.001) z = 0.001 ;
	      z = log(z/0.001)/l10 ;
	    }
	}

      if (slope == -1)
	{
	  if (z < oldz)
	    {
	      oldz = z ; oldx = i ;
	    }
	  else if (z - oldz > delta)
	    {
	      if (oldx > 0) {  nMin++ ; }
	      slope = 1 ; 
	      oldz = z ; oldx = i ;
	      if (nMin >= 5) break ;
	    }
	}
      else if (slope == 1)
	{
	  if (z > oldz)
	    {
	      oldz = z ; oldx = i ;
	    }
	  else if (z - oldz < - delta)
	    {
	      xMax[nMax] = oldx ; nMax++ ;
	      slope = -1 ; 
	      oldz = z ; oldx = i ;
	      if (nMax >= 5) break ;
	    }
	}
    }

  aceOutf (ao, "\t%d", nMax) ;
  for (i = 0 ; i < nMax ; i++)
    aceOutf (ao, "\t%d:%.3f", xMax[i], yMax[i]) ;
  for (; i < 5 ; i++)
    aceOutf (ao, "\t-:-") ;
} /* klExportSmoothedTypology */

/*************************************************************************************/

static void klExportSmoothedProfiles (SX *sx, Array aa, DICT *ids)
{
  int i, ii ;
  double *zp ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".smoothed_profiles", sx->gzo, sx->h) ;
  KLV *up ;

  aceOut (ao, "#\t\t\t\t\t\t\t") ;
  for (i = 0 ; i < KLMX ; i++)
    aceOutf (ao, "\tx%d", i) ; 
  aceOutf (ao, "\n") ;
  for (ii = 0, up = arrp (aa, 0, KLV) ; ii < arrayMax (aa) ; ii++, up++)
    {
      up = arrp (aa, ii, KLV) ; 
      aceOutf (ao, "%s", dictName (ids, up->id)) ;
      klExportSmoothedTypology (ao, up) ;
      aceOutf (ao, "\t%s", dictName (ids, up->id)) ;
      for (i = 0, zp = up->z ; i < KLMX ; i++, zp++)
	aceOutf (ao, "\t%.5f", *zp) ;
      aceOutf (ao, "\n") ;
    }
}  /* klExportSmoothedProfiles */

/*************************************************************************************/

static void klExportMatrix (SX *sx, Array base, Array projections, KLV *energy, int nn)
{
  int i, j, ii ;
  double *z ;
  double w[11] ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".KL.txt", sx->gzo, sx->h) ;

  aceOutf (ao, "# %s\n", timeShowNow ()) ;
  aceOutf (ao, "# Karhunen Loeve first few profiles\nComponent\tEnergy") ;
  aceOutf (ao, "Component") ;
  for (i = 0 ; i < KLMX ; i++)
    aceOutf (ao, "\t%d", i) ;
  aceOutf (ao, "\t\tNumber of profiles with a projection on this component in the range") ;
  for (i = 0 ; i <= 10 ; i++)
    aceOutf (ao, "\t%d%%", i) ;
  for (ii = 0 ; ii < nn && ii < arrayMax (base) ; ii++)
    {
      z = arrp (base, ii, KLV)->z ;
      aceOutf (ao, "\n%d\t%.2f", ii + 1, energy->z[ii]) ;
      for (i = 0 ; i < KLMX ; i++)
	aceOutf (ao, "\t%.2g", z[i]) ;

      /* construct the histogram of the projections on that axis */
      memset (w, 0, sizeof(w)) ;
      for (j = 0 ; j < arrayMax(projections) ; j++)
	{
	  z = arrp (projections, j, KLV)->z ;
	  i = 100 * z[ii] + .5 ; /* find the bin */
	  if (i > 10) i= 10 ;
	  if (i>=0) w[i]++ ;
	}

      /* export the histogram of the projections on that axis */
      aceOutf (ao, "\t\t") ;
      for (i = 0 ; i <= 10 ; i++)
	aceOutf(ao, "\t%.2f", w[i]) ;
    }
  aceOutf (ao, "\n") ;

} /* klExportMatrix */

/*************************************************************************************/

typedef struct kldataStruct { int id, sample, count, cover ;} KLDATA ;

static int kldataOrder (const void *va, const void *vb)
{
  const KLDATA *a = (const KLDATA *)va, *b = (const KLDATA *)vb ;
  int n ;

  n = a->id - b->id ; if (n) return n ;
  n = a->sample - b->sample ; if (n) return n ;

  return 0 ;
} /* kldataOrder */

/*************************************************************************************/
/* 2014: modify this format so that we can parse the 2014 snp format */

static Array klParse (SX *sx, DICT *ids, DICT *samples)
{
  AC_HANDLE h = ac_new_handle () ;
  BitSet ok = bitSetCreate (100000, h) ;
  int id, count, pos, cover, sample, i, ii, jj, nn = 0 ;
  double zz, zz0, zz2 ;
  const char *ccp ;
  Array data = arrayHandleCreate (1000000, KLDATA, h) ;
  Array aa = 0 ;
  KLDATA *up ;
  KLV *v ;
  ACEOUT ao = aceOutCreate (sx->outFileName, ".filter",sx->gzo,h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  Array hhhh = arrayHandleCreate (101, double, h) ;
  int NN = 200 ;
  Array dd = arrayHandleCreate (NN * NN, int, h) ;

  aceOutf (ao, "# %s\n",  timeShowNow ()) ;
  aceOutf (ao, "# Norm of the smoothed histogram and its good part\n") ;
  aceOutf (ao, "# id\tnorm\t0 || >40\n") ;
  /* parse the counts */
  while (aceInCard (sx->ai))
    {
      ccp = aceInWord (sx->ai) ;
      if (! ccp) continue ;

      vtxtClear (txt) ;
      vtxtPrint (txt, ccp) ;
      
      if (1)
	{
	  aceInInt (sx->ai, &pos) ;
	  ccp = aceInWord (sx->ai) ; 
	  vtxtPrintf (txt, ":%d:%s", pos,ccp ? ccp : "-") ;
	  
	  ccp = aceInWord (sx->ai) ; 
	  if (0) vtxtPrintf (txt, ":%s", ccp ? ccp : "-") ;
	  ccp = aceInWord (sx->ai) ; 
	  if (0) vtxtPrintf (txt, ":%s", ccp ? ccp : "-") ;
	}

      dictAdd (ids, vtxtPtr (txt) , &id) ;

      ccp = aceInWord (sx->ai) ;
      if (! ccp) continue ;
      dictAdd (samples, ccp, &sample) ;
      
      if (1)
	{
	  ccp = aceInWord (sx->ai) ; /* genotype, drop it */
	  ccp = aceInWord (sx->ai) ; /* frequency */
	}

      if (! aceInInt (sx->ai, &cover) || ! aceInInt (sx->ai, &count))
	continue ;
      if (count > cover) /* back compatibility with a format use in KL in 2011 */
	{ int dummy = count ; count = cover ; cover = dummy ; }

      up = arrayp (data, nn++, KLDATA) ;
      up->id = id ;
      up->sample = sample ;
      up->count = count ;
      up->cover = cover ;
    }
  
  arraySort (data, kldataOrder) ;  
  
  /* for each id, verify that at least once 
   * it is frequency >= sx->KL, # mut >= 3, count >= 4
   */
  for (ii = id = 0, up = arrp (data, 0, KLDATA) ; ii < arrayMax(data) ; up++, ii++)
    {
      if (up->count >= 5 && up->cover >= 10 && 100 * up->count >= sx->KL * up->cover)
	bitSet (ok, up->id) ;
    }

  /* for each id, create the smooth histo corresponding
   * to the set {count, cover} over all its samples 
   */
  aa = arrayHandleCreate (bitSetCount (ok) + 10, KLV, sx->h) ;

  for (ii = jj = id = 0, v = 0, up = arrp (data, 0, KLDATA) ; ii < arrayMax(data) ; up++, ii++)
    {
      if (! bit (ok, up->id))
	continue ;
      if (up->cover < 10)
	continue ;
      if (up->id != id) 
	v = arrayp (aa, jj++, KLV) ;
      id = up->id ;
      v->id = up->id ;
 
      if (up->count < NN && up->cover < NN)
	arr (dd, up->count * NN + up->cover, int) += 1 ;
      else 
	utSmoothHisto (hhhh,  up->count, up->cover, 1) ;
    }
  {
    int i, k, n ;
    for (k = 0 ; k < NN ; k++)
      for (n = 1 ; n < NN ; n++)
	{
	i = arr (dd, k * NN + n, int) ;
	if (i > 0)
	  utSmoothHisto (hhhh, k, n, i) ;
	}
  }
  

   for (i = 0 ; i < KLMX ; i++)
	v->z[i] += array (hhhh, i, double) ;

  /* normalize the individual histos */
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      v = arrayp (aa, ii, KLV) ;  
      if (! v->id)
	continue ;
      /* normalize */
      for (i = 0, zz = 0 ; i < KLMX ; i++)
	zz += v->z[i] ; /* do not square */
      if (zz > 0)
	{
	  zz = sqrt (zz) ;
	  for (i = 0 ; i < KLMX ; i++)
	    v->z[i] /= zz ;
	}

      /* export the list of good probes 
       * defined as the NOT SQUARED support over the good zone 
       * is 90% of the smoothed histo
       */
       for (i = 0, zz = 0, zz0 = zz2 = 0 ; i < KLMX ; i++)
	 {
	   zz += v->z[i] ;
	   if (i == 0)
	     zz0 += v->z[i] ;
	   if (i >= 40)
	     zz2 += v->z[i] ;
	 }
       if (zz > 0)
	 aceOutf (ao, "%s\t%.3f\t%.3f\n", dictName(ids, v->id), 100.0 * zz0/zz, 100.0 * zz2/zz) ;
    }
  ac_free (h) ;
  
  return aa ;
} /* klParse */

/*************************************************************************************/

static void klAnalyse (SX *sx)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa, base, projections ;
  KLV energy ;
  DICT *ids = dictHandleCreate (10000, sx->h) ;
  DICT *samples = dictHandleCreate (100, sx->h) ;

  /* parse the 4 columns data file */
  aa = klParse (sx, ids, samples) ;

  klExportSmoothedProfiles (sx, aa, ids) ;

  if (0) /* the KL transform is not very useful */
    {
      /* construct the K-L base and project all profile on this basis */
      klTransform (sx, aa, &base, &projections, &energy) ;
      
      /* export the profiles */
      klExportMatrix (sx, base, projections, &energy, 10) ;
    }

  ac_free (h) ;
} /* klAnalyse */

/*************************************************************************************/

static void klParseRunList (SX *sx) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (sx->runFileName, FALSE, h) ;
  
  sx->runDict = dictHandleCreate (100, sx->h) ;
  while (aceInCard (ai))
    {
      int run ;
      char *cq, *cp = aceInWord (ai) ;
      if (! cp || !*cp || *cp == '#')
	continue ;
      cq = strstr (cp, "#@#") ;
      if (cq) *cq = 0 ;
      dictAdd (sx->runDict, cp, &run) ;
      if (! cq || ! cq[3])
	continue ;
      if (! sx->run2title)
	sx->run2title = assHandleCreate (sx->h) ;
      cp = strnew (cq+3, sx->h) ;
      assInsert (sx->run2title, assVoid (run), cp) ;
    }

  ac_free (h) ;
  return ;
} /* klParseRunList */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// histo:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Export a smoothed histogram\n"
	    "// Algorithm\n"
	    "//   -w  [default 100] is the number of point in the histo\n" 
	    "//   -plain Compute a standard hitogram with fixed bin size\n"
	    "//      The titles are taken from the first line of the input if it starts with a #\n"
	    "//   -smooth  smooth the values as a Bayesian binomial\n"
	    "//      In this case input give 2 int per line, observed/total\n"
	    "//         3 6 # gives a large peak around 50%%\n"
	    "//         100 200 # gives a sharp peak around 50%%\n"
	    "//   -columns 1,3,7,...   or   2,3-7,11-15  or 2,3,7+ \n"
	    "//      The optional list of column (default 1) help interpret the input file\n"
	    "//      By default the histo is contructed from column 1, (or 1-2 in -smooth case)\n"
	    "//      1,3,7 will produce three histos from columns 1, 3, 7 (or 1-2, 3-4 and 7-8)\n" 
	    "//      1,3-7 will consider columns 1,3,4,5,6,7\n"
	    "//      2,3,7+ will consider 2,3,7,8,9,10.... to the end of the table\n"
	    "//   -colList <fileName>\n"
	    "//      Sorted list of column names, one per line, control the order of the columns\n"
	    "//      If present, the second column of the file gives associated titles\n"
	    "//   -snp\n"
	    "//      The input is a snp count table, run/cover/count in columns 6/9/10\n"
	    "//      The program exports the corresponding smoothed multi histo\n"
	    "//   -min <number> -max <number>\n"
	    "//      Ignore all values outside [min, max]\n"
	    "// Plot\n"
	    "//   -plot  : requires -o file_name, and the programs \"gnuplot\" and \"gv\"\n"
	    "//          idem but exports the distribution in fil_name.txt and a picture in file_name.ps\n"
	    "//          then directly displays the postscript file file_name.ps using gv\n"

#ifdef JUNK
	    "// Karhunen Loeve\n"
	    "//   -KL <n> : <n> is the min frequency of SNP, say 20 or 5\n"
	    "//      The input has 4 columns: id test count cover\n"
	    "//        For each id in each test report the count/cover as in the option -smooth\n"
	    "//        The system computes the smooth histo for each id integrating over the tests\n"
	    "//        The the Karhunen Loeve transform is applied to this collection of profiles\n"
	    "//        The program finally export the first 10 eigen K-L profiles, their \"energy\"\n"
	    "//        and the histogram of the projection of the id profile on these 10 base profiles\n"
#endif
	    "// Input/Output files\n"
	    "//   -i filename : optional, input file possibly gzipped (default stdin)\n"
	    "//   -o filename : optional, ouput file name prefix (default stdout)\n"
	    "//      In plot case several files are exported with different suffixes\n"
	    "//   -gzi : input file should be unzipped\n"
	    "//      This is the default behaviour if input file is called .gz\n"
	    "//      But it is useful in cases like 'cat *.gz | histo -gzi ...'\n"
	    "//\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' are considered as comments and dropped out\n"
	    "//   Except the eventual first line with a # which is used as a title line\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  const char *ccp = 0 ;
  SX sx ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&sx, 0, sizeof (SX)) ;
  sx.h = h ;


  /* optional arguments */

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  sx.NN = 100 ;
  sx.smooth = getCmdLineOption (&argc, argv, "-smooth", 0) ;
  sx.plain = getCmdLineOption (&argc, argv, "-plain", 0) ;

  sx.gzi = getCmdLineOption (&argc, argv, "-gzi", 0) ;
  sx.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;

  getCmdLineInt (&argc, argv, "-w", &sx.NN) ;
  sx.min = - 1e30 ;
  sx.max =  1e30 ;
  if (getCmdLineFloat (&argc, argv, "-min", &sx.min)) sx.hasMin = TRUE ;
  if (getCmdLineFloat (&argc, argv, "-max", &sx.max)) sx.hasMax = TRUE ;
  if (sx.NN < 1) sx.NN = 1 ; 
  getCmdLineOption (&argc, argv, "-i", &sx.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &sx.outFileName) ;
  getCmdLineOption (&argc, argv, "-title", &sx.title) ;
  getCmdLineOption (&argc, argv, "-columns", &sx.columns) ;
  getCmdLineOption (&argc, argv, "-colList", &sx.runFileName) ;
  getCmdLineInt (&argc, argv, "-KL", &sx.KL) ;
  sx.snp = getCmdLineOption (&argc, argv, "-snp", 0) ;
  sx.plot = getCmdLineOption (&argc, argv, "-plot", 0) ;
  if (sx.plot && ! sx.outFileName)
    usage ("-plot option requires -o option") ;

  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[argc-1])) ;

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  if (sx.runFileName)
    klParseRunList (&sx) ;

  sx.cols = keySetHandleCreate (h) ;
  if (sx.columns)
    {
      int n = 0, jj = 0, isRange = 0 ;

      ccp = sx.columns - 1 ;
      while (*++ccp)
	{
	  if (*ccp >= '0' && *ccp <= '9')
	    n = 10 * n + (*ccp - '0') ;
	  else if (n && *ccp == '+')
	    {
	      for (; jj < 1000 ; jj++)
		keySet (sx.cols, jj) = n++ ;
	      n = 0 ;
	      break ;
	    }
	  else if (*ccp == ',')
	    {
	      if (n) 
		{
		  if (isRange)
		    {
		      int i ;
		      for (i = isRange + 1 ; i < n ; i++)
			keySet (sx.cols, jj++) = i ;
		    }
		  isRange = 0 ;
		  keySet (sx.cols, jj++) = n ;
		}
	      n = 0 ;
	    }
	  else if (*ccp == '-')
	    {
	      if (n) 
		{
		  if (isRange)
		    {
		      int i ;
		      for (i = isRange + 1 ; i < n ; i++)
			keySet (sx.cols, jj++) = i ;
		    }
		  keySet (sx.cols, jj++) = n ;
		}
	      isRange = n ;
	      n = 0 ;
	    }
	}
      if (n) 
	{
	  if (isRange)
	    {
	      int i ;
	      for (i = isRange + 1 ; i < n ; i++)
		keySet (sx.cols, jj++) = i ;
	    }
	  keySet (sx.cols, jj++) = n ;
	}
    }
  else
    keySet (sx.cols, 0) = 1 ; /* default */

  
  if (1)
    {
      sx.ai = aceInCreate (sx.inFileName, sx.gzi, h) ;
      if (!sx.ai)
	usage (messprintf("%s not found", sx.inFileName ? sx.inFileName : "stdin" )) ;
      aceInSpecial (sx.ai,"\t\n") ;
    }

  if (sx.snp)
    {
      histoSnpGetData (&sx) ;
    }
   else if (sx.smooth)
    {
      histoSmoothGetData (&sx) ;
      histoSmoothExport (&sx) ;
    }
 else if (sx.plain)
    {
      histoPlainGetData (&sx) ;
    }
  else if (sx.KL)
    {
      klAnalyse (&sx) ;
    }
  else
    usage ("Please specify -smooyh or -bin") ;
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  
  
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

