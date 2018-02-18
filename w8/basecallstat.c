/*  File: basecallstat.c
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 17 14:20 1996 (ulrich)
 * Created: Tue Jan 16 13:23:59 1996 (ulrich)
 *-------------------------------------------------------------------
 */

/* @(#)basecallstat.c	1.3  11/12/96  */
/*
#define CHRONO
*/
#define ARRAY_CHECK 
#define MALLOC_CHECK

#include "acedb.h"
#include "topology.h"
#include "dnaalign.h"
#include "query.h"
#include "dna.h"
#include "a.h"
#include "pick.h"
#include "plot.h"
#include "lex.h"
#include "bs.h"
#include "tags.h"
#include "systags.h"
#include "classes.h"
#include "display.h"
#include "aligntools.h"
#include "acembly.h"

/***************************************************************/
/**************************** Tests ****************************/
/***************************************************************/

static void statisticsAcetestExport (Array aa, Array bb, Array cc, int nseq, int clone,
				     int fragment, int frequence, int *debut, int *fin)
{ int i, *ip, *ip0, nc = 0, bt = 0, j ;
  BOOL start = FALSE ;

  i = arrayMax (aa) ;
  ip = ip0 = arrp (aa, 0, int) ;
  if (*ip) *debut = 0 ;
  while (i > 0)
    { while (i && *ip)
	{ ip++ ;
	  i-- ;
	  if (!start)
	    start = TRUE ;
	}
      *fin = ip - ip0 ;
      if (!i)
	{ if (start)
	    nc++ ;
	  break ;
	}
      while (i && !(*ip))
	{ ip++ ;
	  i-- ;
	  bt++ ;
	  if (start)
	    { nc++ ;
	      start = FALSE ;
	    }
	}
      if (*debut == -1)
	*debut = ip - ip0 ;
    }
  j = nseq/frequence ;
  array (bb, j, float) += nc ;
  array (cc, j, float) += bt ;
}

/***************************************************************/

static void statisticsMakePseudoShotgun (Array table, int L)
{ int N = 0, e = 0, f = 0, gap, i = 0, j, k, *ip, parite, nbtest, nn, tab, maxN = 0 ;
  int debut = 0, fin = 0, debtot = 0, fintot = 0 ;
  float alpha, beta ;
  Array aa = 0, tablebis = 0, bb = 0, cc = 0 ;

  while (TRUE)
    { maxN = arrayMax (table)/2 ;
      debtot = fintot = 0 ;
      if (!messPrompt("Taille des fragments\n","300","i"))
	return ;
      freeint (&f) ;
      if (!messPrompt(messprintf("Nombre total de sequences (parmis %d)\n", maxN), "1000","i"))
	return ;
      freeint (&N) ;
      if (!messPrompt("Frequence exportation\n","10","i"))
	return ;
      freeint (&e) ;
      if (!messPrompt("Nombre de Test\n","1000","i"))
	return ;
      freeint (&nbtest) ;
      if (!L || !f || !N || !e || !nbtest) return ;
      if (N > maxN) return ;
      nn = nbtest ;
      i = N/e ;
      bb = arrayCreate (i, float) ;
      cc = arrayCreate (i, float) ;
      while (nn--)
	{ aa = arrayReCreate (aa, L, int) ;
	  array (aa, L - 1, int) = 0 ;
	  tablebis = arrayCopy (table) ;
	  for (k = 1 ; k <= N ; k++)
	    { i = -1 ;
	      while(i==-1)
		{ tab = 2*(randint () % maxN) ;
		  i = arr (tablebis, tab, int) ;
		}
	      arr (tablebis, tab, int) = -1 ; /* to prevent double */
	      parite = arr (tablebis, tab+1, int) ;
	      gap = f ;
	      if (i >= L)
		{ gap += L - i - 1 ;
		  i = L - 1 ;
		}
	      if (i < 0)
		{ gap += i ;
		  i = 0 ;
		}
	      if (gap < 0) gap = 0 ;
	      switch (parite)
		{ 
		case 0:
		  ip = arrp (aa, i, int) + 1 ;
		  j = i < gap ? i + 1 : gap ;
		  while (ip--, j--)
		    (*ip)++ ;
		  break ;
		case 1:
		  ip = arrp (aa, i, int) - 1 ;
		  j = L - i > gap ? gap : L - i ;
		  while (ip++, j--)
		    (*ip)++ ;
		  break ;
		default:
		  messerror ("Ca ne marche pas") ;
		  break ;
		}
	      if (!(k%e))
		{ debut = -1 ; fin = 0 ;
		  statisticsAcetestExport (aa, bb, cc, k, L, f, e, &debut, &fin) ;
		}
	    }
	  debtot += debut ;
	  fintot += fin ;
	  arrayDestroy (tablebis) ;
	}
      fprintf (stdout, "Clone : %d, Nombre seq : %d, Fragments : %d, Nombre test : %d, debut : %d, fin : %d\n",
	       L, maxN, f, nbtest, debtot/nbtest, fintot/nbtest) ;
      for (nn = 0 ; nn <= N/e ; nn++)
	{ alpha = array (bb, nn, float)/nbtest ;
	  beta = array (cc, nn, float)/nbtest ;
	  fprintf (stdout, "%d sequences, %f contigs, %f bases dans les trous\n",
		   nn*e, alpha, beta) ;
	}
      arrayDestroy (aa) ;
      arrayDestroy (bb) ;
      arrayDestroy (cc) ;
      if (!messPrompt("Pseudo Shotgun", "1", "i"))
	break ;
    }
}

/***************************************************************/
/*
static void statisticsTruncArrayDna (Array a, int x1, int x2)
{ char *cp, *cq ;

  if (!arrayExists (a) || !arrayMax (a))
    return ;
  if (x1 < 0 || x2 < x1 || x2 > arrayMax (a))
    messcrash 
      ("Bad coordinates x1 = %d, x2 = %d in statisticsTruncArray",
       x1, x2) ;
  cp = arrp (a, 0, char) ;
  cq = arrp (a, x1, char) ;
  i = x2 - x1 ;
  arrayMax (a) = i ;
  while (i--)
    *cp++ = *cq++;
}
*/
/***************************************************************/
/* renvoie le dna correspondant, seulement entre clipvecttop et clipvectend,
   les position des clipping correspondant aux coordonnees bio dans ce dna. */
static Array statisticsGetClipArray (KEY key, int *exce, int *good, int *fair)
{ Array dna = 0, dna2 = 0 ;
  OBJ obj = 0 ;
  KEY dnaKey, seqKey ;
  int x, cvt, cve = 0 ;

  if (class(key) == _VDNA)
    { dnaKey = key ;
      dnaReClass (key, &seqKey) ;
    }
  else if (dnaSubClass (key, &dnaKey))
    { seqKey = key ;
    }
  else goto abort ;
  if (!(obj = bsCreate (seqKey)) || !(dna2 = dnaGet (dnaKey)))
    goto abort ;
  if (!bsFindTag (obj, _Clips))
    { dna = dna2 ;
      dna2 = 0 ;
      *exce = arrayMax (dna) ;
      *good = *fair = *exce ;
      goto abort ;
    }
  if (bsGetData (obj, _Vector_Clipping, _Int, &cvt))
    bsGetData (obj, _bsRight, _Int, &cve) ;
  else cvt = 0 ;
  if (!cvt)
    bsGetData (obj, _Clipping, _Int, &cvt) ;
  if (bsGetData (obj, _Excellent_upto, _Int, &x))
    *exce = x - cvt ;
  else
    *exce = 0 ;
  if (bsGetData (obj, _Good_upto, _Int, &x))
    *good = x - cvt ;
  else
    *good = 0 ;
  if (bsGetData (obj, _Fair_upto, _Int, &x))
    *fair = x - cvt ;
  else
    *fair = 0 ;
  if (cvt) cvt-- ;
  if (!cve || cve > arrayMax(dna2))
    cve = arrayMax (dna2) ;
  if (cve < cvt)
    goto abort ;
  dna = arrayTruncatedCopy (dna2, cvt, cve) ;
  cve -= cvt ;
  if (*exce > cve)
    *exce = cve ;
  if (*good > cve)
    *good = cve ;
  if (*good < *exce)
    *good = *exce ;
  if (*fair > cve)
    *fair = cve ;
  if (*fair < *good)
    *fair = *good ;
 abort:
  arrayDestroy (dna2) ;
  bsDestroy (obj) ;
  return dna ;
}

/***************************************************************/

void statisticsDoMakeErreur (KEY target, KEYSET ks)
{ Array tarDna = 0, erreur = 0, dna = 0, shDna = 0, table = 0, tarDnaR = 0, tarD ;
  KEY key ;
  int i, j, sens, x1, x2, exce = 0, good = 0, fair = 0, itab = 0, pol = 0, nbn = 0, recou = 0 ;
  int nseq = 0, nbad = 0, nbaseE = 0, nbaseG = 0, nbaseF = 0, nbaseA = 0 ;
  int errE[16], errG[16], errF[16], errA[16], max1, max2 ;
  A_ERR *errp ;

  if (!ks || !keySetMax (ks))
    return ;
  if (!target || !(tarDna = dnaGet (target)))
    return ;
  max1 = arrayMax (tarDna) - 1 ;
  tarDnaR = arrayCopy (tarDna) ;
  reverseComplement (tarDnaR) ;
  messStatus ("Aligning fragments" ) ;
  i = 6 ;
  while (i--)
    { errE[i] = 0 ;
      errG[i] = 0 ;
      errF[i] = 0 ;
      errA[i] = 0 ;
    }
  i = keySetMax (ks) ;
  table = arrayCreate (2*i, int) ;
  while (i--)
    { key = keySet (ks, i) ;
      if (!(dna = statisticsGetClipArray (key, &exce, &good, &fair)) ||
	  !exce) /* voir commentaire */
	{ arrayDestroy (dna) ;
	  continue ;
	}
      shDna = arrayTruncatedCopy (dna, 0, exce) ;
      if (!(erreur = dnaAlignCompareDna (tarDna, shDna, &x1, &x2, &sens, FALSE)))
	{ arrayDestroy (shDna) ;
	  arrayDestroy (dna) ;
	  continue ;
	}
      j = arrayMax (erreur) ;
      if (j > arrayMax (shDna) / 6)
	{ printf ("plus de 16.67 %% d'erreur (/zone excellente)") ;
	  nbad++ ;
	  arrayDestroy (shDna) ;
	  arrayDestroy (dna) ;
	  arrayDestroy (erreur) ;
	  continue ;
	}
      if (sens == 1)
	{ array (table, itab++, int) = x1 + 1 ;
	  array (table, itab++, int) = 1 ;
	  tarD = tarDna ;
	  pol = x1 ;
	}
      else
	{ array (table, itab++, int) = x2 + 1 ;
	  array (table, itab++, int) = 0 ;
	  tarD = tarDnaR ;
	  pol = max1 - x2 ;
	}
      arrayDestroy (shDna) ;
      max2 = arrayMax (dna) - 1 ;
      nbn = 1 ;
      newLocalCptErreur (tarD, pol, max1, pol, dna, 0, max2, 0,
			 1, &nbn, &x1, &x2, &recou, erreur) ;
      j = arrayMax (erreur) ;
      nseq++ ;
      nbaseE += exce ;
      nbaseG += good - exce ;
      nbaseF += fair - good ;
      nbaseA += arrayMax (dna) - fair ;
      if (j)
	{ errp = arrp (erreur, 0, A_ERR) - 1 ;
	  while (errp++, j--)
	    {
	      if (errp->iShort > exce)
		{
		  if (errp->iShort > good)
		    {
		      if (errp->iShort > fair)
			errA[errp->type]++ ;
		      else
			errF[errp->type]++ ;
		    }
		  else
		    errG[errp->type]++ ;
		}
	      else
		errE[errp->type]++ ;
	    }
	}
      arrayDestroy (dna) ;
      arrayDestroy (erreur) ;
    }
  fprintf (stdout, "************ Bilan *********\n") ;
  fprintf (stdout, "nb de sequence : %d. (%d rejetees pour trop d'erreur)\n", nseq, nbad) ;
  fprintf (stdout, "Excellent : nb de base : %d\n", nbaseE) ;
  if (nbaseE)
    { fprintf (stdout, "            Ambiguite : %d ; %d o/oo\n", errE[AMBIGUE],
	       1000*errE[AMBIGUE]/nbaseE) ;
      fprintf (stdout, "            Erreur :    %d ; %d o/oo\n", errE[ERREUR],
	       1000*errE[ERREUR]/nbaseE) ;
      i = errE[INSERTION] + errE[INSERTION_DOUBLE] ;
      fprintf (stdout, "            Insertion : %d ; %d o/oo\n", i, 1000*i/nbaseE) ;
      i = errE[TROU] + errE[TROU_DOUBLE] ;
      fprintf (stdout, "            Trou :      %d ; %d o/oo\n", i, 1000*i/nbaseE) ;
    }
  fprintf (stdout, "Good      : nb de base : %d\n", nbaseG) ;
  if (nbaseG)
    { fprintf (stdout, "            Ambiguite : %d ; %d o/oo\n", errG[AMBIGUE],
	       1000*errG[AMBIGUE]/nbaseG) ;
      fprintf (stdout, "            Erreur :    %d ; %d o/oo\n", errG[ERREUR],
	       1000*errG[ERREUR]/nbaseG) ;
      i = errG[INSERTION] + errG[INSERTION_DOUBLE] ;
      fprintf (stdout, "            Insertion : %d ; %d o/oo\n", i, 1000*i/nbaseG) ;
      i = errG[TROU] + errG[TROU_DOUBLE] ;
      fprintf (stdout, "            Trou :      %d ; %d o/oo\n", i, 1000*i/nbaseG) ;
    }
  fprintf (stdout, "Fair      : nb de base : %d\n", nbaseF) ;
  if (nbaseF)
    { fprintf (stdout, "            Ambiguite : %d ; %d o/oo\n", errF[AMBIGUE],
	       1000*errF[AMBIGUE]/nbaseF) ;
      fprintf (stdout, "            Erreur :    %d ; %d o/oo\n", errF[ERREUR],
	       1000*errF[ERREUR]/nbaseF) ;
      i = errF[INSERTION] + errF[INSERTION_DOUBLE] ;
      fprintf (stdout, "            Insertion : %d ; %d o/oo\n", i, 1000*i/nbaseF) ;
      i = errF[TROU] + errF[TROU_DOUBLE] ;
      fprintf (stdout, "            Trou :      %d ; %d o/oo\n", i, 1000*i/nbaseF) ;
    }
  fprintf (stdout, "Bad       : nb de base : %d\n", nbaseA) ;
  if (nbaseA)
    { fprintf (stdout, "            Ambiguite : %d ; %d o/oo\n", errA[AMBIGUE],
	       1000*errA[AMBIGUE]/nbaseA) ;
      fprintf (stdout, "            Erreur :    %d ; %d o/oo\n", errA[ERREUR],
	       1000*errA[ERREUR]/nbaseA) ;
      i = errA[INSERTION] + errA[INSERTION_DOUBLE] ;
      fprintf (stdout, "            Insertion : %d ; %d o/oo\n", i, 1000*i/nbaseA) ;
      i = errA[TROU] + errA[TROU_DOUBLE] ;
      fprintf (stdout, "            Trou :      %d ; %d o/oo\n", i, 1000*i/nbaseA) ;
    }
  fprintf (stdout, "****************************\n") ;

  arrayDestroy (tarDna) ;
  arrayDestroy (tarDnaR) ;
  statisticsMakePseudoShotgun (table, max1+1) ;
  arrayDestroy (table) ;
}

/***************************************************************/

void statisticsMakeErreur (KEYSET ks)
{ KEY target ;
  char *cp ;

  if (!ks || !keySetMax (ks))
    { messout ("Please, load sequences first") ;
      return ;
    }
 again:
  if (messPrompt("Choose the Target", "", "wz"))
    { cp = freeword() ;
      if (!lexword2key(cp, &target, _VSequence))
	{ messout("Unknown Target, should be class Sequence") ;
	  goto again ;
	}
    }
  else return ;
  statisticsDoMakeErreur (target, ks) ;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct { int u1, u2, g ; } FIT ;

void statisticsCountGroup (KEYSET recu)
{ Array dna1 = 0, dna1r = 0, dna1rr, dna2 = 0, mm = 0, fit = 0, temp = 0,
    barreaux = 0, barrnonok = 0, nbok = 0 ;
  KEY key1, key2 ;
  KEYSET ks = 0 ;
  int i, j, paire, k, gg, nn, nn0, nnt, taille, tour, max, zone, paireg, paireb, iii = 3, jjj ;
  int tt[] = {5, 7, 9}, valnn[] = {7, 5, 4}, valmax[] = {10000, 5000, 2000} ;
  FIT *f ;

/*  if (!messPrompt ("Taille\n", "7", "i"))
    return ;
  freeint(&taille) ;
  if (!messPrompt ("nn0\n", "5", "i"))
    return ;
  freeint(&nn0) ;
  if (!messPrompt ("Nb de couples maximum consideres\n", "2000", "i"))
    return ;
  freeint(&max) ;
  if (!messPrompt ("Zone d'etude (nb base a chaque bout)\n", "500", "i"))
    return ;
  freeint(&zone) ;
*/
  if (!recu || (keySetMax (recu) < 2))
    { messout ("Please, select a keySet containing sequences") ;
      return ;
    }
  zone = 500 ;
  while (iii--)
    { taille = tt[iii] ;
      nn0 = valnn[iii] ;
      jjj = 3 ;
      while (jjj--)
	{ max = valmax [jjj] ;
	  paire = 0 ;
	  paireg = 0 ;
	  paireb = 0 ;
	  barreaux = arrayCreate (200, int) ;
	  barrnonok = arrayCreate (200, int) ;
	  nbok = arrayCreate (200, int) ;
	  ks = keySetCopy (recu) ;
	  for (i = 0 ; i < keySetMax (ks) - 1 ; i++)
	    { key1 = keySet (ks, i) ;
	      if (!key1) continue ;
	      if (!(dna1 = dnaGet (key1)))
		{ keySet (ks, i) = 0 ;
		  continue ;
		}
	      if (arrayMax (dna1) < 2*zone)
		{ arrayDestroy (dna1) ;
		  keySet (ks, i) = 0 ;
		  continue ;
		}
	      dna1r = arrayCopy (dna1) ;
	      reverseComplement (dna1r) ;
	      for (j = i + 1 ; j < keySetMax (ks) ; j++)
		{ key2 = keySet (ks, j) ;
		  if (!key2) continue ;
		  if (!(dna2 = dnaGet (key2)))
		    { keySet (ks, j) = 0 ;
		      continue ;
		    }
		  if (arrayMax (dna2) < 2*zone)
		    { arrayDestroy (dna2) ;
		      keySet (ks, j) = 0 ;
		      continue ;
		    }
		  paire++ ;
		  dna1rr = dna1 ;
		  tour = 2 ;
		  nnt = 0 ;
		  while (tour--)
		    { mm = alignToolsMakeShortMatch (dna1rr, dna2, 1, taille, max, zone) ;
		      if (mm && arrayMax(mm))
			{ gg = dnaAlignMakeGroups (mm, &fit, TRUE, taille) ;
			  temp = arrayReCreate (temp, 100, int) ;
			  k = arrayMax(fit) ;
			  f = arrp(fit, 0, FIT) - 1 ;
			  nn = 0 ;
			  while (f++, k--)
			    { array (temp, f->g, int)++ ;
			      if (f->g == gg)
				nn++ ;
			    }
			  if (nn > nn0)
			    { paireg++ ;
			      for (k = 0 ; k < arrayMax (temp) ; k++)
				if ((nn = arr (temp, k, int)))
				  { array (barreaux, nn, int)++ ;
				    if (nn > nn0)
				      nnt++ ;
				  }
			    }
			  else
			    { paireb++ ;
			      for (k = 0 ; k < arrayMax (temp) ; k++)
				if ((nn = arr (temp, k, int)))
				  array (barrnonok, nn, int)++ ;
			    }
			}
		      arrayDestroy (fit) ;
		      arrayDestroy (mm) ;
		      dna1rr = dna1r ;
		    }
		  array (nbok, nnt, int)++ ;
		  arrayDestroy (dna2) ;
		}
	      arrayDestroy (dna1) ;
	      arrayDestroy (dna1r) ;
	    }
	  keySetDestroy (ks) ;
	  k = array (barreaux, arrayMax (barrnonok) - 1, int) ; /* to make same size */
	  k = array (barrnonok, arrayMax (barreaux) - 1, int) ;
	  fprintf (stdout, "************ Bilan *********\n") ;
	  fprintf (stdout, "Taille des oligos : %d\n", taille) ;
	  fprintf (stdout, "Nombre min de barreaux pour accepter : %d\n", nn0) ;
	  fprintf (stdout, "Nombre maximum de couples consides : %d\n", max) ;
	  fprintf (stdout, "Zone d'etude a chaque bout du DNA : %d\n", zone) ;
	  fprintf (stdout, "Nombre de paire etudiees : %d\n", paire) ;
	  fprintf (stdout, "Nombre de paire Ok : %d ; Pas Ok : %d\n", paireg, paireb) ;
	  fprintf (stdout, "Nombre de barreaux trouves dans les groupes\n") ;
	  for (k = 0 ; k < arrayMax (barreaux) ; k++)
	    fprintf (stdout, "      %d   ; nb fois quand OK : %d ; quand pas OK : %d\n", k, 
		     array (barreaux, k, int), array (barrnonok, k, int)) ;
	  fprintf (stdout, "\n\nnombre de groupe superieur a nn0 dans une meme paire\n") ;
	  for (k = 0 ; k < arrayMax (nbok) ; k++)
	    fprintf (stdout, "      %d   ; nb fois : %d\n", k, array (nbok, k, int)) ;
	  fprintf (stdout, "****************************\n\n") ;
/* #ifndef NON_GRAPHIC
	  plotHisto ("Nombre de barreaux si ok", barreaux) ;
	  plotHisto ("Nombre de barreaux sinon", barrnonok) ;
#else */
	  arrayDestroy (barreaux) ;
	  arrayDestroy (barrnonok) ;
/* #endif */
	  arrayDestroy (temp) ;
	  arrayDestroy (nbok) ;
	}
    }
}
  
/***************************************************************/

void statisticsTestOligo (void)
{ 
  DEFCPT look = 0 ;
  KEYSET aa ;
  Associator adna ;
  Array dna = 0 ;

  int nbolig = 0, max, j, ii[] = {2, 5, 10}, i, tt[] = {15, 12, 9, 7, 5}, ntt = 5, taille ;
  KEY key ;
  KEYSET ks = 0 ;
  Array result = 0 ;
  unsigned int test = 0, mask ;
  BOOL depart = FALSE ;
  char *cp, *vp, *vp0 = 0 ;
  char u[9] ;

  u[A_] = 0 ;
  u[T_] = 3 ;
  u[G_] = 1 ;
  u[C_] = 2 ;
#ifndef NON_GRAPHIC
  if (!keySetActive (&aa, 0))
    { 
      messout ("First select a contig") ;
      return ;
    }
#else
  return ;
#endif
/* if (!messPrompt ("Nombre d'Oligos\n", "2", "i"))
     return ;
   freeint (&nbolig) ; */
/* Warning nnt = 5 pour tester 5 et 7 */
  ntt = 4 ;
  while (ntt--)
    { taille = tt[ntt] ;
      i = 1 ; /* i = 3 to test 2, 5 or 10 oligo per sequences */
      mask = (1 << (2*taille)) ;
      while (i--)
	{ nbolig = ii[i] ;
	  adna = assBigCreate (5000) ;
	  look = defCptGetLook (1) ;
	  look->def = queryKey (keySet(aa, 0), ">Assembled_from ; >DNA") ;
	  if (!look->def || !keySetMax (look->def))
	    goto abort ;
	  look->mar = keySetCreate () ;
	  dnaAlignMakeSpecialMotif (look, adna, nbolig, taille) ;
	  ks = keySetCreate () ;
	  dna = dnaGet (keySet (aa, 0)) ;
	  cp = arrp (dna, 0, char) - 1 ;
	  max = arrayMax (dna) ;
	  while (cp++, max--)
	    {
	      switch (*cp)
		{
		case A_: case T_: case G_: case C_:
		  test = ((test << 2) | u[(int)*cp]) ;
		  break ;
		default:
		  depart = FALSE ;
		  test = 1 ;
		  continue ;
		}
	      if (!depart)
		{ if (test & (mask >> 2))
		    depart = TRUE ;
		  continue ;
		}
	      test &= (mask - 1) ;
	      if (!assFind (adna, vp0 + test, &vp))
		continue ;
	      key = vp - vp0 ;
	      if (!keySetFind (look->mar, key, &j))
		messcrash ("Problem in keySetFind") ;
	      keySet (ks, j)++ ;
	    }
	  assDestroy (adna) ;
	  arrayDestroy (dna) ;
	  result = arrayCreate (100, int) ;
	  for (j = 0 ; j < keySetMax (ks) ; j++)
	    array (result, keySet (ks, j), int)++ ;
	  keySetDestroy (ks) ;
/*	  fprintf (stdout, "************ Bilan *********\n") ;
	  fprintf (stdout, "Nombre d'oligos : %d\n", nbolig) ;
	  fprintf (stdout, "Taille des oligos : %d\n", taille) ;
	  fprintf (stdout, "Nombre de sequences : %d\n", keySetMax (look->def)) ;
	  for (j = 0 ; j < arrayMax (result) ; j++)
	    fprintf (stdout, "%d oligo trouves %d fois\n", array (result, j, int), j) ;
	  fprintf (stdout, "****************************\n") ; */
	  fprintf (stdout, "\n\n") ;
	  fprintf (stdout, "Nombre d'oligos\t%d\n", nbolig) ;
	  fprintf (stdout, "Taille des oligos\t%d\n", taille) ;
	  fprintf (stdout, "Nombre de sequences\t%d\n", keySetMax (look->def)) ;
	  fprintf (stdout, "Nb fois\tNb oligo\n") ;
	  for (j = 0 ; j < arrayMax (result) ; j++)
	    fprintf (stdout, "%d\t%d\n", j, array (result, j, int)) ;
	  fprintf (stdout, "****************************\n") ;
/*#ifndef NON_GRAPHIC
          plotHisto ("Nombre d'oligo trouves n fois", result) ;
#else
	  arrayDestroy (result) ;
#endif */
	  arrayDestroy (result) ;
	abort:
	  defCptDestroyLook (look->link) ;
	}
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

