/*  File: boubou.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2004
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the COMPARATOR genome database package, written by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This code works as a client against an acedb server
 *  and should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  The present code and the acedb schema for this project
 *  are available from http://www.aceview.org/aligold
 */

#include "../wac/ac.h"
#include "colours.h"
#include "keyset.h"
#include <errno.h>
#include "bitset.h"
#include "freeout.h"
#include "aceio.h"
#include "vtxt.h"
#include "dict.h"
#include "aceio.h"

#define NFF 128
static int FF[NFF], FF2[NFF] ;
static char*  FFM[NFF] ;
static void usage (void) ;
typedef struct mmStruct {
  AC_DB db ; 
  BOOL bbNmers, G, GC, GA, GT,AT, AC, TC ;
  BOOL manip1, manip2, manip3, green, red, ratio, G1, pAplus, pAmoins, S, expression, single, maqc, histo, histoN, pivot, filter ;
  BOOL getSolexaPeaks, solexaAutoCorrel ;
  char *manip ;
  const char *solexaSitesFileName ;
  DICT *dict ;
  Array mg ;
  BitSet filtered  ;
  DICT *probeDict ;
  int NN ;
} MM ;

typedef struct qqStruct { int nam, nn ; double dz ; } QQ ;

#include "../wabi/boubou.h"
 
/*************************************************************************************/

static int qqOrder (const void *va, const void *vb)
{
  QQ *a = (QQ*)va, *b = (QQ*)vb ;
  double za, zb ;
  za = a->dz/(a->nn ? a->nn : 1) ;
  zb = b->dz/(b->nn ? b->nn : 1) ;
  if (za < zb) return 1 ; /* large first */
  if (za > zb) return -1 ;
  return a->nam - b->nam ;
} 

/*************************************************************************************/

static void bbNmers (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN fi = aceInCreateFromStdin (FALSE, 0, h) ;
  int NN = mm->NN ; /* length of motif beeing analysed */
  char mybuf[200], Motif[NN+1], *cp, *cq, *cp1, *motif ;
  const char *snam1 ;
  int i, ii, jj, nn, ns ;
  double z1, z2, dz, lg2 = log((double)2.0) ;
  float s1, s2 ; 
  Array aa = arrayHandleCreate (10000, QQ, h) ;
  Array Ba = arrayHandleCreate (100, QQ, h) ;
  Array Bt = arrayHandleCreate (100, QQ, h) ;
  Array Bg = arrayHandleCreate (100, QQ, h) ;
  Array Bc = arrayHandleCreate (100, QQ, h) ;
  Array Bcg = arrayHandleCreate (1100, QQ, h) ;
  Array Sa = arrayHandleCreate (100, QQ, h) ;
  Array St = arrayHandleCreate (100, QQ, h) ;
  Array Sg = arrayHandleCreate (100, QQ, h) ;
  Array Sc = arrayHandleCreate (100, QQ, h) ;

  PNS *pns ;
  QQ *qq ;
  int nA, nT, nG, nC ;
  BOOL a2n, t2n, c2n, g2n ;
  DICT *dict = dictHandleCreate (10000, h) ;

  a2n = t2n = g2n = c2n = FALSE ;
  if (mm->G) {a2n = t2n = c2n = TRUE ; } 
  if (mm->GC) {a2n = t2n = TRUE ; } 
  if (mm->GA) {c2n = t2n = TRUE ; }
  if (mm->GT) {c2n = a2n = TRUE ; }
  if (mm->AT) {g2n = c2n = TRUE ; }
  if (mm->AC) {g2n = t2n = TRUE ; }
  if (mm->TC) {g2n = a2n = TRUE ; }
  while (aceInCard (fi))
    {
      const char *ccp = aceInWord(fi) ;
      if (!ccp) continue ; /* first probe name */
      if (mm->manip3 && strncmp(ccp,"BBC",3))
	continue ;
      if (mm->manip2 && strncmp(ccp,"BBB",3))
	continue ;
      if (mm->manip1 && strncmp(ccp,"BBB",3))
	continue ;
      if (!mm->single && !aceInWord(fi)) continue ; /* second probe name */
      if (!(snam1 = aceInWord(fi))) continue ;
      strcpy (mybuf, snam1) ;
      cp1 = strstr (mybuf, ".pair") ;
      if (cp1) *cp1 = 0 ;
      for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
	if (!strcmp (pns->nam, mybuf) || !strcmp (pns->p, mybuf)) break ;
      if (!pns->p)
	continue ;
      snam1 = pns->p ; 
      strcpy (mybuf, snam1) ;
      if (!strstr(pns->nam, "5MS")) continue ;
      if (mm->G1 && !strstr (pns->nam, "G1"))
	continue ;
      if (mm->expression && !pns->isRna)
	continue ;
      if (mm->S && *pns->nam != 'S' && !strstr (pns->nam, "S=") && !strstr (pns->nam, "Samplified="))
	continue ;
      if (mm->pAplus && !strstr (pns->nam, "pA+"))
	continue ;
      if (mm->manip && !strstr (snam1, mm->manip))
	continue ;
      if (mm->pAmoins && !strstr (pns->nam, "NucPA"))
	continue ;
      if (mm->green && !strstr (snam1,"_532")) /* green dye */
	continue ;
      if (mm->red && !strstr (snam1,"_635"))   /* red dye */
	continue ;

      if (mm->manip3) ; /* already tested */
      else
	{
	  if (!mm->manip && !mm->manip3 && !strncmp(snam1,"144",3)) /* manip nouvelle */
	    continue ;
	  if (mm->manip1 && !strstr (pns->nam, "=1"))   /* manip 1 */
	    continue ;
	  if (!mm->manip && mm->manip2 && strstr (pns->nam, "=1"))   /* manip 1 */
	    continue ;
	}

      if (!aceInFloat (fi, &s1)) continue ;
      if (!mm->single)
	{
	  if (!aceInWord(fi)) continue ; /* second signal name == signal name1 */
	  if (!aceInFloat (fi, &s2)) continue ;
	}
      if (!mm->maqc)
	{
	  if (strstr(pns->nam,"AFX"))
	    {
	      if (s1 < 4 || s1 > 30000) continue ;
	    }
	  else
	    {
	      if (s1 < 100 || s1 > 100000) continue ;
	    }
	  z1 = log((double)s1 + pns->damper)/lg2 ; /* 256->too big, 8->too small */
	}
      else
	z1 = s1 ;
      if (!mm->single)
	{
	  if (strstr(pns->nam,"AFX"))
	    {
	      if (s2 < 4 || s2 > 30000) continue ;
	    }
	  else
	    {
	      if (s2 < 100 || s2 > 100000) continue ;
	    }
	  z2 = log((double)s2 + pns->damper)/lg2 ; 
	}
      else
	z2 = 0 ;
      dz = z2 - z1 ;

      for (jj = 0 ; jj < 2 ; jj++)
	{
	  motif = aceInWord(fi) ;
	  if (!motif) continue ;
	  if (strlen(motif) != 45) 
	    continue ;
	  /*
	  nn = strlen (motif) - NN + 1 ;
	  */
	  nn = 1 ;
	  Motif[NN] = 0 ;
	  for (i = 0, cp = motif ; i < nn ; cp++, i++)
	    {
	      strncpy (Motif, cp, NN) ;
	      nA = nT = nG = nC = 0 ;
	      for (cq = Motif ; *cq ; cq++)
		{
		  switch (*cq)
		    {
		    case 'A':
		      nA++ ; 
		      if (a2n) *cq = 'n' ;
		      break ;
		    case 'T':
		      nT++ ; 
		      if (t2n) *cq = 'n' ;
		      break ;
		    case 'G':
		      nG++ ; 
		      if (g2n) *cq = 'n' ;
		      break ;
		    case 'C':
		      nC++ ; 
		      if (c2n) *cq = 'n' ;
		      break ;
		    }
		}
	      if (NN <= 16)
		{
		  dictAdd (dict, Motif, &ii) ;
		  qq = arrayp (aa, ii, QQ) ;
		  qq->nam = ii ;
		  qq->dz += jj ? dz : -dz ;
		  qq->nn += 1 ;
		}
	
	      dictAdd (dict, messprintf("%d-A", nA), &ii) ;
	      qq = arrayp (Ba, nA, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? dz : -dz ;
	      qq->nn += 1 ;

	      qq = arrayp (Sa, nA, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? z2 : z1 ;
	      qq->nn += 1 ;
	      
	      dictAdd (dict, messprintf("%d-T", nT), &ii) ;
	      qq = arrayp (Bt, nT, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? dz : -dz ;
	      qq->nn += 1 ;
	      
	      qq = arrayp (St, nT, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? z2 : z1 ;
	      qq->nn += 1 ;

	      dictAdd (dict, messprintf("%d-G", nG), &ii) ;
	      qq = arrayp (Bg, nG, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? dz : -dz ;
	      qq->nn += 1 ;
	      
	      qq = arrayp (Sg, nG, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? z2 : z1 ;
	      qq->nn += 1 ;

	      dictAdd (dict, messprintf("%d-C", nC), &ii) ;
	      qq = arrayp (Bc, nC, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? dz : -dz ;
	      qq->nn += 1 ;
	      
	      qq = arrayp (Sc, nC, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? z2 : z1 ;
	      qq->nn += 1 ;

	      dictAdd (dict, messprintf("%d", nC - nG), &ii) ;
	      qq = arrayp (Bcg, 1000 + nC - nG, QQ) ;
	      qq->nam = ii ;
	      qq->dz += jj ? dz : -dz ;
	      qq->nn += 1 ;
	    }
	}	
    }

  /* sort */
  arraySort (aa, qqOrder) ;

  arraySort (Ba, qqOrder) ;
  arraySort (Bt, qqOrder) ;
  arraySort (Bg, qqOrder) ;
  arraySort (Bc, qqOrder) ;

  /*
    arraySort (Sa, qqOrder) ;
    arraySort (St, qqOrder) ;
    arraySort (Sg, qqOrder) ;
    arraySort (Sc, qqOrder) ;
  */
  /* report */
  if (! mm->single)
    for (jj = 0 ; jj < 5 ; jj++)
      {
	Array a ;
	switch (jj)
	  {
	  case 0: a=aa ; break ;
	  case 1: 
	    freeOut ("\nProbe/antiprobe ratios (log2)\n") ;
	    a=Ba ; break ;
	  case 2: a=Bt ; break ;
	  case 3: a=Bg ; break ;
	  case 4: a=Bc ; break ;
	  }
	nn = arrayMax (a) ;
	for (ii = 0 ; ii <= nn/2 && ii < 30 ; ii++)
	  {
	    qq = arrp (a, ii, QQ) ; 
	    if (!qq->nn)
	      continue ;
	    freeOutf ("%s\t%d\t%g", dictName (dict, qq->nam),qq->nn, qq->dz/qq->nn) ;
	    qq = arrp (a, nn - ii - 1, QQ) ; 
	    if (qq->nn)
	      freeOutf ("\t\t%s\t%d\t%g", dictName (dict, qq->nam),qq->nn, qq->dz/qq->nn) ;
	    freeOut ("\n") ;
	  } 
	freeOut ("\n") ;
      }

  if (1)
    {
      freeOut ("\nlog2(Probe signal), for each kind of letter in the probe\n") ;
      nn = NN ;
      freeOutf ("#\tA\tT\tG\tC\tnA\tnT\tnG\tnC\n") ;
      for (ii = 0 ; ii <= NN  ; ii++)
	{
	  freeOutf ("%d", ii) ;

	  qq = arrp (Sa, ii, QQ) ; 
	  freeOutf ("\t%g", qq->nn ? qq->dz/qq->nn : 0) ;
	  qq = arrp (St, ii, QQ) ; 
	  freeOutf ("\t%g", qq->nn ? qq->dz/qq->nn : 0) ;
	  qq = arrp (Sg, ii, QQ) ; 
	  freeOutf ("\t%g", qq->nn ? qq->dz/qq->nn : 0) ;
	  qq = arrp (Sc, ii, QQ) ; 
	  freeOutf ("\t%g", qq->nn ? qq->dz/qq->nn : 0) ;

	  qq = arrp (Sa, ii, QQ) ; 
	  freeOutf ("\t%d", qq->nn) ;
	  qq = arrp (St, ii, QQ) ; 
	  freeOutf ("\t%d", qq->nn) ;
	  qq = arrp (Sg, ii, QQ) ; 
	  freeOutf ("\t%d", qq->nn) ;
	  qq = arrp (Sc, ii, QQ) ; 
	  freeOutf ("\t%d", qq->nn) ;

	  freeOutf ("\n") ;
	} 
      freeOut ("\n") ;
    }
  if (arrayMax (Bcg))
    {
      Array a = Bcg ;
      freeOutf("\n\nlog(signal(Probe)/signal(anti probe)) as a function of #C minus #G in the probe\n") ;
      
      for (ii = 0 ; ii < arrayMax(Bcg) ; ii++)
	{
	  qq = arrp (a, ii, QQ) ; 
	  if (qq->nn)
	    break ;
	}
      for (jj = arrayMax (Bcg) ; jj >= 0 ; jj--)
	{
	  qq = arrp (a, jj, QQ) ; 
	  if (qq->nn)
	    break ;
	}
      for (; ii <= jj ; ii++)
	{
	  qq = arrp (a, ii, QQ) ; 
	  if (qq->nn)
	    freeOutf ("%d\t%d\t%g\n", ii - 1000 ,qq->nn, qq->dz/qq->nn) ;
	  else
	    freeOutf ("%d\t%d\t\n", ii - 1000 ,qq->nn) ;
	}
    }

  ac_free (h) ;
} /* bbNmers */

/*************************************************************************************/
/***************************** Histos ************************************************/
/*************************************************************************************/
typedef struct p2sMotifStruct { float tm, zdna ; BOOL filter ; int probe, motif, A,T,G,C,CG,GC,AA,TT,CC, GG, ln ; float ra,rt,rg,rc,rcg ;} P2SM ;
typedef struct p2sStruct { double s ; P2SM m ; int ns ; } P2S ;

/*************************************************************************************/

static int p2sOrder (const void *va, const void *vb)
{
  P2S *za = (P2S*)va, *zb = (P2S*)vb ;
  if (za->s < zb->s) return -1 ;
  if (za->s > zb->s) return 1 ;
  return 0 ;
} 

/*************************************************************************************/

static int htileMotifFilter (PNS *pns, P2S *p2s, char *pNam)
{
  int nnf = 0 ;

  p2s->m.filter = FALSE ;
  if (1) /* filtre recursif Danielle 2007_11_30 */
    {
      /* filtre 1 */
      if (p2s->m.rt + p2s->m.rg <= 18)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r t+g <= 18" ;
      if (p2s->m.rt + p2s->m.rg > 84)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r t+g > 84" ;

      if (p2s->m.ln <= 28)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "len<= 28" ;

      /* filtre 2 */
      if (p2s->m.zdna > 7)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "Z-DNA > 7" ;

      if (p2s->m.G - p2s->m.C <= -24) 
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#G - #C <= 24" ;
      if (p2s->m.CG > 16)
	{ if (!p2s->m.filter) FF2[nnf]++ ;  p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#CpG > 16" ;
      if (p2s->m.GC > 18)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#GpC > 18" ;

      /* filtre 3 */
      if (p2s->m.zdna > 4)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "Z-DNA > 4" ;
      if (p2s->m.G == 0)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#G == 0" ;

      if (p2s->m.ra + p2s->m.rc <= 18 || p2s->m.ra + p2s->m.rc > 78)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a+c > 78" ;
      if (p2s->m.rt + p2s->m.rg <= 24)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r t+g <= 24" ;
      if (p2s->m.CG > 14)
	{ if (!p2s->m.filter) FF2[nnf]++ ;  p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#CpG > 14" ;

      /* filtre 4 */
      if (p2s->m.rg - p2s->m.rc <= -48)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r g-c <= -48" ;
      if (p2s->m.C == 0 || p2s->m.C > 24)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#C== 0 || #C > 24" ;
      if (p2s->m.rc > 54)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r c > 54" ;

      /* filtre 5 */
      if (p2s->m.rg <= 6)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r g <= 6" ;
      if (p2s->m.ra + p2s->m.rc > 72)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a+c > 72" ;

      /* filtre 6 */
      if (p2s->m.tm > 84)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "Tm > 84" ;
      if (p2s->m.ln <= 31)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "ln <= 31" ;
      if (p2s->m.ra + p2s->m.rt <= 18)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a+t <= 18" ;
      if (p2s->m.ra + p2s->m.rc <= 24)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a+c <= 24" ;
      if (p2s->m.rt + p2s->m.rg > 76)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r t+g > 76" ;

      /* filtre 7 */
      if (p2s->m.ra + p2s->m.rg <= 16)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a+g <= 16" ;
      if (p2s->m.ra - p2s->m.rc <= -48)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a-c <= -48" ;
      if (p2s->m.ra - p2s->m.rg > 32)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r a-g > 32" ;
      if (p2s->m.rg + p2s->m.rc > 78)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r g+c > 78" ;
      if (p2s->m.GC > 16)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#GpC > 16" ;
      if (p2s->m.CG > 12)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#CpG > 12" ;

      /* filtre 8 */
      if (p2s->m.G - p2s->m.C <= -18)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#G - #C <= -18" ;
      if (p2s->m.A - p2s->m.T > 25)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "#A - #T > 25" ;
      if (p2s->m.rt + p2s->m.rc > 82)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r t+c > 82" ;
      if (p2s->m.rg - p2s->m.rc <= -30)
	{ if (!p2s->m.filter) FF2[nnf]++ ; p2s->m.filter = TRUE ; FF[nnf]++ ; }
      FFM[nnf++] = "r g-c <= -30" ;
    }
  return p2s->m.filter ;
} /* htileMotifFilter */

/*************************************************************************************/
/* compute the sequence data associated to this motif */
static void p2sMotifData (MM *mm, PNS *pns, P2S *p2s, const char *motif, char *pNam)
{
  const char *ccp ;
  int ln = 0, p=-1 ;
  int i, naa, ntt, ncc, ngg ;
  p2s->m.A = p2s->m.T = p2s->m.G = p2s->m.C = p2s->m.CG = p2s->m.GC =
    p2s->m.AA = p2s->m.TT = p2s->m.GG = p2s->m.CC = p2s->m.ln = p2s->m.zdna = 0 ;
  naa = ntt = ncc = ngg = 0 ;
  p2s->m.ra = p2s->m.rt = p2s->m.rg = p2s->m.rc = p2s->m.rcg = 0 ;
  for (ccp = motif ; *ccp ; ccp++)
    {
      switch (*ccp)
	{
	case 'A':
	  p2s->m.A++ ;
	  if(ccp > motif && *(ccp-1)!='A')
	    {
	      for (i = 1; i<12 ; i++)
		if (*(ccp+i) != 'A')
		  break ;
	      if (i > naa) naa = i ; 
	    }
	  break ;
	case 'T': 
	  p2s->m.T++ ;
	  if(ccp > motif && *(ccp-1)!='T')
	    {
	      for (i = 1; i<12 ; i++)
		if (*(ccp+i) != 'T')
		  break ;
	      if (i > ntt) ntt = i ; 
	    }
	  break ;
	case 'G': 
	  p2s->m.G++ ; 
	  if(ccp > motif && *(ccp-1)=='C') p2s->m.CG += 2 ; 
	  if(ccp > motif && *(ccp-1)!='G')
	    {
	      for (i = 1; i<12 ; i++)
		if (*(ccp+i) != 'G')
		  break ;
	      if (i > ngg) ngg = i ; 
	    }
	  break ;
	case 'C': 
	  p2s->m.C++ ; 
	  if(ccp > motif && *(ccp-1)=='G') p2s->m.GC += 2 ; 
	  if(ccp > motif && *(ccp-1)!='C')
	    {
	      for (i = 1; i<12 ; i++)
		if (*(ccp+i) != 'C')
		  break ;
	      if (i > ncc) ncc = i ; 
	    }
	  break ;
	}
      switch (*ccp)
	{
	case 'A': p = 0 ; break ;
	case 'G': if(p>0 && (p & 0x1)) {p++ ; if (p>8) p2s->m.zdna++;} else p = 0 ; break ;
	case 'T': case 'C': if(p>-1 && !(p & 0x1)) {p++ ; if (p>8) p2s->m.zdna++;} else p = 1 ; break ;
	}
      ln++ ;
    }
  if (ln)
    {
      p2s->m.ra = p2s->m.A * (100.0 / (float)ln) ;
      p2s->m.rt = p2s->m.T * (100.0 / (float)ln) ;
      p2s->m.rg = p2s->m.G * (100.0 / (float)ln) ;
      p2s->m.rc = p2s->m.C * (100.0 / (float)ln) ;
      p2s->m.rcg = p2s->m.CG * (200.0 / (float)ln) ;
      p2s->m.ln = ln ;
      p2s->m.CC = ncc ;
      p2s->m.GG = ngg ;
      p2s->m.AA = naa ;
      p2s->m.TT = ntt ;
    }
  if (1)
    {
      int nn ;
      static int n2=0,n3=0 ;
      
      if (dictAdd (mm->probeDict, motif, &nn))
	{
	  if (0 &&  p2s->m.C > 29)
	    { fprintf (stderr, "%d\t%d\t%s\n",  ++n2, p2s->m.C, motif) ;}
	  if (mm->filter && htileMotifFilter (pns, p2s, pNam))
	    {
	      bitSet (mm->filtered, nn) ;
	      FF[NFF-1]++ ;
	    }
	}
      p2s->m.motif = nn ;
      if (0 && p2s->m.C > 29 && strstr(pns->nam, "G1="))
	{ fprintf (stderr, "\t%d\t%d\t%s\n",  ++n3, p2s->m.C, motif) ;}
      if (bit (mm->filtered, nn))
	p2s->m.filter = TRUE ;	
    }
} /* p2sMotifData */
  
/*************************************************************************************/
/* remove signal <= -900 */
static void p2sAllChipsCompressSignal (MM *mm)
{
  int i, j, n, ns ;
  PNS *pns ;
  P2S *p2s, *p2s2 ;

  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    {
      if (pns->signal && (n = arrayMax(pns->signal)))
	{
	   p2s = p2s2 = arrp (pns->signal, 0, P2S) ;
	   for (i = j = 0 ; i < n ; p2s++, i++)
	    {
	      if (p2s->s < -900)
		continue ;
	      if (j < i) { *p2s2 = *p2s ; }
	      j++ ; p2s2++ ;
	    }
	  arrayMax(pns->signal) = j ;
	}
    }

  return ;
} /* p2sAllChipsCompressSignal */

/*************************************************************************************/
/* Export Histograms of signals for all chips in a single table */
static double p2sDamper (MM *mm, PNS *pns, double *zinflexp, double *bestSlopep)
{
  int i, j, nn ;
  double NN = 20 ;
  double slope, damper, bestDamper = 0, delta, bestDelta = 100000, z1, lg2 = log(2) ;
  P2S *p2s ;
  static Array histo = 0 ;
  int j1 = 0, j2 = 0, h1, h2 ;

  *bestSlopep = 0 ;
  for (damper = 1 ; damper < 4000 ; damper *= sqrt(sqrt(2.0)))
    {
      double zm = 0 ;
      int n, iMax = arrayMax(pns->signal), nMax = -1000 ; 

      histo = arrayReCreate (histo, 40, int) ;
      for (i = nn = 0 ; i < iMax ; i++)
	{
	  p2s = arrp (pns->signal, i, P2S) ;
	  if (0) /* do not normalize */ 
	    { 
	      pns->damper = 1; pns->av = 10 ; 
	    }
	  z1 = exp(lg2 * (p2s->s + pns->av - 10.0)) - pns->damper + damper ;
	  if (z1 <= 0)
	    continue ;
	  z1 = log (z1)/lg2 - pns->av + 10.0 ;
	  if (z1 < 0)
	    continue ;
	  if (damper == pns->damper && (z1 > p2s->s + .1 || z1 < p2s->s - .1))
	    messcrash("error terrible") ;
	  zm += z1 ; nn++ ;
	  j = z1 * NN ;
	  array(histo, j, int)++ ;
	}
      /* locate the inflexion point, measure the slope */
      j1 = j2 = h1 = h2 = 0 ;
      slope = 0 ;
      
      for (j = 7*NN ; j <= 16*NN ; j++)
	{
	  n = array(histo,j,int) ;
	  if (n > nMax) { nMax = n ; }
	}
      for (j = 7*NN ; j <= 16*NN ; j++)
	{
	  n = array(histo,j,int) ;
	  if (!j1 && n > nMax/4) { h1 = n ; j1 = j ;}
	  if (j1 && n > 3*nMax/4) { h2 = n ; j2 = j ; break ; }
	}
      if (j2 > j1) /* slope at inflexion point */
	slope = ((float)h2 - h1)/(j2 - j1) ;	      
      delta = slope - (pns->isRna ? 1200.0 : 600.0) ;
      if (delta < 0) delta = -delta ;
      if (delta < bestDelta)
	{ bestDelta = delta ; bestDamper = damper ; *bestSlopep = slope ; *zinflexp = ((double)j1 + j2)/(2.0 * NN) ; }
    }
  return bestDamper ;
} /* p2sDamper */

/*************************************************************************************/
/* Export Histograms of signals for all chips in a single table */
static void p2sAllChipsSignalHisto (MM *mm)
{
  int i, j, jmax, ns, nn, n, nMax ;
  double NN = 20 ;
  float z01, z5, z20, z25, z40, z60, z75, z80, z95, z99, zmedian, ztop, zInflex, slope ;
  double z1, bestDamper, bestSlope = 0, bestZinflex = 0 ;
  PNS *pns ;
  P2S *p2s ;
  static Array histo = 0 ;
  int j1 = 0, j2 = 0, h1, h2 ;

  freeOutf ("Manip\tmax\t99%%\t95%%\t80%%\t75%%\t60%%\tmean\tmedian\t40%%\t25%%\t20%%\t5%%\t1%%\tmin\ttop\tinflex\t# signal\tManip") ;
  for (j = 7*NN ; j <= 16*NN ; j++)
    freeOutf("\t%.2f", j/NN) ;
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    {
      if (
	  (!mm->ratio || !strstr (pns->nam2, "ratio")) &&
	  pns->signal)
	{
	  double zm = 0, zmin = 10000, zmax = -10000 ;
	  int iMax = arrayMax(pns->signal) ;
	   	
	  bestDamper = p2sDamper (mm, pns, &bestZinflex, &bestSlope) ; 
	  histo = arrayReCreate (histo, 40, int) ;
	  for (i = nn = 0 ; i < iMax ; i++)
	    {
	      p2s = arrp (pns->signal, i, P2S) ;
	      z1 = p2s->s ;
	      /* FORCE centering on inflexion point z1 = z1 + 9.5 - bestZinflex ; */
	      if (z1 < 0)
		continue ;
	      if (0 && p2s->m.tm > 72) continue ;
	      zm += z1 ; nn++ ;
	      if (z1 > zmax) zmax = z1 ;
	      if (z1 < zmin) zmin = z1 ;
	      j = z1 * NN ;
	      array(histo, j, int)++ ;
	    }
	  zm /= nn ;
	  arraySort (pns->signal, p2sOrder) ;
	  i = 1*nn ; i /= 100 ;z01 = arr (pns->signal, i, P2S).s ;
	  i = 5*nn ; i /= 100 ;z5 = arr (pns->signal, i, P2S).s ;
	  i = 20*nn ; i /= 100 ;z20 = arr (pns->signal, i, P2S).s ;
	  i = 25*nn ; i /= 100 ;z25 = arr (pns->signal, i, P2S).s ;
	  i = 40*nn ; i /= 100 ;z40 = arr (pns->signal, i, P2S).s ;
	  i = 50*nn ; i /= 100 ;zmedian = arr (pns->signal, i, P2S).s ;
	  i = 60*nn ; i /= 100 ;z60 = arr (pns->signal, i, P2S).s ;
	  i = 75*nn ; i /= 100 ;z75 = arr (pns->signal, i, P2S).s ;
	  i = 80*nn ; i /= 100 ;z80 = arr (pns->signal, i, P2S).s ;
	  i = 95*nn ; i /= 100 ;z95 = arr (pns->signal, i, P2S).s ;
	  i = 99*nn ; i /= 100 ;z99 = arr (pns->signal, i, P2S).s ;
 
  	  zInflex = ztop = slope = 0 ; j1 = j2 = jmax = nMax = 0 ;
	  if (1) /* find the max */
	    {
	      for (j = 7*NN ; j <= 16*NN ; j++)
		{
		  n = array(histo,j,int) ;
		  if (n > nMax) { nMax = n ; jmax = j ; ztop = (float)jmax/NN ; }
		}
	    }
	  h1 = h2 = 0 ; slope = zInflex = 0 ;
	  if (1) /* locate the inflexion point, measure the slope */
	    {
	      for (j = 7*NN ; j <= 16*NN ; j++)
		{
		  n = array(histo,j,int) ;
		  if (!j1 && n > nMax/4) { h1 = n ; j1 = j ;}
		  if (j1 && n > 3*nMax/4) { h2 = n ; j2 = j ; break ; }
		}
	      if (j2 > j1) 
		{
		  slope = ((float)h2 - h1)/(j2 - j1) ;	      
		  zInflex = (float)(j1+j2)/(2.0 * NN) ;
		}
	    }
	  freeOutf ("\n%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%s"
	            , pns->nam, zmax, z99, z95, z80, z75, z60, zm, zmedian, z40, z25, z20, z5, z01, zmin, ztop, zInflex, nn, pns->nam) ;
	  if (1)	      /* export the histo */
	    {
	      for (j = 7*NN ; j <= 16*NN ; j++)
		{
		  n = array(histo,j,int) ;
		  freeOutf ("\t%d", n) ;
		}
	    }
	  freeOutf("\t%.2f\tSlope\t%g\t%g\tDAMPER\t%g\t%g\tInflexion\t%g\t%g\tAv\t%g\t%g"
		   , jmax/NN,slope,bestSlope, pns->damper, bestDamper, zInflex, bestZinflex, pns->av, pns->av + bestZinflex - zInflex) ;
	}
    }
} /* p2sAllChipsSignalHisto */

/*************************************************************************************/

/* Export Histograms of signals for all chips in a single table */
static void p2sSingleChipManyHistos (MM *mm, PNS *pns, int ns)
{
  int NN = 13 ; /* number of points in each histo */
  int i, j, nf, np, nHisto ;
  int nMax = arrayMax(pns->signal) ;
  double x, x1, step, z1, z2, sigma ;
  static Array histo = 0, val1 = 0, val2 = 0 ;
  P2S *p2s ;
  typedef enum { TM, LN, GpC, GmC, ApT, AmT, ZDNA, RApT, RApG, RApC, RTpG, RTpC, RGpC, 
		 RAmT, RAmG, RAmC, RTmG, RTmC, RGmC, 
		 NCG, NGC, NGG, NCC, NAA, NTT,
		 RA, RT, RG, RC, NA, NT, NG, NC, ZERO } HTP ;
  typedef struct typeStruct { HTP type ; const char *nam ; int min, max ;} HTYPE ;
  HTYPE *tt ;
  HTYPE histoTypes [] = {
    { NC, "#C", 0, 36},
    { NG, "#G", 0, 36},
    { NA, "#A", 0, 48},
    { NT, "#T", 0, 48},
    { ApT, "#A+T", 0, 72 },
    { GpC, "#G+C", 16, 40 },
    { AmT, "#A-T", -36, 36 },
    { GmC, "#G-C", -36, 36 },
    { NCG, "CpG", 0, 24 },
    { NGC, "GpC", 0, 24 },
    { NCC, "polyC", 0, 12 },
    { NGG, "polyG", 0, 12 },
    { NAA, "polyA", 0, 12 },
    { NTT, "polyT", 0, 12 },

    { ZDNA, "ZDNA", 0, 24 },
    { TM, "TM", 66, 90 },
    { LN, "Length", 22, 94 },

    { RA, "rate_A", 0, 72},
    { RT, "rate_T", 0, 72},
    { RG, "rate_G", 0, 72},
    { RC, "rate_C", 0, 72},

    { RApG, "rate_A+G", 16, 88 },
    { RTpC, "rate_T+C", 16, 88 },
    { RApC, "rate_A+C", 18, 90 },
    { RApT, "rate_A+T", 18, 90 },
    { RTpG, "rate_T+G", 18, 90 },
    { RGpC, "rate_G+C", 34, 82 },
 
    { RAmC, "rate_A-C", -48, 48 },
    { RAmG, "rate_A-G", -48, 48 },
    { RAmT, "rate_A-T", -48, 48 },
    { RGmC, "rate_G-C", -48, 48 },
    { RTmC, "rate_T-C", -48, 48 },
    { RTmG, "rate_T-G", -48, 48 },
    /*  */
    { ZERO, 0, 0, 0 }} ;

  for (nHisto=1, tt = histoTypes ; tt->nam ; nHisto++, tt++)
    {
      histo = arrayReCreate (histo, 40, int) ;
      val1 = arrayReCreate (val1, 40, double) ;
      val2 = arrayReCreate (val2, 40, double) ;
      step = ((double)tt->max - tt->min)/(NN-1) ;
      freeOutf ("\nG1=00\t%03d.%s", nHisto, tt->nam) ;
      for (j = 0 ; j <= NN ; j++)
	freeOutf ("\t%g", tt->min + j*step) ; /* bin value */  
      freeOutf ("\n") ;
      /* accumulate the data */
      for (i = nf = np = 0, p2s = arrp (pns->signal, i, P2S) ; i < nMax ; i++, p2s++)
	{  
	  np++ ;
	  if (mm->filter && p2s->m.filter)
	    { nf++; continue ; }
	  
	  /* choose the value x to be plotted */
	 x = 0 ;
	 switch (tt->type)
	    {
	    case TM: x = p2s->m.tm ; break ; 
	    case ZDNA: x = p2s->m.zdna ; break ; 
	    case LN: x = p2s->m.ln ; break ; 
	    case NA: x = p2s->m.A ; break ;
	    case NT: x = p2s->m.T ; break ;
	    case NG: x = p2s->m.G ; break ;
	    case NC: x = p2s->m.C ; break ;
	    case GpC: x = p2s->m.G + p2s->m.C ; break ; 
	    case GmC: x = p2s->m.G - p2s->m.C ; break ; 
	    case ApT: x = p2s->m.A + p2s->m.T ; break ; 
	    case AmT: x = p2s->m.A - p2s->m.T ; break ; 
	    case RA: x = p2s->m.ra ; break ;
	    case RT: x = p2s->m.rt ; break ;
	    case RG: x = p2s->m.rg ; break ;
	    case RC: x = p2s->m.rc ; break ;

	    case RApT: x = p2s->m.ra + p2s->m.rt ; break ;
	    case RApG: x = p2s->m.ra + p2s->m.rg ; break ;
	    case RApC: x = p2s->m.ra + p2s->m.rc ; break ;
	    case RTpG: x = p2s->m.rt + p2s->m.rg ; break ;
	    case RTpC: x = p2s->m.rt + p2s->m.rc ; break ;
	    case RGpC: x = p2s->m.rg + p2s->m.rc ; break ;

	    case RAmT: x = p2s->m.ra - p2s->m.rt ; break ;
	    case RAmG: x = p2s->m.ra - p2s->m.rg ; break ;
	    case RAmC: x = p2s->m.ra - p2s->m.rc ; break ;
	    case RTmG: x = p2s->m.rt - p2s->m.rg ; break ;
	    case RTmC: x = p2s->m.rt - p2s->m.rc ; break ;
	    case RGmC: x = p2s->m.rg - p2s->m.rc ; break ;

	    case NCG: x = p2s->m.CG ; break ;
	    case NGC: x = p2s->m.GC ; break ;
	    case NCC: x = p2s->m.CC ; break ;
	    case NGG: x = p2s->m.GG ; break ;
	    case NAA: x = p2s->m.AA ; break ;
	    case NTT: x = p2s->m.TT ; break ;
	    default: break ;
	    }
	  /* find the value j of the correct bin */
	  for (j = 0, x1 = tt->min ; x1 < tt->max ; j++, x1 += step)
	    if (x <= x1)
	      break ;
	  array (histo, j, int)++ ;
	  if (0 && j >= 10)
	    {
	      static int nnmine = 0 ;
	      fprintf (stderr, "########### %d\tx=%g\tj=%d\thh=%d\t%s\t%s\t%s\t%g\n"
		       , ++nnmine, x, j, array (histo, j, int), pns->nam
		       ,dictName(mm->probeDict, p2s->m.probe)
		       ,dictName(mm->probeDict, p2s->m.motif)
		       ,p2s->s) ;
	    }
	  array (val1, j, double) += p2s->s ;
	  array (val2, j, double) += p2s->s * p2s->s ;
	}
      if (0 && tt == histoTypes && mm->filter)
	fprintf (stderr, "%s\tkilled %d/%d\t%.1f%%\n", pns->nam, nf, np, (100.0 * nf)/np) ;
      /* regroup the last bins */
      if (0)
	{
	  for (j=10;j<24;j++)
	    fprintf(stderr, " %d",  array (histo, j, int)) ;
	  fprintf(stderr, "\n") ;
	}
      for (j = NN + 1; j > 0 ; j--)
	{
	  i = array (histo, j, int) ;
	  if (i < 50)
	    {
	      array (histo, j-1, int) += i ;
	      array (histo, j, int) =  0 ;
	      z1 = array (val1, j, double) ; 
	      z2 = array (val2, j, double) ;
	      array (val1, j, double) = 0 ;
	      array (val2, j, double) = 0 ;
	      array (val1, j-1, double) += z1 ;
	      array (val2, j-1, double) += z2 ; 
	    }
	}
     if (0)
       {
	 for (j=10;j<24;j++)
	   fprintf(stderr, " %d",  array (histo, j, int)) ;
	 fprintf(stderr, "\n") ;
       }
      /* regroup the first bins */
      for (j = 0 ; j < NN ; j++)
	{
	  i = array (histo, j, int) ;
	  if (i < 50)
	    {
	      array (histo, j+1, int) += i ;
	      array (histo, j, int) =  0 ;
	      z1 = array (val1, j, double) ; 
	      z2 = array (val2, j, double) ;
	      array (val1, j, double) = 0 ;
	      array (val2, j, double) = 0 ;
	      array (val1, j+1, double) += z1 ;
	      array (val2, j+1, double) += z2 ; 
	    }
	}
     if (0)
       {
	 for (j=10;j<24;j++)
	   fprintf(stderr, " %d",  array (histo, j, int)) ;
	 fprintf(stderr, "\n") ;
       }
      /* compute the average and the variance */
      for (j = 0 ; j <= NN ; j++)
	{
	  i = array (histo, j, int) ;
	  if (i)
	    {
	      z1 = array (val1, j, double) / i ;
	      z2 = array (val2, j, double) / i ;
	      sigma  = sqrt (z2 - z1*z1) ;
	      array (val1, j, double) = z1 ;
	      array (val2, j, double) = sigma ;
	    }
	}
      /* export the plot */
      if (!mm->pivot)
	{
	  freeOutf ("\n%s#%03d\t%03d.%s", pns->nam, ns, nHisto, tt->nam) ;
	  if (mm->histoN)
	    {
	      for (j = 0 ; j <= NN ; j++)
		{
		  z1 = array (histo, j, int) ;
		  if (z1 > 0)
		    freeOutf ("\t%g", z1) ;
		  else
		    freeOutf ("\t") ;
		}
	    }
	  else
	    {
	      for (i = 0 ; i < 2 ; i++)
		{
		  if (i == 1)
		    freeOutf ("\n%s#%03d\t%03d.sigma_%s", pns->nam,  ns, nHisto, tt->nam) ;
		  for (j = 0 ; j <= NN ; j++)
		    {
		      if (i == 0)
			z1 = array (val1, j, double) ;
		      else
			z1 = array (val2, j, double) ;
		      if (z1 > 0)
			freeOutf ("\t%g", z1) ;
		      else
			freeOutf ("\t") ;
		    }
		}
	    }
	}
      else  /* pivot table method */
	{
	  for (j = 0 ; j <= NN ; j++)
	    {
	      freeOutf ("%d", j+1) ; /* bin */
	      step = (tt->max - tt->min)/(NN-1) ;
	      freeOutf ("\t%g", tt->min + j*step) ; /* bin value */  
	      freeOutf("\t%s", pns->nam) ;
	      freeOutf ("\t%s", tt->nam) ; /* histo type */
	      /* type de cellule */
	      if (strstr(pns->nam, "ES"))
		freeOutf("\tES") ;
	      else if (strstr(pns->nam, "MS"))
		freeOutf("\tMS") ;
	      else
		freeOutf("\tX") ;
	      /* type d'echantillon */
	      if (strstr(pns->nam, "G1="))
		freeOutf("\tG1") ;
	      else if (strstr(pns->nam, "S="))
		freeOutf("\tS") ;
	      else if (strstr(pns->nam, "Samplified="))
		freeOutf("\tSamplified") ;
	      else if (strstr(pns->nam, "G1amplified="))
		freeOutf("\tG1amplified") ;
	      else if (strstr(pns->nam, "total="))
		freeOutf("\ttotal=") ;
	      else if (strstr(pns->nam, "NucPA-="))
		freeOutf("\tNucPA-=") ;
	      else if (strstr(pns->nam, "pA"))
		freeOutf("\tpA") ;
	      else
		freeOutf("\tX") ;
	      /* manip */
	      if (strstr(pns->nam, "=1"))
		freeOutf("\tmanip1") ;
	      else if (strstr(pns->nam, "=2"))
		freeOutf("\tmanip2") ;
	      else 
		freeOutf("\tmanip3") ;
	      /* chip */
	      if (strstr(pns->nam, "=1"))
		freeOutf("\tchip1") ;
	      else if (strstr(pns->nam, "=2"))
		freeOutf("\tchip1") ;
	      else 
		freeOutf("\tchip2") ;
	      
	      /* N  */
	      i = array (histo, j, int) ;
	      freeOutf ("\t%d", i) ;
	      /* valeur  */
	      z1 = array (val1, j, double) ;
	      freeOutf ("\t%g", z1) ;
	      /* valeur */
	      z1 = array (val2, j, double) ;
	      freeOutf ("\t%g", z1) ;

	      freeOutf ("\n") ;
	    }
	}
      if ((mm->expression && strstr (pns->nam, "=35ES")) ||
	  (mm->S && strstr (pns->nam, "S=15MS")) ||	  
	  (mm->S && strstr (pns->nam, "S=25MS")) ||	  
	  strstr (pns->nam, "G1=35ES") ||	  
	  !strncmp (pns->nam2, "zz", 2)
	  )
	{
	  freeOutf ("\nzz=zz00#999\t%03d.%s", nHisto, tt->nam) ;
	  for (j = 0 ; j <= NN ; j++)
	    freeOutf ("\t%g", tt->min + j*step) ; /* bin value */  
	  if (!strcmp (pns->nam2, "zzA") || !strstr (pns->nam, "G1=22MS"))
	    freeOutf ("\nzz=zz%s#996\t%03d.%s", "A", nHisto, tt->nam) ;
	  else if (!strcmp (pns->nam2, "zzB") || !strcmp (pns->nam, "G1=22MS"))
	    freeOutf ("\nzz=zz%s#997\t%03d.%s", "B", nHisto, tt->nam) ;
	  else
	    freeOutf ("\nzz=zz%s#998\t%03d.%s", "SecondChip", nHisto, tt->nam) ;
	  for (j = 0 ; j <= NN ; j++)
	    {
	      i = array (histo, j, int) ;
	      freeOutf ("\t%d", i) ;
	    }
	  freeOutf ("\n") ;
	}

    }
  ac_free (histo) ;
  ac_free (val1) ;
  ac_free (val2) ;
} /* p2sSingleChipManyHistos */
  
/*************************************************************************************/
/* Export per chip a bunch of histos for various Sequence parameters */
static void bbHisto (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN fi = aceInCreateFromStdin (FALSE, 0, h) ;
  int i, nn = 0, ns, probe ;
  const char *snam1 ;
  char motif[300], mybuf[200] ;
  float s1, tm ;
  double z1, lg2 = log((double)2.0) ;
  PNS *pns ;
  P2S *p2s = 0, *p2sold ;
  char pNam[400], *cp1 ;

  for (i = 0 ; i < NFF ; i++) 
    { FF[i] = 0 ; FF2[i] = 0 ; FFM[i] = "" ; }
  FFM[NFF-1] = "total" ;
  if (!mm->single)
    messcrash ("-histo needs an input file of type -single") ;

  while (nn++ > -1 && aceInCard (fi))
    {
      const char *ccp = aceInWord(fi) ;
      if (!ccp) continue ; /* first probe name */
      strcpy (pNam, ccp) ;

      if (!(snam1 = aceInWord(fi))) continue ;
      strcpy (mybuf, snam1) ;
      if (0 && !strstr(mybuf, "5MS")) continue ;
      cp1 = strstr (mybuf, ".pair") ;
      if (cp1) *cp1 = 0 ;
      for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
	if (!strcmp (pns->nam, mybuf) || !strcmp (pns->p, mybuf)) 
	  break ;
      if (!pns->p)
	continue ;
      snam1=pns->p ; 
      strcpy (mybuf, snam1) ;
#ifdef BOU_DEF
      if (mm->G1 && ! (pns->flag & PGG_G1))
	continue ;
      if (mm->expression && !pns->isRna)
	continue ; 
      if (mm->S && ! (pns->flag & PGG_S))
	continue ;  
      if (mm->pAplus && ! (pns->flag & PGG_PA))
	continue ;  
      if (mm->manip && !strstr (snam1, mm->manip))
	continue ;
      if (mm->pAmoins && ! (pns->flag & PGG_NUC))
	continue ;
      if (mm->green && !strstr (pns->p,"_532")) /* green dye */
	continue ;
      if (mm->red && !strstr (pns->p,"_635"))   /* red dye */
	continue ;
      if (mm->ratio && !strstr (pns->p,"_532") && !strstr (pns->p,"_635"))
	continue ;
#endif
      dictAdd (mm->probeDict, pNam, &probe) ;

      if (!aceInFloat (fi, &s1)) continue ;
      p2sold = p2s ;
	{
	  ccp = aceInWord(fi) ;
	  strncpy (motif, ccp, 300) ;
	  aceInFloat (fi, &tm) ;
	  if (!pns->signal)
	    pns->signal = arrayCreate (1000000, P2S) ;
	  p2s = arrayp (pns->signal, arrayMax(pns->signal), P2S) ;
	  p2s->m.probe = probe ;
	  z1 = s1 ;
	  if (0) /* do not normalize */ 
	    {
	      pns->damper = 1; pns->av = 10 ;
	    }
	  if (strstr(pns->p,"AFX"))
	    {
	      s1 = log((double)z1 + pns->damper)/lg2 ; 
	    }
	  else
	    {
	      s1 = log((double)z1 + pns->damper)/lg2 ; 
	    }
	  p2s->s = s1  - pns->av + 10.0 ;
	  p2s->ns = ns ;
	  p2s->m.tm = tm ;
	  p2sMotifData (mm, pns, p2s, motif, pNam) ;
	}
      if (p2s && mm->ratio && p2sold && strstr (pns->nam2,"ratio"))
	{
	  if (p2s->m.probe == p2sold->m.probe &&
	      p2sold->ns == ns + 1
	      )
	    {  p2s->s += 10 - p2sold->s ;  p2sold->s = -1000 ; }
	  else
	    p2s->s = -1000 ;
	}
    }
  for (i = 0 ; i < NFF ; i++) 
    if (FF[i] > 0) fprintf(stderr, "Filter number\t%2d\tkills\t%d\t%d\tprobes\t%s\n", i, FF[i],FF2[i],FFM[i]) ;

  p2sAllChipsCompressSignal (mm) ;
  /* Export Histograms of signals for all chips in a single table */
  if (!mm->pivot)
    {
      p2sAllChipsSignalHisto (mm) ;
      freeOutf ("\n\nHISTOS\n") ;
    }
  /* Export per chip a bunch of histos for various Sequence parameters */
  for (ns = 0, pns = pnsAll ; pns->p ; ns++, pns++)
    if (pns->signal &&
	(!mm->ratio || strstr (pns->nam2, "ratio"))
	)
      p2sSingleChipManyHistos (mm, pns, ns) ; 
  freeOutf ("\n") ;
  ac_free (h) ;
} /* bbHisto */

/*************************************************************************************/
/*************************** Find the peaks  *****************************************/
/*************************************************************************************/
/*
  Solexa STAT3_IL21 t7_Mm7_39453_37_0.STAT3_IL21.slx     RED
  Solexa STAT3_noIL21 t7_Mm7_39453_37_0.STAT3_noIL21.slx

  Solexa STAT3_IL6_CD4 t7_Mm7_39453_37_0.STAT3_IL6_CD4.slx    MAGENTA
  Solexa STAT3_noIL6_CD4 t7_Mm7_39453_37_0.STAT3_noIL6_CD4.slx

  Solexa IRF4_IL21 t7_Mm7_39453_37_0.IRF4_IL21.slx       GREEN
  Solexa IRF4_noIL21 t7_Mm7_39453_37_0.IRF4_noIL21.slx   

  Solexa mIgG1_IL6_CD4 t7_Mm7_39453_37_0.mIgG1_IL6_CD4.slx

*/
#define slxMax  2048
typedef struct solexaDataStruct { double g ; int n, x1, a1 ; } SLXD ;
typedef struct solexaStruct { Array aa[slxMax] ; } SLX ;

static DICT *peakDict = 0 ;

static Array bbGetSolexaArray (const char *fnam, int g1, AC_HANDLE h)
{
  static Array ww = 0 ;
  Array aa = 0 ;
  ACEIN ai ;
  int step = 10, width = 100 ;
  int ii = 0, jj, j1, j2, nn, x1, a1, iMax = 0 ;
  int dx = 0 ;
  double w, dw, dxsq, sigma2, wdummy ;
  SLXD *slxd ;
  BOOL debug = FALSE ;

  ai = aceInCreateFromFile (messprintf("SOLEXA/%s.u", fnam), "r", 0, h) ;
  if (!ai) return 0 ;
  if (!ww) /* lazy store the exponentials */ 
    ww = arrayCreate (100, double) ;
  
  /* get the data */
  while (aceInCard(ai))
    {
      if (aceInInt (ai, &x1) &&  /* coordinate on cosmid */
	  aceInInt (ai, &nn) &&  /* tag count */
	  aceInWord (ai) &&      /* chromosome */
	  aceInInt (ai, &a1)     /* intmap */
	  )
	{
	  if (!aa)
	    {
	      aa = arrayHandleCreate (100000, SLXD, h) ;
	    }
	  if (debug && (a1 < 44110000 || a1 > 44180000))
	    continue ;
	  ii = x1/step ;
	  slxd = arrayp (aa, ii, SLXD) ;
	  slxd->n = nn ;
	  slxd->x1 = x1 ;
	  slxd->a1 = a1 ;
	  if (ii > iMax) iMax = ii ;
	}
    }
  ac_free (ai) ;
  if (iMax)
    {
      /* create the coordinates also between tags */
      for (ii = 0, slxd = arrp (aa, ii, SLXD) ; ii < iMax ; slxd++, ii++)
	{
	  slxd->x1 = ii * step ;
	  slxd->a1 = g1 + ii * step ;
	}

      /* gauss average over width bp */
      sigma2 = (double) (width * width * 2)/(step *step) ;
      for (ii = 0 ; ii <= iMax ; ii++)
	{
	  slxd = arrp (aa, ii, SLXD) ;
	  if (0 && slxd->a1 != 6894685)
	    slxd->n = 0 ;
	  nn = slxd->n ;
	  if (!nn)
	    continue ;
	  /* spread the data over the neighbouring kb */
	  j1 = ii - 5* width/step ;
	  j2 = ii + 5* width/step ;
 	  if (j1 < 0) j1 = 0 ;      
	  if (j2 > iMax-1) j2 = iMax -1 ;
	  
	  /* measure the total relevant weight */
	  for (w = 0, jj = j1 ; jj <= j2 ; jj++)
	    {
	      dx = ii - jj ;
	      if (dx < 0) dx = -dx ;
	      dw = (dx >= arrayMax (ww)) ? 0 : arr (ww, dx, double) ;
	      if (dw == 0)
		{
		  dxsq = dx * dx/sigma2 ;
		  dw = array (ww, dx, double) = exp (- dxsq) ;
		}
	      w += dw ;
	    }
	  /* spread the data */
	  wdummy = 0 ;
	  for (jj = j1 ; jj <= j2 ; jj++)
	    {
	      dx = ii - jj ;
	      if (dx < 0) dx = -dx ;
	      slxd = arrayp (aa, jj, SLXD) ;
	      slxd->g += nn * arr (ww, dx, double)/w ;
	      wdummy +=  nn * arr (ww, dx, double)/w ;
	    }
	  jj++ ;
	}
      if (debug)
	for (jj = 0 ; jj < arrayMax(ww) ; jj++)
	  fprintf (stderr, "ww %d %g\n", jj, arr (ww, jj, double)) ;
    }
  return aa ;
}

/*************************************************************************************/

typedef struct dfeStruct {int tag, run, control, ln, area, total, yes, color, type, oldnam, title, contrib[64] ;} DFE2 ;

static int bbGetSolexaData (AC_OBJ seq, SLX *slx, Array dfes, AC_HANDLE h)
{
  int ir, ii, jj, a1, nn = 0 ;
  AC_HANDLE h1 = ac_new_handle () ;
  AC_TABLE iMap = ac_tag_table (seq, "IntMap", h1) ;
  AC_TABLE tbl = ac_tag_table (seq, "Solexa", h1) ;
  const char *ccp, *ccq ;
  Array aa = 0 ;
  DFE2 *dfe ;
  BOOL debug = FALSE ;
  
  if (iMap && iMap->rows && iMap->cols >2)
    {
      a1 = ac_table_int (iMap, 0, 1, 0) ;
      
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, 0) ;
	  if (!ccp) continue ;
	  dictAdd (peakDict, ccp, &jj) ;
	  if (jj >= slxMax) 
	    messcrash ("Please increase slxMax") ;
	  ccq = ac_table_printable (tbl, ir, 1, 0) ;
	  if (!ccq) continue ;
	  for (ii = 0 ; ii < arrayMax(dfes) ; ii++)
	    {
	      dfe = arrp (dfes, ii, DFE2) ;
	      if (dfe->tag == jj)
		break ;
	    }
	  if (ii >= arrayMax(dfes) || dfe->tag != jj)
	    {
	      if (0) messout ("Cannot find the definition of track %s %s\n", dictName (peakDict, jj), ccq) ;
	      continue ;
	    }
	  aa = bbGetSolexaArray (ccq, a1, h) ;
	  if (aa)
	    nn++ ;
	  slx->aa[jj] = aa ;
	  if (debug && aa && jj == 0)
	    {
	      SLXD *slxd ;
	      int i ;
	      for (i = 0 ; i < arrayMax(aa) ; i++)
		{
		  slxd = arrp (aa, i, SLXD) ;
		  if (slxd->n || slxd->g > 0)
		    fprintf(stderr, "%d\t%d\t%d\t%g\n",slxd->x1,slxd->a1,slxd->n,slxd->g) ;
		}
	    }
	}
    }
  ac_free (h1) ;
  return nn ;
} /* bbGetSolexaData */

/*************************************************************************************/

static int bbExportMotifs (const char *dna, int ww)
{
  int nn = 0, i2, delta ;
  const char **cm, *cp1, *cp2 ; 
  const char *motifs[] = {   "cttcct", "aggaag"
			     , "ttcctg", "caggaa"
			     , "acttcc", "ggaagt"
			     , "gggcgg", "ccgccc"
			     , "gactca", "tgagtc"
			     , "ggcggg", "cccgcc"
			     , "tgactc", "gagtca"
			     , "cggaag", "cttccg"
			     , "ccccgc", "gcgggg"
			     , "cgcccc", "ggggcg"
			     , "ccggaa", "ttccgg"
			     , "gccgcc", "ggcggc"
			     , "ggcgga", "tccgcc"
			     , "atgact", "agtcat"
			     , "gtttca", "tgaaac"
			     , "actcat", "atgagt"
			     , "gaaacc", "ggtttc"
			     , "agtttc", "gaaact"
			     , "gcgcgc"
			     , "gaaac", "gtttc"
			     , "agtca", "tgact"
			     , "tttct", "agaaa"
			     , "ggaag", "cttcc"
			     , "ggcgg", "ccgcc"
                             , "gcgcgc"
			     , "tgactca", "tgagtca"
			   , 0} ;

  const char *motifs2[] = { 
      "ttc", "gaa"
    , "tac", "gaa"
    , "ttc", "gta"
    , "ctc", "gaa"
    , "ttc", "gag"
    , "tcc", "gaa"
    , "ttc", "gga"
    , "tgc", "gaa"
    , "ttc", "gca"
    , "gtc", "gaa"
    , "ttc", "gac"
    , "tta", "gaa"
    , "ttc", "taa"
    , "ttg", "gaa"
    , "ttc", "caa"
    , "atc", "gaa"
    , "ttc", "gat"
    , "ttt", "gaa"
    , "ttc", "aaa"
    , "agg", "gaa"
    , "ttc", "cct"
    , "gaa", "agg"
    , "cct", "ttc"
    , "gga", "aag"
    , "ctt", "tcc"
    , "aag", "gga"
    , "tcc", "ctt"
    , "cca", "cct"
    , "agg", "tgg"
    , "cac", "ctg"
    , "cag", "gtg"
    , "ttc", "ttt"
    , "aaa", "gaa"
    , "ccg", "ccg"
    , "cgg", "cgg"
    , "ccg", "cct"
    , "agg", "cgg"
    , "tcc", "aag"
    , "ctt", "gga"
    , "agg", "ccg"
    , "cgg", "cct"
    , "cag", "ccg"
    , "cgg", "ctg"
    , "ccc", "ccg"
    , "cgg", "ggg"
    , "gcc", "ccc"
    , "ggg", "ggc"
    , 0, 0
  } ;
 

  const char *motifs3[] = { 
        "ttc", "a", "gaa"
      , "ttc", "t", "gaa"
      , "ttc", "g", "gaa"
      , "ttc", "c", "gaa"
      , 0, 0, 0
  } ;
 
  const char *motifsIRF[] = { 
     "tttc", "tt"
    , "aa", "gaaa"
    , 0, 0
  } ;

  for (cm = motifs ; *cm ; cm++)
    {
      cp1 = dna ;
      while ((cp1 = strstr (cp1, *cm))) 
	{
	  i2 = cp1 - dna ;
	  if (1) /* here, we could limit the geometrical distance */
	    {
	      nn++ ;
	      /* export the centre of the motif relative to the centre of the peak */
	      freeOutf("%s %d\n", *cm, i2 - ww + strlen(*cm)/2) ;
	    }
	  cp1++ ;
	}
    }
  for (cm = motifs2 ; *cm ; cm+=2)
    {
      cp1 = dna ;
      while ((cp1 = strstr (cp1, *cm))) 
	{
	  i2 = cp1 - dna ;
	  if (1)
	    if ((cp2 = strstr(cp1+6,*(cm+1))) && cp2 - cp1 == 6)
	      {
		nn++ ;
		freeOutf("%s___%s %d\n", *cm, *(cm+1), i2 - ww + 5) ;
	      }
	  cp1++ ;
	}
    }
  if (0) /* not informative */
    for (cm = motifs3 ; *cm ; cm+=3)
      {
	cp1 = dna ;
	while ((cp1 = strstr (cp1, *cm))) 
	  {
	    i2 = cp1 - dna ;
	    if (1)
	      if ((cp2 = strstr(cp1+6,*(cm+2))) && cp2 - cp1 == 6 &&
		  *(cp1 + 4) == **(cm+1)
		  )
		{
		  nn++ ;
		  freeOutf("%s_%s_%s %d\n", *cm, *(cm+1), *(cm+2), i2 - ww + 5) ;
		}
	    cp1++ ;
	  }
      }
  for (cm = motifsIRF, delta=6 ; *cm ; delta -= 2, cm += 2)
    {
      cp1 = dna ;
      while ((cp1 = strstr (cp1, *cm))) 
	{
	  i2 = cp1 - dna ;
	  if (1)
	    if ((cp2 = strstr(cp1+delta,*(cm+1))) && cp2 - cp1 == delta)
	      {
		nn++ ;
		freeOutf("%s__%s %d\n", *cm, *(cm+1), i2 - ww + 5) ;
	      }
	  cp1++ ;
	}
    }
  return nn  ;
} /* bbExportMotifs */

/*************************************************************************************/

/* create a DFE structure from the Solexanames.txt file */
static Array bbGetSolexaSeqPeaksMake (DICT *dict, AC_HANDLE h)
{
  int i, ii = 0, n ;
  const char *ccp ;
  DFE2 dd, *dfe ;
  Array dfes = arrayHandleCreate (256, DFE2, h) ;
  ACEIN ai = 0 ;

#ifdef WARREN
  ai = aceInCreateFromFile ("/home/mieg/mm/LeonardInfo/SolexaNames.txt", "r", 0, h) ;
#else
  ai = aceInCreateFromFile ("/home/mieg/36a/SeqcInfo/SolexaNames.txt", "r", 0, h) ;
#endif

  if (! ai) messcrash("Cannot find the config file SolexaNames.txt") ;
  aceInSpecial (ai, "\n") ; /* so \t will not be expanded */
  while (aceInCard (ai))
    {
      /* 1 tag name
	 2 run
	 3 run igg
	 4 effective clone length
	 5 area threshold
	 6 nb tags
	 7 ys
	 8 color
	 9 type
	 10 olf tag name
	 11 titre
	 12-31 les positifs
	 32-37 les negatis
      */
      memset (&dd, 0, sizeof(dd)) ;
      /* get tag name */
      ccp = aceInWordCut (ai, "\t", 0) ; 
      if (!ccp || *ccp == '#') continue ; 
      dictAdd (dict, ccp, &n) ; 
      dd.tag = n ;
      
      /* get run name */
      ccp = aceInWordCut (ai, "\t", 0) ; 
      if (ccp && *ccp != '-')
	{
	  dictAdd (dict, ccp, &n) ; 
	  dd.run = n ;
	}
      /* get control name */
      ccp = aceInWordCut (ai, "\t", 0) ; 
      if (ccp && *ccp != '-')
	{
	  dictAdd (dict, ccp, &n) ; 
	  dd.control = n ;
	}
      /* jump 4 and 5, do not look at the threshold */
      for (i = 4 ; i <= 5 ; i++)
	ccp = aceInWordCut (ai, "\t", 0) ; 
      /* get total number of uniquely aligned tags */
      ccp = aceInWordCut (ai, "\t", 0) ; 
      if (ccp && *ccp && (n=atoi(ccp))>0)
	dd.total = n ;
      else
	dd.total = 3000000 ;
      /* jump to yes integrate */
      for (i = 7 ; i <= 7 ; i++)
	ccp = aceInWordCut (ai, "\t", 0) ; 
      if (ccp && strstr (ccp, "y"))
	dd.yes = 1 ;
      if (ccp && strstr (ccp, "y"))
	dd.yes |= 0x2 ;
      /* jump to title */
      for (i = 8 ; i <= 11 ; i++)
	ccp = aceInWordCut (ai, "\t", 0) ; 
      /* find the positive runs */
      {
	int j = 0 ;
	for (j=0 ; j < 64 ; j++)
	  dd.contrib[j] = 0 ;
	for (j=0, i = 12 ; i <= 31 ; i++)
	  {
	    ccp = aceInWordCut (ai, "\t", 0) ; 
	    if (ccp && *ccp != '-')
	      {
		dictAdd (dict, ccp, &n) ; 
		dd.contrib[j++] = n ;
	      }
	  }
	/* find the negative runs */
	for (i = 32 ; i <= 37 ; i++)
	  {
	    ccp = aceInWordCut (ai, "\t", 0) ; 
	    if (ccp && *ccp != '-')
	      {
		dictAdd (dict, ccp, &n) ; 
		dd.contrib[j++] = - n ;
	      }
	  }
      }
      if (dd.run)
	{
	  dd.contrib[0] = dd.run ;
	  dd.contrib[1] = 0 ;
	}
      dfe = arrayp (dfes, ii++, DFE2) ;
      *dfe = dd ;
    }
  ac_free (ai) ;
  /* edit the contrib table */
  {
    KEYSET ks = keySetCreate () ;
    int  jj, sign, nc ;
    /* register */
    for (ii = 0 ; ii < arrayMax(dfes) ; ii++)
      {
	dfe = arrayp (dfes, ii, DFE2) ;
	keySet (ks, dfe->run) = ii+1 ;
      }
    /* replace all contrib run name by their offset in dfes */
    for (ii = 0 ; ii < arrayMax(dfes) ; ii++)
      {
	dfe = arrayp (dfes, ii, DFE2) ;
	if (dfe->control)
	  {
	    nc = keySet (ks, dfe->control) ;
#ifdef WARREN_OLD
	    if (!nc) 
	      messcrash ("Definition of line %s has no associated control\n"
			 , dictName (dict, dfe->tag)
			 ) ;
#endif
	  }
	for (jj=0;jj<64;jj++)
	  {
	    sign = 1 ; n = dfe->contrib[jj] ;
	    if (n<0) { sign = -1 ; n = -n ; }
	    if(n==0) break ;
	    nc = keySet (ks, n) ;
#ifdef WARREN_OLD 
	    if (!nc) 
	      messcrash ("Definition of line %s has no associated control\n"
			 , dictName (dict, dfe->tag)
			 ) ;
#endif
	    dfe->contrib[jj] = sign * nc ;
	  }
      }

    /* add the controls as offset in dfes */
    for (ii = 0 ; ii < arrayMax(dfes) ; ii++)
      {
	int j, j0, jj, nc ;
	DFE2 *dfe2 ;
	
	dfe = arrayp (dfes, ii, DFE2) ;
	/* find the first empty spot */
	for (jj=0;jj<64;jj++)
	  if (dfe->contrib[jj] == 0)
	    break ;
	j0 = jj ;
	/* add the controls */      
	for (j=0; j < j0 ; j++)
	  {
	    sign = 1 ; n = dfe->contrib[j] ;
	    if (n<0) { sign = -1 ; n = -n ; }

	    if (n<1 || n > arrayMax(dfes))
	      messcrash ("Confusion while adding controls") ;
	    dfe2 = arrayp (dfes, n-1, DFE2) ;
	    nc = dfe2->control ;
	    if (!nc) continue ;
	    continue ;             /* DO NOT REMOVE the IgG controls */
	    nc = keySet (ks, nc) ;
	    dfe->contrib[jj++] = - sign * (1000 + nc) ;
	  }
	dfe->contrib[jj++] = 0 ;
      }
    }
  return dfes ;
} /* bbGetSolexaSeqPeaksMake */ 

/*************************************************************************************/

static int bbGetSolexaSeqPeaks (AC_OBJ seq, SLX *slx, Array dfes, AC_HANDLE h)
{
  int iMax = 0, idfe, ii, jj, kk, n, nc, sign, nn = 0, a1=0, p1=0, np = 0, pmax, x1 = 0, total ;
  double zoom[slxMax], ddz, ddzControl, ddzmax, cumul = 0, cumulControl = 0 ;
  double zoom3 = 25.50 ;     /* effective zoom at gaussian smoothing 200bp */
  double zmin=1, min = 5.0/zoom3 ;   /* peak starts/ends when crossing 5 if seen at gaussian smoothing 200bp */
  double minCumul = 4 ; /* true tag count in the minimal peak, per 3M tags */
  KEY map ;
  BOOL debug = FALSE ;
  SLXD *slxd ; 
  char buf100[100] ;
  DFE2 *dfe, *dfeb ;
  Array aa = 0 ;
  
  /* adjust the zoom */
  memset (zoom, 0, sizeof(zoom)) ;
  for (ii = 0 ; ii < arrayMax(dfes) ; ii++)
    {
      dfe = arrayp (dfes, ii, DFE2) ;
      zoom[ii] = dfe->total > 0 ? (double)3000000/dfe->total : 1 ;
    }
  freeOutf ("//\t#%s\n",  timeShow (timeNow(), buf100, 100)) ;
  freeOutf ("//\t#Experiment\tChromosome\tCenter\tTop\tSequence\tCenter\tTop\tWidth\tThreshold\tArea\n") ;
  map = ac_tag_key (seq, "IntMap", 0) ;
  for (jj = 0 ; jj < slxMax ; jj++)
    {
      if (slx->aa[jj] && arrayMax (slx->aa[jj]) &&
	  (iMax == 0 || iMax < arrayMax (slx->aa[jj]))
	  )
	iMax = arrayMax (slx->aa[jj]) ;
    }
  for (idfe = 0 ; idfe < arrayMax (dfes) ; idfe++) /* loop on all peak types */
    {
      dfe = arrp (dfes, idfe, DFE2) ;
      if (!dfe->yes & 0x1) continue ;
      p1 = np = pmax = 0 ;
      cumul = cumulControl = ddzmax = 0 ;
      /* loop on all positions */
      for (ii = 0 ; ii < iMax ; ii++)
	{
	  /* get the positive data */
	  ddz = ddzControl = 0 ; zmin = 0 ;
	  total = 0 ;
	  for (kk = 0 ; kk < 64 ; kk++)
	    {
	      n = dfe->contrib[kk] ;

	      nc = 0 ;
	      if (n < 0) { sign = -1 ; n = -n ; }
	      else if (n>0) { sign =  1; }
	      else break ;
	      if (n >= 1000) { nc = 1 ; n -= 1000 ; }
	      if (n==0) continue ;
	      if (n < 1 || n > arrayMax (dfes))
		messcrash ("Confusion 1 in contribs of line %s", dictName (peakDict, dfe->tag)) ;
	      dfeb = arrp (dfes, n-1, DFE2) ;  /* contributing track */
	      aa = slx->aa[dfeb->tag] ;
	      if (!dfeb->tag)
		messcrash ("Confusion 2 in contribs of line %s", dictName (peakDict, dfe->tag)) ;

	      if (arrayExists(aa) && ii < arrayMax (aa))
		{		  
		  slxd = arrp (aa, ii, SLXD) ;
		  ddz += sign * slxd->g * zoom[n-1] ;
		  if (nc && sign > 0) ddzControl += slxd->g ;
		  a1 = slxd->a1 ;
		  x1 = slxd->x1 ;
		  if (sign > 0)
		    {
		      zmin++ ;
		      total += dfe->total ;
		    }
		}
	    }

	  if (!zmin) continue ;
	  /* rescale the base threshold and take the difference */

	  zmin = sqrt (total/((double)3000000.0)) ;
	  zmin = zmin * sqrt(zmin) ;
	  
	  ddz -= zmin*min ;
	  
	  if (ddz > 0) /* found a peak */
	    {
	      if (!p1) 
		{
		  cumul = cumulControl = ddzmax = 0 ;
		  if (debug) fprintf (stderr, "peak starts at %d ddz %g cumul %g\n", a1, ddz, cumul) ;
		  p1 = pmax = a1 ; 
		}
	      else   if (debug) fprintf (stderr, "peak inside at %d ddz %g cumul %g\n", a1, ddz, cumul) ;
	      
	      if (ddz == ddzmax && (randint() & 0x1)) { ddzmax = ddz ; pmax = a1 ; }
	      else if (ddz > ddzmax) { ddzmax = ddz ; pmax = a1 ; }
	      cumul += ddz ;
	      cumulControl += ddzControl ;
	      np++ ;
	    }
	  else
	    {
	      cumul += zmin * min * np ; 
	      if (p1 && cumul > zmin * minCumul && cumul > 3 * cumulControl)
		{
		  int xc =  (a1+p1)/2 - (a1 - x1) ;
		  int ww = (a1 - p1 + 2)/2; /* half width */
		  /* 		  int xmax = pmax - (a1 - x1) ; */
		  const char *dna = 0 ;
		  const char *dnac = 0 ;
#ifdef WARREN
		  dna = ac_zone_dna (seq, xmax - 300, xmax + 300, h) ;
		  dnac = ac_zone_dna (seq, xc - ww, xc + ww, h) ;		  
#endif
		  nn++ ;
		  freeOutf ("//\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%g\t%g\t%s\n"
			    , dictName (peakDict, dfe->tag)
			    , ac_key_name(map), (a1+p1)/2, pmax
			    , ac_name(seq), xc, pmax - (a1 + p1)/2
			    , a1 - p1, min, cumul, dna ? dna : "-") ;
		  freeOutf("Sequence %s\nPBS %s_%s_%d %d %d\n\n",  ac_name(seq), dictName (peakDict, dfe->tag), ac_key_name(map), (a1+p1)/2, xc - 300, xc + 300) ;
		  freeOutf("PBS %s_%s_%d\n",  dictName (peakDict, dfe->tag),  ac_key_name(map), (a1+p1)/2) ;
		  freeOutf("%s\n", dictName (peakDict, dfe->tag)) ;
		  freeOutf("IntMap %s %d %d\n",  ac_key_name(map), (a1+p1)/2-300,  (a1+p1)/2+300) ;
		  freeOutf("Max %g %d\nArea %g\nWidth %d\n", zoom3*ddzmax, 10* ((int)((pmax - (a1 + p1)/2 + 5)/10)), cumul, 2*ww+1) ;
		  if (dnac) bbExportMotifs (dnac, ww) ;
		  freeOutf("\n") ;
		}
	      np = 0 ;
	      cumul = cumulControl = ddzmax = 0 ;
	      p1 = 0 ;
	    }
	}
      
      cumul += zmin * min * np ; 
      if (p1 && cumul  > zmin * minCumul && cumul > 3 * cumulControl)
	{
	  int xc =  (a1+p1)/2 - (a1 - x1) ;
	  int ww = (a1 - p1 + 2)/2; /* half width */
	  /* 	  int xmax = pmax - (a1 - x1) ; */
	  const char *dna = 0 ;
	  const char *dnac = 0 ;
#ifdef WARREN
		  dna = ac_zone_dna (seq, xmax - 300, xmax + 300, h) ;
		  dnac = ac_zone_dna (seq, xc - ww, xc + ww, h) ;		  
#endif
	  nn++ ;
	  freeOutf ("//\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%g\t%g\t%s\n"
		    , dictName (peakDict, dfe->tag)
		    , ac_key_name(map), (a1+p1)/2, pmax
		    , ac_name(seq), xc, pmax - (a1 + p1)/2
		    , a1 - p1, min, cumul, dna ? dna : "-") ;
	  freeOutf("Sequence %s\nPBS %s_%s_%d\n\n",  ac_name(seq), dictName (peakDict, dfe->tag), ac_key_name(map), (a1+p1)/2) ;
	  freeOutf("PBS %s_%s_%d\n",  dictName (peakDict, dfe->tag),  ac_key_name(map), (a1+p1)/2) ;
	  freeOutf("%s\n", dictName (peakDict, dfe->tag)) ;
	  freeOutf("IntMap %s %d %d\n",  ac_key_name(map), (a1+p1)/2-300,  (a1+p1)/2+300) ;
	  freeOutf("Max %g %d\nArea %g\nWidth %d\n", zoom3*ddzmax, pmax - (a1 + p1)/2, cumul, 2*ww+1) ;
	  if (dnac) bbExportMotifs (dnac, ww) ;
	  freeOutf("\n") ;
	}
    } /* kk loop on all peak types */
  
  return nn ;
} /* bbGetSolexaSeqPeaks */

/*************************************************************************************/

static int bbGetSolexaPeaks (MM *mm)
{
  int n1, nn = 0 ;
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ seq = 0 ;
  SLX slx ;
  Array dfes = 0 ;

  memset (&slx, 0, sizeof(slx)) ;
  if (!mm->db) 
    messcrash ("-db needed to find the names of the sections and files \n") ;
  if (!peakDict)
    peakDict = dictHandleCreate (64, h) ;
  dfes = bbGetSolexaSeqPeaksMake (peakDict, h) ;
  iter = ac_dbquery_iter (mm->db, "Find sequence solexa  ", h) ; /* IS t10_Mm10_39532_37_62 */
  while (ac_free (seq), ac_free (h1), seq = ac_iter_obj (iter))
    {
      fprintf (stderr, "%s\n", ac_name(seq)) ;
      h1 = ac_new_handle () ;
      n1 = bbGetSolexaData (seq, &slx, dfes, h1) ;
      nn += bbGetSolexaSeqPeaks (seq, &slx, dfes, h1) ;
    }
 ac_free (h) ;
 freeOutf ("\n// Done\n") ;
  return nn ;
} /* bbGetSolexaPeaks */

/*************************************************************************************/
/*************************************************************************************/
/* spread over a gauss curve of 100 bp */
static void bbGauss (Array aa, int ii)
{
  double w, dw, *wp, dxsq, sigma2 ;
  int j1, j2, jj, step = 1, width = 30, dx ;
  static Array ww = 0 ;

  if (! ww) 
    ww = arrayCreate (200, double) ;
  sigma2 = (double) (width * width * 2)/(step *step) ;
  
  /* spread the data over the neighbouring kb */
  j1 = ii - 5* width/step ;
  j2 = ii + 5* width/step ;
  if (j1 < 0) j1 = 0 ;      
  w = array (aa, j2, double) ; /* make space */
	  
  /* measure the total relevant weight */
  for (w = 0, jj = j1 ; jj <= j2 ; jj++)
    {
      dx = ii - jj ;
      if (dx < 0) dx = -dx ;
      dw = (dx >= arrayMax (ww)) ? 0 : arr (ww, dx, double) ;
      if (dw == 0)
	{
	  dxsq = dx * dx/sigma2 ;
	  dw = array (ww, dx, double) = exp (- dxsq) ;
	}
      w += dw ;
    }
  /* spread the data */
  for (jj = j1 ; jj <= j2 ; jj++)
    {
      dx = ii - jj ;
      if (dx < 0) dx = -dx ;
      wp = arrp (aa, jj, double) ;
      *wp += arr (ww, dx, double)/w ;
    }
  return ;
}

/*************************************************************************************/
/* autocorrelation des signaux forward et reverse */
typedef struct siteStruct { int site, map, a1, a2 ; } SITE ;
static int bbSolexaAutocorrel (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  Array dd = arrayHandleCreate (100000000, double,h) ;
  Array rr = arrayHandleCreate (100000000, double,h) ;
  int a1, a2, map, dx, x1, x2, ii, jj, iMax ;
  double uu, vv, uv[500], u2[500], z, *z1p, *z2p ;
  Array sites[26] ;
  SITE *up ;
  ACEIN ai = 0 ;
  
  memset (sites, 0, sizeof(sites)) ;
  if (mm->solexaSitesFileName)
    {
      ai = aceInCreateFromFile (mm->solexaSitesFileName, "r", 0, h) ;
      if (ai)
	{ 
	  while (aceInCard (ai))
	    {
	      if (aceInInt (ai, &map) &&
		  map >0 && map < 25 &&
		  aceInInt (ai, &a1) &&
		  aceInInt (ai, &a2)
		  )
		{
		  if (!sites[map])
		    sites[map] = arrayHandleCreate (32, SITE, h) ;
		  up = arrayp (sites[map], arrayMax(sites[map]), SITE) ;
		  up->map = map ;
		  up->a1 = a1 ;
		  up->a2 = a2 ;
		}
	    }
	  ac_free (ai) ;
	}
    }

  ai = aceInCreateFromStdin (0, 0, h) ;
  while (aceInCard (ai))
    {
      if (aceInInt (ai, &map) && aceInInt (ai, &x1) && aceInInt (ai, &x2))
	{
	  if (mm->solexaSitesFileName)
	    {
	      BOOL ok = FALSE ;
	      int iSite ;

	      if (map > 0 && map < 25 && sites[map])
		{
		  for (iSite = 0 ; iSite < arrayMax(sites[map]) ; iSite++)
		    {
		      up = arrp (sites[map], iSite, SITE) ;
		      if (up->a1 <= x1 && up->a2 >= x2)
			{ ok = TRUE ; break ; }
		    }
		}
	      if (!ok) continue ;
	    }
	  if (x1 < x2) bbGauss (dd, x1) ;
	  else bbGauss (rr, x1) ;
	}
    }
  ii = arrayMax (dd) ;
  jj = arrayMax (rr) ;
  iMax = ii < jj ? ii : jj ; iMax -= 1000 ;
  uu = vv = 0 ; 
  for (dx = 0 ; dx < 500 ; dx++)
    uv[dx] = u2[dx] = 0 ;
  for (ii = 0, z1p = arrp(dd,0,double), z2p = arrp(rr,0,double) ; ii < iMax ; z1p++, z2p++, ii++)
    { 
      z = *z1p ; uu += z*z ;  
      for (dx = 0 ; dx < 500 ; dx++)
	{ u2[dx] += z * (z1p[dx]) ; uv[dx] += z * (z2p[dx]) ; }
      z = *z2p ; vv += z*z ; 
    }
  for (dx = 0 ; dx < 500 ; dx++)
    {
      uv[dx] /= sqrt (uu * vv) ;
      u2[dx] /= uu ;
    }
  
  for (dx = 0 ; dx < 500 ; dx++)
    freeOutf ("%d\t%g\t%g\n", dx, u2[dx], uv[dx]) ;
  ac_free (h) ;
  return 0 ;
} /* bbSolexaAutocorrel */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: boubou ACEDB [-o outfile] [-mapProbe] [-mrna] [-cds] [-intron] [-t template] [-k keyset]\n") ;
  fprintf (stderr, "// Example:  boubou locusid 384D8-2 \n") ;
  fprintf (stderr, "//   -single: Expect in input file: probe manip signal motif, no antiprobe \n") ;
  fprintf (stderr, "//   -histo: gives the value > 95/100 all probes per chip \n") ;
  fprintf (stderr, "//   -histoN: gives the counts for the above histos\n") ;
  fprintf (stderr, "//   -maqc: expect a p2s_maqc.txt input file\n") ;
  fprintf (stderr, "//   -bbNmers: probe anti probe 6mers sensitivity\n") ;
  fprintf (stderr, "//   -manip1: all chips denoted G1*=1\n") ;
  fprintf (stderr, "//   -manip2: all chips denoted G1*=x x!= 1 except manip 3\n") ;
  fprintf (stderr, "//   -manip3: all chips denoted G1 with plate name 144*\n") ;
  fprintf (stderr, "//   -manip plateName: all chips with given *plateName*\n") ;
  fprintf (stderr, "//   -filter remove all probes with bad properties\n") ;
  fprintf (stderr, "//   -pivot export as a big table for excel pivot method \n") ;
  fprintf (stderr, "//   -green: only chips denoted _532\n") ;
  fprintf (stderr, "//   -red: only chips denoted _635\n") ;
  fprintf (stderr, "//   -pAmoins: only chips denoted pAmoins\n") ;
  fprintf (stderr, "//   -pAplus: only chips denoted pAplus\n") ;
  fprintf (stderr, "//   -expression: only expression chips\n") ;
  fprintf (stderr, "//   -S: only chips denoted S*\n") ;
  fprintf (stderr, "//   -G1: only chips denoted G1\n") ;
  fprintf (stderr, "//   -ratio: ration xx/G1\n") ;
  fprintf (stderr, "//   -getSolexaPeaks hard coded for now\n") ;
  fprintf (stderr, "//   -solexaAutoCorrel [-solexaSites file_name] measure the effective clone length\n") ;

  fprintf (stderr, "//   -out filename : redirects the output in that file, useful for batch jobs\n") ;

  exit (1) ;
}

/*************************************************************************************/

static int bbTest (void)
{
  char buf[7] ;
  char **cpp, *bp = "atgc" ;
  int ii, i, n[256], nbp[256] ;
  char *words[] = {
      "ttcgaa", "ttctaa", "ttagaa"
    , "cttcct", "aggaag"
    , "caggaa", "ttcctg"
    , "aggaaa", "tttcct"
    , "acttcc", "ggaagt" 
    , "tcctgt", "acagga"
    , "gcttcc", "ggaagc"
    , "gggcgg", "ccgccc"
    , 0,  0
  } ;
  
  strcpy (buf, "aaaaaa") ; 

  memset (n,0,sizeof(n)) ;
  memset (nbp,0,sizeof(nbp)) ;

  for (ii = 0 ; ii < 10000000 ; ii++)
    {
      for (i=0;i<5;i++) buf[i]=buf[i+1] ;
      i = randint () % 4 ;
      buf[5] =  bp[i] ;
      nbp[i]++ ;
      for (i=0, cpp = words ; *cpp; i++, cpp++)
	if (!strncmp (buf,*cpp, 6)) n[i]++ ;
    }
  for (i=0 ; i < 4 ; i++)
    printf ("\t%c=%d", bp[i], nbp[i]) ;
  for (i=0, cpp = words ; *cpp; i++, cpp++)
    printf("\n%s\t%d", *cpp, n[i]) ;
  printf("\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  const char *outfilename = 0 ;
  const char *dbName = 0, *ccp ;
  char buf100[100] ;
  int outlevel = 0 ;
  MM mm ;
  
  if(0) bbTest () ;
  memset (&mm, 0, sizeof (MM)) ;
  mm.filtered = bitSetCreate (1000000, 0) ;
  mm.probeDict = dictHandleCreate (1000000, 0) ;

  /* consume optional args */
  getCmdLineOption (&argc, argv, "-out", &outfilename ) ;
  if (getCmdLineOption (&argc, argv, "-bbNmers", &ccp) &&
      sscanf(ccp, "%d", &(mm.NN)) == 1 &&
      mm.NN > 0)
    mm.bbNmers = TRUE ;
  mm.getSolexaPeaks = getCmdLineOption (&argc, argv, "-getSolexaPeaks", 0) ;
  mm.solexaAutoCorrel = getCmdLineOption (&argc, argv, "-solexaAutoCorrel", 0) ;
  getCmdLineOption (&argc, argv, "-solexaSites",  &mm.solexaSitesFileName) ;
  mm.single = getCmdLineOption (&argc, argv, "-single", 0) ;
  mm.maqc = getCmdLineOption (&argc, argv, "-maqc", 0) ;
  mm.histo = getCmdLineOption (&argc, argv, "-histo", 0) ;
  mm.histoN = getCmdLineOption (&argc, argv, "-histoN", 0) ;
  if (mm.histoN)
    mm.histo = TRUE ;
  if(mm.histo)
    mm.single = TRUE ;
  mm.filter = getCmdLineOption (&argc, argv, "-filter", 0) ;
  mm.pivot = getCmdLineOption (&argc, argv, "-pivot", 0) ;
  mm.G = getCmdLineOption (&argc, argv, "-G", 0) ;
  mm.GC = getCmdLineOption (&argc, argv, "-GC", 0) ;
  mm.GA = getCmdLineOption (&argc, argv, "-GA", 0) ;
  mm.GT = getCmdLineOption (&argc, argv, "-GT", 0) ;
  mm.AT = getCmdLineOption (&argc, argv, "-AT", 0) ;
  mm.AC = getCmdLineOption (&argc, argv, "-AC", 0) ;
  mm.TC = getCmdLineOption (&argc, argv, "-TC", 0) ;

  mm.manip1 = getCmdLineOption (&argc, argv, "-manip1", 0) ;
  mm.manip2 = getCmdLineOption (&argc, argv, "-manip2", 0) ;
  mm.manip3 = getCmdLineOption (&argc, argv, "-manip3", 0) ;
  if (getCmdLineOption (&argc, argv, "-manip", &ccp))
    mm.manip = strnew (ccp, 0) ;

  mm.green = getCmdLineOption (&argc, argv, "-green", 0) ;
  mm.red = getCmdLineOption (&argc, argv, "-red", 0) ;

  mm.G1 = getCmdLineOption (&argc, argv, "-G1", 0) ;
  mm.pAplus = getCmdLineOption (&argc, argv, "-pAplus", 0) ;
  mm.S = getCmdLineOption (&argc, argv, "-S", 0) ;
  mm.ratio = getCmdLineOption (&argc, argv, "-ratio", 0) ;
  mm.expression = getCmdLineOption (&argc, argv, "-expression", 0) ;
  mm.pAmoins = getCmdLineOption (&argc, argv, "-pAmoins", 0) ;

  /* read absolute args */
  getCmdLineOption (&argc, argv, "-db", &dbName) ;

  if (mm.getSolexaPeaks)
    {
      const char *s ;
      if (!dbName) 
	usage () ;
      mm.db = ac_open_db (dbName, &s);
      if (!mm.db)
	messcrash ("Failed to open db %s, error %s", dbName, s) ;
    }

  if (outfilename)
    {
      f = filopen (outfilename, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  freeOutf ("# %s ", timeShow (timeNow(), buf100, 100)) ;

  if (mm.histoN) freeOutf (" histoN") ;
  else if (mm.histo) freeOutf (" histo") ;
  if (mm.ratio) freeOutf (" Ratio") ;
  if (mm.G1) freeOutf (" G1") ;
  if (mm.pAplus) freeOutf (" pA+") ;
  if (mm.expression) freeOutf (" Expression") ;
  if (mm.S) freeOutf (" S") ;
  if (mm.pAmoins) freeOutf (" NucPA-") ;

  if (mm.single) freeOutf (" Single probe") ;

  if (mm.manip1) freeOutf (" manip1 ") ;
  if (mm.manip2) freeOutf (" manip2 ") ;
  if (mm.manip3) freeOutf (" manip3 ") ;
  if (mm.manip) freeOutf (" manip %s ", mm.manip) ;

  if (mm.green) freeOutf (" just green:532 ") ;
  if (mm.red) freeOutf (" just red:635 ") ;

  if (!mm.single)
    {
      if (mm.G) freeOutf ("ATC mapped to n ") ;
      if (mm.GC) freeOutf ("AT mapped to n ") ;
      if (mm.GA) freeOutf ("CT mapped to n ") ;
      if (mm.GT) freeOutf ("CA mapped to n ") ;
      if (mm.AT) freeOutf ("GC mapped to n ") ;
      if (mm.AC) freeOutf ("GT mapped to n ") ;
      if (mm.TC) freeOutf ("GA mapped to n ") ;
    }

  freeOutf ("\n") ;

  if (mm.histo)
    bbHisto (&mm) ;
  else if (mm.bbNmers)
    bbNmers (&mm) ;
  else if (mm.getSolexaPeaks)
    bbGetSolexaPeaks (&mm) ;
  else if (mm.solexaAutoCorrel)
    bbSolexaAutocorrel (&mm) ;
  else
    usage () ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  ac_db_close (mm.db) ;

  sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

