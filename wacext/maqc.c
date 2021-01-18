/*  File: maqc.c
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
 * This file is part of the MAQC project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *       and Damir Herman, oct-nov 2005
 *
 *  This code works as a client against an acedb server
 *  and should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  The present code and the acedb schema for this project
 *  are available from mieg@ncbi.nlm.nih.gov
 */

#include "../wac/ac.h"
#include "keyset.h"
#include <errno.h>
#include <math.h>
#include "bitset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"

#define  LIMIT_PLATFORM 3
#define  LIMIT_AFX 4
#define NLAB 3
#define NTISSUE 4
#define NREPLICA 10

typedef enum {FLAGBOF=0, FLAG, NOFLAG} FLAGS ;
typedef enum {CENTREBOF=0, MEDIANE, MOYENNE, NOCENTRE} CENTERS ;
typedef enum {BASICBOF=0, BASIC, NOBASIC} BASICS ;
static void usage (void) ;
typedef struct maqcStruct {
  AC_HANDLE h ;
  AC_DB db ;
  DICT *dict ;
  DICT *dictNM ;
  DICT *dictGene ;
  DICT *dictGeneId ; 
  DICT *dictBasicGene ;
  KEYSET nm2geneId, titratingGeneAb, titratingGeneaB, titratingGeneAny ;
  KEYSET testedBy3 ;
  KEYSET testedGeneId, signalGeneId ;
  KEYSET testedGene, signalGene ;
  const char *template ;
  const char *keysetName ;
  const char *labName ;
  BOOL useNm, random ;
  BOOL c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, d1, allAgree ;
  int damir ;
  FLAGS flag ;
  CENTERS centre ;
  BASICS basic ;
  Array aa ;
  double s0, m0, nSigma ;
  FILE *aceSignalOut ;
} MAQC ;

typedef struct ttStruct {
  int A_1n, A_2n, A_3n, 
    B_1n, B_2n, B_3n, 
    C_1n, C_2n, C_3n, 
    D_1n, D_2n, D_3n ;
  int A_1nn[NREPLICA], A_2nn[NREPLICA], A_3nn[NREPLICA], 
    B_1nn[NREPLICA], B_2nn[NREPLICA], B_3nn[NREPLICA], 
    C_1nn[NREPLICA], C_2nn[NREPLICA], C_3nn[NREPLICA], 
    D_1nn[NREPLICA], D_2nn[NREPLICA], D_3nn[NREPLICA] ;
  double A_1mm[NREPLICA], A_2mm[NREPLICA], A_3mm[NREPLICA], 
    B_1mm[NREPLICA], B_2mm[NREPLICA], B_3mm[NREPLICA], 
    C_1mm[NREPLICA], C_2mm[NREPLICA], C_3mm[NREPLICA], 
    D_1mm[NREPLICA], D_2mm[NREPLICA], D_3mm[NREPLICA] ;
  /* ATTENTION : we use a union-like struct starting at A_1m in the parser */
  double A_1m, A_1[NREPLICA] ; /* tissue A, lab 1, mean then replica 0-4 */
  double A_2m, A_2[NREPLICA] ; /* tissue A, lab 2, mean then replica 0-4 */
  double A_3m, A_3[NREPLICA] ; /* tissue A, lab 3, mean then replica 0-4 */
  double B_1m, B_1[NREPLICA] ; /* tissue B, lab 1, mean then replica 0-4 */
  double B_2m, B_2[NREPLICA] ; /* tissue B, lab 2, mean then replica 0-4 */
  double B_3m, B_3[NREPLICA] ; /* tissue B, lab 3, mean then replica 0-4 */
  double C_1m, C_1[NREPLICA] ; /* tissue C, lab 1, mean then replica 0-4 */
  double C_2m, C_2[NREPLICA] ; /* tissue C, lab 2, mean then replica 0-4 */
  double C_3m, C_3[NREPLICA] ; /* tissue C, lab 3, mean then replica 0-4 */
  double D_1m, D_1[NREPLICA] ; /* tissue D, lab 1, mean then replica 0-4 */
  double D_2m, D_2[NREPLICA] ; /* tissue D, lab 2, mean then replica 0-4 */
  double D_3m, D_3[NREPLICA] ; /* tissue D, lab 3, mean then replica 0-4 */
} TT ;

typedef struct mmStruct {
  int probe ; /* value in the dict */
  int nm ;    /* value in the dictNM */
  int nm2 ;   /* value in the dictNM */
  int geneId ; /* value in the dictGeneId , zero if grade > 2 */
  int gene ;  /* value in the dictGene, zero if grade > 2  */
  BOOL flag ; /* true if flagged in Damir's file */
  BOOL valid ; /* true for validated titrating probes */
  BOOL basic ; /* true if a gene used to establish the base line of the chip */
  int probe_class ; /* 1 excellent, 2 ok, 3 multi, 4 antisense, 5 genome */
  KEY id ;
  char lab[4], bench[4] ;
  double tm, score ;
  double signal, foldChange ;
  BOOL titrating, titratable, cotitrating, nontitrating, antititrating ;
  int A_1n, A_2n, A_3n, 
    B_1n, B_2n, B_3n, 
    C_1n, C_2n, C_3n, 
    D_1n, D_2n, D_3n ;
  /* ATTENTION : we use a union-like struct starting at A_1m in the parser */
  double A_1m,  A_2m, A_3m ;
  double B_1m,  B_2m, B_3m ;
  double C_1m,  C_2m, C_3m ;
  double D_1m,  D_2m, D_3m ;

  double A_1[NREPLICA] ; /* tissue A, lab 1, mean then replica 0-4 */
  double A_2[NREPLICA] ; /* tissue A, lab 2, mean then replica 0-4 */
  double A_3[NREPLICA] ; /* tissue A, lab 3, mean then replica 0-4 */
  double B_1[NREPLICA] ; /* tissue B, lab 1, mean then replica 0-4 */
  double B_2[NREPLICA] ; /* tissue B, lab 2, mean then replica 0-4 */
  double B_3[NREPLICA] ; /* tissue B, lab 3, mean then replica 0-4 */
  double C_1[NREPLICA] ; /* tissue C, lab 1, mean then replica 0-4 */
  double C_2[NREPLICA] ; /* tissue C, lab 2, mean then replica 0-4 */
  double C_3[NREPLICA] ; /* tissue C, lab 3, mean then replica 0-4 */
  double D_1[NREPLICA] ; /* tissue D, lab 1, mean then replica 0-4 */
  double D_2[NREPLICA] ; /* tissue D, lab 2, mean then replica 0-4 */
  double D_3[NREPLICA] ; /* tissue D, lab 3, mean then replica 0-4 */
} MM ;

#define ABI     0x001
#define AFX_123 0x002
#define AFX_456 0x004
#define AGL     0x008
#define AGL_cy3 0x008
#define AGL_cy5 0x010
#define GEH     0x020 
#define ILM     0x040
#define NCI     0x080
#define GEX     0x100
#define QGN     0x200
#define TAQ     0x400
#define EPP     0x800
#define AFX_5  0x1000
#define AFX_6  0x2000
#define AFX_523  0x4000

typedef struct labostruct { 
  KEY id ;
  char lab[4], bench[4] ; 
  int nReplica, nLab, nTissue, jumpCol, pattern, nn ;
  int nSignal ;
  double meanSignal, sigma ;
  FLAGS flag ;
  CENTERS centre ;
  BASICS basic ;
  KEYSET gene2mmAb ;
  KEYSET gene2mmaB ;
  int nGeneTestedValidatedAny ;
} LAB ;

/*       use -damir (DAMIR) switch to use the DAMIR processed (non raw) data
 * ABI, raw, moyenne , 22     (previousy we used method=12 jumpcol=0)
 * AFX, raw, moyenne , 20
 * AGL, raw, mediane , 21 ( previously 1 in 2 colors, 2 in 1 color, jumpcol=0
 * GEH, raw, moyenne , 20 (previousy we used method=15 jumpcol=2)
 * ILM, raw, moyenne , 20      (was mediane before jun 20, previousy we used method=14 jumpcol=0)
 * NCI, processed, mediane,  (they used lowess ratio jumpcol=2)

 * EPP, no-flag, moyenne
 * TAQ, no-flag, no-centre
 * GEX, no-flag, no-centre
 * QGN, FLAG, no-centre
 */
/* there are 36M data points in AFX, 12M in others, total 48M data points */
LAB labos [] = {
  { AFX_456, "AFX", "456", 5,3,4,2,27, 604258, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC,0,0,0}, /* just parse afx 456 for the moment */

  { AFX_123, "AFX", "123", 5,3,4,2,26, 604258, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0}, /* just parse afx 123 for the moment */

  { AFX_523, "AFX", "523", 5,3,4,2,235, 604258, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0}, /* just parse afx 123 for the moment */

  { ABI, "ABI", "", 5,3,4, /* 0,12 */ 2, 22, 32878, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0} , /* flag < 3 => non detectable , A_2[5], B_2[3], D_3[2] un peu moins bonnes  et les 2 sites  D_1[1] et B_1[5] sont horribles */
  { AGL, "AGL", "", 5,3,4, /* 0, 2 */ 2,21 , 41000, 0, 0.0, 0.0, NOFLAG, MEDIANE, NOBASIC, 0,0,0} ,/* flag P=path, F=failed; ignorer  A_2[1]; C_2[4] */
  { GEH, "GEH", "", 5,3,4, /* 2,15 */ 2,24, 55038, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0} ,  /* flag G=good, L=absent; SMIC->ignore; ignorer le site 2 */
  { ILM, "ILM", "", 5,3,4, /* 0,14 */ 2,25, 47293, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0} ,    /* signal > 0.99 = present ; < 0.99 = absent */
  { NCI, "NCI", "", 5,3,4,2,3 /* was 3 */, 35237, 0, 0.0, 0.0, NOFLAG, MEDIANE, NOBASIC, 0,0,0} ,    /* ignore A_2[3], A_3[3], B_3[3] ; ignorer tout NCI */

  { EPP, "EPP", "", 5,3,4,0,16,   294, 0, 0.0, 0.0, NOFLAG, MOYENNE, NOBASIC, 0,0,0} , /* EPP-NOFLAG has more validated 132 versus 101 and a higer rate 98.7 versus 98.1 than with flags */ 
  { GEX, "GEX", "", 3,1,4,0,13,   207, 0, 0.0, 0.0, NOFLAG, NOCENTRE, NOBASIC, 0,0,0} , /* GEX has no flags */
  { QGN, "QGN", "", 3,1,4,1,11,   245, 0, 0.0, 0.0, FLAG, NOCENTRE, NOBASIC, 0,0,0} ,  /* FLAG is better, 9 more confirmed genes */
  { TAQ, "TAQ", "", 4,1,4,3,10,  1004, 0, 0.0, 0.0, NOFLAG, NOCENTRE, NOBASIC, 0,0,0} ,  /* 3 jump both col 1 and col 3, flag detruit des bons genes */

  { 0x000, "", "",    0,0,0,0,0,    0, 0, 0.0, 0.0, NOFLAG, NOCENTRE, NOBASIC, 0,0,0} ,
  { AFX_6, "AFX", "456", 5,3,4,0,26, 54675, 0, 0.0, 0.0, NOFLAG, MOYENNE, BASIC, 0,0,0}, /* just parse afx 6 for the moment */

  { AGL_cy3, "AGL", "cy3", 10,3,2,2,1, 41675, 0, 0.0, 0.0, NOFLAG, NOCENTRE, NOBASIC, 0,0,0} ,/* flag P=path, F=failed; ignorer  A_2[1]; C_2[4] */
  { AGL_cy5, "AGL", "cy5", 10,3,2,2,2, 41675, 0, 0.0, 0.0, NOFLAG, NOCENTRE, NOBASIC, 0,0,0} , /* flag P=path, F=failed; ignorer  A_2[1]; C_2[4] */
  { AFX_123, "AFX", "123", 5,3,4,0,0, 54675, 0, 0.0, 0.0, NOFLAG, NOCENTRE, BASIC, 0,0,0} ,
  { AFX_456, "AFX", "456", 5,3,4,0,5, 54675, 0, 0.0, 0.0, NOFLAG, NOCENTRE, BASIC, 0,0,0} , /* ignorer le site 6 */
} ;
/*
  GEH comparison of the different files  
                 #probes A&B      titrating     tit/acv/grade12   geneid          follow geneid     
  raw flagged     41117           25543           21732           14175           12316
  raw no flag     52356           24261           20712           12985           11314
  leming flagged  52419           23979           20548           12872           11196
*/

/* abi afx agl geh ilm nci*/
/* ABI | AFX_123 | AFX_456 | AGL_cy3 | AGL_cy5 | GEH | ILM | NCI */
static KEY big6 = ABI | AFX_123 | AGL  | GEH | ILM | NCI ;

/*************************************************************************************/

static BOOL maqcGetTxtDataLab (MAQC *maqc, LAB *lab)
{
  char cutter, *cp ;
  int i, level, line = 0, nGeneTested = 0 ;
  int nTissue, nLab, nReplica ;
  int probe, nn = 0 ; /* number of probe results in the table */
  MM *mm ;
  double *zp ;
  BOOL ok ;
  FILE *f=0 ;
  int DAMIR ;
  BOOL flag = FALSE, JEAN_FLAG, DAMIR_FLAG ;
  KEYSET geneTested = keySetCreate () ;

  if (maqc->flag == NOFLAG) ;
  else if (maqc->flag == FLAG || lab->flag  == FLAG)
    flag = TRUE ;
  
  JEAN_FLAG = DAMIR_FLAG = flag ;
  DAMIR = maqc->damir ;
  
  if (DAMIR == 0)
    f = filopen (messprintf ("DATA2/%s%s", lab->lab, lab->bench), "txt", "r") ;
  else  if (DAMIR == 1)
    {
      switch (lab->pattern)
	{
	case 12: lab->pattern = 22 ; break ; /* ABI */
	case  2: lab->pattern = 21 ; break ; /* AGL */
	case 15: lab->pattern = 24 ; break ; /* GEH */
	case 14: lab->pattern = 25 ; break ; /* ILM */
	case  3: lab->pattern = 23 ; break ; /* NCI */

	}
      lab->pattern = lab->pattern == 3 ? 23 : 20 ; 
      lab->jumpCol = 2 ;
      f = filopen (messprintf ("DATA3/%s%s_data", lab->lab, lab->bench), "txt", "r") ;
    }
  else if (DAMIR == 2)
    {
      if (lab->pattern == 3)
	{
	  JEAN_FLAG = DAMIR_FLAG = FALSE ;
	  lab->pattern = 23 ;
	}
      else
	lab->pattern = 20 ; 
      lab->jumpCol = 2 ;
      f = filopen (messprintf ("DATA4/%s%s_raw_data", lab->lab, lab->bench), "txt", "r") ;
    }
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  freespecial ("\n") ; /* so we can correctly scan the \t */
  ok = TRUE ; nn = 0 ;
  freecard (level) ; /* jump caption line */
  while (freecard (level))
    {
      ok = FALSE ;
      line++ ;

      /* parse the probe name */
      if (lab->jumpCol == 3 || lab->jumpCol == 1)
	{
	  cp = freewordcut ("\t",&cutter) ;
	  if (!cp) continue ;
	}
      
      cp = freewordcut ("\t",&cutter) ;
      if (!cp) continue ;   
      if (0 && strcmp (cp, "H200001924"))
	continue ;
      if (strncmp(lab->lab,cp,3))
	dictAdd (maqc->dict, messprintf ("%s_%s|%s", lab->lab, cp, lab->bench),  &probe) ;
      else
	dictAdd (maqc->dict, messprintf ("%s|%s", cp, lab->bench),  &probe) ;
      mm = arrayp (maqc->aa, probe, MM) ;
      mm->id = lab->id ;
      strcpy (mm->lab, lab->lab) ;
      strcpy (mm->bench, lab->bench) ;
      mm->tm = mm->score = 0 ;
      /* parse the probe name */
      if (lab->jumpCol == 2 || lab->jumpCol == 3) 
	{
	  cp = freewordcut ("\t",&cutter) ;
	  if (!cp && !cutter) continue ;
	}
       
      /* parse and discard the spot id */
      cp = freewordcut ("\t",&cutter) ;
      if (!cp && !cutter) continue ;
      
      zp = &(mm->A_1[0]) ; ok =  TRUE ;
      switch (lab->pattern)
	{
	case 10: /* TAQ's format */
	  {
	    double b1[120] ; /* 4 replica, 4 tissue, 6 labs */
	    int i ;

	    for (i=0 ; ok && i < 16 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter)) /* negative log is ok */
		b1[i] = UT_NON_DOUBLE ;

	    for (i=0 ; JEAN_FLAG && ok && i < 16 ; i++)
	      if ((cp = freewordcut ("\t",&cutter)))
		{
		  if (*cp == 'A')
		    {
		      if (JEAN_FLAG)
			b1[i] = UT_NON_DOUBLE ;
		    }
		}

	    for (i=0 ; ok && i < 4 ; i++)
	      {
		mm->A_1[i] = b1 [0+i] ;     
		mm->B_1[i] = b1 [4+i] ;     
		mm->C_1[i] = b1 [8+i] ;     
		mm->D_1[i] = b1 [12+i] ;     
	      }	 
	  }
	  break ;
	case 11: /* QGN's format */
	  {
	    double b1[120] ; /* 4 replica, 4 tissue, 6 labs */
	    int i ;

	    for (i=0 ; ok && i < lab->nLab * lab->nReplica * lab->nTissue  ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter)) /* negative log is ok */
		b1[i] = UT_NON_DOUBLE ;
	      else ; /* we are already in logs */

	    for (i=0 ; JEAN_FLAG && ok &&  i < lab->nLab * lab->nReplica * lab->nTissue  ; i++)
	      if ((cp = freewordcut ("\t",&cutter)))
		{
		  if (*cp == 'A' && b1[i] == UT_NON_DOUBLE)
		    b1[i] = -5.0 ; 
		}
	      else
		b1[i] = UT_NON_DOUBLE ;	/* no flag is a non measured spot */

	    for (i=0 ; ok && i < lab->nReplica ; i++)
	      {
		mm->A_1[i] = b1 [0 * lab->nReplica + i] ;
		mm->B_1[i] = b1 [1 * lab->nReplica + i] ;     
		mm->C_1[i] = b1 [2 * lab->nReplica + i] ;
		mm->D_1[i] = b1 [3 * lab->nReplica + i] ;
	      }	 
	  }
	  break ;
	case 12: /* ABI's format */
	case 14: /* ILM's format */
	  {
	    double b1[120] ; /* 4 replica, 4 tissue, 6 labs */
	    double flags[120] ;
	    int i ;
	    int iReplica, iLab, iTissue ;
	    
	    for (i=0 ; ok && i < lab->nLab * lab->nReplica * lab->nTissue  ; i++)
	      {
		if (lab->pattern == 14 && i == 2 + 2 * lab->nLab * lab->nReplica)
		  b1[i] = UT_NON_DOUBLE ; /* 1C3 :: is missing in ILM */
		else if (! freedoublecut (&b1[i], "\t", &cutter)) /* negative log is ok */
		  b1[i] = UT_NON_DOUBLE ;
	      }
	    for (iLab = 0 ; JEAN_FLAG && iLab < lab->nLab ; iLab++)
	      for (iTissue = 0 ;  iTissue < lab->nTissue  ; iTissue++)
		for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
		  switch (lab->pattern)
		    {
		    case 12:
		      if (! freedoublecut (&flags[i], "\t", &cutter) ||
			  flags[i] < 3.0)
			b1[iReplica + iLab * lab->nReplica + iTissue * lab->nLab * lab->nReplica ] 
			  = UT_NON_DOUBLE ;	/* no flag is a non measured spot */
		      break ;
		    case 14:
		      if (iLab == 0 && iTissue == 2 && iReplica  == 2)               /* ILM C_1[3] */
			continue ;
		      if (! freedoublecut (&flags[i], "\t", &cutter) ||
			  flags[i] < .99)
			b1[iReplica + iLab * lab->nReplica + iTissue * lab->nLab * lab->nReplica ] 
			  = UT_NON_DOUBLE ;	/* no flag is a non measured spot */
		      break ;
		    }

	    for (i=0 ; ok && i < lab->nReplica ; i++)
	      {
		mm->A_1[i] = b1 [0 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->B_1[i] = b1 [1 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;     
		mm->C_1[i] = b1 [2 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->D_1[i] = b1 [3 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;

		mm->A_2[i] = b1 [0 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;
		mm->B_2[i] = b1 [1 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;     
		mm->C_2[i] = b1 [2 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;
		mm->D_2[i] = b1 [3 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;

		mm->A_3[i] = b1 [0 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->B_3[i] = b1 [1 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;     
		mm->C_3[i] = b1 [2 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->D_3[i] = b1 [3 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
	      }	 
	    if (lab->pattern == 12) /* ABI A_2[5] D_3[2] */
	      {
		mm->A_2[4] = mm->D_3[1] = UT_NON_DOUBLE ;
		if (1) /* also remove according to danielle */
		  mm->B_2[2] = mm->B_3[2]  = mm->B_3[4] = UT_NON_DOUBLE ;
	      }
	    if (lab->pattern == 14) /* ILM */
	      { 
		mm->B_1[2] = UT_NON_DOUBLE ; /* also remove according to danielle */
	      }
	  }
	  break ;
	case 13: /* GEX's format */
	  {
	    double b1[120] ; /* 4 replica, 4 tissue, 6 labs */
	    int i ;

	    for (i=0 ; ok && i < lab->nLab * lab->nReplica * lab->nTissue  ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter)) /* negative log is ok */
		b1[i] = UT_NON_DOUBLE ;
	      else ; /* we are already in logs */

	      /* this is GOD talking, there are no flags */

	    for (i=0 ; ok && i < lab->nReplica ; i++)
	      {
		mm->A_1[i] = b1 [0 * lab->nReplica + i] ;
		mm->B_1[i] = b1 [1 * lab->nReplica + i] ;     
		mm->C_1[i] = b1 [2 * lab->nReplica + i] ;
		mm->D_1[i] = b1 [3 * lab->nReplica + i] ;
	      }	 
	  }
	  break ;
	case 27: /* afx6 damir's format */
	case 26: /* afx123 damir's format */
	case 25: /* ILM damir raw format */
	case 24: /* GEH damir raw format */
	case 23: /* NCI damir raw format */
	case 22: /* ABI damir's format */
	case 21: /* AGL damir's format */
	case 20: /* damir's format */
	  {
	    double flag, b1[60] ; /* 5 replica, 4 tissue, 63labs */
	    BOOL f[60] ;
	    int i ;
	    /* freewordcut ("\t",&cutter) ;  jump the position on the slide */
	    for (i=0 ; ok && i < 60 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter) || b1[i] <= 0)
		b1[i] = UT_NON_DOUBLE ;
	      else
		b1[i] =  lab->pattern == 23 ? b1[i] : log2 (b1[i]) ;

	    for (i=0 ; ok && i < 60 ; i++)
	      {
		if (! freedoublecut (&flag, "\t", &cutter) || flag != 1)
		  {
		    f[i] = TRUE ;
		    if (DAMIR_FLAG) b1[i] = UT_NON_DOUBLE ;
		  }
		else
		  f[i] = FALSE ;
	      }	
	    for (i=0 ; ok && i < 5 ; i++)
	      {
		if (lab->pattern == 27)
		  {
		    mm->A_1[i] = UT_NON_DOUBLE ;
		    mm->B_1[i] = UT_NON_DOUBLE ;
		    mm->C_1[i] = UT_NON_DOUBLE ;
		    mm->D_1[i] = UT_NON_DOUBLE ;
		    
		    mm->A_2[i] = UT_NON_DOUBLE ;
		    mm->B_2[i] = UT_NON_DOUBLE ;
		    mm->C_2[i] = UT_NON_DOUBLE ;
		    mm->D_2[i] = UT_NON_DOUBLE ;

		    mm->A_3[i] = b1 [40+i] ;     
		    mm->B_3[i] = b1 [45+i] ;     
		    mm->C_3[i] = b1 [50+i] ;     
		    mm->D_3[i] = b1 [55+i] ;     
		  }
		else
		  {
		    mm->A_1[i] = b1 [0+i] ;     
		    mm->B_1[i] = b1 [5+i] ;     
		    mm->C_1[i] = b1 [10+i] ;     
		    mm->D_1[i] = b1 [15+i] ;     
		    
		    mm->A_2[i] = b1 [20+i] ;     
		    mm->B_2[i] = b1 [25+i] ;     
		    mm->C_2[i] = b1 [30+i] ;     
		    mm->D_2[i] = b1 [35+i] ;     
		    
		    mm->A_3[i] = b1 [40+i] ;     
		    mm->B_3[i] = b1 [45+i] ;     
		    mm->C_3[i] = b1 [50+i] ;     
		    mm->D_3[i] = b1 [55+i] ;     
		  }
	      }
	    /* 2006/04/20: AGL asks that we remove 
	     * A_1[1]  A_2[3] D_2[2] B_3[3] but actually we only remove 2 of them
	     */ 
	    if (lab->pattern == 21) /* AGL */
	     {
	       if (0)
		 mm->A_1[0] = mm->A_2[2] = mm->D_2[1] = mm->B_3[2] = UT_NON_DOUBLE ;
	       else /* on remet A_1[1]  B_3[3]  jun 19 */
		 mm->A_2[2] = mm->D_2[1] = UT_NON_DOUBLE ;
	       mm->A_1[0] = UT_NON_DOUBLE ;
	     }
	    if (lab->pattern == 23) /* NCI */
	      mm->A_2[2] = mm->A_3[2] = mm->B_3[2] =
		mm->A_2[4] = mm->B_2[4] = mm->C_2[3] = 
		mm->D_2[3] = mm->C_2[4] = mm->D_2[4] = UT_NON_DOUBLE ;
	    if (lab->pattern == 22) /* ABI */
	      mm->A_2[4] = mm->D_3[1] = UT_NON_DOUBLE ;

	    if (1) /* also remove according to danielle */
	      {
		if (lab->pattern == 21) /* AGL  b2.* c2.1345 D2.1 */ 
		  mm->B_2[0] = mm->B_2[1] = mm->B_2[2]  = mm->B_2[3] = mm->B_2[4] = 
		  mm->C_2[0] = mm->C_2[2]  = mm->C_2[3] = mm->C_2[4] = 
		  mm->D_2[0] = UT_NON_DOUBLE ;
		if (lab->pattern == 23) /* NCI  b2.3  b3.1245  c3.125 d3.* */
		  mm->B_2[2] = 
		    mm->B_3[0] = mm->B_3[1] = mm->B_3[3] = mm->B_3[4] = 
		    mm->C_3[0] = mm->C_3[1] = mm->C_3[4] = 
		    mm->D_3[0] = mm->D_3[1] = mm->D_3[2]  = mm->D_3[3] = mm->D_3[4] = UT_NON_DOUBLE ;
		if (lab->pattern == 22) /* ABI b2.3 b3.3 b3.5 */
		  mm->B_2[2] = mm->B_3[2]  = mm->B_3[4] = UT_NON_DOUBLE ;
		if (lab->pattern == 24) /* GEH ok */
		  ;
		if (lab->pattern == 25) /* ILM b1.3 */
		  mm->B_1[2] = UT_NON_DOUBLE ;
		if (lab->pattern == 26) /* AFX a1.123 */
		  mm->A_1[0] = mm->A_1[1] = mm->A_1[2] = UT_NON_DOUBLE ;
	      }

	    /* a probe is flagged if all its values are flagged */
	    mm->flag = TRUE ;
	    for (i = 0 ; mm->flag && i < 60 ; i++)
	      if (b1[i] != UT_NON_DOUBLE && !f[i])
		mm->flag = FALSE ;
	  }
	  break ;
	case 523: /* joined damir's format */
	  {
	    double flag, b1[160 + NREPLICA] ; /* 5 replica, 4 tissue, 6 labs */
	    int i ;
	    /* freewordcut ("\t",&cutter) ;  jump the position on the slide */
	    for (i=0 ; ok && i < 60 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter) || b1[i] <= 0)
		b1[i] = UT_NON_DOUBLE ;
	      else
		b1[i] = log2 (b1[i]) ;
	    for (i=0 ; ok && i < 60 ; i++)
	      if (! freedoublecut (&flag, "\t", &cutter) || flag != 1)
		if (DAMIR_FLAG) b1[i] = UT_NON_DOUBLE ;
	    freewordcut ("\t",&cutter) ;
	    freewordcut ("\t",&cutter) ;
	    for (i=60 ; ok && i < 120 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter) || b1[i] <= 0)
		b1[i] = UT_NON_DOUBLE ;
	      else
		b1[i] = log2 (b1[i]) ;
	    for (i=60 ; ok && i < 120 ; i++)
	      if (! freedoublecut (&flag, "\t", &cutter) || flag != 1)
		if (DAMIR_FLAG) b1[i] = UT_NON_DOUBLE ;

	    for (i=0 ; ok && i < 5 && i < NREPLICA ; i++)
	      {
		mm->A_1[i] = b1 [140+0+i] ;      /* really lab 5 */ 
		mm->B_1[i] = b1 [140+5+i] ;     
		mm->C_1[i] = b1 [140+10+i] ;     
		mm->D_1[i] = b1 [140+15+i] ;     
		
		mm->A_2[i] = b1 [20+i] ;     
		mm->B_2[i] = b1 [25+i] ;     
		mm->C_2[i] = b1 [30+i] ;     
		mm->D_2[i] = b1 [35+i] ;     
		
		mm->A_3[i] = b1 [40+i] ;     
		mm->B_3[i] = b1 [45+i] ;     
		mm->C_3[i] = b1 [50+i] ;     
		mm->D_3[i] = b1 [55+i] ;     
	      }
	  }
	  /* mm->A_1[0] = UT_NON_DOUBLE ;  artificial test for AGL */
	  break ;
	case 0:
	  for (i = nTissue = 0 ; nTissue < lab->nTissue ; nTissue++)
	    for (nLab = 0 ; nLab < lab->nLab ; nLab++)
	      for (nReplica = 0; ok && nReplica < lab->nReplica ; i++, zp++, nReplica++)
		{ /* initial zp++ jumps the X_xm */
		  *zp = 0 ;
		  if (nReplica < lab->nReplica)
		    if (! freedoublecut (zp, "\t",&cutter))
		      ok = FALSE ;
		}
	  break ;
	case 5: /* afx 456 eliminate site 6 */
	  for (i = nTissue = 0 ; nTissue < lab->nTissue ; nTissue++)
	    for (nLab = 0 ; nLab < lab->nLab ; nLab++)
	      for (nReplica = 0 ; nReplica < NREPLICA ; i++, zp++, nReplica++)
		{ /* initial zp++ jumps the X_xm */
		  *zp = 0 ;
		  if (nReplica < lab->nReplica)
		    if (! freedoublecut (zp, "\t",&cutter))
		      ok = FALSE ; 
		  /* having read all the data we eliminate site 6 */
		  if (nLab == 2)
		    *zp = UT_NON_DOUBLE ;
		}
	  break ;
	case 1:
	  /* Agilent */
	  /* experiments A and B are tissueA/tissueA, then tissueB/tissueB, dye swapped */
	  /* experiment C is [A]cy3 then [B]cy5, experiment D is [A]cy5 then [B]cy3 */
	  /* so in fact we have 4 measure of tissue [A]: Aa Ab Ca Db (respectivelly cy3 cy5 cy3 cy5 */
	  /* so in fact we have 4 measure of tissue [B]: Ba Bb Cb Da (respectivelly cy3 cy5 cy5 cy3 */
	  /* case 1: we consider AGL/colorant Cy3 */
	  /* case 2: we consider AGL/colorant Cy5 */
	  /* use col 1,3,5,7,9,21,23,25,27,29,41,43,45,47,49 for tissue A, replica 12345 lab 1 2 3*/
	  /* use col 61... for tissue B*/
	  /* ignore A_2[0] and C_2[3] */
	  {
	    double b1[120] ;
	    int i ;
	    
	    for (i=0 ; ok && i < 120 ; i++)
	      if (! freedoublecut (&b1[i], "\t",&cutter))
		ok = FALSE ;
	    for (i=0 ; ok && i < 5 ; i++)
	      {
		mm->A_1[i] = b1 [2*i] ;         /* experiment Aa, col  1, 3, 5, 7, 9, */
		mm->A_1[5+i] = b1 [10 + i] ;    /* experiment Ca, col 11,12,13,14,15 */
		mm->A_2[i] = b1 [20+2*i] ;      /* experiment Aa, col 21,23,25,27,29, */
		mm->A_2[5+i] = b1 [30 + i] ;    /* experiment Ca, col 31,32,33,34,35 */
		mm->A_3[i] = b1 [40+2*i] ;      /* experiment Aa, col 41,43,45,47,49 */
		mm->A_3[5+i] = b1 [50 + i] ;    /* experiment Ca, col 51,52,53,54,55 */
		
		mm->B_1[i] = b1 [60+2*i] ;      /* experiment Ba, col 61,63,65,67,69, */
		mm->B_1[5+i] = b1 [70 + i] ;    /* experiment Da, col 71,72,73,74,75 */
		mm->B_2[i] = b1 [80+2*i] ;      /* experiment Ba, col 81,83,85,87,89, */
		mm->B_2[5+i] = b1 [90 + i] ;    /* experiment Da, col 91,92,93,94,95 */
		mm->B_3[i] = b1 [100+2*i] ;     /* experiment Ba, col 101,103,105,107,109 */
		mm->B_3[5+i] = b1 [110 + i] ;   /* experiment Da, col 111,112,113,114,115 */
	      }
	    if (1) { mm->A_2[0] = mm->A_2[8] =  mm->B_2[8] = UT_NON_DOUBLE ; }
	  }
	  break ;
	case 2:
	  /* case 2: we consider AGL/colorant Cy5 */ 
	  /* ignore A_2[0] and C_2[3] */
	  {
	    double b1[120 + NREPLICA] ;
	    int i ;
	    
	    for (i=0 ; i < 120 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter))
		b1[i] = UT_NON_DOUBLE ;
	    for (i=0 ; ok && i < 5 && i < NREPLICA ; i++)
	      {
		mm->A_1[i] = b1 [1+2*i] ;       /* experiment Ab, col  2, 4, 6, 8,10 */
		mm->A_1[5+i] = b1 [15 + i] ;    /* experiment Db, col 16,17,18,19,20 */
		mm->A_2[i] = b1 [21+2*i] ;      /* experiment Ab, col 22,24,26,28,30 */
		mm->A_2[5+i] = b1 [35 + i] ;    /* experiment Db, col 36.37.38.39.40 */
		mm->A_3[i] = b1 [41+2*i] ;      /* experiment Ab, col 42,44,46,48,50 */
		mm->A_3[5+i] = b1 [55 + i] ;    /* experiment Db, col 56,57,58,59,60 */

		mm->B_1[i] = b1 [61+2*i] ;      /* experiment Bb, col 62,64,66,68,70 */
		mm->B_1[5+i] = b1 [75 + i] ;    /* experiment Cb, col 76,77,78,79,80 */
		mm->B_2[i] = b1 [81+2*i] ;      /* experiment Bb, col 82,84,86,88,90 */
		mm->B_2[5+i] = b1 [95 + i] ;    /* experiment Cb, col 96,97,98,99,100 */
		mm->B_3[i] = b1 [101+2*i] ;     /* experiment Bb, col 102,104,106,108,110 */
		mm->B_3[5+i] = b1 [116 + i] ;   /* experiment Cb, col 116,117,118,119,120 */
	      }
	    if (1) 
	      { 
	       if (0)
		 mm->A_1[0] = mm->A_2[2] = mm->D_2[1] = mm->B_3[2] = UT_NON_DOUBLE ;
	       else /* on remet A_1[1]  B_3[3]  jun 19 */
		 mm->A_2[2] = mm->D_2[1] = UT_NON_DOUBLE ;

	       if (1) /* also remove according to danielle */
		 /* AGL  b2.* c2.1345 D2.1 */ 
		  mm->B_2[0] = mm->B_2[1] = mm->B_2[2]  = mm->B_2[3] = mm->B_2[4] = 
		  mm->C_2[0] = mm->C_2[2]  = mm->C_2[3] = mm->C_2[4] = 
		  mm->D_2[0] = UT_NON_DOUBLE ;
	      }
	  }
	  break ;
	case 15: /* GEH */
	  /*  GEH's original flag information (CODELINKFLAG) is also provided:
	      G = "Good"
	      L = "Absent" (below noise since the signal mean is less than background mean plus 
                  1.5 standard deviations of local background)
	      S = "Saturated"
	      M = "Manufacturing Defect"
	      I = "Irregular Shape"
	      C = "Contaminated Spot"
	      Only use spots which are flagged as "G" or "L". 
	  */
	  {
	    double b1[120] ;
	    int i, iLab , iTissue, iReplica ;
	    
	    for (i=0 ; i < lab->nLab * lab->nReplica * lab->nTissue ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter))
		b1[i] = UT_NON_DOUBLE ;
	    for (iLab = 0 ; JEAN_FLAG && iLab < lab->nLab ; iLab++)
	      for (iTissue = 0 ;  iTissue < lab->nTissue  ; iTissue++)
		for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
		  if ((cp = freewordcut ("\t",&cutter)))
		    {
		      switch (*cp)
			{
			case 'G':
			case 'L':
			  break ;
			case 'S':
			case 'M':
			case 'I':
			case 'C':
			default:
			  b1[iReplica + iLab * lab->nReplica + iTissue * lab->nLab * lab->nReplica ] 
			    = UT_NON_DOUBLE ;
			  break ;
			}
		    }
		  else
		    b1[iReplica + iLab * lab->nReplica + iTissue * lab->nLab * lab->nReplica ] 
		      = UT_NON_DOUBLE ;	/* no flag is a non measured spot */
	    for (i=0 ; ok && i < lab->nReplica ; i++)
	      {
		mm->A_1[i] = b1 [0 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->B_1[i] = b1 [1 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;     
		mm->C_1[i] = b1 [2 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->D_1[i] = b1 [3 * lab->nLab * lab->nReplica + 0 * lab->nReplica + i] ;

		mm->A_2[i] = b1 [0 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;
		mm->B_2[i] = b1 [1 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;     
		mm->C_2[i] = b1 [2 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;
		mm->D_2[i] = b1 [3 * lab->nLab * lab->nReplica + 1 * lab->nReplica + i] ;

		mm->A_3[i] = b1 [0 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->B_3[i] = b1 [1 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;     
		mm->C_3[i] = b1 [2 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->D_3[i] = b1 [3 * lab->nLab * lab->nReplica + 2 * lab->nReplica + i] ;
	      }	 
	  }
	  break ;
	case 16: /* EPP */
	  /*   
	      P = Present = acceptable
	      M = Marginal = saturated
	      A = Absent = low intensity
	  */
	  {
	    double b1[120] ;
	    int i ;
	    for (i=0 ; i < 7 ; i++) /* jump 7 junk fields */
	      freewordcut ("\t",&cutter) ;
	    for (i=0 ; i < lab->nLab * lab->nReplica * lab->nTissue ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter) ||
		  b1[i] <= 0.00001)
		b1[i] = UT_NON_DOUBLE ;
	      else
		b1[i] = log2(b1[i]) ;
	    for (i=0 ; JEAN_FLAG && ok &&  i < lab->nLab * lab->nReplica * lab->nTissue  ; i++)
	      if ((cp = freewordcut ("\t",&cutter)))
		{
		  switch (*cp)
		    {
		    case 'P':
		      break ;
		    case 'M':
		    case 'A':
		      b1[i] = UT_NON_DOUBLE ;
		      break ;
		    }
		}
	      else
		b1[i] = UT_NON_DOUBLE ;	/* no flag is a non measured spot */
	    for (i=0 ; ok && i < lab->nReplica ; i++)
	      {
		mm->A_1[i] = b1 [0 * lab->nTissue * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->B_1[i] = b1 [0 * lab->nTissue * lab->nReplica + 1 * lab->nReplica + i] ;     
		mm->C_1[i] = b1 [0 * lab->nTissue * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->D_1[i] = b1 [0 * lab->nTissue * lab->nReplica + 3 * lab->nReplica + i] ;

		mm->A_2[i] = b1 [1 * lab->nTissue * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->B_2[i] = b1 [1 * lab->nTissue * lab->nReplica + 1 * lab->nReplica + i] ;     
		mm->C_2[i] = b1 [1 * lab->nTissue * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->D_2[i] = b1 [1 * lab->nTissue * lab->nReplica + 3 * lab->nReplica + i] ;

		mm->A_3[i] = b1 [2 * lab->nTissue * lab->nReplica + 0 * lab->nReplica + i] ;
		mm->B_3[i] = b1 [2 * lab->nTissue * lab->nReplica + 1 * lab->nReplica + i] ;     
		mm->C_3[i] = b1 [2 * lab->nTissue * lab->nReplica + 2 * lab->nReplica + i] ;
		mm->D_3[i] = b1 [2 * lab->nTissue * lab->nReplica + 3 * lab->nReplica + i] ;
	      }	 
	  }
	  break ;
	case 3:
	  /* NCI */
	  /* array A_2[3], A_3[3] B_3[3] must be rejected */
	  /* use col 1,3,5,7,9,26,23,25,27,29,45,47,49 for tissue A, replica 12345 lab 1 2 3*/
	  /* use col 61... for tissue B*/

	  /* site 1 : 20 manip */
	  /* site 2 : 14 manip a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 d1 d2 d3 BUT a3 is rubbish */
	  /* site 3 : 20 manip a3 et b3 = rubbish */
	  /* site 4 : name_again + 8 values log(b/a), replacing site 3 */

	  /* flags: 0 good, -50 not-found , -75 absent, -100 bad */
	  {
	    double b1[120], flag ;
	    int iTissue, iLab, iReplica ;
	    int i = 0 ;
	    
	    for (i=0 ; i < 54 ; i++)
	      if (! freedoublecut (&b1[i], "\t", &cutter))
		b1[i] = UT_NON_DOUBLE ;
	    for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	      for (iTissue = 0 ;  iTissue < lab->nTissue  ; iTissue++)
		{
		  for (iReplica = 0 ; iReplica < 5 ; iReplica++)
		    {
		      if (
			  (iLab == 2  && iTissue == 0 && iReplica == 4) ||
			  (iLab == 2  && iTissue == 1 && iReplica == 4) ||
			  (iLab == 2  && iTissue == 2 && iReplica == 3) ||
			  (iLab == 2  && iTissue == 2 && iReplica == 4) ||
			  (iLab == 2  && iTissue == 3 && iReplica == 3) ||
			  (iLab == 2  && iTissue == 3 && iReplica == 4)
			  )
			continue ;
		      if (freedoublecut (&flag, "\t", &cutter) && JEAN_FLAG &&
			  flag <= -70)
		      /* removing -100=bad is favorable, 
		       * removing -75 has no effect, 
		       * removing -50, we loose good weak signals
		       */
			b1[iReplica + iLab * lab->nReplica + iTissue * lab->nLab * lab->nReplica ] 
			  = UT_NON_DOUBLE ;
		    }
		  for (; iReplica < lab->nReplica;  iReplica++)
		    {
		      mm->A_1[iReplica] = mm->A_2[iReplica] = mm->A_3[iReplica] = 
		      mm->B_1[iReplica] = mm->B_2[iReplica] = mm->B_3[iReplica] = 
		      mm->C_1[iReplica] = mm->C_2[iReplica] = mm->C_3[iReplica] =
		      mm->D_1[iReplica] = mm->D_2[iReplica] = mm->D_3[iReplica] = UT_NON_DOUBLE ;
		    }
		}
	    for (i=0 ; i < 5 ; i++) 
	      {
		mm->A_1[i] = b1 [i] ;
		mm->A_2[i] = b1 [5+i] ; /* but only 4 available */
		mm->A_3[i] = b1 [9+i] ;
		
		mm->B_1[i] = b1 [14+i] ;
		mm->B_2[i] = b1 [19+i] ; /* but only 4 available */
		mm->B_3[i] = b1 [23+i] ;

		mm->C_1[i] = b1 [28+i] ;
		mm->C_2[i] = b1 [33+i] ; /* but only 3 available */
		mm->C_3[i] = b1 [36+i] ;

		mm->D_1[i] = b1 [41+i] ;
		mm->D_2[i] = b1 [46+i] ; /* but only 3 available */
		mm->D_3[i] = b1 [49+i] ;
	      }
	    mm->A_2[2] = mm->A_3[2] = mm->B_3[2] = 
	      mm->A_2[4] = mm->B_2[4] = mm->C_2[3] = mm->D_2[3] = mm->C_2[4] = mm->D_2[4] = UT_NON_DOUBLE ;
	    if (1) /* also remove according to danielle */
	      { /* NCI  b2.3  b3.1245  c3.125 d3.* */
		mm->B_2[2] = 
		  mm->B_3[0] = mm->B_3[1] = mm->B_3[3] = mm->B_3[4] = 
		  mm->C_3[0] = mm->C_3[1] = mm->C_3[4] = 
		  mm->D_3[0] = mm->D_3[1] = mm->D_3[2]  = mm->D_3[3] = mm->D_3[4] = UT_NON_DOUBLE ;
	      }
	    if (0) /* NEW NCI data jun 17, number of probe drop,s sensitivity drops 7 points */
	      { /* we now replace site 3 */
		for (i=0 ; i < 5 ; i++) 
		  {
		    if (1) mm->A_3[i] = mm->B_3[i] = UT_NON_DOUBLE ;
		    if (0) mm->C_3[i] = mm->D_3[i] = UT_NON_DOUBLE ;
		  }
		if ((cp = freeword ()) && *cp) /* repetition of the name of the probe */
		  {
		    double z ;
		    Array zz = arrayCreate (8, double) ;
		    int j ;
		    
		    for (i=j=0 ; i < 8 ; i++)  /* read the 8 values */
		      if (freedoublecut (&z, "\t", &cutter))
			array (zz, j++, double) = z ;
		      else
			array (zz, j++, double) = UT_NON_DOUBLE ;
		    if (1)
		      {
			/* note that we may distort the data by sorting, since later we shift the whole chip */
			arraySort (zz, doubleOrder) ; 
			j = arrayMax(zz) - 5 ;
			if (j<=0) j = 0 ; else j = (j+1)/2 ;
			for (i=0 ; i < 5 && i+j < arrayMax(zz) ; i++) 
			  {  mm->A_3[i] = 0 ; mm->B_3[i] = array (zz, i+j, double) ; }
		      }
		    else
		      {
			mm->A_2[4] = 0 ; 
			mm->B_2[4] = array (zz, 0, double) ;
			for (i=0 ; i < 5 ; i++) 
			  { mm->A_3[i] = 0 ; mm->B_3[i] = array (zz, i, double) ; }
		      }
		    arrayDestroy (zz) ;
		  }
	      }
	  }
	  break ;
	}
      if (ok) /* keep this line */
	mm->probe = probe ;  /* keep this line */
      else
	printf ("ERROR parsing file DATA2/%s%s.txt line %d\n", lab->lab, lab->bench, line) ;

      nn++ ;
    }
  freeclose (level) ; /* will close f */

  if ((f = filopen (messprintf ("OSP/osp.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      double tm = 0 ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  if (! freedouble (&tm)) continue ;
	  
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->tm = tm ;
	}
      freeclose (level) ; /* will close f */
    }
  if ((f = filopen (messprintf ("OSP/ospscore.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      double score = 0 ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  if (! freedouble (&score)) continue ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->score = score ;
	}
      freeclose (level) ; /* will close f */
    }
  if (maqc->useNm && (f = filopen (messprintf ("P2G/probe2nm.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      int nm ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictNM, cp, &nm) ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  if (mm->nm && mm->nm != nm) mm->nm2 = nm ;
	  else mm->nm = nm ;
	  mm->geneId = keySet (maqc->nm2geneId, nm) ; /* the geneid of the nm */
	}
      freeclose (level) ; /* will close f */
    }

  if (maqc->useNm && (f = filopen (messprintf ("P2G/probe2gene.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      int nmGene ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictGene, cp, &nmGene) ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->gene = nmGene ; 
	}
      freeclose (level) ; /* will close f */
    }

   if ((f = filopen ("P2G/basicGene", "txt", "r")))
    {
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the gene name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictBasicGene, cp, 0) ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2G/probe2Confirmed_gene.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      int acvGene ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictGene, cp, &acvGene) ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->gene = acvGene ;
	  if (dictFind (maqc->dictBasicGene, cp, 0))
	    mm->basic = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }
  if (!maqc->useNm && (f = filopen (messprintf ("P2G/probe2Confirmed_geneId.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      int acvGeneId ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictGene, cp, &acvGeneId) ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->geneId = acvGeneId ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2G/pcl.%s", lab->lab), "txt", "r")))
    {
      int probe, pcl ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ;  
	  if (! freeint (&pcl)) continue ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->probe_class = (pcl % 100) ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2T/cotitrating.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ; 
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->valid = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2T/titrating.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ; 
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->titrating = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2T/cotitrating.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ; 
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->cotitrating = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2T/antititrating.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ; 
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->antititrating = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }

  if (!maqc->useNm && (f = filopen (messprintf ("P2T/titratable.%s", lab->lab), "txt", "r")))
    {
      int probe ;
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  if (!dictFind (maqc->dict, messprintf ("%s|%s", cp, lab->bench), &probe))
	    continue ; 
	  mm = arrayp (maqc->aa, probe, MM) ;
	  mm->titratable = TRUE ;
	}
      freeclose (level) ; /* will close f */
    }

  /* kill genes not probe_class 1 or 2 */
  for (probe = 0 ; probe < arrayMax(maqc->aa) ; probe++)
    {
      mm = arrayp (maqc->aa, probe, MM) ;
      if (strcmp (mm->lab,lab->lab))
	continue ; 
      if (mm->gene)
	switch (mm->probe_class/10)
	  {
	  case 1: 
	  case 2:
	    keySet (geneTested, nGeneTested++) = mm->gene ;
	    break ;
	  default:
	    mm->gene =  mm->geneId = 0 ;
	    break ;
	  }
    }

  if (! keySetMax(maqc->titratingGeneaB) &&
      (f = filopen (maqc->useNm ? "TitratingNmGenes_all" : "TitratingAcvGenes_all","txt"
		    , "r")))
    {
      int gene, ng, type ;
      KEYSET ks, ks2, ks60 ;
      KEYSET ksaB60 = keySetCreate () ;
      KEYSET ksAb60 = keySetCreate () ;
      float foldChange = 0, signal = 0 ;
      
      level = freesetfile (f, 0) ;
      
      while (freecard (level))
	{
	  /* parse the probe name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictGene, cp, &gene) ;
	  cp = freeword () ; /* platform name */ 
	  cp = freeword () ; /* A or B */
	  if (*cp == 'A') { ks60 = ksAb60 ;}
	  else if (*cp == 'B') { ks60 = ksaB60 ;}
	  else continue ;
	  freefloat (&signal)  ; 
	  freefloat (&foldChange)  ; 
	  type = 0 ;
	  freeint (&type) ;
	  if (type & 0x1) /* all60 */
	    {
	      ng = keySet (ks60, gene) ;
	      keySet (ks60, gene) = ng + 1 ;
	    }
	}
      freeclose (level) ; /* will close f */

      ks = ksaB60 ;
      keySet (maqc->titratingGeneaB, 0) = 0 ;
      for (ng = gene = 0 ; gene < keySetMax(ks) ; gene++)
	if (keySet (ks, gene) >= LIMIT_PLATFORM)
	  keySet (maqc->titratingGeneaB, ng++) = gene ;
      keySetSort (maqc->titratingGeneaB) ;
      keySetCompress (maqc->titratingGeneaB) ;
      keySetDestroy (ks) ;

      ks = ksAb60 ;
      keySet (maqc->titratingGeneAb, 0) = 0 ;
      for (ng = gene = 0 ; gene < keySetMax(ks) ; gene++)
	if (keySet (ks, gene) >= LIMIT_PLATFORM)
	  keySet (maqc->titratingGeneAb, ng++) = gene ;
      keySetSort (maqc->titratingGeneAb) ;
      keySetCompress (maqc->titratingGeneAb) ;
      keySetDestroy (ks) ;

      ks2 = keySetOR (maqc->titratingGeneAb, maqc->titratingGeneaB) ;
      maqc->titratingGeneAny = ks2 ;
   }
  keySetSort (geneTested) ;
  keySetCompress (geneTested) ;
  {
    KEYSET ks = 0 ;

    ks = keySetAND (maqc->titratingGeneAny, geneTested) ;
    lab->nGeneTestedValidatedAny = keySetMax(ks) ? keySetMax(ks) : 1 ;
    keySetDestroy (ks) ;
  }
  keySetDestroy (geneTested) ;
  return nn > 0 ? TRUE : FALSE ;
} /* maqcGetTxtData */

/*************************************************************************************/

static BOOL maqcGetNm2GeneId (MAQC *maqc)
{
  FILE *f ;
  BOOL ok = FALSE ;
  int level, nm, geneId ;
  const char *ccp ;

  if ((f = filopen ("P2G/nm2geneid", "txt", "r")))
    {
      level = freesetfile (f, 0) ;
      while (freecard (level))
	{
	  /* parse the NM */
	  ccp = freeword () ;
	  if (!ccp) continue ;
	  dictAdd (maqc->dictNM, ccp, &nm) ;
	  ccp = freeword () ;
	  if (!ccp) continue ;
	  dictAdd (maqc->dictGene, ccp, &geneId) ;
	  ok = TRUE ;
	  keySet (maqc->nm2geneId, nm) = geneId ; 
	}
      freeclose (level) ; /* will close f */
    }
  return ok ;
} /* maqcGetNm2GeneId */

/*************************************************************************************/

static BOOL maqcGetTxtData (MAQC *maqc)
{
  LAB *lab ;
  BOOL ok = FALSE ;

  maqc->dict = dictHandleCreate (10000, maqc->h) ;
  maqc->dictNM = dictHandleCreate (10000, maqc->h) ;
  maqc->dictGene = dictHandleCreate (10000, maqc->h) ;
  maqc->dictGeneId = dictHandleCreate (10000, maqc->h) ;
  maqc->dictBasicGene = dictHandleCreate (10000, maqc->h) ;
  maqc->nm2geneId = keySetHandleCreate (maqc->h) ;
  maqc->testedBy3 = keySetHandleCreate (maqc->h) ;
  dictAdd (maqc->dict, "toto", 0) ; /* avoid 0 */
  dictAdd (maqc->dictNM, "toto", 0) ; /* avoid 0 */
  dictAdd (maqc->dictGene, "toto", 0) ; /* avoid 0 */
  dictAdd (maqc->dictGeneId, "toto", 0) ; /* avoid 0 */
  maqc->aa = arrayHandleCreate (10000, MM, maqc->h) ;

  maqcGetNm2GeneId (maqc) ;

  for (lab = &labos [0] ; *lab->lab ; lab++)
    {
      if (maqc->labName && strcmp (maqc->labName,lab->lab) &&
	  strcmp (maqc->labName, messprintf("%s_%s", lab->lab, lab->bench)))			    
	continue ;
      ok = maqcGetTxtDataLab (maqc, lab) ;
    }
  /* read the probe 2 gene mapping */
  {
    FILE *f = filopen ("P2G/P2G_any","txt","r") ;
    if (f)
      {
	int gene, n3 = 0, level = freesetfile (f, 0) ;
	char *cp ;

	freecard (level) ; /* jump caption line */
	while (freecard (level))
	  {
	    cp = freeword () ;
	    dictAdd (maqc->dictGene, cp, &gene) ;
	    keySet (maqc->testedBy3, n3++) = gene ;
	  }
	freeclose (level) ; /* will close f */
      }
  }
  keySetSort (maqc->testedBy3) ;
  keySetNCompress (maqc->testedBy3, 3) ;

  return ok ;
} /* maqcGetTxtDataLab */

/*************************************************************************************/

static BOOL maqcGetBinData (int version, MAQC *maqc)
{
  int fd = -1 ;
  int nn, n ;
  char *buf, *cp ;
  BOOL ok = FALSE ;

  fd = open (messprintf ("maqc.%d.binary_data", version), O_RDWR, 0644) ;
  if (fd >= 0 && (n = read (fd, &nn, sizeof (int))) > 0 )
    {
      ok = TRUE ;
      maqc->aa = arrayHandleCreate (nn+1, MM, maqc->h) ;
      array (maqc->aa, nn - 1, MM).probe = 0 ; /* make room */
      if ((n = read (fd, arrp (maqc->aa, 0, MM), nn * sizeof (MM)))  != nn * sizeof (MM))
	ok = FALSE ;

      nn = read (fd, &nn, sizeof (int)) ;
      cp = buf = messalloc (nn) ; 
      if ((n = read (fd, buf, nn))  != nn)
	ok = FALSE ;
      maqc->dict = dictHandleCreate (10000, maqc->h) ;
      while (*cp)
	{
	  dictAdd (maqc->dict, cp, 0) ;
	  n = strlen (cp) + 1 ;
	  cp += n ;
	}
      messfree (buf) ;

      nn = read (fd, &nn, sizeof (int)) ;
      cp = buf = messalloc (nn) ; 
      if ((n = read (fd, buf, nn))  != nn)
	ok = FALSE ;
      maqc->dictNM = dictHandleCreate (10000, maqc->h) ;
      while (*cp)
	{
	  dictAdd (maqc->dictNM, cp, 0) ;
	  n = strlen (cp) + 1 ;
	  cp += n ;
	}
      messfree (buf) ;

      nn = read (fd, &nn, sizeof (int)) ;
      cp = buf = messalloc (nn) ; 
      if ((n = read (fd, buf, nn))  != nn)
	ok = FALSE ;
      maqc->dictGeneId = dictHandleCreate (10000, maqc->h) ;
      while (*cp)
	{
	  dictAdd (maqc->dictGeneId, cp, 0) ;
	  n = strlen (cp) + 1 ;
	  cp += n ;
	}
      messfree (buf) ;

      nn = read (fd, &nn, sizeof (int)) ;
      cp = buf = messalloc (nn) ; 
      if ((n = read (fd, buf, nn))  != nn)
	ok = FALSE ;
      maqc->dictGene = dictHandleCreate (10000, maqc->h) ;
      while (*cp)
	{
	  dictAdd (maqc->dictGene, cp, 0) ;
	  n = strlen (cp) + 1 ;
	  cp += n ;
	}
      messfree (buf) ;
    }
  if (fd >= 0) close (fd) ;
  if (0)
    maqc->dict = dictHandleCreate (10000, maqc->h) ;
  return ok ;
}

/*************************************************************************************/

static BOOL maqcWriteBinData (int version, MAQC *maqc)
{
  int i, fd = -1 ;
  int nn, n ;
  char *cp, *buf ;
  const char *cq ;
  BOOL ok = FALSE ;

  nn = maqc->aa ? arrayMax (maqc->aa) : 0 ;
  if (! nn)
    return FALSE ;

  fd = open (messprintf ("maqc.%d.binary_data", version), O_RDWR | O_CREAT, 0644) ;
  if (fd >= 0)
    {
      ok = TRUE ;
      nn = arrayMax (maqc->aa) ;
      write (fd, &nn, sizeof (int)) ;
      if ((n = write (fd, arrp (maqc->aa, 0, MM), nn * sizeof (MM)))  != nn * sizeof (MM))
	ok = FALSE ;

      for (nn = 0, i = 1 ; i <= dictMax (maqc->dict) ; i++)
	{
	  cq = dictName (maqc->dict, i) ;
	  nn += strlen (cq) + 1 ;
	}
      nn++ ; /* double zero to close the buffer */
      cp = buf = messalloc (nn) ;
      for (i = 1 ; i <= dictMax (maqc->dict) ; i++)
	{
	  cq = dictName (maqc->dict, i) ;
	  while ((*cp++ = *cq++)) ;
	}
      cp++ ;
      write (fd, &nn, sizeof (int)) ;
      if ((n = write (fd, buf, nn))  != nn)
	ok = TRUE ;
      messfree (buf) ;

      for (nn = 0, i = 1 ; i <= dictMax (maqc->dictNM) ; i++)
	{
	  cq = dictName (maqc->dictNM, i) ;
	  nn += strlen (cq) + 1 ;
	}
      nn++ ; /* double zero to close the buffer */
      cp = buf = messalloc (nn) ;
      for (i = 1 ; i <= dictMax (maqc->dictNM) ; i++)
	{
	  cq = dictName (maqc->dictNM, i) ;
	  while ((*cp++ = *cq++)) ;
	}
      cp++ ;
      write (fd, &nn, sizeof (int)) ;
      if ((n = write (fd, buf, nn))  != nn)
	ok = TRUE ;
      messfree (buf) ;


      for (nn = 0, i = 1 ; i <= dictMax (maqc->dictGene) ; i++)
	{
	  cq = dictName (maqc->dictGene, i) ;
	  nn += strlen (cq) + 1 ;
	}
      nn++ ; /* double zero to close the buffer */
      cp = buf = messalloc (nn) ;
      for (i = 1 ; i <= dictMax (maqc->dictGene) ; i++)
	{
	  cq = dictName (maqc->dictGene, i) ;
	  while ((*cp++ = *cq++)) ;
	}
      cp++ ;
      write (fd, &nn, sizeof (int)) ;
      if ((n = write (fd, buf, nn))  != nn)
	ok = TRUE ;
      messfree (buf) ;


      for (nn = 0, i = 1 ; i <= dictMax (maqc->dictGeneId) ; i++)
	{
	  cq = dictName (maqc->dictGeneId, i) ;
	  nn += strlen (cq) + 1 ;
	}
      nn++ ; /* double zero to close the buffer */
      cp = buf = messalloc (nn) ;
      for (i = 1 ; i <= dictMax (maqc->dictGeneId) ; i++)
	{
	  cq = dictName (maqc->dictGeneId, i) ;
	  while ((*cp++ = *cq++)) ;
	}
      cp++ ;
      write (fd, &nn, sizeof (int)) ;
      if ((n = write (fd, buf, nn))  != nn)
	ok = TRUE ;
      messfree (buf) ;
    }
  if (fd >= 0) close (fd) ;

  return ok ;
}

/*************************************************************************************/

static BOOL maqcGetData (MAQC *maqc)
{
  int version = 2 ;

  maqc->titratingGeneAb = keySetCreate () ;
  maqc->titratingGeneaB = keySetCreate () ;

  if (0 && maqcGetBinData (version, maqc))
    return TRUE ;
  if (maqcGetTxtData (maqc))
    {
      if (0) maqcWriteBinData (version, maqc) ;
      return TRUE ;
    }
  return FALSE ;
}

/*************************************************************************************/
/* force mean = m0, sigma = s0 */
static double maqcRenormalize (MAQC *maqc, double x, double xm, double s)
{
  double s0 = maqc->s0, m0 = maqc->m0 ? maqc->m0 : 10.0 ;
  if (m0 != -1)
    {
      x -= xm ;
      if (s0 > 0 && s > 0)
	x = x * s0 / sqrt(s) ;
      x += m0 ;
    }
  return x ;
}

/*************************************************************************************/
/* force the average value of each replica to 10, the current average is in tt */
static void  maqcForceAverage (MAQC *maqc, LAB *lab, TT *tt, TT *tt2)
{
  int i, ii ;
  MM *mm ;

  for (ii = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (mm->id != lab->id)
	continue ;

      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  if (lab->nTissue > 0)
	    {
	      if ( mm->A_1[i] != UT_NON_DOUBLE)
		mm->A_1 [i] = maqcRenormalize (maqc, mm->A_1 [i], tt->A_1[i], tt2->A_1[i]) ;
	      if ( mm->A_2[i] != UT_NON_DOUBLE)
		mm->A_2 [i] = maqcRenormalize (maqc, mm->A_2 [i], tt->A_2[i], tt2->A_2[i]) ;
	      if ( mm->A_3[i] != UT_NON_DOUBLE)
		mm->A_3 [i] = maqcRenormalize (maqc, mm->A_3 [i], tt->A_3[i], tt2->A_3[i]) ;
	    }

	  if (lab->nTissue > 1)
	    {
	      if ( mm->B_1[i] != UT_NON_DOUBLE)
		mm->B_1 [i] = maqcRenormalize (maqc, mm->B_1 [i], tt->B_1[i], tt2->B_1[i]) ;
	      if ( mm->B_2[i] != UT_NON_DOUBLE)
		mm->B_2 [i] = maqcRenormalize (maqc, mm->B_2 [i], tt->B_2[i], tt2->B_2[i]) ;
	      if ( mm->B_3[i] != UT_NON_DOUBLE)
		mm->B_3 [i] = maqcRenormalize (maqc, mm->B_3 [i], tt->B_3[i], tt2->B_3[i]) ;
	    }

	  if (lab->nTissue > 2)
	    {
	      if ( mm->C_1[i] != UT_NON_DOUBLE)
		mm->C_1 [i] = maqcRenormalize (maqc, mm->C_1 [i], tt->C_1[i], tt2->C_1[i]) ;
	      if ( mm->C_2[i] != UT_NON_DOUBLE)
		mm->C_2 [i] = maqcRenormalize (maqc, mm->C_2 [i], tt->C_2[i], tt2->C_2[i]) ;
	      if ( mm->C_3[i] != UT_NON_DOUBLE)
		mm->C_3 [i] = maqcRenormalize (maqc, mm->C_3 [i], tt->C_3[i], tt2->C_3[i]) ;
	    }

	  if (lab->nTissue > 3)
	    {
	      if ( mm->D_1[i] != UT_NON_DOUBLE)
		mm->D_1 [i] = maqcRenormalize (maqc, mm->D_1 [i], tt->D_1[i], tt2->D_1[i]) ;
	      if ( mm->D_2[i] != UT_NON_DOUBLE)
		mm->D_2 [i] = maqcRenormalize (maqc, mm->D_2 [i], tt->D_2[i], tt2->D_2[i]) ;
	      if ( mm->D_3[i] != UT_NON_DOUBLE)
		mm->D_3 [i] = maqcRenormalize (maqc, mm->D_3 [i], tt->D_3[i], tt2->D_3[i]) ;
	    }

	}
    }
} /* maqcForceAverage */

/*************************************************************************************/
/* force the average value of each replica to 10, the current average is in tt 
 * centre = 0 : just measure the median, 1: force, 2: force the GLOBAL median
 */
static void  maqcForceMedian (MAQC *maqc, LAB *lab, MM *mm3, int centre)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, jj, iLab, iTissue, iReplica, k, nd2, nn ;
  MM *mm ;
  int NNN = NLAB * NTISSUE * NREPLICA ;
  Array CC[NNN] ;
  double deltaMedian[NNN], *zp, d2 ; ;
  BOOL isBasic = centre >= 10 ? TRUE : FALSE ;

  centre = centre % 10 ;
  memset (deltaMedian, 0, sizeof(deltaMedian)) ;
  memset (CC, 0, sizeof(CC)) ;

  for (iLab = 0 ; iLab < lab->nLab ; iLab++)
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
      for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	{
	  k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	  CC[k] = arrayHandleCreate (30000, double, h) ;
	}

  /* gather the signals */
  for (ii = jj = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (isBasic && ! mm->basic)
	continue ;
      if (mm->id != lab->id)
	continue ;
      zp = &(mm->A_1[0]) ;
      for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
	  for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	    {
	      k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	      nn = arrayMax (CC[k]) ;
	      if (zp[k] != UT_NON_DOUBLE)
		array (CC[k], nn, double) = zp[k] ;
	      else if (0) /* stupid idea */
		array (CC[k], nn, double) = - 10000 ;
	    }
    }

  
  /* sort and find the median store the result in mm3 */
  zp = &(mm3->A_1[0]) ;
  d2 = 0 ; nd2 = 0 ;
  for (iLab = 0 ; iLab < lab->nLab ; iLab++)
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
      for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	{
	  k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	  if ((nn=arrayMax(CC[k])))
	    {
	      arraySort (CC[k], doubleOrder) ;
	      zp[k] = array(CC[k], nn/2, double) ;
	      deltaMedian[k] = 10.0 - zp[k] ;
	      nd2++ ; d2 += deltaMedian[k] ;
	    }
	}
  d2 /= (nd2 ? nd2 : 1) ;
  /* recentre */
  if (centre)
    for (ii = jj = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
      {
	if (! mm->probe)
	  continue ;
	if (mm->id != lab->id)
	  continue ;
	zp = &(mm->A_1[0]) ;
	for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	  for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
	    for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	      {
		k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
		if (zp[k] != UT_NON_DOUBLE)
		  zp[k] += (centre == 2 ? d2 : deltaMedian[k]) ;
	      }
      }
  ac_free (h) ;
} /* maqcForceMedian  */

/*************************************************************************************/
/*  -exportP2Nm : export something */
static void maqcC1 (MAQC *maqc, LAB *lab)
{
  int nn, ii, jj, i, pass = 0, nsA, nsB ;
  TT tt, tt2, tb ;
  MM mm3, mmb3, *mm ;
  double total, nTotal, s, sigma, nSigma, sA, sB ;
  BOOL gotolao = FALSE ;

 lao:
  lab->gene2mmAb = keySetReCreate (lab->gene2mmAb) ;
  lab->gene2mmaB = keySetReCreate (lab->gene2mmaB) ;
  total = 0, nTotal = 0 ;
  memset (&tt, 0, sizeof(TT)) ;
  memset (&tb, 0, sizeof(TT)) ;
  memset (&tt2, 0, sizeof(TT)) ;
  memset (&mm3, 0, sizeof(MM)) ;
  memset (&mmb3, 0, sizeof(MM)) ;

  if (!maqc->aa || ! arrayMax (maqc->aa))
    {
      printf ("C1: no data files, sorry\n") ;
      return ;
    }
  for (ii = nn = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (mm->id != lab->id)
	continue ;
      if (mm->gene)
	keySet (maqc->testedGene, mm->gene) |= lab->id ;
      if (mm->geneId)
	keySet (maqc->testedGeneId, mm->geneId) |= lab->id ;
      nn++ ;
      nsA = nsB = 0 ; sA = sB = 0 ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  if (lab->nTissue > 0)
	    {
	      if ( mm->A_1[i] != UT_NON_DOUBLE)
		{ 
		  sA += mm->A_1[i] ; nsA++ ;
		  tt.A_1[i] += mm->A_1[i] ; tt.A_1nn[i]++ ; nTotal++ ; total += mm->A_1[i] ; 
		  if (mm->basic) {tb.A_1[i] += mm->A_1[i] ; tb.A_1nn[i]++ ;}
		}
	      if (lab->nLab > 1 && mm->A_2[i] != UT_NON_DOUBLE)
		{ 
		  sA += mm->A_2[i] ; nsA++ ;
		  tt.A_2[i] += mm->A_2[i] ; tt.A_2nn[i]++ ; nTotal++ ; total += mm->A_2[i] ; 
		  if (mm->basic) {tb.A_2[i] += mm->A_2[i] ; tb.A_2nn[i]++ ;}
		}
	      if (lab->nLab > 2 &&  mm->A_3[i] != UT_NON_DOUBLE)
		{ 
		  sA += mm->A_3[i] ; nsA++ ;
		  tt.A_3[i] += mm->A_3[i] ; tt.A_3nn[i]++ ; nTotal++ ; total += mm->A_3[i] ; 
		  if (mm->basic) {tb.A_3[i] += mm->A_3[i] ; tb.A_3nn[i]++ ;}
		}
	    }

	  if (lab->nTissue > 1)
	    {
	      if ( mm->B_1[i] != UT_NON_DOUBLE)
		{ 
		  sB += mm->B_1[i] ; nsB++ ;
		  tt.B_1[i] += mm->B_1[i] ; tt.B_1nn[i]++ ; nTotal++ ; total += mm->B_1[i] ;
		  if (mm->basic) {tb.B_1[i] += mm->B_1[i] ; tb.B_1nn[i]++ ;}
		}
	      if (lab->nLab > 1 &&  mm->B_2[i] != UT_NON_DOUBLE)
		{ 
		  sB += mm->B_2[i] ; nsB++ ;
		  tt.B_2[i] += mm->B_2[i] ; tt.B_2nn[i]++ ; nTotal++ ; total += mm->B_2[i] ;
		  if (mm->basic) {tb.B_2[i] += mm->B_2[i] ; tb.B_2nn[i]++ ;}
		}
	      if (lab->nLab > 2 &&   mm->B_3[i] != UT_NON_DOUBLE)
		{ 
		  sB += mm->B_3[i] ; nsB++ ;
		  tt.B_3[i] += mm->B_3[i] ; tt.B_3nn[i]++ ; nTotal++ ; total += mm->B_3[i] ;
		  if (mm->basic) {tb.B_3[i] += mm->B_3[i] ; tb.B_3nn[i]++ ;}
		}
	    }

	  if (lab->nTissue > 2)
	    {
	      if ( mm->C_1[i] != UT_NON_DOUBLE)
		{ 
		  tt.C_1[i] += mm->C_1[i] ; tt.C_1nn[i]++ ; nTotal++ ; total += mm->C_1[i] ; 
		  if (mm->basic) {tb.C_1[i] += mm->C_1[i] ; tb.C_1nn[i]++ ;}
		}
	      if (lab->nLab > 1 &&  mm->C_2[i] != UT_NON_DOUBLE)
		{ 
		  tt.C_2[i] += mm->C_2[i] ; tt.C_2nn[i]++ ; nTotal++ ; total += mm->C_2[i] ; 
		  if (mm->basic) {tb.C_2[i] += mm->C_2[i] ; tb.C_2nn[i]++ ;}
		}
	      if (lab->nLab > 2 &&   mm->C_3[i] != UT_NON_DOUBLE)
		{ 
		  tt.C_3[i] += mm->C_3[i] ; tt.C_3nn[i]++ ; nTotal++ ; total += mm->C_3[i] ; 
		  if (mm->basic) {tb.C_3[i] += mm->C_3[i] ; tb.C_3nn[i]++ ;}
		}
	    }

	  if (lab->nTissue > 3)
	    {
	      if ( mm->D_1[i] != UT_NON_DOUBLE)
		{ 
		  tt.D_1[i] += mm->D_1[i] ; tt.D_1nn[i]++ ; nTotal++ ; total += mm->D_1[i] ;
		  if (mm->basic) {tb.D_1[i] += mm->D_1[i] ; tb.D_1nn[i]++ ;}
		}
	      if (lab->nLab > 1 &&  mm->D_2[i] != UT_NON_DOUBLE)
		{ 
		  tt.D_2[i] += mm->D_2[i] ; tt.D_2nn[i]++ ; nTotal++ ; total += mm->D_2[i] ;
		  if (mm->basic) {tb.D_2[i] += mm->D_2[i] ; tb.D_2nn[i]++ ;}
		}
	      if (lab->nLab > 2 &&  mm->D_3[i] != UT_NON_DOUBLE)
		{ 
		  tt.D_3[i] += mm->D_3[i] ; tt.D_3nn[i]++ ; nTotal++ ; total += mm->D_3[i] ;
		  if (mm->basic) {tb.D_3[i] += mm->D_3[i] ; tb.D_3nn[i]++ ;}
		}
	    }

	   if (lab->nTissue > 0)
	    {
	      if ( mm->A_1[i] != UT_NON_DOUBLE)
		{ tt2.A_1[i] += mm->A_1[i] * mm->A_1[i] ; tt2.A_1nn[i]++ ;}
	      if (lab->nLab > 1 &&  mm->A_2[i] != UT_NON_DOUBLE)
		{ tt2.A_2[i] += mm->A_2[i] * mm->A_2[i] ; tt2.A_2nn[i]++ ;}
	      if (lab->nLab > 2 &&   mm->A_3[i] != UT_NON_DOUBLE)
		{ tt2.A_3[i] += mm->A_3[i] * mm->A_3[i] ; tt2.A_3nn[i]++ ;}
	    }

	   if (lab->nTissue > 1)
	    {
	      if ( mm->B_1[i] != UT_NON_DOUBLE)
		{ tt2.B_1[i] += mm->B_1[i] * mm->B_1[i] ; tt2.B_1nn[i]++ ;}
	      if (lab->nLab > 1 &&  mm->B_2[i] != UT_NON_DOUBLE)
		{ tt2.B_2[i] += mm->B_2[i] * mm->B_2[i] ; tt2.B_2nn[i]++ ;}
	      if (lab->nLab > 2 &&   mm->B_3[i] != UT_NON_DOUBLE)
		{ tt2.B_3[i] += mm->B_3[i] * mm->B_3[i] ; tt2.B_3nn[i]++ ;}
	    }

	   if (lab->nTissue > 2)
	    {
	      if ( mm->C_1[i] != UT_NON_DOUBLE)
		{ tt2.C_1[i] += mm->C_1[i] * mm->C_1[i] ; tt2.C_1nn[i]++ ;}
	      if (lab->nLab > 1 &&  mm->C_2[i] != UT_NON_DOUBLE)
		{ tt2.C_2[i] += mm->C_2[i] * mm->C_2[i] ; tt2.C_2nn[i]++ ;}
	      if (lab->nLab > 2 &&   mm->C_3[i] != UT_NON_DOUBLE)
		{ tt2.C_3[i] += mm->C_3[i] * mm->C_3[i] ; tt2.C_3nn[i]++ ;}
	    }

	   if (lab->nTissue > 3)
	     {
	       if ( mm->D_1[i] != UT_NON_DOUBLE)
		 { tt2.D_1[i] += mm->D_1[i] * mm->D_1[i] ; tt2.D_1nn[i]++ ;}
	       if (lab->nLab > 1 &&  mm->D_2[i] != UT_NON_DOUBLE)
		 { tt2.D_2[i] += mm->D_2[i] * mm->D_2[i] ; tt2.D_2nn[i]++ ;}
	       if (lab->nLab > 2 &&   mm->D_3[i] != UT_NON_DOUBLE)
		 { tt2.D_3[i] += mm->D_3[i] * mm->D_3[i] ; tt2.D_3nn[i]++ ;}
	     }
	}
      if (nsA >0 &&  nsB > 0)
	{
	  sA /= nsA ; sB /= nsB ;
	  mm->signal =  (sA + sB) / 2.0 ;
	  mm->foldChange =  sA - sB ;

	  if (mm->gene && sA > sB)
	    {
	      double old = 0 ;
	      jj = keySet (lab->gene2mmAb, mm->gene) ;
	      if (jj)
		old = (arrp (maqc->aa, jj, MM))->foldChange ;
	      if (!jj || mm->foldChange > old)
		keySet (lab->gene2mmAb, mm->gene) = ii ;
	    }
	  if (mm->gene && sA < sB)
	    {
	      double old = 0 ;
	      jj = keySet (lab->gene2mmaB, mm->gene) ;
	      if (jj)
		old = (arrp (maqc->aa, jj, MM))->foldChange ;
	      if (!jj || mm->foldChange < old)
		keySet (lab->gene2mmaB, mm->gene) = ii ;
	    }
	}
    }
  if (!pass)
    printf ("Found %d lines in file %s%s until d22, is i used some basic, the median shown was that of the basic genes\n", nn, lab->lab,lab->bench) ;
  if (!nn) nn = 1 ;

  /* divide by n, and compute sigma2 = sum_x^2/nn - (sum_x/nn)^2
   */
  for (i = 0 ; i < lab->nReplica ; i++)
    {
      tt.A_1[i] /= tt.A_1nn[i] ? tt.A_1nn[i] : 1 ;
      tt.A_2[i] /= tt.A_2nn[i] ? tt.A_2nn[i] : 1 ;
      tt.A_3[i] /= tt.A_3nn[i] ? tt.A_3nn[i] : 1 ;
      
      tt.B_1[i] /= tt.B_1nn[i] ? tt.B_1nn[i] : 1 ;
      tt.B_2[i] /= tt.B_2nn[i] ? tt.B_2nn[i] : 1 ;
      tt.B_3[i] /= tt.B_3nn[i] ? tt.B_3nn[i] : 1 ;
      
      tt.C_1[i] /= tt.C_1nn[i] ? tt.C_1nn[i] : 1 ;
      tt.C_2[i] /= tt.C_2nn[i] ? tt.C_2nn[i] : 1 ;
      tt.C_3[i] /= tt.C_3nn[i] ? tt.C_3nn[i] : 1 ;
      
      tt.D_1[i] /= tt.D_1nn[i] ? tt.D_1nn[i] : 1 ;
      tt.D_2[i] /= tt.D_2nn[i] ? tt.D_2nn[i] : 1 ;
      tt.D_3[i] /= tt.D_3nn[i] ? tt.D_3nn[i] : 1 ;
      
      tb.A_1[i] /= tb.A_1nn[i] ? tb.A_1nn[i] : 1 ;
      tb.A_2[i] /= tb.A_2nn[i] ? tb.A_2nn[i] : 1 ;
      tb.A_3[i] /= tb.A_3nn[i] ? tb.A_3nn[i] : 1 ;
      
      tb.B_1[i] /= tb.B_1nn[i] ? tb.B_1nn[i] : 1 ;
      tb.B_2[i] /= tb.B_2nn[i] ? tb.B_2nn[i] : 1 ;
      tb.B_3[i] /= tb.B_3nn[i] ? tb.B_3nn[i] : 1 ;
      
      tb.C_1[i] /= tb.C_1nn[i] ? tb.C_1nn[i] : 1 ;
      tb.C_2[i] /= tb.C_2nn[i] ? tb.C_2nn[i] : 1 ;
      tb.C_3[i] /= tb.C_3nn[i] ? tb.C_3nn[i] : 1 ;
      
      tb.D_1[i] /= tb.D_1nn[i] ? tb.D_1nn[i] : 1 ;
      tb.D_2[i] /= tb.D_2nn[i] ? tb.D_2nn[i] : 1 ;
      tb.D_3[i] /= tb.D_3nn[i] ? tb.D_3nn[i] : 1 ;
      
      tt2.A_1[i] /= tt2.A_1nn[i] ? tt2.A_1nn[i] : 1 ;
      tt2.A_2[i] /= tt2.A_2nn[i] ? tt2.A_2nn[i] : 1 ;
      tt2.A_3[i] /= tt2.A_3nn[i] ? tt2.A_3nn[i] : 1 ;
      
      tt2.B_1[i] /= tt2.B_1nn[i] ? tt2.B_1nn[i] : 1 ;
      tt2.B_2[i] /= tt2.B_2nn[i] ? tt2.B_2nn[i] : 1 ;
      tt2.B_3[i] /= tt2.B_3nn[i] ? tt2.B_3nn[i] : 1 ;
      
      tt2.C_1[i] /= tt2.C_1nn[i] ? tt2.C_1nn[i] : 1 ;
      tt2.C_2[i] /= tt2.C_2nn[i] ? tt2.C_2nn[i] : 1 ;
      tt2.C_3[i] /= tt2.C_3nn[i] ? tt2.C_3nn[i] : 1 ;
      
      tt2.D_1[i] /= tt2.D_1nn[i] ? tt2.D_1nn[i] : 1 ;
      tt2.D_2[i] /= tt2.D_2nn[i] ? tt2.D_2nn[i] : 1 ;
      tt2.D_3[i] /= tt2.D_3nn[i] ? tt2.D_3nn[i] : 1 ;
      
      tt2.A_1[i] -= tt.A_1[i] * tt.A_1[i] ;
      tt2.A_2[i] -= tt.A_2[i] * tt.A_2[i] ;
      tt2.A_3[i] -= tt.A_3[i] * tt.A_3[i] ;
      
      tt2.B_1[i] -= tt.B_1[i] * tt.B_1[i] ;
      tt2.B_2[i] -= tt.B_2[i] * tt.B_2[i] ;
      tt2.B_3[i] -= tt.B_3[i] * tt.B_3[i] ;
      
      tt2.C_1[i] -= tt.C_1[i] * tt.C_1[i] ;
      tt2.C_2[i] -= tt.C_2[i] * tt.C_2[i] ;
      tt2.C_3[i] -= tt.C_3[i] * tt.C_3[i] ;
      
      tt2.D_1[i] -= tt.D_1[i] * tt.D_1[i] ;
      tt2.D_2[i] -= tt.D_2[i] * tt.D_2[i] ;
      tt2.D_3[i] -= tt.D_3[i] * tt.D_3[i] ;
    }
  
  lab->nSignal = nTotal ;  
  total = total / (nTotal ? nTotal : 1) ;  
  lab->meanSignal = total ;

  gotolao = FALSE ;
  if ((maqc->centre == MEDIANE || (lab->centre == MEDIANE && maqc->centre ==  CENTREBOF))
      && ! pass++) /* force average */
    { 
      gotolao = TRUE ; 
      if ((maqc->basic == BASIC || (maqc->basic == BASICBOF && lab->basic == BASIC)))
	{
	  maqcForceMedian (maqc, lab, &mm3, 0) ;    /* just measure the median of the whole chip  */
	  maqcForceMedian (maqc, lab, &mmb3, 11) ;  /* force the median of the basic genes */
	} 
      else
	{
	  maqcForceMedian (maqc, lab, &mmb3, 10) ;  /* just measure the median of the basic genes */
	  maqcForceMedian (maqc, lab, &mm3, 1) ;    /* force the median of the whole chip  */
	}
    }
  else if ((maqc->centre == MOYENNE || (lab->centre == MOYENNE && maqc->centre ==  CENTREBOF))
	   && ! pass++) /* force average */
    {
      gotolao = TRUE ; 
      maqcForceMedian (maqc, lab, &mmb3, 10) ; /* just measure the median of the basic genes */
      maqcForceMedian (maqc, lab, &mm3, 0) ;   /* just measure the median of the whole chip  */
      if ((maqc->basic == BASIC || (maqc->basic == BASICBOF && lab->basic == BASIC)))
	maqcForceAverage (maqc, lab, &tb, &tt2) ; 
      else
	maqcForceAverage (maqc, lab, &tt, &tt2) ; 
    }
  else if ((maqc->centre == NOCENTRE || (lab->centre == NOCENTRE && maqc->centre ==  CENTREBOF))
	   && ! pass++) /* force average */
    {  
      maqcForceMedian (maqc, lab, &mmb3, 10) ; /* just measure the median */
      maqcForceMedian (maqc, lab, &mm3, 2) ;   /* force the global median the whole lab */
      gotolao = TRUE ; 
    }
  else
    {
      maqcForceMedian (maqc, lab, &mm3, 0) ;   /* just measure the median of the whole chip  */
      maqcForceMedian (maqc, lab, &mmb3, 10) ; /* just measure the median of the basic genes */
    }

  printf ("%s\t%s\t%s\t%s\t%s\t\t", "Tissue", "A", "B", "C", "D") ;
  printf ("%s\t%s\t%s\t%s\t%s\n", "Tissue", "A", "B", "C", "D") ;
  printf ("%s\n","Lab/Replica Moyenne all probes                         basic genes") ;
  if (lab->nLab > 0)
    {
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("1_%d",i+1) ;
	  printf ("\t%.2f", tt.A_1[i]) ;
	  printf ("\t%.2f", tt.B_1[i]) ;
	  printf ("\t%.2f", tt.C_1[i]) ;
	  printf ("\t%.2f", tt.D_1[i]) ;
	  printf ("\t\t1_%d",i+1) ;
	  printf ("\t%.2f", tb.A_1[i]) ;
	  printf ("\t%.2f", tb.B_1[i]) ;
	  printf ("\t%.2f", tb.C_1[i]) ;
	  printf ("\t%.2f", tb.D_1[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ; 
    }

  if (lab->nLab > 1)
    {
      printf ("%s\n","Lab/Replica Moyenne all probes                         basic genes") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("2_%d",i+1) ;
	  printf ("\t%.2f", tt.A_2[i]) ;
	  printf ("\t%.2f", tt.B_2[i]) ;
	  printf ("\t%.2f", tt.C_2[i]) ;
	  printf ("\t%.2f", tt.D_2[i]) ;
	  printf ("\t\t2_%d",i+1) ;
	  printf ("\t%.2f", tb.A_2[i]) ;
	  printf ("\t%.2f", tb.B_2[i]) ;
	  printf ("\t%.2f", tb.C_2[i]) ;
	  printf ("\t%.2f", tb.D_2[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  if (lab->nLab > 2)
    {
      printf ("%s\n","Lab/Replica Moyenne all probes                         basic genes") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("3_%d",i+1) ;
	  printf ("\t%.2f", tt.A_3[i]) ;
	  printf ("\t%.2f", tt.B_3[i]) ;
	  printf ("\t%.2f", tt.C_3[i]) ;
	  printf ("\t%.2f", tt.D_3[i]) ;
	  printf ("\t\t3_%d",i+1) ;
	  printf ("\t%.2f", tb.A_3[i]) ;
	  printf ("\t%.2f", tb.B_3[i]) ;
	  printf ("\t%.2f", tb.C_3[i]) ;
	  printf ("\t%.2f", tb.D_3[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  printf ("%s\t%s\t%s\t%s\t%s\n", "Tissue", "A", "B", "C", "D") ;
  printf ("%s\n","Lab/Replica Mediane all probes                         basic genes") ;
  if (lab->nLab > 0)
    {
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("1_%d",i+1) ;
	  printf ("\t%.2f", mm3.A_1[i]) ;
	  printf ("\t%.2f", mm3.B_1[i]) ;
	  printf ("\t%.2f", mm3.C_1[i]) ;
	  printf ("\t%.2f", mm3.D_1[i]) ;
	  printf ("\t\t3_%d",i+1) ;
	  printf ("\t%.2f", mmb3.A_1[i]) ;
	  printf ("\t%.2f", mmb3.B_1[i]) ;
	  printf ("\t%.2f", mmb3.C_1[i]) ;
	  printf ("\t%.2f", mmb3.D_1[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ; 
    }

  if (lab->nLab > 1)
    {
      printf ("%s\n","Lab/Replica Mediane all probes                         basic genes") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("2_%d",i+1) ;
	  printf ("\t%.2f", mm3.A_2[i]) ;
	  printf ("\t%.2f", mm3.B_2[i]) ;
	  printf ("\t%.2f", mm3.C_2[i]) ;
	  printf ("\t%.2f", mm3.D_2[i]) ;
	  printf ("\t\t3_%d",i+1) ;
	  printf ("\t%.2f", mmb3.A_2[i]) ;
	  printf ("\t%.2f", mmb3.B_2[i]) ;
	  printf ("\t%.2f", mmb3.C_2[i]) ;
	  printf ("\t%.2f", mmb3.D_2[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  if (lab->nLab > 2)
    {
      printf ("%s\n","Lab/Replica Mediane all probes                         basic genes") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("3_%d",i+1) ;
	  printf ("\t%.2f", mm3.A_3[i]) ;
	  printf ("\t%.2f", mm3.B_3[i]) ;
	  printf ("\t%.2f", mm3.C_3[i]) ;
	  printf ("\t%.2f", mm3.D_3[i]) ;
	  printf ("\t\t3_%d",i+1) ;
	  printf ("\t%.2f", mmb3.A_3[i]) ;
	  printf ("\t%.2f", mmb3.B_3[i]) ;
	  printf ("\t%.2f", mmb3.C_3[i]) ;
	  printf ("\t%.2f", mmb3.D_3[i]) ;
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  printf ("%s\t%s\t%s\t%s\t%s\n", "Tissue", "A", "B", "C", "D") ;
  printf ("%s\n","Lab/Replica sigma") ;
  sigma = 0 ; nSigma = 0 ;
  if (lab->nLab > 0)
    {
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("1_%d",i+1) ;
	  s = sqrt (tt2.A_1[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.B_1[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.C_1[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.D_1[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  if (lab->nLab > 1)
    {
      printf ("%s\n","Lab/Replica sigma") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("2_%d",i+1) ;
	  s = sqrt (tt2.A_2[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.B_2[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.C_2[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.D_2[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  if (lab->nLab > 2)
    {
      printf ("%s\n","Lab/Replica sigma") ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  printf ("3_%d",i+1) ;
	  s = sqrt (tt2.A_3[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.B_3[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.C_3[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  s = sqrt (tt2.D_3[i]) ; printf ("\t%.2f", s) ; if (s > 0) { sigma += s ; nSigma++ ; }
	  printf ("\n") ;
	}
    }
  lab->sigma = sigma / (nSigma ? nSigma : 1) ;
  printf ("lab->sigma : %.2f\n", lab->sigma) ;
  printf ("\n") ;
  if (gotolao)
    goto lao ;
} /* exportP2Nm */

/*************************************************************************************/
/* if dx = -1 a < B, if dx = 1 A > b, if all replica ok, add flag in keyset */
static BOOL maqcIsAllab (MAQC *maqc, LAB *lab, MM *mm, int ii, KEYSET ksAllab, int dx)
{
  int i ;
  BOOL ok = TRUE ;
  
  if (! maqc->allAgree)
    return TRUE ;
  for (i = 0 ; ok && i < lab->nReplica ; i++)
    if (
	(lab->nLab > 0 &&
	 mm->A_1[i] != UT_NON_DOUBLE && 
	 mm->B_1[i] != UT_NON_DOUBLE && 
	 mm->C_1[i] != UT_NON_DOUBLE && 
	 mm->D_1[i] != UT_NON_DOUBLE &&
	 (
	  (dx ==  1 &&  (mm->A_1[i] < mm->C_1[i] || mm->C_1[i] < mm->D_1[i] || mm->D_1[i] < mm->B_1[i])) ||
	  (dx == -1 &&  (mm->A_1[i] > mm->C_1[i] || mm->C_1[i] > mm->D_1[i] || mm->D_1[i] > mm->B_1[i]))
	  )
	 ) ||
	(lab->nLab > 0 &&
	 mm->A_2[i] != UT_NON_DOUBLE && 
	 mm->B_2[i] != UT_NON_DOUBLE && 
	 mm->C_2[i] != UT_NON_DOUBLE && 
	 mm->D_2[i] != UT_NON_DOUBLE &&
	 (
	  (dx ==  1 &&  (mm->A_2[i] < mm->C_2[i] || mm->C_2[i] < mm->D_2[i] || mm->D_2[i] < mm->B_2[i])) ||
	  (dx == -1 &&  (mm->A_2[i] > mm->C_2[i] || mm->C_2[i] > mm->D_2[i] || mm->D_2[i] > mm->B_2[i]))
	  )
	 ) ||
	(lab->nLab > 0 &&
	 mm->A_3[i] != UT_NON_DOUBLE && 
	 mm->B_3[i] != UT_NON_DOUBLE && 
	 mm->C_3[i] != UT_NON_DOUBLE && 
	 mm->D_3[i] != UT_NON_DOUBLE &&
	 (
	  (dx ==  1 &&  (mm->A_3[i] < mm->C_3[i] || mm->C_3[i] < mm->D_3[i] || mm->D_3[i] < mm->B_3[i])) ||
	  (dx == -1 &&  (mm->A_3[i] > mm->C_3[i] || mm->C_3[i] > mm->D_3[i] || mm->D_3[i] > mm->B_3[i]))
	  )
	 )
	)
      ok = FALSE ;
    
  if (ok) 
    keySet (ksAllab, ii) |= lab->id << (dx == 1 ? 16 : 0) ;
  return ok ;
} /* maqcIsAllab */

/*************************************************************************************/
/*  -c2 c2: export number of tissue sensitive genes */
static void maqcC2 (MAQC *maqc, LAB *lab
		    , KEYSET ks2, KEYSET ks4
		    , KEYSET ksg2, KEYSET ksg4 
		    , KEYSET ksAllab
		    )
{
  double dSignal = .1 ; /* used to be dx = 1, with abi->sigma = 2.8 */
  int nn, ii, i, pcl ;
  double z, dx = .01 ;
  int NP = ((int)(50/dSignal)) | 0x1 ; /* number of points */
  int NC = NP / 2 ;                    /* centre point where fold change == 0 */
  int nAb [NP] ;
  int nACDB [NP] ;
  int nBDCA [NP] ;
  int nTitPcl [NP*10] ;
  int nPcl [NP*10] ;
  int nvACDB [NP] ;
  int nvBDCA [NP] ;
  int ngACDB [NP] ;
  int ngBDCA [NP] ;
  int ngidACDB [NP] ;
  int ngidBDCA [NP] ;
  int titratable [NP] ;
  int titrating [NP] ;
  int cotitrating [NP] ;
  int antititrating [NP] ;
  int nontitrating [NP] ;
  int titratable3 [NP] ;
  int titrating3 [NP] ;
  int cotitrating3 [NP] ;
  int antititrating3 [NP] ;
  int nontitrating3 [NP] ;
  double sAb [NP] ;
  double sACDB [NP] ;
  double sBDCA [NP] ;
  double svACDB [NP] ;
  double svBDCA [NP] ;
  double sgACDB [NP] ;
  double sgBDCA [NP] ;
  double sgidACDB [NP] ;
  double sgidBDCA [NP] ;
  TT tt, tt2 ;
  MM *mm ;

  /*
  dSignal = lab->sigma * maqc->nSigma  ;
  */
  memset (&tt, 0, sizeof(TT)) ;
  memset (&tt2, 0, sizeof(TT)) ;
  memset (nPcl, 0, sizeof (nPcl)) ;
  memset (nTitPcl, 0, sizeof (nTitPcl)) ;

  maqc->aceSignalOut = filopen (messprintf ("aceSignalOut.%s%s", lab->lab, lab->bench), "ace","w") ;

  for (i = 0 ; i < NP ; i++)
    {
      nACDB[i] = nBDCA[i] = nvACDB[i] = nvBDCA[i] = nAb[i] = ngACDB[i] = ngBDCA[i] = ngidACDB[i] = ngidBDCA[i] = 0 ;
      sACDB[i] = sBDCA[i] = sAb[i] = svACDB[i] = svBDCA[i] = sgACDB[i] = sgBDCA[i] = sgidACDB[i] = sgidBDCA[i] = 0 ;
      titratable[i] = titrating[i] = cotitrating[i] = antititrating[i] = nontitrating[i] = 0 ;
      titratable3[i] = titrating3[i] = cotitrating3[i] = antititrating3[i] = nontitrating3[i] = 0 ;
    }
  if (!maqc->aa || ! arrayMax (maqc->aa))
    {
      printf ("C1: no data files, sorry\n") ;
      return ;
    }
  ii = dictMax (maqc->dictNM) ;
  keySet (ks2, ii) = 0 ;
  keySet (ks4, ii) = 0 ;
  ii = dictMax (maqc->dictGene) ;
  keySet (ksg2, ii) = 0 ;
  keySet (ksg4, ii) = 0 ;

  for (ii = nn = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (mm->id != lab->id)
	continue ;
      nn++ ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  if ( mm->A_1[i] != UT_NON_DOUBLE)
	    { mm->A_1m += mm->A_1[i] ; mm->A_1n++ ; }
	  if (lab->nLab > 1 &&  mm->A_2[i] != UT_NON_DOUBLE)
	    { mm->A_1m += mm->A_2[i] ; mm->A_1n++ ; }
	  if (lab->nLab > 2 &&  mm->A_3[i] != UT_NON_DOUBLE)
	    { mm->A_1m += mm->A_3[i] ; mm->A_1n++ ; }

	  if ( mm->B_1[i] != UT_NON_DOUBLE)
	    { mm->B_1m += mm->B_1[i] ; mm->B_1n++ ; }
	  if (lab->nLab > 1 &&  mm->B_2[i] != UT_NON_DOUBLE)
	    { mm->B_1m += mm->B_2[i] ; mm->B_1n++ ; }
	  if (lab->nLab > 2 &&  mm->B_3[i] != UT_NON_DOUBLE)
	    { mm->B_1m += mm->B_3[i] ; mm->B_1n++ ; }

	  if ( mm->C_1[i] != UT_NON_DOUBLE)
	    { mm->C_1m += mm->C_1[i] ; mm->C_1n++ ; }
	  if (lab->nLab > 1 &&  mm->C_2[i] != UT_NON_DOUBLE)
	    { mm->C_1m += mm->C_2[i] ; mm->C_1n++ ; }
	  if (lab->nLab > 2 &&  mm->C_3[i] != UT_NON_DOUBLE)
	    { mm->C_1m += mm->C_3[i] ; mm->C_1n++ ; }

	  if ( mm->D_1[i] != UT_NON_DOUBLE)
	    { mm->D_1m += mm->D_1[i] ; mm->D_1n++ ; }
	  if (lab->nLab > 1 &&  mm->D_2[i] != UT_NON_DOUBLE)
	    { mm->D_1m += mm->D_2[i] ; mm->D_1n++ ; }
	  if (lab->nLab > 2 &&  mm->D_3[i] != UT_NON_DOUBLE)
	    { mm->D_1m += mm->D_3[i] ; mm->D_1n++ ; }
	}

      if (mm->A_1n)
	mm->A_1m /= mm->A_1n ;
      if (mm->B_1n)
	mm->B_1m /= mm->B_1n ;
      if (mm->C_1n)
	mm->C_1m /= mm->C_1n ;
      if (mm->D_1n)
	mm->D_1m /= mm->D_1n ;
      
      if (maqc->aceSignalOut)
	{
	  char *cp, buff[1000] ;
	  strcpy (buff, dictName (maqc->dict, mm->probe)) ;
	  cp = strstr (buff, "|") ;
	  if (cp) *cp = 0 ;
	  fprintf (maqc->aceSignalOut, "Probe \"%s\"\n", buff) ;
	  fprintf (maqc->aceSignalOut, "-D Signal_A\n-D Signal_C\n-D Signal_D\n-D Signal_B\n") ;
	  if (mm->A_1n)
	    fprintf (maqc->aceSignalOut, "Signal_A %.2f\n", mm->A_1m) ;
	  if (mm->C_1n)
	    fprintf (maqc->aceSignalOut, "Signal_C %.2f\n", mm->C_1m) ;
	  if (mm->D_1n)
	    fprintf (maqc->aceSignalOut, "Signal_D %.2f\n", mm->D_1m) ;
	  if (mm->B_1n)
	     fprintf (maqc->aceSignalOut, "Signal_B %.2f\n", mm->B_1m) ;
	  fprintf (maqc->aceSignalOut, "\n") ;
	}
      if (mm->gene) /* register the largest signal for each gene */
	{
	  double umax = 0 ;
	  int uu ;
	  if (mm->A_1n) umax = mm->A_1m ;
	  if (mm->B_1n && umax < mm->B_1m) umax = mm->B_1m ;
	  if (mm->C_1n && umax < mm->C_1m) umax = mm->C_1m ;
	  if (mm->D_1n && umax < mm->D_1m) umax = mm->D_1m ;
	  
	  umax = 10 + (umax - 10)/dSignal ;
	  uu = umax ; if (uu < 0) uu = 0 ;
	  if (uu > 31) uu = 31 ;
	  
	  keySet (maqc->signalGene, mm->gene) |= 0x1 << uu ;
      }

      for (i = 1 ; i < NP - 1 ; i++)  /* iterate on fold change */
	if (mm->A_1n && mm->B_1n &&
	    mm->B_1m < mm->A_1m - (NC - .5) * dSignal + i * dSignal &&
	    mm->B_1m > mm->A_1m - (NC + .5) * dSignal + i * dSignal)
	  {
	    nAb [i]++ ;
	    sAb[i] += (mm->A_1m + mm->B_1m)/2 ;
	    pcl = mm->probe_class/10 + 1 ;
	    nPcl [i * 10 + pcl]++ ;
	    {
	      int nnn = 0 ;
	      for (pcl = 0 ; pcl < 9 ; pcl++)
		nnn += nPcl [i * 10 + pcl] ;
	      if (nnn > nAb [i])
		invokeDebugger () ;
	    }
	    if (i>= NP || pcl > 9) messcrash ("bad i, pcl") ;

	    if (mm->titratable) 
	      {
		if (mm->probe_class < 30)
		  {
		    titratable[i]++ ;
		    if (mm->titrating) 
		      {
			titrating[i]++ ;
			if (mm->antititrating) 
			  antititrating[i]++ ;
			if (mm->cotitrating)
			  cotitrating[i]++ ;
		      }
		    else
		      nontitrating[i]++ ;
		  }
		else
		  {
		    titratable3[i]++ ;
		    if (mm->titrating) 
		      {
			titrating3[i]++ ;
			if (mm->antititrating) 
			  antititrating3[i]++ ;
			if (mm->cotitrating)
			  cotitrating3[i]++ ;
		      }
		    else
		      nontitrating3[i]++ ;
		  }
	      }

	    if (*lab->bench && ! strcmp (lab->lab, "AGL") && lab->nTissue == 2)
	      { /* if just AGL cy3/cy5 */
		if (i < NC && mm->geneId)
		  {
		    keySet (ks2, mm->geneId) |= lab->id ;
		    if (i < NC - 1)
		      keySet (ks4, mm->geneId) |= lab->id ;
		  }
		if (i < NC && mm->gene)
		  keySet (ksg2, mm->gene) |= lab->id ;
		if (i < (NC - 1)  && mm->gene)
		  keySet (ksg4, mm->gene) |= lab->id ;

		if (i > NC && mm->geneId)
		  {
		    keySet (ks2, mm->geneId) |= lab->id << 16 ;
		    if (i > NC + 1)
		      keySet (ks4, mm->geneId) |= lab->id << 16 ;
		  } 
		if (i > NC && mm->gene)
		  keySet (ksg2, mm->gene) |= lab->id << 16 ;
		if (i > (NC + 1) && mm->gene)
		  keySet (ksg4, mm->gene) |= lab->id << 16 ;
	      }
	    else if (lab->nTissue == 4 && mm->C_1n && mm->D_1n &&
		     mm->A_1m < mm->C_1m - dx &&
		     mm->C_1m < mm->D_1m - dx &&
		     mm->D_1m < mm->B_1m - dx &&
		     maqcIsAllab (maqc, lab, mm, ii, ksAllab, -1)
		     )
	      {
		nACDB[i]++ ;
		sACDB[i] += (mm->A_1m + mm->B_1m)/2 ;
		if (mm->valid) { nvACDB[i]++ ; svACDB[i] += (mm->A_1m + mm->B_1m)/2 ; }
		pcl = mm->probe_class/10 + 1 ;
		nTitPcl [i *  10 + pcl]++ ;
		if (i>= NP || pcl > 9) messcrash ("bad i, pcl") ;
		if (mm->gene) 
		  {
		    ngACDB[i]++ ;
		    sgACDB[i] += (mm->A_1m + mm->B_1m)/2 ;
		  }
		if (mm->geneId) 
		  {
		    ngidACDB[i]++ ;
		    sgidACDB[i] += (mm->A_1m + mm->B_1m)/2 ;
		  }
		if (mm->geneId)
		  {
		    if (i < NC)
		      keySet (ks4, mm->geneId) |= lab->id ;
		    if (i <= NC)
		      keySet (ks2, mm->geneId) |= lab->id ;
		  }
		if (mm->gene)
		  {
		    if (i < NC)
		      keySet (ksg4, mm->gene) |= lab->id ;
		    if (i <= NC)
		      keySet (ksg2, mm->gene) |= lab->id ;
		  }

	      }
	    else if (lab->nTissue == 4 && mm->C_1n && mm->D_1n &&
		     mm->A_1m > mm->C_1m  + dx &&
		     mm->C_1m > mm->D_1m  + dx &&
		     mm->D_1m > mm->B_1m  + dx &&
		     maqcIsAllab (maqc, lab, mm, ii, ksAllab, 1) 
		     )
	      {
		nBDCA[i]++ ; 
		sBDCA[i] += (mm->A_1m + mm->B_1m)/2 ;
		if (mm->valid) { nvBDCA[i]++ ; svBDCA[i] += (mm->A_1m + mm->B_1m)/2 ; }
		pcl = mm->probe_class/10 + 1 ;
		nTitPcl [i *  10 + pcl]++ ;
		if (i>= NP || pcl > 9) messcrash ("bad i, pcl") ;
		if (mm->gene) 
		  {
		    ngBDCA[i]++ ;
		    sgBDCA[i] += (mm->A_1m + mm->B_1m)/2 ;
		  }
		if (mm->geneId) 
		  {
		    ngidBDCA[i]++ ;
		    sgidBDCA[i] += (mm->A_1m + mm->B_1m)/2 ;
		  }
		if (mm->geneId)
		  {
		    if (i > NC)
		      keySet (ks4, mm->geneId) |= (lab->id << 16) ;
		    if (i >= NC)
		      keySet (ks2, mm->geneId) |= (lab->id << 16) ;
		  }
		if (mm->gene)
		  {
		    if (i > NC)
		      keySet (ksg4, mm->gene) |= (lab->id << 16) ;
		    if (i >= NC)
		      keySet (ksg2, mm->gene) |= (lab->id << 16) ;
		  }
	      }
	  }
    
      if (0) for (i = 24 ; i < 25; i++)
	if (maqc->dict && mm->A_1m > mm->B_1m - (NC + .5) + i)
	  {
	    printf ("%s %.2f %.2f\n", dictName (maqc->dict, mm->probe), mm->A_1m, mm->B_1m) ;
	  }
    }  
  printf ("Found %d lines in file %s%s\n", nn, lab->lab,lab->bench) ;
  printf ("Fold change log2(b/a)\tAll probes\tTitrating probes\tTitrating validated probes\tTitrating gene specific probe\tsame for gene with geneId\tProbe grade 1\tProbe grade 2\tProbe grade 3\tProbe grade 4\tProbe grade 5\tProbe grade 6\tTitrating probe grade 1\tTitrating probe grade 2\tTitrating probe grade 3\tTitrating probe grade 4\tTitrating probe grade 5\tTitrating probe grade 6\tFold change log2(b/a)\tTitratable probes\t%% co-titrating\t%% non-titrating\t%% anti-titrating\t%% titrating\tTitratable grade 3 probes\t%% co-titratinggrade 3 \t%% non-titratinggrade 3 \t%% anti-titratinggrade 3 \t%% titrating grade 3 \n") ;
  /* 
     col 1: Fold change log2(b/a), gene is mostly expresses in B tissue at end of table
     col 2: number of probes with corresponding fold change
     col 3: titrating probes
     col 4: titrating probes testing an aceview gene
     col 5: titrating probes testing a gene with geneId
     col 6: follow geneId
  */  
  if (0)
    {
      int ntt, nt, nct, nat, nnt ;
      
      ntt = nt = nct = nat = nnt = 0 ;
      for (i = NC ; i < NP ; i++)
	{
	  ntt += titratable[i] ; titratable[i] = ntt ;
	  nt += titrating[i] ; titrating[i] = nt ;
	  nct += cotitrating[i] ; cotitrating[i] = nct ;
	  nat += antititrating[i] ; antititrating[i] = nat ;
	  nnt += nontitrating[i] ; nontitrating[i] = nnt ;
	}
      ntt = nt = nct = nat = nnt = 0 ;
      for (i = NC ; i >= 0 ; i--)
	{
	  ntt += titratable[i] ; titratable[i] = ntt ;
	  nt += titrating[i] ; titrating[i] = nt ;
	  nct += cotitrating[i] ; cotitrating[i] = nct ;
	  nat += antititrating[i] ; antititrating[i] = nat ;
	  nnt += nontitrating[i] ; nontitrating[i] = nnt ;
	}
      ntt = nt = nct = nat = nnt = 0 ;
      for (i = NC ; i < NP ; i++)
	{
	  ntt += titratable3[i] ; titratable3[i] = ntt ;
	  nt += titrating3[i] ; titrating3[i] = nt ;
	  nct += cotitrating3[i] ; cotitrating3[i] = nct ;
	  nat += antititrating3[i] ; antititrating3[i] = nat ;
	  nnt += nontitrating3[i] ; nontitrating3[i] = nnt ;
	}
      ntt = nt = nct = nat = nnt = 0 ;
      for (i = NC ; i >= 0 ; i--)
	{
	  ntt += titratable3[i] ; titratable3[i] = ntt ;
	  nt += titrating3[i] ; titrating3[i] = nt ;
	  nct += cotitrating3[i] ; cotitrating3[i] = nct ;
	  nat += antititrating3[i] ; antititrating3[i] = nat ;
	  nnt += nontitrating3[i] ; nontitrating3[i] = nnt ;
	}
    }
  
  {
    int i1, i2, n, nn, tot[100] ;
    for (i = 0 ; i < 100 ; i++)
      tot[i] = 0 ;
    
    /* on coupe 1% de chaque cote */
    for (i = nn = 0 ; i < NP ; i++)   /* nn = total number of probes */
      nn += nAb [i] ;
    
    for (i1 = n = 0 ; 100 * n < nn  && i1 < NP ; i1++)
      n += nAb [i1] ;
    for (i2 = NP - 1, n = 0 ; 100 * n < nn  && i2 >= 0 ; i2--)
      n += nAb [i2] ;
    i1-- ; i2++ ;
    
    for (i = 0 ; i < NP ; i++)
      {
	if (0 || (i >= i1 && i <= i2))
	  printf ("%.1f\t%8d\t%8d\t%8d\t%8d\t%8d"
		  , (-NC+i) * dSignal
		  , nAb [i]
		  , nACDB [i] + nBDCA [i]
		  , nvACDB [i] + nvBDCA [i]
		  , ngACDB [i] + ngBDCA [i]
		  , ngidACDB [i] + ngidBDCA [i]
		  ) ;
	tot[0] += nAb [i] ;
	tot[1] += nACDB [i] + nBDCA [i] ;
	tot[2] += nvACDB [i] + nvBDCA [i] ;
	tot[3] += ngACDB [i] + ngBDCA [i] ;
	tot[4] += ngidACDB [i] + ngidBDCA [i] ;
	
	for (pcl = 1 ; pcl < 7 ; pcl++)
	  {
	    if (0 ||  (i >= i1 && i <= i2))
	      printf ("\t%8d", nPcl [i * 10 + (pcl+1)]) ;
	    tot[4 + pcl] += nPcl [i * 10 + (pcl+1)] ;
	    if (0 && i >= i1 && i <= i2)
	      printf ("\ttot6=%d", tot[6]) ;
	  }
	for (pcl = 1 ; pcl < 7 ; pcl++)
	  {
	    if (i >= i1 && i <= i2)
	      printf ("\t%8d", nTitPcl [i *  10 + (pcl+1)]) ;
	    tot[10 + pcl] += nTitPcl [i *  10 + (pcl+1)] ;
	  }
	z = 100.0/(titratable[i] ? titratable[i] : 1) ;
	if (i >= i1 && i <= i2)
	  printf ("\t%.1f\t%d\t%.2f\t%.2f\t%.2f\t%.2f"
		  , (-NC+i) * dSignal
		  , titratable[i]
		  , z*cotitrating[i]
		  , z*nontitrating[i]
		  , z*antititrating[i]
		  , z*titrating[i]
		  ) ;
	tot[17] += titratable[i] ;
	tot[18] += cotitrating[i] ;
	tot[29] += nontitrating[i] ;
	tot[20] += antititrating[i] ;
	tot[21] += titrating[i] ;
	z = 100.0/(titratable3[i] ? titratable3[i] : 1) ;
	if (i >= i1 && i <= i2)
	  printf ("\t%d\t%.2f\t%.2f\t%.2f\t%.2f"
		  , titratable3[i]
		  , z*cotitrating3[i]
		  , z*nontitrating3[i]
		  , z*antititrating3[i]
		  , z*titrating3[i]
		  ) ;
	tot[23] += titratable3[i] ;
	tot[24] += cotitrating3[i] ;
	tot[25] += nontitrating3[i] ;
	tot[26] += antititrating3[i] ;
	tot[27] += titrating3[i] ;
	if (i >= i1 && i <= i2)
	  printf ("\n") ;
      }
    printf ("Total") ;
    for (i = 0 ; i < 28 ; i++)
      printf ("\t%8d", tot[i]) ;
    printf ("\n") ;
    
  }

  printf ("\n") ;
  printf ("Now we give the average signal in each category\n") ;
  printf ("Fold change log2(b/a)\tAll probes\tTitrating probes\tTitrating validated probes\tTitrating gene specific probe\tsame for gene with geneId\n") ;
  
  /* 
     col 1: Fold change log2(b/a), gene is mostly expresses in B tissue at end of table
     col 2: number of probes with corresponding fold change
     col 3: titrating probes
     col 4: titrating probes testing an aceview gene
     col 5: titrating probes testing a gene with geneId
  */  
  {
    int i1, i2, n, nn, tot[8] ;
    double stot[8] ;
    for (i = 0 ; i < 8 ; i++)
      { tot[i] = 0 ; stot[i] = 0 ; }
    
    /* on coupe 1% de chaque cote */
    for (i = nn = 0 ; i < NP ; i++)   /* nn = total number of probes */
      nn += nAb [i] ;

    for (i1 = n = 0 ; 100 * n < nn  && i1 < NP ; i1++)
      n += nAb [i1] ;
    for (i2 = NP - 1, n = 0 ; 100 * n < nn  && i2 >= 0 ; i2--)
      n += nAb [i2] ;
    i1-- ; i2++ ;

    for (i = i1 ; i <= i2 ; i++)
      {
	printf ("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
		, (-NC+i) * dSignal
		, sAb [i]/ (nAb [i] ? nAb [i] : 1)
		, (sACDB [i] + sBDCA [i]) / ( nACDB [i] + nBDCA [i]? nACDB [i] + nBDCA [i] : 1)
		, (svACDB [i] + svBDCA [i]) / ( nvACDB [i] + nvBDCA [i]? nvACDB [i] + nvBDCA [i] : 1)
		, (sgACDB [i] +  sgBDCA [i]) / (ngACDB [i] + ngBDCA [i] ? ngACDB [i] + ngBDCA [i]  : 1) 	
		, (sgidACDB [i] +  sgidBDCA [i]) / (ngidACDB [i] + ngidBDCA [i] ? ngidACDB [i] + ngidBDCA [i]  : 1) 	
		) ;
	stot[0] += sAb [i] ;
	stot[1] += sACDB [i] + sBDCA [i] ;
	stot[2] += svACDB [i] + svBDCA [i] ;
	stot[3] += sgACDB [i] + sgBDCA [i] ;
	stot[4] += sgidACDB [i] + sgidBDCA [i] ;

	tot[0] += nAb [i] ;
	tot[1] += nACDB [i] +  nBDCA [i] ;
	tot[2] += nvACDB [i] +  nvBDCA [i] ;
	tot[3] += ngACDB [i] + ngBDCA [i] ;
	tot[4] += ngidACDB [i] + ngidBDCA [i] ;
      }
    for (i = 0 ; i < 7 ; i++)
      if (tot[i] == 0) tot[i] = 1 ;
    printf ("Average\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
	    ,  stot[0]/tot[0], stot[1]/tot[1], stot[2]/tot[2], stot[3]/tot[3], stot[4]/tot[4]
	    ) ;
  }

  if (maqc->aceSignalOut)
    filclose(maqc->aceSignalOut) ;
  maqc->aceSignalOut = 0 ;

  printf ("\n") ;
} /* maqC2 */

/*************************************************************************************/
/*  -c3 c3: export dispersion interlab as a function of TM */
static void maqcC3 (MAQC *maqc, LAB *lab)
{
  int n0, nn, i, ii ;
  int Xn0, Xn1, Xn2 ;
  double Yn1, Yn2 ;
  int tmNn [300], scoreNn [300], nnscore = 0 ;
  int titrating [300], titratable [300], cotitrating [300], nontitrating [300], antititrating [300] ;
  int stitrating [300], stitratable [300], scotitrating [300], snontitrating [300], santititrating [300] ;
  double disp, scoreDisp [300], scoreAA[300], tmDisp [300], AA[300] ;
  TT tt, tt2 ;
  MM *mm ;

  memset (&tt, 0, sizeof(TT)) ;
  memset (&tt2, 0, sizeof(TT)) ;

  for (i = 0 ; i < 300 ; i++)
    {
      tmDisp [i] = 0 ; tmNn [i] = 0 ; AA[i] = 0 ; 
      scoreDisp [i] = 0 ; scoreNn [i] = 0 ; scoreAA[i] = 0 ;  
      titrating [i] = titratable [i] = cotitrating [i] = nontitrating [i] = antititrating [i] = 0 ;
      stitrating [i] = stitratable [i] = scotitrating [i] = snontitrating [i] = santititrating [i] = 0 ;
    }

  if (!maqc->aa || ! arrayMax (maqc->aa))
    {
      printf ("C1: no data files, sorry\n") ;
      return ;
    }
  Xn0 = Xn1 = Xn2 = 0 ;
  Yn1 = Yn2 = 0 ;
  for (ii = nn = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (mm->id != lab->id)
	continue ; 
      i = 1 ;
      if ( mm->A_1[i] != UT_NON_DOUBLE)
	{ Xn1++ ; Yn1 +=  mm->A_1[i] ;}
      else
	Xn0++ ;
      if (mm->tm < 1 || ! mm->score)
	continue ;
      nn++ ;
      nnscore++ ;
      memset (&tt, 0, sizeof(TT)) ;
      memset (&tt2, 0, sizeof(TT)) ;
      for (i = 0 ; i < lab->nReplica ; i++)
	{
	  if (lab->nTissue > 0)
	    {
	      if ( mm->A_1[i] != UT_NON_DOUBLE)
		{ tt.A_1m += mm->A_1[i] ; tt.A_1n++ ;}
	      if ( mm->A_2[i] != UT_NON_DOUBLE)
		{ tt.A_2m += mm->A_2[i] ; tt.A_2n++ ;}
	      if ( mm->A_3[i] != UT_NON_DOUBLE)
		{ tt.A_3m += mm->A_3[i] ; tt.A_3n++ ;}
	    }
	
	  if (lab->nTissue > 1)
	    {
	      if ( mm->B_1[i] != UT_NON_DOUBLE)
		{ tt.B_1m += mm->B_1[i] ; tt.B_1n++ ;}
	      if ( mm->B_2[i] != UT_NON_DOUBLE)
		{ tt.B_2m += mm->B_2[i] ; tt.B_2n++ ;}
	      if ( mm->B_3[i] != UT_NON_DOUBLE)
		{ tt.B_3m += mm->B_3[i] ; tt.B_3n++ ;}
	    }
	  
	  if (lab->nTissue > 2)
	    {
	      if ( mm->C_1[i] != UT_NON_DOUBLE)
		{ tt.C_1m += mm->C_1[i] ; tt.C_1n++ ;}
	      if ( mm->C_2[i] != UT_NON_DOUBLE)
		{ tt.C_2m += mm->C_2[i] ; tt.C_2n++ ;}
	      if ( mm->C_3[i] != UT_NON_DOUBLE)
		{ tt.C_3m += mm->C_3[i] ; tt.C_3n++ ;}
	    }
	  
	  if (lab->nTissue > 3)
	    {
	      if ( mm->D_1[i] != UT_NON_DOUBLE)
		{ tt.D_1m += mm->D_1[i] ; tt.D_1n++ ;}
	      if ( mm->D_2[i] != UT_NON_DOUBLE)
		{ tt.D_2m += mm->D_2[i] ; tt.D_2n++ ;}
	      if ( mm->D_3[i] != UT_NON_DOUBLE)
		{ tt.D_3m += mm->D_3[i] ; tt.D_3n++ ;}
	    }
	  
	  if (lab->nTissue > 0)
	    {
	      if ( mm->A_1[i] != UT_NON_DOUBLE)
		{ tt2.A_1m += mm->A_1[i] * mm->A_1[i] ; tt2.A_1n++ ;}
	      if ( mm->A_2[i] != UT_NON_DOUBLE)
		{ tt2.A_2m += mm->A_2[i] * mm->A_2[i] ; tt2.A_2n++ ;}
	      if ( mm->A_3[i] != UT_NON_DOUBLE)
		{ tt2.A_3m += mm->A_3[i] * mm->A_3[i] ; tt2.A_3n++ ;}
	    }
	  
	  if (lab->nTissue > 1)
	    {
	      if ( mm->B_1[i] != UT_NON_DOUBLE)
		{ tt2.B_1m += mm->B_1[i] * mm->B_1[i] ; tt2.B_1n++ ;}
	      if ( mm->B_2[i] != UT_NON_DOUBLE)
		{ tt2.B_2m += mm->B_2[i] * mm->B_2[i] ; tt2.B_2n++ ;}
	      if ( mm->B_3[i] != UT_NON_DOUBLE)
		{ tt2.B_3m += mm->B_3[i] * mm->B_3[i] ; tt2.B_3n++ ;}
	    }
	  
	  if (lab->nTissue > 2)
	    {
	      if ( mm->C_1[i] != UT_NON_DOUBLE)
		{ tt2.C_1m += mm->C_1[i] * mm->C_1[i] ; tt2.C_1n++ ;}
	      if ( mm->C_2[i] != UT_NON_DOUBLE)
		{ tt2.C_2m += mm->C_2[i] * mm->C_2[i] ; tt2.C_2n++ ;}
	      if ( mm->C_3[i] != UT_NON_DOUBLE)
		{ tt2.C_3m += mm->C_3[i] * mm->C_3[i] ; tt2.C_3n++ ;}
	    }
	  
	  if (lab->nTissue > 3)
	    {
	      if ( mm->D_1[i] != UT_NON_DOUBLE)
		{ tt2.D_1m += mm->D_1[i] * mm->D_1[i] ; tt2.D_1n++ ;}
	      if ( mm->D_2[i] != UT_NON_DOUBLE)
		{ tt2.D_2m += mm->D_2[i] * mm->D_2[i] ; tt2.D_2n++ ;}
	      if ( mm->D_3[i] != UT_NON_DOUBLE)
		{ tt2.D_3m += mm->D_3[i] * mm->D_3[i] ; tt2.D_3n++ ;}
	    }
	}
      disp = 0 ;

      i = tt.A_1n ?  tt.A_1n : 1 ;
      disp += tt2.A_1m/i - (tt.A_1m/i) * (tt.A_1m/i) ;
      i = tt.A_2n ?  tt.A_2n : 1 ;
      disp += tt2.A_2m/i - (tt.A_2m/i) * (tt.A_2m/i) ;
      i = tt.A_3n ?  tt.A_3n : 1 ;
      disp += tt2.A_3m/i - (tt.A_3m/i) * (tt.A_3m/i) ;

      i = tt.B_1n ?  tt.B_1n : 1 ;
      disp += tt2.B_1m/i - (tt.B_1m/i) * (tt.B_1m/i) ;
      i = tt.B_2n ?  tt.B_2n : 1 ;
      disp += tt2.B_2m/i - (tt.B_2m/i) * (tt.B_2m/i) ;
      i = tt.B_3n ?  tt.B_3n : 1 ;
      disp += tt2.B_3m/i - (tt.B_3m/i) * (tt.B_3m/i) ;

      i = tt.C_1n ?  tt.C_1n : 1 ;
      disp += tt2.C_1m/i - (tt.C_1m/i) * (tt.C_1m/i) ;
      i = tt.C_2n ?  tt.C_2n : 1 ;
      disp += tt2.C_2m/i - (tt.C_2m/i) * (tt.C_2m/i) ;
      i = tt.C_3n ?  tt.C_3n : 1 ;
      disp += tt2.C_3m/i - (tt.C_3m/i) * (tt.C_3m/i) ;

      i = tt.D_1n ?  tt.D_1n : 1 ;
      disp += tt2.D_1m/i - (tt.D_1m/i) * (tt.D_1m/i) ;
      i = tt.D_2n ?  tt.D_2n : 1 ;
      disp += tt2.D_2m/i - (tt.D_2m/i) * (tt.D_2m/i) ;
      i = tt.D_3n ?  tt.D_3n : 1 ;
      disp += tt2.D_3m/i - (tt.D_3m/i) * (tt.D_3m/i) ;

      disp /=  (lab->nLab * lab->nTissue) ;

      n0 = tt.A_1n +  tt.A_2n + tt.A_3n + 
	tt.B_1n +  tt.B_2n + tt.B_3n + 
	tt.C_1n +  tt.C_2n + tt.C_3n +
	tt.D_1n +  tt.D_2n + tt.D_3n ;
      if (n0)
	{
	  Xn2 += n0 ; 
	  Yn2 +=  tt.A_1m +  tt.A_2m + tt.A_3m + 
	    tt.B_1m +  tt.B_2m + tt.B_3m + 
	    tt.C_1m +  tt.C_2m + tt.C_3m +
	    tt.D_1m +  tt.D_2m + tt.D_3m ; 
	}
      for (i = 1 ; i < 300 ; i++)
	if (n0 && mm->tm >= i && mm->tm < i + 1)
	  {
	    tmDisp [i] += disp ; tmNn[i]++ ;
	    if (mm->titratable)
	      {
		titratable[i]++ ;
		if (mm->titrating) 
		  {
		    titrating[i]++ ;
		    if (mm->antititrating) 
		      antititrating[i]++ ;
		    if (mm->cotitrating)
		      cotitrating[i]++ ;
		  }
		else
		  nontitrating[i]++ ;
	      }
	    AA[i] += 
	      (
	       tt.A_1m +  tt.A_2m + tt.A_3m + 
	       tt.B_1m +  tt.B_2m + tt.B_3m + 
	       tt.C_1m +  tt.C_2m + tt.C_3m +
	       tt.D_1m +  tt.D_2m + tt.D_3m
	       ) / n0 ; 
	    break ;
	  }
      for (i = 1 ; i < 300 ; i++)
	if (n0 && mm->score >= i && mm->score < i + 1)
	  {
	    if (mm->titratable)
	      {
		stitratable[i]++ ;
		if (mm->titrating) 
		  {
		    stitrating[i]++ ;
		    if (mm->antititrating) 
		      santititrating[i]++ ;
		    if (mm->cotitrating)
		      scotitrating[i]++ ;
		  }
		else
		  snontitrating[i]++ ;
	      }
	    scoreDisp [i] += disp ; scoreNn[i]++ ;
	    scoreAA[i] += 
	      (
	       tt.A_1m +  tt.A_2m + tt.A_3m + 
	       tt.B_1m +  tt.B_2m + tt.B_3m + 
	       tt.C_1m +  tt.C_2m + tt.C_3m +
	       tt.D_1m +  tt.D_2m + tt.D_3m
	       ) / n0 ;
	    break ;
	  }
    }

  if (1)
    {
      float avs = 0 ; 
      int n, ns = 0, i1, i2, j ;
      /* group into bins of 2 */
      if (1)
	{
	  for (i = j = 0 ; i < 300 ; i += 2, j++)
	    {
	      n = tmNn[i] + tmNn[i+1] ;
	      tmNn[j] = n ;
	      n = AA[i] + AA[i+1] ;
	      AA[j] = n ;
	      n = tmDisp[i] + tmDisp[i+1] ;
	      tmDisp[j] = n ;
	      n = titrating [i] + titrating [i+1] ;
	      stitrating [j]  = n ;
	      n = titratable [i] + titratable [i+1] ;
	      titratable [j]  = n ;
	      n = nontitrating [i] + nontitrating [i+1] ;
	      nontitrating [j]  = n ;
	      n = cotitrating [i] + cotitrating [i+1] ;
	      cotitrating [j]  = n ;
	      n = antititrating [i] + antititrating [i+1] ;
	      antititrating [j]  = n ;
	    } 
	  for (;j<300; j++)
	    tmNn[j] = AA[j] = tmDisp[j] = 0 ;
	}
      for (i = 0 ; i < 300 ; i++)
	{
	  ns += tmNn[i] ;
	  if (tmNn[i]) avs += AA[i] ; 
	}
      avs /= ns ? ns : 1 ;
      printf ("Verif Xn0 = %d Xn1 = %d Yn1 = %.2f  yn1/xn1= %.2f  Xn2 = %d Yn2 = %.2f yn2/xn2= %.2f \n"
	      , Xn0, Xn1, Yn1, Yn1/(Xn1>0?Xn1:1)  ,  Xn2, Yn2, Yn2/(Xn2>0?Xn2:1) ) ;
      printf ("Found %d probes with TM in file %s%s, %d with signal average %.2f\n"
	      , nn, lab->lab,lab->bench
	      , ns, avs) ;
      printf ("For each TM, average over all probes with this TM of sum(lab * tissue)(sigma2 inside replica), signal\n") ;

      for (i = i1 = 0 ; i1 == 0 && i < 300 ; i++)
	if (tmDisp[i] > 0 || AA[i] > 0)
	  i1 = i ;
      for (i = i2 = 299 ; i2 == 299 && i >= 0 ; i--)
	if (tmDisp[i] > 0 || AA[i] > 0)
	  i2 = i  ;
      i1 = 25 ; i2 = 50 ; 
      printf ("TM\tlog2(signal)\tNumber of probes\tSignal dispersion") ;
      printf ("\tTitratable\tCotitrating\tNon titrating\tAnti-titrating") ;
      printf ("\tTM\t%% Cotitrating\t%% Non titrating\t%% Anti-titrating") ;
      for (i = i1 ; i <= i2 ; i++)
	{
	  n = tmNn[i] ? tmNn[i] : 1 ;
	  printf ("\n%3d\t%.2f\t%5d\t%.2f", 2*i, AA[i]/n, tmNn[i], tmDisp[i]/n) ;
	  printf ("\t%d\t%d\t%d\t%d"
		  , titratable [i], cotitrating [i]
		  , nontitrating [i], antititrating [i]) ;
	  printf ("\t%d\t%.2f\t%.2f\t%.2f"
		  , 2*i
		  , (100.0*cotitrating [i])/(titratable [i] ?titratable [i] : 1)
		  , (100.0*nontitrating [i])/(titratable [i] ?titratable [i] : 1)
		  , (100.0*antititrating [i])/(titratable [i] ?titratable [i] : 1)
		  ) ;
	}
      printf ("\n\n") ;
    }

  if (1)
    {
      int n, j, i1, i2 ;

      printf ("Found %d probes with score in file %s%s\n", nnscore, lab->lab,lab->bench) ;
      printf ("For each score, average over all probes with this score of sum(lab * tissue)(sigma2 inside replica), signal\n") 
;
      for (i = j = 0 ; i < 300 ; i += 4, j++)
	{
	  n = scoreNn[i] + scoreNn[i+1] + scoreNn[i+2] + scoreNn[i+3] ;
	  scoreNn[j] = n ;
	  n = scoreDisp[i] + scoreDisp[i+1] + scoreDisp[i+2] + scoreDisp[i+3] ;
	  scoreDisp[j] = n ;
	  n = scoreAA[i] + scoreAA[i+1] + scoreAA[i+2] + scoreAA[i+3] ;
	  scoreAA[j] = n ;
	  n = stitrating [i] + stitrating [i+1] + stitrating [i+2] + stitrating [i+3] ;
	  stitrating [j]  = n ;
	  n = stitratable [i] + stitratable [i+1] + stitratable [i+2] + stitratable [i+3] ;
	  stitratable [j]  = n ;
	  n = snontitrating [i] + snontitrating [i+1] + snontitrating [i+2] + snontitrating [i+3] ;
	  snontitrating [j]  = n ;
	  n = scotitrating [i] + scotitrating [i+1] + scotitrating [i+2] + scotitrating [i+3] ;
	  scotitrating [j]  = n ;
	  n = santititrating [i] + santititrating [i+1] + santititrating [i+2] + santititrating [i+3] ;
	  santititrating [j]  = n ;
	}
      for (;j<300; j++)
	scoreNn[j] = scoreDisp[j] = scoreAA[j] = 0 ;
      for (i1 = 0 ; i1 < 300 ; i1++)
	if (scoreNn[i1]) break ;
      for (i2 = 299 ; i2 >= 0 ; i2--)
	if (scoreNn[i2]) break ;
      i1 = 0 ; i2 = 25 ; 
      printf ("Secondary structure\tlog2(signal)\tNumber of probes\tSignal dispersion") ;
      printf ("\tTitratable\tCotitrating\tNon titrating\tAnti-titrating") ;
      printf ("\tSecondary structure\t%% Cotitrating\t%% Non titrating\t%% Anti-titrating") ;
      for (i = i1 ; i <= i2 ; i++)
	{
	  n = scoreNn[i] ? scoreNn[i] : 1 ;
	  printf ("\n%3d\t%.2f\t%5d\t%.2f", 4*i, scoreAA[i]/n, scoreNn[i], scoreDisp[i]/n) ;
	  printf ("\t%d\t%d\t%d\t%d"
		  , stitratable [i], scotitrating [i]
		  , snontitrating [i], santititrating [i]) ;
	  printf ("\t%d\t%.2f\t%.2f\t%.2f"
		  , 4*i
		  , (100.0*scotitrating [i])/(stitratable [i] ?stitratable [i] : 1)
		  , (100.0*snontitrating [i])/(stitratable [i] ?stitratable [i] : 1)
		  , (100.0*santititrating [i])/(stitratable [i] ?stitratable [i] : 1)
		  ) ;
	}
      printf ("\n\n") ;
    }
} /* maqC3 */

/*************************************************************************************/
/*  -c4: number of NM/AceView/geneId hit */
static void maqcC4 (MAQC *maqc, LAB *lab)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn, ii ;
  MM *mm ;
  BitSet bbA, bbNm, bbNmGid ;
 
  nn = 100000 ;
  bbA = bitSetCreate (nn, h) ;
  bbNm = bitSetCreate (nn, h) ;
  bbNmGid = bitSetCreate (nn, h) ;

  for (ii = nn = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      if (! mm->probe)
	continue ;
      if (mm->id != lab->id)
	continue ;
      if (mm->nm) bitSet (bbNm, mm->nm) ;  
      if (mm->geneId) bitSet (bbNmGid, mm->geneId) ;
      if (mm->gene) bitSet (bbA, mm->gene) ;
    }
  printf ("Lab %s%s Hits %lu NM -> %lu GeneId // %lu AceView genes\n\n"
	  , lab->lab, lab->bench
	  , bitSetCount (bbNm)
	  , bitSetCount (bbNmGid)
	  , bitSetCount (bbA)
	  ) ;
  ac_free (h) ;
} /* maqC4 */

/*************************************************************************************/
typedef struct titratingStruct {
  KEYSET ksAll60ChipsGenesAb ; /* NM geneId with A>C>D>B */
  KEYSET ksAll60ChipsGenesaB ; /* NM geneId with A<C<D<B */ 
  KEYSET ksPoolGenesAb, ksPoolGenesaB ;
  int nProbeAb, nProbeaB ;
  int nGeneAb, nGeneaB, nGeneAny, nGeneBoth ;
  int nGeneAbValidatedAb, nGeneaBValidatedaB, nGeneAnyValidatedAny, nGeneValidatedBoth ;
  int nGeneVal60Ab, nGeneVal60aB, nGeneVal60Any, nGeneVal60Both ;
  int nGeneAnyTestedBy3 ;
} TTG ;

static void maqcCountTitratingGenes  (MAQC *maqc, LAB *lab
				      , TTG *ttg
				      , KEYSET discardedChips
				      , BOOL isAll
				      )
{
  AC_HANDLE h = ac_new_handle () ;
  int ns, nt ;
  double s, ss[4], *zp ;
  double dx = .002 ;
  int i, ii, iLab, iTissue, iReplica, k ;
  MM *mm ;
  static Array a[4] = { 0, 0, 0, 0} ;
  int n[4] ;
  KEYSET ksPoolAb, ksPoolaB ;
  int nksPoolAb = 0, nksPoolaB = 0 ;

  ksPoolAb = keySetHandleCreate (h) ;
  ksPoolaB = keySetHandleCreate (h) ;

  if (! a[0])
    for (i = 0 ; i < 4 ; i++)
      a[i] = arrayCreate (60, double) ;

  ttg->nProbeAb = 0 ;
  ttg->nProbeaB = 0 ;
  ttg->nGeneAb = 0 ;
  ttg->nGeneaB = 0 ;
  ttg->nGeneAny = 0 ;
  ttg->nGeneBoth = 0 ;

  ttg->nGeneVal60Ab = ttg->nGeneVal60aB = ttg->nGeneVal60Any = ttg->nGeneVal60Both = 0 ;
  ttg->nGeneAbValidatedAb = ttg->nGeneaBValidatedaB = ttg->nGeneAnyValidatedAny = ttg->nGeneValidatedBoth = 0 ;

  for (ii = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
    {
      /* compute the average signal for this probe using remaining chips */
      memset (ss, 0, sizeof(ss)) ;
      for (i = 0 ; i < 4 ; i++)
	{ n[i] = 0 ; a[i] = arrayReCreate (a[i], 60, double) ; }
      zp = &(mm->A_1[0]) ; 
      if (!strcmp (dictName(maqc->dict,  mm->probe), "NCI_H200001924"))
	invokeDebugger () ;	
      nt = 0 ; /* number of measured tissues */
      for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
	{
	  s = 0 ; ns = 0 ;
	  for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	    for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	      {
		k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
		if (keySet (discardedChips, k)) /* avoid discarded chips */
		  continue ;
		if (zp[k] != UT_NON_DOUBLE)
		  {
		    s += zp[k] ; ns++ ;
		    array(a[iTissue],n[iTissue],double) = zp[k] ;
		    n[iTissue]++ ;
		  }
	      }
	  if (ns) 
	    {
	      ss [iTissue] = s/ns ; nt++ ;
	    } 
	}
      if (nt < 4)
	continue ;
      for (i = 0 ; i < 4 ; i++)
	arraySort (a[i], doubleOrder) ;

      /* count the titrating genes */
      if (ss[0] < ss[2] - dx && ss[2] < ss[3] - dx && ss[3] < ss[1] - dx)
	{ (ttg->nProbeaB)++ ; keySet (ksPoolaB, nksPoolaB++) = mm->gene ; }
      if (ss[0] > ss[2] + dx && ss[2] > ss[3] + dx && ss[3] > ss[1] + dx)
	{ (ttg->nProbeAb)++ ; keySet (ksPoolAb, nksPoolAb++) = mm->gene ; }
    }

  keySetSort (ksPoolAb) ;
  keySetSort (ksPoolaB) ;

  if (strcmp (lab->lab, "AFX"))
    {
      keySetCompress (ksPoolAb) ;
      keySetCompress (ksPoolaB) ;
    }
  else  /* keep only genes seen by at least LIMIT_AFX probes */
    {
      keySetNCompress (ksPoolAb, LIMIT_AFX) ;
      keySetNCompress (ksPoolaB, LIMIT_AFX) ;
    }

  /* record titrating genes */
  ttg->nGeneAb = keySetMax (ksPoolAb) ;
  ttg->nGeneaB = keySetMax (ksPoolaB) ;

  {
    KEYSET ks = 0 ;

    ks = keySetOR (ksPoolAb, ksPoolaB) ;
    ttg->nGeneAny =  keySetMax (ks) ;
    keySetDestroy (ks) ;

    ks = keySetAND (ksPoolAb, ksPoolaB) ;
    ttg->nGeneBoth =  keySetMax (ks) ;
    keySetDestroy (ks) ;
  }

  /* record validated genes */
  {
    KEYSET ks1, ks3, ksAny, ksAnyValidatedAny, ks ;

    ks1 = keySetAND (maqc->titratingGeneAb, ksPoolAb) ;
    ttg->nGeneAbValidatedAb = keySetMax(ks1) ;

    ks3 = keySetAND (maqc->titratingGeneaB, ksPoolaB) ;
    ttg->nGeneaBValidatedaB = keySetMax(ks3) ;

    ksAny = keySetOR (ksPoolAb, ksPoolaB) ;
    ksAnyValidatedAny = keySetOR (ks1, ks3) ;

    ks = keySetAND (maqc->testedBy3, ksAny) ;
    ttg->nGeneAnyTestedBy3 = keySetMax(ks) ? keySetMax(ks) : 1 ; ;
    keySetDestroy (ks) ;

    ttg->nGeneAnyValidatedAny = keySetMax(ksAnyValidatedAny) ;
    keySetDestroy (ks) ;

    ks = keySetAND (ks1, ks3) ;
    ttg->nGeneValidatedBoth = keySetMax(ks) ;
    keySetDestroy (ks) ;
    
    keySetDestroy (ks1) ;
    keySetDestroy (ks3) ;
    keySetDestroy (ksAny) ;
    keySetDestroy (ksAnyValidatedAny) ;
  }

  /* store for later export the case with maximal number of chips */
  if (isAll)
    {
      KEYSET ks ;

      ttg->nGeneVal60Ab = ttg->nGeneAbValidatedAb ;
      ttg->nGeneVal60aB = ttg->nGeneaBValidatedaB ;
      ttg->nGeneVal60Any = ttg->nGeneAnyValidatedAny ;
      ttg->nGeneVal60Both = ttg->nGeneValidatedBoth ;

      ks = ksPoolaB ;
      for (i = 0 ; i < keySetMax(ks) ; i++)
	keySet (ttg->ksAll60ChipsGenesaB, i) = keySet (ks, i) ;
      ks = ksPoolAb ;
      for (i = 0 ; i < keySetMax(ks) ; i++)
	keySet (ttg->ksAll60ChipsGenesAb, i) = keySet (ks, i) ;
    }

  ac_free (h) ;
} /* maqcCountTitratingGenes */

/*************************************************************************************/
/*  -correlation: compute an equivalent of the color plot of leming */
/* Let Vi be the vector with N-probe component of (replica/lab) i
 * we compute barycentre B[j] = <Vi[j]>
 * the centralized vector  Wi = Vi - B
 * and the cosinus Cij =  Wi.Wj / sqrt(Wi.Wi  Wj.Wj
 * we then compute the global score of the line as
 *       Li  = sum(j) (1 - Cij^2) / (sum (j) 1)
 */

static void maqcCorrelation (MAQC *maqc, LAB *lab)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, iLab, iTissue, iReplica, k, zn, znn ;
  int jLab, jTissue, jReplica, k2 ;
  int NNN = NLAB * NTISSUE * NREPLICA ;
  MM *mm ;
  double *zp, zz, zz1, zz2, zz3, B[NNN], Cij[NNN * NNN], CC[NNN] ;
  int firstScore [NNN] ;
  char *tissue[4] = {"A","B","C","D" };
  BOOL debug = FALSE, linear = FALSE, isAll ;
  Array aa = arrayHandleCreate (30000, double, h) ;
  TTG ttg ;
  BOOL isRandom = maqc->random ;
  int nTest = 0 ;
  double snAll = 0, spAll = 0, sn , sp ;



  
  memset (&ttg, 0, sizeof(TTG)) ;
  ttg.ksPoolGenesAb = keySetHandleCreate (h) ;
  ttg.ksPoolGenesaB = keySetHandleCreate (h) ;
  ttg.ksAll60ChipsGenesAb = keySetHandleCreate (h) ;
  ttg.ksAll60ChipsGenesaB = keySetHandleCreate (h) ;

  memset (firstScore, 0, sizeof(firstScore)) ;
  memset (B, 0, sizeof(B)) ;
  memset (CC, 0, sizeof(CC)) ;
  memset (Cij, 0, sizeof(CC)) ;
  if (debug)
    {
      arrayMax (maqc->aa) = 6 ;
      lab->nReplica = 2 ;
      lab->nTissue = 2 ;
    }
  memset (firstScore, 0, sizeof(firstScore)) ;
  memset (B, 0, sizeof(B)) ;
  for (iLab = 0 ; !isRandom && iLab < lab->nLab ; iLab++)
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
      for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	{
	  aa = arrayReCreate (aa, 30000, double) ;
	  k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	  zn = 0 ; zz = 0 ; 
	  if (debug)  printf("\n%s:%d:%d",tissue[iTissue],iLab+1,iReplica+1) ;
	  for (ii = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
	    {
	      if (! mm->probe)
		continue ;
	      if (mm->id != lab->id)
		continue ;
	      zp = &(mm->A_1[0]) ;
	      if (zp[k] != UT_NON_DOUBLE)
		{ 
		  array (aa, zn, double) = zp[k] ; 
		  zn++ ; zz += zp[k] ;
		  if (debug)  printf("\t%.2f",zp[k]) ;
		}
	    }
	  
	  if (zn)
	    {
	      if (linear)
		B[k] = zz/zn ; 
	      else
		{
		  arraySort (aa, doubleOrder) ;
		  zn = arrayMax(aa)/2 ;
		  B[k] = array (aa, zn, double) ;
		}
	    }
	  else
	    B[k] = 0 ;
	  if (debug)  printf("\tn=%d B+=%.2f\n", zn, zz) ;
	}

  for (iLab = 0 ; ! isRandom && iLab < lab->nLab ; iLab++)
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
      for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	{
	  znn = 0 ;
	  k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	  if (debug)  printf("\n%s:%d:%d",tissue[iTissue],iLab+1,iReplica+1) ;
	  for (jLab = 0 ; jLab < lab->nLab ; jLab++)
	    for (jTissue = 0 ; jTissue < lab->nTissue ; jTissue++)
	      for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
		{
		  k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
		  if (debug)  printf("\n\t%s:%d:%d",tissue[jTissue],jLab+1,jReplica+1) ;
		  zz1 = 0 ; zz2 = 0 ; zz3 = 0 ; zn = 0 ;
		  aa = arrayReCreate (aa, 30000, double) ;
		  for (ii = 0, mm = arrp (maqc->aa, 0, MM) ; ii < arrayMax (maqc->aa) ; ii++, mm++)
		    {
		      if (! mm->probe)
			continue ;
		      if (mm->id != lab->id)
			continue ;
		      zp = &(mm->A_1[0]) ;
		      if (zp[k] != UT_NON_DOUBLE &&  zp[k2] != UT_NON_DOUBLE)
			{
			  if (linear) /* linear case, classic cosinus theta */
			    { /* this would be correct in limear data */
			      zz1 += (zp[k] - B[k]) * (zp[k] - B[k]) ;
			      zz2 += (zp[k2] - B[k2]) * (zp[k2] - B[k2]) ;
			      zz3 += (zp[k] - B[k]) * (zp[k2] - B[k2]) ;
			    }
			  else /* log case:  sum of absolute values of fold change */
			    { 
			      zn++ ;
			      zz1 = zp[k] - zp[k2] ;
			      zz3 += zz1 > 0 ? zz1 : -zz1 ;
			      array (aa, zn, double) = zz1 ;
			    } 
			  if (debug)  printf("\n\t\tzp=%.2f zp-B=%.2f zp2=%.2f zp2-B=%.2f zz1=%.2f zz2=%.2f zz2=%.2f"
					     , zp[k], zp[k]- B[k], zp[k2], zp[k2] - B[k2], zz1, zz2, zz3) ;
			}
		    }
		  if (zz1 == 0) zz1 = 1 ;
		  if (zz2 == 0) zz2 = 1 ;
		  /* here we optimize the abs(fold change) around its median */
		  if (!linear && zn)
		    {
		      int i ; double z ;

		      arraySort (aa, doubleOrder) ;
		      i = arrayMax(aa)/2 ;
		      z = arr (aa, i, double) ;
		      zz3 = 0 ;
		      for (i = 0 ; i < zn ; i++)
			{
			  zz1 = arr (aa, i, double) - z ;
			  zz3 += (zz1 > 0 ? zz1 : -zz1) ;
			}
		    }

		  if (zn && iTissue == jTissue)
		    {
		      znn++ ;
		      if (linear)		      
			CC[k] +=  (zz3 * zz3) / (zz1 * zz2) ;
		      else
			CC[k] +=  zz3 / zn ;
		    }
		  if (!zn) zn = 1 ;
		  if (linear)
		    Cij [k * NNN + k2] = zz3/sqrt(zz1*zz2) ;  /* cosinus of experiment k1 and k2 */ 
		  else
		    Cij [k * NNN + k2] = zz3 / zn ;  /* average abs(fold change) */
		  if (debug)  printf("\t==== cos=%.4f", Cij [k * NNN + k2]) ;
		}
	  if (znn > 1)
	    CC[k] /= (znn - 1) ; /* self does not count */
	}
  
  
  printf ("Lab %s%s Hits  COEFFICIENTS de correlation\n"
	  , lab->lab, lab->bench
	  ) ;
  printf("\nTest") ;
  if (0) /* show all the 60 titles */
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
    {
      printf("\t.") ;
      for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	  printf("\t%s:%d:%d",tissue[iTissue],iLab+1,iReplica+1) ; 
    }

  if (0) /* show full correlation table */
    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++)
    { 
      printf("\n") ;
      for (iLab = 0 ; iLab < lab->nLab ; iLab++)
	for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
	  {
	    printf("\n%s:%d:%d",tissue[iTissue],iLab+1,iReplica+1) ; 
	    k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
	    for (jTissue = 0 ; jTissue < lab->nTissue ; jTissue++)
	      {
		printf("\t.") ;
		for (jLab = 0 ; jLab < lab->nLab ; jLab++)
		  for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
		    {
		      k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
		      printf("\t%d", (int)(1000*Cij [k * NNN + k2])) ;
		    } 
	      }
	  }
    }
  printf("\n\n") ;
  
  if (0) /* show all the 60 cumulated columns */
    for (jTissue = 0 ; jTissue < lab->nTissue ; jTissue++)
      {
	printf("\t.") ;
	for (jLab = 0 ; jLab < lab->nLab ; jLab++)
	  for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
	    {
	      k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
	      printf("\t%d", (int)(1000*CC [k2])) ;
	    }
      }
  printf("\n\n") ;
  
  if (1) /* recursivelly remove a column */
    /*  for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++) */
      {
	int wk = 0, nRemoved = 0 ;
	double w = 0 ;
	KEYSET discardedChips = 0 ;
	BOOL showDetails = FALSE ;
	int wTissue = 0,  wLab = 0,  wReplica = 0 ; 
	discardedChips = keySetReCreate (discardedChips) ;

	printf("\n\nnn\tArray\tOriginal score\tScore when removed") ;
	if (showDetails)
	  {
	    for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++) 
	      for (iLab = 0 ; iLab < lab->nLab ; iLab++)
		for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
		  printf("\t%s:%d:%d",tissue[iTissue],iLab+1,iReplica+1) ; 
	  }
	printf("\tTitrating probe A>B\tTitrating probe B>A\tAll titrating probes") ;

	printf("\tTitrating gene A>B\tTitrating gene B>A\tAll titrating genes\tGene titrating in both directions") ;

	printf("\tConsensual titrating gene A>B\tConsensual titrating gene B>A\tAll consensual titrating genes\tGene consensually titrating in both directions") ;

	/* true sensitivity: (nGeneAb seen && validated) / (lab->tested && validated) */
	printf("\tSensitivity : %% seen of tested consensually titrating genes A>B") ; 
	printf("\tSensitivity : %% seen of tested consensually  titrating genes B>A") ; 

	/* real specificity: (nGeneAb seen && validated)  / (nGeneAb && maqc->testedBy3labs) */
	printf("\tSpecificity A>b\tSpecificity a<B") ; 



	/* iterate, print the cost, find worst, eliminate and recompute costs */
	for (nRemoved = 0, isAll = TRUE ; nRemoved < lab->nTissue * lab->nLab * lab->nReplica ;  nRemoved++)
	  {
	    if (nRemoved)
	      printf("\n%d\t%s:%d:%d\t%d\t%d"
		     , nRemoved, tissue[wTissue],wLab+1,wReplica+1
		     , firstScore[wk]
		     , (int) (1000 * w)) ;
	    else
	      {
		printf("\n%d\t \t \t%d", nRemoved, (int) (1000 * w)) ;
		for (k2 = 0 ; k2 < lab->nTissue * NLAB * NREPLICA ; k2++)
		  firstScore[k2] = (int)(1000*CC [k2]) ;
	      }
	    w = -1 ; wk = -1 ;
	    /* non ramdom : select wk = worst platform for removal afetr ptrintout */
	    if (!isRandom)
	      {
		for (jTissue = 0 ; ! isRandom && jTissue < lab->nTissue ; jTissue++) 
		  {
		    for (jLab = 0 ; jLab < lab->nLab ; jLab++)
		      for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
			{
			  k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
			  if (showDetails)
			    printf("\t%d", (int)(1000*CC [k2])) ;
			  if ((w == -1 || CC [k2] > w) &&  CC [k2] > 0 && !keySet (discardedChips, k2))
			    { 
			      w = CC [k2] ; wk = k2 ;
			      wTissue = jTissue ; wLab = jLab ; wReplica = jReplica ; 
			    }
			}
		  }
	      }
	    /* random: remove all but 4 platforms */
	    else
	      {
		for (jTissue = 0 ; jTissue < lab->nTissue ; jTissue++) 
		  for (jLab = 0 ; jLab < lab->nLab ; jLab++)
		    for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
		      {
			k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
			keySet (discardedChips, k2) = 1 ;
		      }
		jLab = randint () % lab->nLab ;
		jReplica = randint ()  % lab->nReplica ;
		for (jTissue = 0 ; jTissue < lab->nTissue ; jTissue++)
		  {
		    k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
		    keySet (discardedChips, k2) = 0 ;
		  }
	      }
	    maqcCountTitratingGenes (maqc, lab, &ttg, discardedChips, isAll) ;
	    isAll = FALSE ;

	    printf("\t%d\t%d\t%d", ttg.nProbeAb, ttg.nProbeaB, ttg.nProbeAb + ttg.nProbeaB) ;
	    printf("\t%d\t%d\t%d\t%d", ttg.nGeneAb, ttg.nGeneaB, ttg.nGeneAny, ttg.nGeneBoth) ;

	    /* validate genes */
	    printf("\t%d\t%d\t%d\t%d"
		   , ttg.nGeneAbValidatedAb, ttg.nGeneaBValidatedaB
		   , ttg.nGeneAnyValidatedAny, ttg.nGeneValidatedBoth) ;

	    /* true sensitivity: (nGeneAb seen && validated) / (lab->tested && validated) */
	    sn = (100.0*ttg.nGeneAnyValidatedAny)/lab->nGeneTestedValidatedAny ;
	    snAll += sn ;
	    printf("\t%.1f", sn) ;

	    /* real specificity: (nGeneAb seen && validated)  / (nGeneAb && maqc->testedBy3labs) */
	    sp = (100.0*ttg.nGeneAnyValidatedAny)/ttg.nGeneAnyTestedBy3 ;
	    spAll += sp ;
	    printf("\t%.1f", sp) ;

	    if (lab->nGeneTestedValidatedAny > 1)
	      nTest++ ;

	    if (! isRandom)
	      {
		if (0 || wk < 0)
		  break ;
		keySet (discardedChips, wk) = 1 ;
		
		for (iTissue = 0 ; iTissue < lab->nTissue ; iTissue++) 
		  for (iLab = 0 ; iLab < lab->nLab ; iLab++)
		    for (iReplica = 0 ; iReplica < lab->nReplica ; iReplica++)
		      {
			k = iReplica + iLab * NREPLICA + iTissue * NLAB * NREPLICA ;
			CC [k] = 0 ; znn = 0 ; 
			jTissue = iTissue ;
			for (jLab = 0 ; jLab < lab->nLab ; jLab++)
			  for (jReplica = 0 ; jReplica < lab->nReplica ; jReplica++)
			    {
			      k2 = jReplica + jLab * NREPLICA + jTissue * NLAB * NREPLICA ; 
			      if (keySet (discardedChips, k) == 0 && keySet (discardedChips, k2) == 0)
				{
				  CC [k] += Cij [k * NNN + k2] ;
				  znn++ ;
				}
			    }
			if (znn > 1)
			  CC [k] /= znn - 1 ;
		      }
	      }	  
	  }
      }
  printf("\n") ;
  if (!nTest) nTest = 1 ;
  printf ("Global sensitivity\t%.2f\tSpecificity\t%.2f\n",
	  snAll/nTest,  spAll/nTest
	  ) ;

  if (1)
    {
      FILE *f ;
      int gene, n1 ;
      KEYSET ksAll60 ;
      int gg ;

      dictFind (maqc->dictGene, "DDEF1", &gg) ;
      if ((f = filopen (messprintf("Titrating%sGenes/Ab_%s", maqc->useNm ? "Nm" : "Acv", lab->lab), "txt", "w")))
	{
	  ksAll60 = ttg.ksAll60ChipsGenesAb ;
	  for (ii = 0 ; ii < keySetMax (ksAll60) ; ii++)
	    {
	      gene =  keySet (ksAll60, ii) ;
	      n1 = keySet (lab->gene2mmAb, gene) ;
	      if (n1 < 0 || n1 >= arrayMax (maqc->aa)) n1 = 0 ;
	      mm = arrp (maqc->aa, n1, MM) ;
	      fprintf (f, "%s\t%s\tA\t%f\t%f\t1\tall60\t\n"
		       , dictName(maqc->dictGene, gene)
		       , lab->lab
		       , mm->signal
		       , mm->foldChange
		       ) ;
	    }
	  filclose (f) ;
	}
      if ((f = filopen (messprintf("Titrating%sGenes/aB_%s", maqc->useNm ? "Nm" : "Acv", lab->lab), "txt", "w")))
	{
	  ksAll60 = ttg.ksAll60ChipsGenesaB ;
	  for (ii = 0 ; ii < keySetMax (ksAll60) ; ii++)
	    {
	      gene =  keySet (ksAll60, ii) ;
	      n1 = keySet (lab->gene2mmaB, gene) ;
	      if (n1 < 0 || n1 >= arrayMax (maqc->aa)) n1 = 0 ;
	      mm = arrp (maqc->aa, n1, MM) ;
	      fprintf (f, "%s\t%s\tB\t%f\t%f\t1\tall60\t\n"
		       , dictName(maqc->dictGene, gene)
		       , lab->lab
		       , mm->signal
		       , mm->foldChange
		       ) ;
	    }
	  filclose (f) ;
	}
    }

  {
    KEYSET ks = keySetOR (ttg.ksAll60ChipsGenesaB, ttg.ksAll60ChipsGenesaB) ;
    if (maqc->useNm)
      printf ("exported %d A> b NM genes, %d a<B NM gene, %d a<B OR A>b\n"
	      , keySetMax (ttg.ksAll60ChipsGenesAb),  keySetMax (ttg.ksAll60ChipsGenesaB)
	      , keySetMax(ks) 
	      ) ;
    else
      printf ("exported %d A> b AceView genes, %d a<B AceView genes, %d a<B OR A>b\n"
	      , keySetMax (ttg.ksAll60ChipsGenesAb),  keySetMax (ttg.ksAll60ChipsGenesaB)
	      , keySetMax(ks) 
	      ) ;
    keySetDestroy (ks) ;
  }

  printf("\n") ;
  ac_free (h) ;
} /* maqcCorrelation */

/*************************************************************************************/
/*  analyse of influence of TM on dye swap in AGL */
static void maqcD1 (MAQC *maqc)
{
  MM *m1, *m2, *mm, *mm2 ;
  LAB *lab1, *lab2 ;
  int i, j, n ;
  double za, za2 ; 
  double zb, zb2 ;
  double t3, t5, correl ;
  int na, nb, tm, tm_min, tm_max ;
  double u, v, uu, vv, uv ;

  for (i = 0, lab1 = labos ; lab1->lab[0] ; lab1++)
    if (*lab1->bench && !strcmp (lab1->lab, "AGL"))
      break ;
  for (lab2 = lab1+1 ; lab2->lab[0] ; lab2++)
    if (*lab2->bench && !strcmp (lab2->lab, "AGL"))
      break ;
  if (!lab1->lab[0] || ! lab2->lab[0]) /* i did not find the 2 colorant for AGL */
    return ;
  tm_min = 1000 ; tm_max = -9 ;
  for (i = 0, mm = arrp (maqc->aa, 0, MM) ; i < arrayMax (maqc->aa) ; i++, mm++)
    {
      if (! mm->lab[0] ||  (*lab1->bench && strcmp (mm->lab, "AGL")))
	continue ;
      if (mm->tm < tm_min) tm_min = mm->tm ;
      if (mm->tm > tm_max) tm_max = mm->tm ;
      /* compute the average of all replicas */
      za = za2 = zb = zb2 = 0 ;
      for (na = nb = j = 0 ; j < lab1->nReplica ; j++)
	{
	  if (mm->A_1[j] != UT_NON_DOUBLE) { za += mm->A_1[j] ; na++ ; }
	  if (mm->A_2[j] != UT_NON_DOUBLE) { za += mm->A_2[j] ; na++ ; }
	  if (mm->A_3[j] != UT_NON_DOUBLE) { za += mm->A_3[j] ; na++ ; }

	  if (mm->B_1[j] != UT_NON_DOUBLE) { zb += mm->B_1[j] ; nb++ ; }
	  if (mm->B_2[j] != UT_NON_DOUBLE) { zb += mm->B_2[j] ; nb++ ; }
	  if (mm->B_3[j] != UT_NON_DOUBLE) { zb += mm->B_3[j] ; nb++ ; }
	}
      za /= na ? na : 1 ;
      zb /= nb ? nb : 1 ;
      mm->A_1m = za ;  mm->B_1m = zb ;
    }
  /* position on the first record of each dye */
   for (i = 0, m1 = 0, m2 = 0, mm = arrp (maqc->aa, 0, MM) ; i < arrayMax (maqc->aa) && (! m1 || ! m2) ; i++, mm++)
     {
       if (! m1 && mm->id == lab1->id) m1 = mm ;
       if (! m2 && mm->id == lab2->id) m2 = mm ;
     }
   /* we now want to compute the correlation matrix between tissue A/B as seen with the 2 dyes, per TM */
   /* compute the 2 averages accross all genes */
  printf ("\n###############################################################\n") ;
  printf ("## TM\t#probes\tCorrel\n") ;
   for (tm = tm_min ; tm <= tm_max ; tm++)
     {
       for (n = 0, t3 = 0, mm = m1 ; mm < m2 ; mm++) 
	 if (mm->tm - tm < 1 && mm->tm >= tm)
	   { n++ ; t3 += mm->A_1m - mm->B_1m ; }
       if (!n) n = 1 ;
       t3 /= n ;
       
       for (n = 0, t5 = 0, mm = m2 ; mm->id == lab1->id ; mm++) 
	 if (mm->tm - tm < 1 && mm->tm >= tm)
	   { n++ ; t5 += mm->A_1m - mm->B_1m ; }
       if (!n) n = 1 ;
       t5 /= n ;

       uu = vv = uv = 0 ;
       for (mm = m1, mm2 = m2 ; mm < m2 ; mm++, mm2++)
	 {
	   u = v = 0 ;
	   if (mm->tm - tm < 1 && mm->tm >= tm)
	     {
	       u += mm->A_1m - mm->B_1m - t3 ;
	       v += mm2->A_1m - mm2->B_1m - t5 ;
	       uu += u * u ; vv += v * v ; uv += u * v ;
	     }
	 }
       if (tm > 0 && n > 1)
	 {
	   correl = uv / sqrt (uu * vv) ;
	   printf ("%3d\t%5d\t%.2f\n", tm, n, correl) ;
	 }
     }

} /* maqcD1 */

/*************************************************************************************/

static void maqcVenn (MAQC *maqc, KEYSET ks2, KEYSET ks4, KEY big, BOOL isUp, BOOL isAceView)
{
  AC_HANDLE h = ac_new_handle () ;
  LAB *lab, *mylab[3] ;
  BitSet bs[3], bb ;
  KEY *kp2, *kp4, *kg ;
  KEYSET tested = isAceView ? maqc->testedGene : maqc->testedGeneId ;
  int i, ii, n, n1, n2 ;

  for (ii = 0, lab = labos ; ii < 3 && lab->lab[0] ; lab++)
    if (lab->id & big) mylab[ii++] = lab ;

  printf ("Venn Diagram %s between: ", isUp ? "A<B" : "A>B") ;
  for (ii = 0 ; ii < 3 ; ii++)
    printf (" %s%s", mylab[ii]->lab, mylab[ii]->bench) ; 
  printf ("\n") ;

  n = keySetMax (ks2) ;
  for (ii = 0 ; ii < 3 ; ii++)
    bs[ii] = bitSetCreate (n, h) ;

  for (i = 0, kp2 = arrp (ks2, 0, KEY), kp4 = arrp (ks4, 0, KEY), kg = arrp (tested, 0, KEY) ;
       i < keySetMax (ks2) ; kg++, kp2++, kp4++, i++)
    if ((big & *kg) == big)
      for (ii = 0 ; ii < 3 ; ii++)
	if ((*kp4 >> (isUp ? 16 : 0)) &&                  /* somebody is very strong */
	    ((*kp2 >> (isUp ? 16 : 0)) & mylab[ii]->id))  /* then assess even weak guys */
	  bitSet (bs[ii], i) ;
  for (ii = 0 ; ii < 3 ; ii++)
    printf ("%s:%lu ", mylab[ii]->lab, bitSetCount (bs[ii])) ; 
  printf ("\n") ; 
  for (ii = 0 ; ii < 3 ; ii++) 
    {
      bb = bitSetCreate (n, h) ;
      bitSetOR (bb, bs[(ii + 1) % 3]) ;
      bitSetOR (bb, bs[(ii + 2) % 3]) ;
      printf ("%s OR %s:%lu ", mylab[(ii+1) % 3]->lab, mylab[(ii+2) % 3]->lab, bitSetCount (bb)) ;
    }
  printf ("\n") ; 
  for (ii = 0 ; ii < 3 ; ii++) 
    {
      bb = bitSetCreate (n, h) ;
      bitSetOR (bb, bs[(ii + 1) % 3]) ;
      bitSetAND (bb, bs[(ii + 2) % 3]) ;  
      bitSetMINUS (bb, bs[ii]) ;
      printf ("%s AND %s NOT %s:%lu ", mylab[(ii+1) % 3]->lab, mylab[(ii+2) % 3]->lab,  mylab[ii]->lab, bitSetCount (bb)) ;
    }
  printf ("\n") ; 
  for (ii = 0 ; ii < 3 ; ii++) 
    {
      bb = bitSetCreate (n, h) ;
      bitSetOR (bb, bs[ii]) ;
      bitSetMINUS (bb, bs[(ii + 1) % 3]) ;
      bitSetMINUS (bb, bs[(ii + 2) % 3]) ;  
      printf ("%s NOT %s NOT %s:%lu ", mylab[ii]->lab, mylab[(ii+1) % 3]->lab, mylab[(ii+2) % 3]->lab, bitSetCount (bb)) ;
    }
  printf ("\n") ; 
  bb = bitSetCreate (n, h) ;
  bitSetOR (bb, bs[0]) ;
  bitSetOR (bb, bs[1]) ;
  bitSetOR (bb, bs[2]) ;
  n1 = bitSetCount (bb) ;

  printf ("union of all3\t%d\n", n1) ;

  bb = bitSetCreate (n, h) ;
  bitSetOR (bb, bs[0]) ;
  bitSetAND (bb, bs[1]) ;
  bitSetAND (bb, bs[2]) ;
  n2 = bitSetCount (bb) ;

  printf ("triple intersect\t%d\t%.2f%%\n\n", n2, (100.0 * n2) / (n1 ? n1 : 1)) ;
  ac_free (h) ;
} /* maqcVenn */

/*************************************************************************************/
/*  analyse concistency between platforms for genes differentially expressed */
/* Venn diagrams of platform triplets and #differentially expressed genes */
/* common to 7 6 ..1 lab , export the found_in_1.list gene list files */
 
static void maqcD2 (char *title, MAQC *maqc, KEYSET ks2, KEYSET ks4, BOOL isAceView)
{
  int n, i, ii, naB, nAb, n1, n2, nbad, nbadinterne, nGene1 = 0, nGene2 = 0 ;
  int over1[32], over2[32], overCommon1[32], overCommon2[32] ;
  KEY key, k1, k2, k0, aB = 0xffff, Ab = 0xffff0000 ;
  KEYSET tested = isAceView ? maqc->testedGene : maqc->testedGeneId ;
  FILE *outf1 = NULL ;
  KEYSET ksG[300] ;
  KEYSET ksGdiff[300] ;
  
  for (n1 = 0 ; n1 < 300 ; n1++)
    {
      ksG[n1] = keySetCreate () ;
      ksGdiff[n1] = keySetCreate () ;
    }
  if (isAceView && dictMax (maqc->dictGene) > 10)
    outf1 = filopen ("found_in_1", "ace", "w") ;
  n = keySetMax (ks2) ;
  if (keySetMax (ks4) < n)
    keySet (ks4,n-1) = 0 ; /* adjust sizes */
  n = keySetMax (ks4) ;
  if (keySetMax (ks2) < n)
    keySet (ks2,n-1) = 0 ; /* adjust sizes */
			 
  for (i = 0 ; i < 32 ; i++)
    overCommon1[i] = overCommon2[i] = over1[i] = over2[i] = 0 ;
  nAb = naB = nbad = nbadinterne = 0 ;

  printf ("\n###############################################################\n") ;
  printf ("## Discordant genes at %.2f sigma \n", maqc->nSigma) ;
  for (ii = 0 ; ii < keySetMax (tested); ii++)
    {
      key = keySet (ks4, ii) ; /* should be really different at least once */
      if (0 && !key)
	continue ;
      {
	KEY tt = keySet (tested, ii) ;
	for (n1 = i = 0, k0 = 0x1; i < 15 ; k0 <<= 1, i++)
	  if (k0 & tt & big6)
	    n1++ ;
	if (outf1 && n1 >= 1)
	  keySet (ksG[n1],keySetMax(ksG[n1])) = ii ;
      }
      nGene1++ ;
      if ((big6 & keySet (tested, ii)) == big6)
	nGene2++ ;
      key = keySet (ks2, ii) ;
      k1 = ((key & Ab) >> 16) & aB ;
      k2 = (key & aB) ;
      k1 &= big6 ;       k2 &= big6 ; 
      if (k1) nAb++ ;
      if (k2) naB++ ;
      if (1 && k1 && k2)      /* put 0 to authorize conflicts */
	{
	  {
	    nbad++ ;
	    if (0 && dictMax(maqc->dictGeneId) >= ii)
	      printf ("  %s\n", dictName (maqc->dictGeneId, ii)) ;
	  }
	  if (k1 & k2) 
	    {
	      nbadinterne++ ; 
	      if (! isAceView && dictMax(maqc->dictGeneId) >= ii)
		printf ("  %s\n", dictName (maqc->dictGeneId, ii)) ;
	    }
	}
      else
	{
	  /* count the libs where we have overexpression */
	  for (n1 = n2 = i = 0, k0 = 0x1; i < 15 ; k0 <<= 1, i++)
	    if (k0 & big6)
	      {
		if (k0 & k1) n1++ ;
		if (k0 & k2) n2++ ;
	      }
	  over1[n1]++ ;
	  over2[n2]++ ;
	  if ((big6 & keySet (tested, ii)) == big6)
	    {
	      overCommon1[n1]++ ;
	      overCommon2[n2]++ ;
	    }
	  if (outf1 && (n1 >= 1 || n2 >= 1))
	    keySet (ksGdiff[n1],keySetMax(ksGdiff[n1])) = ii ;
	}
    }
  printf ("\n") ;
  if (outf1) 
    {
      unsigned int xs, xs0, xs1 ;
      for (n1 = 1 ; n1 < 32 ; n1++)
	{
	  if (keySetMax (ksG[n1]))
	    {
	      fprintf (outf1, "KeySet GeneTouchedBy%d\n", n1) ;
	      for (ii = 0 ; ii < keySetMax (ksG[n1]) ; ii++)
		fprintf (outf1, "Gene \"%s\"\n", dictName (maqc->dictGene, keySet (ksG[n1],ii))) ;
	      fprintf (outf1, "\n") ;
	    }
	}
      for (n1 = 1 ; n1 < 32 ; n1++)
	{
	  if (keySetMax (ksG[n1]))
	    {
	      for (ii = 0 ; ii < keySetMax (ksG[n1]) ; ii++)
		fprintf (outf1, "Gene \"%s\"\nMAQC_N_test %d\n\n"
			 , dictName (maqc->dictGene, keySet (ksG[n1],ii))
			 , n1) ;
	    }
	  keySetDestroy (ksG[n1]) ;
	}
      for (n1 = 1 ; n1 < 32 ; n1++)
	{
	  if (keySetMax (ksGdiff[n1]))
	    {
	      for (ii = 0 ; ii < keySetMax (ksGdiff[n1]) ; ii++)
		fprintf (outf1, "Gene \"%s\"\nMAQC_N_acedb %d\n\n"
			 , dictName (maqc->dictGene, keySet (ksGdiff[n1],ii))
			 , n1) ;
	    }
	  keySetDestroy (ksGdiff[n1]) ;
	}

      for (n1 = 0 ; n1 < 32 ; n1++)
	{
	  xs0 = xs1 = 0x1 << n1 ;
	  if (n1 < 31) xs1 <<= 1 ;
	  fprintf (outf1, "KeySet GeneMaqcLevel%d\n", n1) ;
	  for (ii = 0 ; ii < keySetMax (maqc->signalGene) ; ii++)
	    {
	      xs = keySet (maqc->signalGene, ii) ;
	      if ((xs & xs0) && (n1 == 31 || xs < xs1))
		fprintf (outf1, "Gene \"%s\"\n", dictName (maqc->dictGene, ii)) ;
	    }
	  fprintf (outf1, "\n") ;
	}
      for (ii = 0 ; ii < keySetMax (maqc->signalGene) ; ii++)
	{
	  xs = keySet (maqc->signalGene, ii) ;
	  for (n1 = 0 ; n1 < 32 ; n1++)
	    {
	      xs0 = xs1 = 0x1 << n1 ;
	      if (n1 < 31) xs1 <<= 1 ;
	      if ((xs & xs0) && (n1 == 31 || xs < xs1))
		fprintf (outf1, "Gene \"%s\"\nMAQC_max_signal %d\n\n"
			 , dictName (maqc->dictGene, ii)
			 , n1) ;
	    }
	}
      fprintf (outf1, "\n") ;

      filclose (outf1) ;
    }

  printf ("We count the geneId where in at least one platform is ACeDB with A/B > 2 -> %d genes or B/A > 2 -> %d genes\n"
	  , nAb, naB) ;
  printf ("There are %d genes where one platform is also ACeDB but gives the opposite order\n", nbad) ;
  printf ("There are %d genes where another probe of the same platform is also ACeDB but gives the opposite order\n", nbadinterne) ;
  printf ("Number of genes where A is underexpressed, then overexpressed in at least n of the big platforms\n"
	  "not counting the genes where there is at least one contradiction\n"
	  "\n"
	  "Col 1/2 = number of genes tested by at least 1 of the platfors, 3/4 by all\n"
	  ) ;
  /* cumulate the data */
  for (n1 = n2 = 0, i = 31 ; i >= 0 ; i--)
    {
      n1 += over1[i] ; n2 += over2[i] ;
      over1[i] = n1 ; over2[i] = n2 ;
    }
  for (n1 = n2 = 0, i = 31 ; i >= 0 ; i--)
    {
      n1 += overCommon1[i] ; n2 += overCommon2[i] ;
      overCommon1[i] = n1 ; overCommon2[i] = n2 ;
    }
  over1[0] = over2[0] = nGene1 ;
  overCommon1[0] = overCommon2[0] = nGene2 ;
  printf ("\n%s\t%s\t%s\n", title, "a < B", "A > b") ;
  for (i = 0 ; i < 9 ; i++)
    printf ("%d\t%d\t%d\t%d\t%d\n"
	    , i
	    , over1[i], over2[i]
	    , overCommon1[i], overCommon2[i]
	    ) ;
  printf ("\n") ;
  overCommon2[0] = 0 ; over2[0] = 0 ;
  printf ("\n%s\t%s\t%s\n", title, "a <> B", "common") ;
  for (i = 0 ; i < 9 ; i++)
    printf ("%d\t%d\t%d\n"
	    , i
	    , over1[i] + over2[i]
	    , overCommon1[i] + overCommon2[i]
	    ) ;
  printf ("\n") ;
  /* various Venn */
  /* big6 = ABI | AFX_123 | AGL_cy3 | AGL_cy5 | GEH | ILM */

  maqcVenn (maqc, ks2, ks4, ABI | AGL_cy3 | AGL_cy5, 0, isAceView) ;
  maqcVenn (maqc, ks2, ks4, ABI | AGL_cy3 | AGL_cy5, 1, isAceView) ;

  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | AFX_456, 0, isAceView) ;
  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | AFX_456, 1, isAceView) ;

  maqcVenn (maqc, ks2, ks4, ABI | GEH | ILM, 0, isAceView) ;
  maqcVenn (maqc, ks2, ks4, ABI | GEH | ILM, 1, isAceView) ;

  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | ILM, 0 , isAceView) ; 
  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | ILM, 1 , isAceView) ; 

  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, ABI | AFX_123 | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, ABI | NCI | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, ABI | NCI | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, AFX_123 | ILM | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, AFX_123 | ILM | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, AFX_456 | ILM | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, AFX_456 | ILM | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, AFX_123 | ILM | NCI, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, AFX_123 | ILM | NCI, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, AFX_123 | AGL_cy3 | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, AFX_123 | AGL_cy3 | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, ILM | AGL_cy3 | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, ILM | AGL_cy3 | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, AFX_123 | AGL_cy5 | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, AFX_123 | AGL_cy5 | GEH, 1 , isAceView) ;

  maqcVenn (maqc, ks2, ks4, ILM | AGL_cy5 | GEH, 0 , isAceView) ;
  maqcVenn (maqc, ks2, ks4, ILM | AGL_cy5 | GEH, 1 , isAceView) ;
} /* maqcD2 */

/*************************************************************************************/
/* square correlation table of all the big6 methods */
#define NLAB1 32
#define NLAB2 1024  /* NLAB1 * NLAB1 */
#define cc(_i1,_i2) correl[(_i1)*NLAB1 + (_i2)]
static void maqcD3 (char *title, MAQC *maqc, KEYSET ks2, KEYSET ks4, BOOL dx, BOOL isAceView)
{
  double correl [NLAB2] ;
  LAB *l1, *l2 ;
  int ii, ii1, ii2 ;
  KEY key, key4, test, mask1, mask2 ;
  KEYSET tested = isAceView ? maqc->testedGene : maqc->testedGeneId ;

  for (ii = 0 ; ii < NLAB2 ; ii++)
    correl [ii] = 0 ;
  /* compute in correl the number of genes flagged both in lab1-lab2 */
  for (ii = 0 ; ii < keySetMax (ks2); ii++)
    {
      test = keySet (tested, ii) ;
      key = keySet (ks2, ii) ;
      key4 = keySet (ks4, ii) ;
      if (! key)
	continue ;
      for (ii1 = 0, l1 = labos ; l1->id ; ii1++, l1++)
	{
	  if (!(big6 & l1->id))
	    continue ;
	  mask1 = l1->id | (l1->id << 16) ;
	  if (! (test & mask1))
	    continue ; 
	  if ((dx && (key4 & mask1)) || (!dx && (key & mask1)))
	    cc(ii1,ii1)++ ;

	  for (ii2 = ii1+1, l2 = l1 + 1 ; l2->id ; ii2++, l2++)
	    {
	      if (!(big6 & l2->id))
		continue ;
	      mask2 = l2->id | (l2->id << 16) ;
	      if (!(test & mask2))
		continue ;
	      if ( (!dx || (dx && (key4 & (mask1 | mask2)))) && (key & (mask1 | mask2)))
		{
		  cc(ii2,ii1)++ ; /* union tested by both */
		  if ((key & mask1) && (key & mask2))
		    cc(ii1,ii2)++ ; /* intersect */
		}
	    }
	}
    } 
  /* print-out */
  printf ("%s", title) ;
  for (ii1 = 0, l1 = labos ; l1->id ; ii1++, l1++)
    if (big6 & l1->id)
      printf ("\t%s%s", l1->lab, l1->bench) ;
  for (ii1 = 0, l1 = labos ; l1->id ; ii1++, l1++)
    if (big6 & l1->id)
      {
	printf ("\n%s%s", l1->lab, l1->bench) ;
	for (ii2 = 0, l2 = labos ; l2->id ; ii2++, l2++)
	  if (big6 & l2->id)
	    {
	      if (ii1 < ii2) printf ("\t%d"
				     , (int)cc(ii1,ii2)
				    ) ;
	      else printf ("\t%d", (int)cc(ii1,ii2)) ;
	    }
      }
  printf ("\n\n") ;   
  printf ("%s", title) ;
  for (ii1 = 0, l1 = labos ; l1->id ; ii1++, l1++)
    if (big6 & l1->id)
      printf ("\t%s%s", l1->lab, l1->bench) ;
  for (ii1 = 0, l1 = labos ; l1->id ; ii1++, l1++)
    if (big6 & l1->id)
      {
	printf ("\n%s%s", l1->lab, l1->bench) ;
	for (ii2 = 0, l2 = labos ; l2->id ; ii2++, l2++)
	  if (big6 & l2->id)
	    {
	      if (ii1 < ii2) printf ("\t%.1f%%"
				     , 100.0 *  cc(ii1,ii2) / (cc(ii2,ii1) > 0 ? cc(ii2,ii1) : 0)) ; 
	      else printf ("\t%d", (int)cc(ii1,ii2)) ;
	    }
      }
  printf ("\n\n") ;   
} /* maqcD3 */

/*************************************************************************************/
/*************************************************************************************/

static BOOL maqcCountConfirmedGenesTitratingProbes_c9 (MAQC *maqc)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp, **lab, **lab2
    , *mylabs[] = {  "ABI",       "AFX", "AGL", "GEH", "ILM", "NCI","EPP", "GEX", "QGN","TAQ","*", 0 } 
    , *mylabs2[] = { "ABI","AFX-set-4p", "AGL", "GEH", "ILM", "NCI","EPP", "GEX", "QGN","TAQ","*", 0 }
  ;
  const char *qType[] = {
      "Find Gene maqc_acdb_plus && ! maqc_acdb_minus"
    , "Find Gene ! maqc_acdb_plus &&  maqc_acdb_minus"
    , "Find Gene maqc_acdb_plus && maqc_acdb_minus"
    , 0 } ;

  /* limit the study to genes tested by a given platform (taq or qgn)
    const char *qTypeZ[] = { 
    "Find Probe qgn_* ; >Confirmed_gene  maqc_acdb_plus && ! maqc_acdb_minus"
    , "Find Probe qgn_* ; >Confirmed_gene  ! maqc_acdb_plus &&  maqc_acdb_minus"
    , "Find Probe qgn_* ; >Confirmed_gene  maqc_acdb_plus && maqc_acdb_minus"
    , 0 } ;
  */
  const char *qSubType[] = {
    "acdb_plus"
    , "acdb_minus"
    , "! acdb_plus  && ! acdb_minus"
    , 0 } ;

  AC_KEYSET ks, ks1, ksInter = 0, ksInterExact = 0, gks[100] ;
  int n, n0, iLab, gType, subType, nn[8], nnn[100] ;

  printf ("\n\n") ;
  printf ("Table per platform of the consistency of the gene specific titrating probes\n") ;
  printf ("Platform"
	  "\tGenes over-expressed in A tested in the platform"
	  "\tProbes"
	  "\tCo-titrating probes"
	  "\tDiscordant probes"
	  "\tNon titrating probes"
	  "\tPlatform"
	  "\t%% co-titrating probes in genes over-expressed in A "
	  "\t%% discordant probes"
	  "\t%% non titrating probes"
	  "\tPlatform"
	  "\tGenes over-expressed in B tested in the platform"
	  "\tProbes"
	  "\tCo-titrating probes"
	  "\tDiscordant probes"
	  "\tNon titrating probes"
	  "\tPlatform"
	  "\t%% co-titrating probes in genes over-expressed in B "
	  "\t%% discordant probes"
	  "\t%% non titrating probes"
	  "\tPlatform"
	  "\tGenes over-expressed both ways tested in the platform"
	  "\tProbes"
	  "\tCo-titrating probes"
	  "\tDiscordant probes" 
	  "\tNon titrating probes"
	  "\tPlatform"
	  "\t%% co-titrating probes in genes over-expressed both way"
	  "\t%% discordant probes"
	  "\t%% non titrating probes"
	  "\tPlatform"
	  "\tGenes differentially expressed tested in the platform"
	  "\tProbes"
	  "\tCo-titrating probes"
	  "\tDiscordant probes" 
	  "\tNon titrating probes"
	  "\tPlatform"
	  "\t%% co-titrating probes in genes differentially expressed"
	  "\t%% discordant probes"
	  "\t%% non titrating probes"
	  ) ;

  for (gType = 0 ; gType < 3 ; gType++)
    {
      gks[10*gType + 0] = ac_dbquery_keyset (maqc->db, qType[gType], h) ;
      gks[10*gType + 1] = ac_ksquery_keyset (gks[10*gType + 0]
					     , "{>maqc_probe_set ; IS AFX*} SETOR {>maqc_probe ;NOT IS AFX* && NOT IS MLT_*}"
					     , h) ;
      for (subType = 2 ; subType < 5 ; subType++)
	gks[10*gType + subType] = ac_ksquery_keyset (gks[10*gType + 1], qSubType[subType-2], h) ; 
    }

  if (1)
    { 
      /* flip the A<B A>B probes in A<B genes */
      ks = gks[10*1+2]; gks[10*1+2]= gks[10*1+3] ; gks[10*1+3] = ks ; 
      /* fuse the A<B and A>B probes in the double titrating genes */
      ac_keyset_or (gks[10*2+2], gks[10*2+3]) ; gks[10*2+3] = ac_dbquery_keyset (maqc->db, "Find probe zztoto", h) ;
    }

  for (lab = mylabs, lab2 = mylabs2 ; *lab ; lab++, lab2++)
    {
      for (gType = 0 ; gType < 3 ; gType++)
	{
	  printf ("%s%s", gType ? "\t" : "\n", *lab2) ;
	  ccp = *lab ;
	  if (*ccp == '*')
	    ks = ac_ksquery_keyset (gks[10*gType + 0], "Maqc_probe != MLT_*", h) ;
	  else 
	    ks = ac_ksquery_keyset (gks[10*gType + 0], messprintf ("Maqc_probe = %s_*", *lab), h) ;
	  n = ac_keyset_count (ks) ;
	  ac_free (ks) ;

	  printf ("\t%d", n) ;
	  nnn[10*gType] = n ;
	  for (n0 = 0, subType = 1 ; subType < 5 ; subType++)
	    {
	      ccp = *lab ;
	      if (*ccp == '*')
		ks = ac_ksquery_keyset (gks[10*gType + subType], "! IS MLT_*", h) ;
	      else 
		ks = ac_ksquery_keyset (gks[10*gType + subType], messprintf ("IS %s_*", *lab), h) ;
	      n = ac_keyset_count (ks) ;
	      ac_free (ks) ;
	      printf ("\t%d", n) ;
	      nn[subType] = n ;
	      nnn[10*gType+subType] = n ;
	      if (subType == 1) n0 = n ? n : 1 ;
	    }

	  printf ("\t%s", *lab2) ;
	  for (subType = 2 ; subType < 5 ; subType++)
	    printf ("\t%.2f", (100.0 * nn[subType])/n0) ;
	} 
      printf ("\t%s", *lab) ;
      n = nnn[0] + nnn[10] + nnn[20] ; printf ("\t%d", n) ;
      n = nnn[1] + nnn[11] + nnn[21] ; printf ("\t%d", n) ;
      n = nnn[2] + nnn[12] + nnn[22] ; printf ("\t%d", n) ;
      n = nnn[3] + nnn[13] + nnn[23] ; printf ("\t%d", n) ;
      n = nnn[4] + nnn[14] + nnn[24] ; printf ("\t%d", n) ;
 
      printf ("\t%s", *lab2) ;
      n0 = nnn[1] + nnn[11] + nnn[21] ;
      n = nnn[2] + nnn[12] + nnn[22] ;
      printf ("\t%.2f", (100.0 * n)/n0) ;
      n = nnn[3] + nnn[13] + nnn[23] ;
      printf ("\t%.2f", (100.0 * n)/n0) ;
      n = nnn[4] + nnn[14] + nnn[24] ;
      printf ("\t%.2f", (100.0 * n)/n0) ;

    }
  printf ("\n\n") ;

  printf ("\nPlatform"
	  "\tEntrez gene with exact match"
	  "\tNovel gene with exact match"
	  "\tEntrez gene with gene specific probe"
	  "\tNovel gene with gene specific probe"
	  "\n"
	  ) ;

  for (iLab = 0, lab = mylabs ; 1 && *lab ; iLab++, lab++)
    {
      /* exact match gene with geneid or novel */
      ccp = *lab ;
      if (*ccp == '*')
	ks = ac_dbquery_keyset (maqc->db, "Find probe ! IS MLT_* && (NM_exact_hit || mrna_exact_hit) && probe_grade1 ; >confirmed_gene" , h) ;
      else if (strcmp (ccp, "AFX"))
	ks = ac_dbquery_keyset (maqc->db,
				messprintf ("Find probe IS %s_* && (NM_exact_hit || mrna_exact_hit) && probe_grade1 ; >confirmed_gene"
					    , *lab) , h) ;
      else
	ks = ac_dbquery_keyset (maqc->db, "Find probe_set AFX* confirmed_gene && COUNT {>probe  NM_exact_hit || mrna_exact_hit} >= 9; >confirmed_gene", h) ;
      ks1 = ac_ksquery_keyset (ks, "geneid", h) ;
      n0 = ac_keyset_count (ks) ;
      n = ac_keyset_count (ks1) ;
      ccp = *lab ;
      if (*ccp == '*') ccp = "Union" ;
      else if (!strcmp (ccp, "AFX")) ccp = "AFX probe set with >= 9 exact hit" ;
      printf ("\t%s\t%d\t%d", ccp, n, n0 - n) ; 
      if (!iLab)
	ksInterExact = ac_copy_keyset (ks, h) ;
      else if (iLab < 6)
	ac_keyset_and (ksInterExact, ks) ;	
       ac_free (ks) ;  ac_free (ks1) ;
      /* gene with gene specific probe */
      if (!strcmp (ccp, "AFX"))
	ks = ac_dbquery_keyset (maqc->db, "Find gene ; maqc_probe_set == AFX", h) ;
      else if (*ccp == '*')
	ks = ac_dbquery_keyset (maqc->db, "Find gene ; (maqc_probe != mlt_*) || maqc_probe_set", h) ;
      else
	ks = ac_dbquery_keyset (maqc->db, messprintf ("Find gene ; maqc_probe == %s*", *lab), h) ;
      if (!iLab)
	ksInter = ac_copy_keyset (ks, h) ;
      else if (iLab < 6)
	ac_keyset_and (ksInter, ks) ;	
      ks1 = ac_ksquery_keyset (ks, "geneid", h) ;
      n0 = ac_keyset_count (ks) ;
      n = ac_keyset_count (ks1) ;
      printf ("\t%d\t%d", n, n0 - n) ; 
      ac_free (ks) ;  ac_free (ks1) ;
      printf ("\n") ;
    }
  if (ksInterExact)
    {
      n0 = ac_keyset_count (ksInterExact) ; 
      ks1 = ac_ksquery_keyset (ksInterExact, "geneid", h) ; 
      n = ks1 ? ac_keyset_count (ks1) : 0 ;
      ac_free (ks1) ;
      printf ("\nIntersection\t%d\t%d", n, n0 - n) ; 
      n0 = ac_keyset_count (ksInter) ; 
      ks1 = ac_ksquery_keyset (ksInter, "geneid", h) ; 
      n = ac_keyset_count (ks1) ;
      ac_free (ks1) ;
      ac_free (ksInterExact) ; ac_free (ksInter) ;
      printf ("\t%d\t%d", n, n0 - n) ; 
    }
  printf ("\n\n") ;

  ac_free (h) ;

  return TRUE ;
} /*  maqcCountGenesTitratingProbes_c9 */

/*************************************************************************************/
#ifdef JUNK
static BOOL maqcCountAllGenesTitratingProbes_c10 (MAQC *maqc)
{
  /*
    AC_HANDLE h = ac_new_handle () ;
  const char *ccp, **lab, **lab2
    , *mylabs[] = {  "ABI",       "AFX", "AGL", "GEH", "ILM", "NCI","EPP", "GEX", "QGN","TAQ","*", 0 } 
    , *mylabs2[] = { "ABI","AFX-set-4p", "AGL", "GEH", "ILM", "NCI","EPP", "GEX", "QGN","TAQ","*", 0 }
  ;
  const char *qType[] = {
      "Find Gene maqc_acdb_plus && ! maqc_acdb_minus"
    , "Find Gene ! maqc_acdb_plus &&  maqc_acdb_minus"
    , "Find Gene maqc_acdb_plus && maqc_acdb_minus"
    , 0 } ;
  const char *qSubType[] = {
    "acdb_plus"
    , "acdb_minus"
    , "! acdb_plus  && ! acdb_minus"
    , 0 } ;
  */

  return TRUE ;
} /*  maqcCountGenesTitratingProbes_c9 */
#endif
/*************************************************************************************/
/* export a general table
 *   genecategory category
 *     p = probed    by 1 to 14 platforms
 *     e = expressed in 1 to 14 platforms
 *     t = titrating in 1 to 14 platforms
 *     g = category
 *       0: has NM
 *       1: has GeneId
 *       2: has neither
 *       3: any
 */
   
#define CG(_p,_e,_t,_g) _CG[_g + 5 * _t + 100 * _e + 20000 * _p]   
 
typedef struct ggStruct { 
  int gene, geneId ;
  unsigned int np ;   /* number of platform testing this gene */
  unsigned int nt     ;  /* number of A titrating platform, their list */
  unsigned int nta, A ;  /* number of A titrating platform, their list */
  unsigned int ntb, B ;  /* number of B titrating platform, their list */
  BOOL hasNm ;
} GG ;

/*************************************************************************************/

static BOOL maqcCountGenes_c2 (MAQC *maqc)
{
  AC_HANDLE h  = ac_new_handle () ;
  GG *gg ;
  int _CG[400000], zt[20] ;
  float z, proba[1000] ;
  FILE *f = 0 ;
  int level, x, y, ii, jj, lab, gene, geneId, t, nt ;
  const char *cp ; 
  DICT *dictLab ;
  Array aa = arrayCreate (40000, GG) ;
  char *title[] = {"NM gene", "AceView gene with geneId"
		   ,  "All genes with geneId", "Novel AceView gene", "Any gene"} ;  

  char *methodTitle[] = {"Titrating", "Co-titrating", "Strictly cotitrating"} ;  
  int method = 0 ;

  maqc->dictGene = dictCreate (40000) ;
  maqc->dictGeneId = dictCreate (40000) ;
  dictAdd (maqc->dictGene, "toto", 0) ;
  dictAdd (maqc->dictGeneId, "toto", 0) ;
  dictLab = dictCreate (15) ;
 

  /* read the probe 2 gene mapping */
  f = filopen ("P2G/P2G_any","txt","r") ;
  if (!f) exit(1) ;
  level = freesetfile (f, 0) ;
  freecard (level) ; /* jump caption line */
  while (freecard (level))
    {
      cp = freeword () ;
      dictAdd (maqc->dictGene, cp, &gene) ;
      gg = arrayp (aa, gene, GG) ;
      gg->gene = gene ;
      gg->np++ ;
    }
  freeclose (level) ; /* will close f */
  if (0) printf ("Found %d genes arrayMax(aa)=%d\n", dictMax (maqc->dictGene),arrayMax(aa)) ;
  /* read the gene2geneid relation */
  f = filopen ("P2G/gene2geneid","txt","r") ;
  if (!f) exit(1) ;
  level = freesetfile (f, 0) ;
  freecard (level) ; /* jump caption line */
  while (freecard (level))
    {
      cp = freeword () ;
      dictAdd (maqc->dictGene, cp, &gene) ;
      gg = arrayp (aa, gene, GG) ;
      gg->gene = gene ;
      cp = freeword () ;
      dictAdd (maqc->dictGeneId, cp, &geneId) ;
      gg->geneId = geneId ;
    }
  freeclose (level) ; /* will close f */
  if (0) printf ("Found %d genes arrayMax(aa)=%d\n", dictMax (maqc->dictGene),arrayMax(aa)) ;
  
  /* read the geneHasNm relation */
  f = filopen ("P2G/gene2nm2tg","txt","r") ;
  if (!f) exit(1) ;
  level = freesetfile (f, 0) ;
  freecard (level) ; /* jump caption line */
  while (freecard (level))
    {
      cp = freeword () ;
      dictAdd (maqc->dictGene, cp, &gene) ;
      gg = arrayp (aa, gene, GG) ;
      gg->gene = gene ;
      gg->hasNm = TRUE ;
    }
  freeclose (level) ; /* will close f */
  if (0) printf ("Found %d genes arrayMax(aa)=%d\n", dictMax (maqc->dictGene),arrayMax(aa)) ;
  
  /* read the probe 2 titration file */
  f = filopen ("TitratingAcvGenes_all","txt","r") ;
  if (!f) exit(1) ;
  level = freesetfile (f, 0) ;
  freecard (level) ; /* jump caption line */
  while (freecard (level))
    {
      cp = freepos () ;
      if (! strstr (cp, "all60"))
	continue ;
      cp = freeword () ;
      dictAdd (maqc->dictGene, cp, &gene) ;  
      cp = freeword () ; /* platform */
      dictAdd (dictLab, cp, &lab) ;
      gg = arrayp (aa, gene, GG) ;
      gg->gene = gene ;
      cp = freeword () ; /* A or B */
      if (! (gg->A & (0x1 << lab)) && ! (gg->B & (0x1 << lab))) gg->nt++ ;
      if (*cp == 'A') { gg->nta++ ; gg->A |= (0x1 << lab) ; }
      if (*cp == 'B') { gg->ntb++ ; gg->B |= (0x1 << lab) ; }
    }
  freeclose (level) ; /* will close f */
  if (0) printf ("Found %d genes arrayMax(aa)=%d\n", dictMax (maqc->dictGene),arrayMax(aa)) ;
  
  memset (proba, 0, sizeof (proba)) ;
  for (method = 0 ; method < 3 ; method++)
    {
      memset (_CG, 0, sizeof (_CG)) ;
      for (gene = 0, gg = arrp (aa, gene, GG) ; gene < arrayMax (aa) ; gene++, gg++) 
	{
	  
	  switch (method)
	    {
	    case 0: /* count genes with titrating probes in any direction */
	      break ;
	    case 1: /* just count consistent titration */
	      gg->nt = gg->nta > gg->ntb ? gg->nta : gg->ntb  ;
	      break ;
	    case 2:  /*  exclude genes with anti probe */
	      gg->nt = gg->nta > gg->ntb ? gg->nta : gg->ntb  ;
	      if (gg->nta && gg->ntb) gg->nt = 0 ;
	      break ;
	    }
	  if (gg->gene && gg->nt)
	    {
	      t = gg->hasNm ? 0 : (gg->geneId ? 1 :3) ;
	      CG(gg->np,0,gg->nt, t)++ ;
	      CG(gg->np,0,gg->nt, 4)++ ; 
	      if (t < 2) 
		CG(gg->np,0,gg->nt, 2)++ ;
	      if (gg->nt > gg->np || gg->nt > 14 || gg->np > 14)
		printf ("ERROR  %s np=%d nt=%d\n", dictName(maqc->dictGene,gg->gene), gg->np, gg->nt) ;
	    }
	  if (gg->gene && gg->np)
	    {
	      t = gg->hasNm ? 0 : (gg->geneId ? 1 : 3) ;
	      CG(gg->np, 0, 0, t)++ ;
	      CG(gg->np, 0, 0, 4)++ ;
	      if (t < 2) CG(gg->np, 0, 0, 2)++ ;
	    }
	}
      
      for (t = 0 ; t < 5 ; t++)
	{
	  for (jj = 1 ; jj < 11 ; jj++)
	    zt[jj] = 0 ;
	  printf ("%s %s\nTested by n platforms", methodTitle[method], title[t]) ; 
	  for (ii = 1 ; ii < 8 ; ii++)
	    printf ("\tTitrating in %d", ii) ;
	  printf ("\t more than %d\tNon titrating\tTotal\t %% %s titrating", ii - 1, methodTitle[method]) ;
	  
	  for (ii = 1 ; ii < 11 ; ii++)
	    {
	      nt = CG(ii, 0, 0, t) ;
	      printf ("\n%d", ii) ;
	      for (y=0, jj = 1 ; jj < 8 ; jj++)
		{
		  x = CG(ii,0,jj,t) ;
		  printf ("\t%d", x) ;
		  nt -= x ;
		  y += x ;
		  zt[jj] += x ;
		}
	      for (x = 0; jj < 11 ; jj++) /* titrating in more */
		x += CG(ii,0,jj,t) ; 
	      zt[8] += x ; zt[9] += nt - x ; zt[10] += CG(ii, 0, 0, t) ;
	      y += x ;
	      z = (100.0 * y)/(CG(ii, 0, 0, t) ? CG(ii, 0, 0, t) : 1) ;
	      printf ("\t%d\t%d\t%d", x, nt - x,CG(ii, 0, 0, t)) ;
	      printf ("\t%.2f", z) ;
	      proba [20*method+ii] = z ; 
	    }  
	  printf ("\nTotal") ;
	  for (jj = 1 ; jj < 11 ; jj++)
	    printf ("\t%d", zt[jj]) ;
	  printf ("\n\n") ;
	}
    }
  
  printf ("Tested in") ;
  for (method = 0 ; method < 3 ; method++)
    printf ("\t%s", methodTitle[method])  ;
  for (ii = 1 ; ii < 11 ; ii++)
    {
      printf ("\n%d", ii) ;
      for (method = 0 ; method < 3 ; method++)
	{
	  printf ("\t%.2f", proba [20*method+ii] ) ;
	}
    }
  
  
  if ((f = filopen ( "TitratingAcvGenes_all", "txt", "r")))
    {
      int gene, i, iLab, isA = 0 ;
      int nA, nB, nAB, nvA, nvB, nvAB ; 
      KEYSET ksLabA[10], ksLabB[10], ksA, ksB ;
      char *myLab[] = {"ABI","AFX", "AGL", "GEH", "ILM", "NCI", "EPP", "GEX", "QGN", "TAQ", 0 } ;
      level = freesetfile (f, 0) ;
      
      for (iLab = 0 ; iLab < 10 ; iLab++)
	{ ksLabA[iLab] = keySetHandleCreate (h) ; ksLabB[iLab] = keySetHandleCreate (h) ; }
      ksA = keySetHandleCreate (h) ; ksB = keySetHandleCreate (h) ; 
       
      while (freecard (level))
	{
	  /* parse the gene name */
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dictGene, cp, &gene) ;
	  cp = freeword () ; /* platform name */ 
	  for (iLab = 0 ; iLab < 10 ; iLab++)
	    if (!strncmp (cp, myLab[iLab], 3))
	      break ;
	  if (iLab >= 10)
	    continue ;
	  
	  cp = freeword () ; /* A or B */
	  if (cp && *cp == 'A') isA = 1 ;
	  else if (cp && *cp == 'B')  isA = 2 ;
	  else continue ;
	  cp = freepos () ;
	  if (!strstr(cp, "all60"))
	    continue ;
	  if (isA == 1)
	    {
	      if (keySetInsert (ksLabA[iLab], gene))
		keySet (ksA, gene) ++ ;
	    }
	  if (isA == 2)
	    {
	      if (keySetInsert (ksLabB[iLab], gene))
		keySet (ksB, gene) ++ ;
	    }
	}
      freeclose (level) ; /* will close f */
      
      /* count gene A and B */
      nA = keySetMax(ksA) ;  
      nB = keySetMax(ksB) ;
      if (nA > nB) keySet (ksB, nA-1) = 0 ;
      else keySet (ksA, nB-1) = 0 ;
      /* count validated gene A and B */
      nA = nB = nAB = nvA = nvB = nvAB = 0 ;
      for (i = 0 ; i < keySetMax (ksA) ; i++)
	{
	  if (keySet (ksA, i) >= LIMIT_PLATFORM)
	    {
	      if (keySet (ksB, i) >= LIMIT_PLATFORM)
		nvAB++ ;
	      else
		nvA++ ;
	    }
	  else if (keySet (ksB, i) >= LIMIT_PLATFORM)
	    nvB++ ;
	}
      
      /* count gene A and B */
      for (i = 0 ; i < keySetMax (ksA) ; i++)
	{
	  if (keySet (ksA, i) > 0)
	    {
	      if (keySet (ksB, i) > 0)
		nAB++ ;
	      else
		nA++ ;
	    }
	  else if (keySet (ksB, i) > 0)
	    nB++ ;
	}
      
      printf ("\n\nPlatform\tGenes\tTitrating\tCo-titrating\t%%\n") ;
      printf ("All\tA\t%d\t%d\t%.2f\n", nA, nvA, (100.0 * nvA) / (nA ? nA : 1)) ;
      printf ("All\tB\t%d\t%d\t%.2f\n", nB, nvB, (100.0 * nvB) / (nB ? nB : 1)) ;
      printf ("All\tAB\t%d\t%d\t%.2f\n", nAB, nvAB, (100.0 * nvAB) / (nAB ? nAB : 1)) ;
      printf ("All\tAny\t%d\t%d\t%.2f\n", nA+nB+nAB, nvA+nvB+nvAB, (100.0 * (nvA+nvB+nvAB)) / (nA+nB+nAB ? nA+nB+nAB : 1)) ;
      for (iLab = 0 ; iLab < 10 ; iLab++)
	{
	  /* count gene A and B */
	  nA = nB = nAB = nvA = nvB = nvAB = 0 ;
	  for (i = 0 ; i < keySetMax (ksLabA[iLab]) ; i++)
	    {
	      nA++ ;
	      gene =  keySet (ksLabA[iLab], i) ;
	      if (keySetFind (ksLabB[iLab], gene, 0))
		nAB++ ;
	      if (keySet (ksA, gene) >= LIMIT_PLATFORM)
		nvA++ ;
	    }
	  for (i = 0 ; i < keySetMax (ksLabB[iLab]) ; i++)
	    {
	      nB++ ;
	      gene =  keySet (ksLabB[iLab], i) ;
	      if (keySet (ksB, gene) >= LIMIT_PLATFORM)
		nvB++ ;
	    }
	  printf ("%s\tA\t%d\t%d\t%.2f\n", myLab[iLab], nA, nvA, (100.0 * nvA) / (nA ? nA : 1)) ;
	  printf ("%s\tB\t%d\t%d\t%.2f\n", myLab[iLab], nB, nvB, (100.0 * nvB) / (nB ? nB : 1)) ;
	  printf ("%s\tAB\t%d\t%d\t%.2f\n", myLab[iLab], nAB, nvAB, (100.0 * nvAB) / (nAB ? nAB : 1)) ;
	  printf ("%s\tany\t%d\t%d\t%.2f\n", myLab[iLab], nA+nB+nAB, nvA+nvB+nvAB, (100.0 * (nvA+nvB+nvAB)) / (nA+nB+nAB ? nA+nB+nAB : 1)) ;	}
    }
  
  ac_free (h) ;
  return TRUE ;
} /* maqcCountGenes_c2 */

/*************************************************************************************/
/*************************************************************************************/
static void maqcSignalPerProbeClassRegister (int *CCP0, int *CCP, int pcl, float s)
{
  int x ;
  x = s/.2 - 1 ;
  if (x > 100) x = 100 ;
  if (x < 0) x = 0 ;
  CCP[101*pcl + x]++ ;
  CCP0[101*pcl + x]++ ;
}

static BOOL maqcSignalPerProbeClass_c4 (MAQC *maqc)
{
  FILE *f ;
  float s, S ;
  int level, ii, x, probe, pcl, ns, nLab, iType ;
  LAB *lab, *oldLab ;
  MM *mm = 0 ;
  const char *cp ;
  int cumul[21], maxClasse = 0 ;
  int *CCP, *CCP0, CC[40200] ; /* 20 labos, 20 classes(class 10 assumed empty), 101 points */
  char *title[] = { "Gene specific", "Titrating", "Titrating validated", "Ambiguous", "Genome only", "Antisense", "Unmapped", "", "", "", "In titrator", "Co-titrating", "Non titrating", "Discordant"} ;
  char *type[] = { "titrating", "titratable", "cotitrating", "antititrating", 0} ;
  maqc->dict = dictCreate (1000000) ;
  maqc->aa = arrayCreate (1000000, MM) ;

  memset (CC, 0, sizeof(CC)) ;
  CCP0 = CC ; CCP = CC + 2000 ;
  for (oldLab = 0, lab = labos ; lab->id ; lab++, CCP += 2000)
    {
      if (oldLab && !strcmp (lab->lab, oldLab->lab))
	continue ;
      else if (lab->id == QGN || lab->id == EPP || lab->id == GEX  || lab->id == TAQ  || lab->id == NCI)
	continue ;

      oldLab = lab ;
      f = filopen (messprintf ("P2G/pcl.%s", lab->lab), "txt", "r") ;
      if (!f) continue ;
      level = freesetfile (f, 0) ;
      while (freecard (level))
	{
	  cp = freeword () ;
	  if (!cp) continue ;
	  dictAdd (maqc->dict, cp, &probe) ;
	  if (! freeint (&pcl)) continue ;
	  mm = arrayp (maqc->aa, probe, MM) ;
	  pcl = (pcl%100)/10 ;
	  if (pcl == 1) pcl = 2 ; /* coalesce class 1 and 2 */
	  mm->probe_class = pcl ;
	  if (pcl > maxClasse) maxClasse = pcl ;
	}
      freeclose (level) ; /* will close f */ 



      for (iType = 0 ; type[iType] ; iType++) 
	if ((f = filopen (messprintf ("P2T/%s.%s", type[iType], lab->lab), "txt", "r")))
	  {
	    int probe ;
	    level = freesetfile (f, 0) ;
	    
	    while (freecard (level))
	      {
		/* parse the probe name */
		cp = freeword () ;
		if (!cp) continue ;
		if (!dictFind (maqc->dict, messprintf ("%s", cp), &probe))
		continue ; 
		mm = arrayp (maqc->aa, probe, MM) ;
		switch (iType)
		  {
		  case 0: mm->titrating = TRUE ; break ;
		  case 1: mm->titratable = TRUE ; break ;
		  case 2: mm->cotitrating = TRUE ; break ;
		  case 3: mm->antititrating = TRUE ; break ;
		  }
	      }
	    freeclose (level) ; /* will close f */
	  }  

      f = filopen (messprintf ("aceSignalOut.%s", lab->lab), "ace", "r") ;
      if (!f) continue ;
      level = freesetfile (f, 0) ; S = 0 ; ns = 0 ; pcl = 0 ;
      /* sa = sc = sd = sb = 0 ; */
      while (freecard (level))
	{
	  cp = freeword () ;
	  if (cp && !strcmp (cp, "Probe"))
	    {
	      if (ns && pcl) 
		{
		  /* switch pcl=2 and the titrating/tit-validated columns */
		  maqcSignalPerProbeClassRegister (CCP0, CCP, pcl==2 ? 0 : pcl, S/ns) ;
		  /*
		    if (ns == 4 && (pcl == 1 || pcl == 2) &&
		    (
		    (sa < sc && sc < sd && sd < sb ) ||
		    (sa > sc && sc > sd && sd > sb )
		    )
		    )
		    {
		    maqcSignalPerProbeClassRegister (CCP0, CCP, 1, S/ns) ;
		    if (valid)
		    maqcSignalPerProbeClassRegister (CCP0, CCP, 2, S/ns) ;		    
		    }
		  */
		  if (mm->titrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 1, S/ns) ;
		  if (mm->cotitrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 2, S/ns) ;
		  if (mm->titratable)
		    { 
		      maqcSignalPerProbeClassRegister (CCP0, CCP, 10, S/ns) ;
		      if (mm->cotitrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 11, S/ns) ;
		      if (!mm->titrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 12, S/ns) ; 
		      if (mm->antititrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 13, S/ns) ;
		    }
		}
	      probe = 0 ; S = 0 ; ns = 0 ; pcl = 6 ;
	      if ((cp = freeword ()) &&
		  dictFind (maqc->dict, cp, &probe) &&
		  probe < arrayMax(maqc->aa))
		{
		  mm = arrayp (maqc->aa, probe, MM) ;
		  pcl = mm->probe_class ;
		}
	    }
	  else if (cp && probe && !strcmp (cp, "Signal_A") && freefloat (&s))
	    { S += s ; ns++ ;/* sa = s ; */  }
	  else if (cp && probe && !strcmp (cp, "Signal_C") && freefloat (&s))
	    { S += s ; ns++ ; /* sc = s ; */ }
	  else if (cp && probe && !strcmp (cp, "Signal_D") && freefloat (&s))
	    { S += s ; ns++ ; /* sd = s ; */ }
	  else if (cp && probe && !strcmp (cp, "Signal_B") && freefloat (&s))
	  { S += s ; ns++ ; /* sb = s ; */ }
	} 
      if (ns) /* last gene must also be registered */
	{
	  maqcSignalPerProbeClassRegister (CCP0, CCP, pcl==2 ? 0 : pcl, S/ns) ;
	  if (mm->titrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 1, S/ns) ;
	  if (mm->cotitrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 2, S/ns) ;
	  if (mm->titratable)
	    { 
	      maqcSignalPerProbeClassRegister (CCP0, CCP, 10, S/ns) ;
	      if (mm->cotitrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 11, S/ns) ;
	      if (!mm->titrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 12, S/ns) ; 
	      if (mm->antititrating)  maqcSignalPerProbeClassRegister (CCP0, CCP, 13, S/ns) ;
	    }
	}
      
      freeclose (level) ; /* will close f */ 
    }
  maxClasse = 13 ;
  for (oldLab = 0, nLab = 1, lab = labos ; lab->id ; nLab++, lab++)
    {
      int i0, i1, i2, i3, i4, x1, x2 ;

      if (oldLab && !strcmp (lab->lab, oldLab->lab))
	continue ;
      else if (lab->id == QGN || lab->id == EPP || lab->id == GEX  || lab->id == TAQ  || lab->id == NCI)
	continue ;
      if (nLab) oldLab = lab ;
      memset (cumul, 0, sizeof(cumul)) ;
      
      if (!nLab) lab-- ;
      printf ("\n%s: Histogram of the Log2(signal) of each probe, averaged over all RNA samples and all replicas \n"
	      "sorted by probe class\n"
	      , nLab ? lab->lab : "All platforms together") ;
      /* find the first and last useful signal */
      for (x = pcl = i0 = 0 ; x <= 100 ; x++)
	i0 += CC[2000*nLab+101*pcl+x] ;
      /* cut left 1 % */
      for (x1 = pcl = i1 = 0 ; 100 * i1 < i0 &&  x1 <= 100 ; x1++)
	i1 += CC[2000*nLab+101*pcl+x1] ;
      if (x1 > 0) x1-- ;      if (x1 > 0) x1-- ;
      for (x2 = 100,  pcl = i2 = 0 ; 100 * i2 < i0 &&  x2 >= 0 ; x2--)
	i2 += CC[2000*nLab+101*pcl+x2] ;
      if (x2 < 100) x2++ ;
      
      for (pcl = 0 ; pcl <= maxClasse ; pcl++)
	if (pcl != 2) printf ("\t%s", title[pcl]) ;
      for (x = x1 ; x <= x2 ; x++)
	{
	  printf ("\n%.1f",x/5.0) ;
	  for (pcl = 0 ; pcl <= maxClasse ; pcl++)
	    {
	      if (pcl == 9)
		{ printf ("\t%.1f", x/5.0) ; continue ; }
	      if (pcl == 2) continue ;
	      i0 = CC[2000*nLab+101*pcl+x-2] ;
	      i1 = CC[2000*nLab+101*pcl+x-1] ;
	      i2 = CC[2000*nLab+101*pcl+x] ;
	      i3 = CC[2000*nLab+101*pcl+x+1] ;
	      i4 = CC[2000*nLab+101*pcl+x+2] ;
	      ii = (i0 + 4*i1 + 6*i2 + 4*i3 + i4) / 16 ;
	      ii = (i1 + 2*i2 + i3) / 4 ;
	      printf ("\t%d", ii) ;
	      cumul[pcl] += ii ;
	    }
	}
      printf ("\nCumul") ;
      for (pcl = 0 ; pcl <= maxClasse ; pcl++)
	printf ("\t%d", cumul[pcl]) ;
      printf ("\n\n") ;
    }
  
  return TRUE ;
} /* maqcSignalPerProbeClass_c4  */

/*************************************************************************************/
/*************************************************************************************/
/* GenesTouchedPerLab , data for MS-3, table 2 */
static BOOL  maqcGenesTouchedPerLab_c5 (MAQC *maqc)
{
  AC_HANDLE h ;
  char *myLabs[] = { "ABI", "AGL", "GEH", "ILM", "NCI", "AFX", 0 } ;
  /*
    char *myLabs2[] = { "EPP", "GEX", "TAQ", "QGN",  0 } ;
    char *myLabs3[] = { "EPP", "GEX", "TAQ", "QGN", "ABI", "AGL", "GEH", "ILM", "NCI", "AFX", 0 } ;
  */
  char **lab ;
  AC_KEYSET ksp, ksm1, ksm11, ksm2, ksm3 ;
  AC_KEYSET ksv1, ksv2, ksv3, ksv4 ;
  int nProbes, nm1, nm11, nm3, nv1, nv3, nv4 ;
  int nnp = 0, nnm1 = 0, nnm11 = 0, nnv1 = 0 ;
  AC_KEYSET ksum3 = 0, ksuv3 = 0, ksuv4 = 0 ;
  AC_KEYSET ksim3 = 0, ksiv3 = 0, ksiv4 = 0 ;

  if (! maqc->db)
    return FALSE ;

  printf ("\tPlatform\tProbes with known sequence\t"
	  "\tGene specific probe with NM match"
	  "\t%%Gene specific probe with NM match"
	  "\tProbe with NM exact match"
	  "\t%%Probe with NM exact match"
	  "\tGene specific probe with AceView match"
	  "\t%%Gene specific probe with AceView match"
	  "\tGeneId seen through NM  exact hit"
	  "\tGeneId seen through AceView gene specific hit"
	  "\tNovelGene seen through AceView gene specific hit"
	  "\n"
	  ) ;
  for (lab = myLabs ; *lab ; lab++)
    { 
      h = ac_new_handle () ;
      ksp = ac_dbquery_keyset (maqc->db
			      , messprintf ("Find probe IS %s* AND Motif", *lab)
			      , h) ;
      nProbes = ac_keyset_count (ksp) ;
      nnp +=  nProbes ;
      ksm1 = ac_ksquery_keyset (ksp, "Confirmed_gene && (NM_exact_hit || NM_hit)", h) ;
      nm1 = ac_keyset_count (ksm1) ;
      nnm1 += nm1 ;
      ksm11 = ac_ksquery_keyset (ksp, "NM_exact_hit", h) ;
      nm11 = ac_keyset_count (ksm11) ;
      nnm11 += nm11 ;
      ksm2 = ac_ksquery_keyset (ksm11, ">NM_exact_hit", h) ;

      ksm3 = ac_ksquery_keyset (ksm2, ">geneid", h) ;
      nm3 = ac_keyset_count (ksm3) ;
      if (ksum3) 
	ac_keyset_or (ksum3, ksm3) ;
      else
	ksum3 = ac_copy_keyset (ksm3, 0) ;
      if (ksim3) 
	ac_keyset_and (ksim3, ksm3) ;
      else
	ksim3 = ac_copy_keyset (ksm3, 0) ;

      ksv1 = ac_ksquery_keyset (ksp, "Confirmed_gene", h) ;
      nv1 = ac_keyset_or (ksv1, ksm1) ;
      nnv1 += nv1 ;
      ksv1 = ac_ksquery_keyset (ksp, "Confirmed_gene", h) ;
      ksv2 = ac_ksquery_keyset (ksv1, ">Confirmed_gene", h) ;

      ksv3 = ac_ksquery_keyset (ksv2, ">geneid", h) ;
      nv3 = ac_keyset_or (ksv3, ksm3) ;
      if (ksuv3) 
	ac_keyset_or (ksuv3, ksv3) ;
      else
	ksuv3 = ac_copy_keyset (ksv3, 0) ;
      if (ksiv3) 
	ac_keyset_and (ksiv3, ksv3) ;
      else
	ksiv3 = ac_copy_keyset (ksv3, 0) ;
      ksv4 = ac_ksquery_keyset (ksv2, "NOT geneid", h) ;
      nv4 = ac_keyset_count (ksv4) ;
      if (ksuv4) 
	ac_keyset_or (ksuv4, ksv4) ;
      else
	ksuv4 = ac_copy_keyset (ksv4, 0) ;
      if (ksiv4) 
	ac_keyset_and (ksiv4, ksv4) ;
      else
	ksiv4 = ac_copy_keyset (ksv4, 0) ;
 
      printf ("\t%s\t%d", *lab, nProbes) ;
      if (!nProbes) nProbes = 1 ;
      printf ("\t%d\t%.2f", nm1, 100.0 * nm1/nProbes) ;
      printf ("\t%d\t%.2f", nm11, 100.0 * nm11/nProbes) ;
      printf ("\t%d\t%.2f", nv1, 100.0 * nv1/nProbes) ;
      printf ("\t%d\t%d", nm3, nv3) ;
      printf ("\t%d", nv4) ;
      printf ("\n") ;

      fflush (stdout) ;
      ac_free (h) ;
    }

  nProbes = nnp ;
  nm3 = ac_keyset_count (ksum3) ;
  nv3 = ac_keyset_count (ksuv3) ;
  nv4 = ac_keyset_count (ksuv4) ;

  printf ("\t%s\t%d", "Union", nProbes) ;
  if (!nProbes) nProbes = 1 ;
  printf ("\t%d\t%.2f", nnm1, 100.0 * nnm1/nProbes) ;
  printf ("\t%d\t%.2f", nnm11, 100.0 * nnm11/nProbes) ;
  printf ("\t%d\t%.2f", nnv1, 100.0 * nnv1/nProbes) ;
  printf ("\t%d\t%d", nm3, nv3) ;
  printf ("\t%d", nv4) ;
  printf ("\n") ;

  nm3 = ac_keyset_count (ksim3) ;
  nv3 = ac_keyset_count (ksiv3) ;
  nv4 = ac_keyset_count (ksiv4) ;

  printf ("\t%s\tNA", "Intersection") ;
  printf ("\tNA\tNA") ;
  printf ("\tNA\tNA") ;
  printf ("\t%d\t%d", nm3, nv3) ;
  printf ("\t%d", nv4) ;
  printf ("\n") ;

   return TRUE ;
} /* maqcGenesTouchedPerLab_c5 */

/*************************************************************************************/
/*************************************************************************************/
/* ProbeMappingPerGradePerLab , data for MS-5, table 2 */
static BOOL  maqcGenesTouchedPerGrade_c7 (MAQC *maqc)
{
  AC_HANDLE h = ac_new_handle () ;
  /*
  char *myLabs3[] = { "EPP", "GEX", "TAQ", "QGN", "ABI", "AGL", "GEH", "ILM", "NCI", "AFX", 0 } ;
  char *myLabs2[] = { "ABI", "AGL", "GEH", "ILM", "NCI", "AFX", 0 } ;
  */
  char *myLabs[] = { "EPP", "GEX", "TAQ", "QGN",  0 } ;
  char **lab ;
  AC_ITER iter ;
  AC_OBJ gene = 0, probe = 0 ;
  printf ("Platform\tProbes with known sequence\t"
	  "\tGene specific NM match\tGene specific AceView match"
	  "\tGeneId seen through NM  gene specific hit"
	  "\tGeneId seen through AceView gene specific hit"
	  "\tNovelGene seen through AceView gene specific hit"
	  "\n"
	  ) ;

  const char **title, *titles [] = {  "inside a RefSeq"
				    , "Specific of Entrez gene"
				    , "Specific of a discovery gene"
				    , "Ambiguous"
				    , "Genome or predicted gene"
				    , "Antisense"
				    , "Unmapped"
				    , 0 } ;
  char **cqq, *qq[] = { "confirmed_gene && (NM_hit || NM_exact_hit)"
			, "confirmed_gene && !(NM_hit || NM_exact_hit) && Confirmed_geneid"
			, "confirmed_gene && !(NM_hit || NM_exact_hit) && !Confirmed_geneid"
			, "probe_grade3"
			, "probe_grade4"
			, "probe_grade5"
			, "probe_grade6"
			, 0 } ;
  AC_KEYSET ks ;
  int nn, iLab, ii, CC[200] ;
 
  memset (CC, 0, sizeof(CC)) ;
  if (! maqc->db)
    return FALSE ;

  iter = ac_dbquery_iter (maqc->db, "Find gene maqc_probe", h) ;
  while (ac_free (gene), (gene = ac_iter_obj (iter)))
    { 
      for (iLab = 0, lab = myLabs ; *lab ; lab++, iLab++)
	{ 
	  printf ("\n%s", *lab) ;
	  for (ii =  0, cqq = qq ; *cqq ; cqq++, ii++)
	    {
	      ks = ac_objquery_keyset (probe
				       , messprintf (*cqq, *lab)
				       , 0) ; 
	      nn = ac_keyset_count (ks) ;
	      if (nn > 0)
		CC [iLab * 20 + ii]++ ;

	      ac_free (ks) ;     
	    }
	}
    }
  printf ("\nPlatform") ;
  for (title = titles ; *title ; title++)
    printf ("\t%s", *title) ;
  for (iLab = 0, lab = myLabs ; *lab ; lab++, iLab++)
    { 
      printf ("\n%s", *lab) ;
      for (ii = 0, title = titles ; *title; ii++, title++)
	printf ("\t%d", CC [iLab * 20 + ii]) ;
    }
  printf ("\n") ;

  ac_free (h) ;
  return TRUE ;
} /* maqcGenesTouchedPerGrade_c7 */

/*************************************************************************************/
/*************************************************************************************/
/* ProbeMappingPerGradePerLab , data for MS-5, table 2 */
static BOOL maqcProbeMappingPerGradePerLab_c6 (MAQC *maqc)
{
  /*
    char *myLabs2[] = { "EPP", "GEX", "TAQ", "QGN",  0 } ;
    char *myLabs3[] = { "EPP", "GEX", "TAQ", "QGN", "ABI", "AGL", "GEH", "ILM", "NCI", "AFX", 0 } ;
  */

  char *myLabs[] = { "ABI", "AFX", "AGL", "GEH", "ILM", "NCI", 0 } ;

  char **lab ;
  AC_KEYSET ks ;
  int nn ;
  const char **title, *titles [] = {  "inside a RefSeq"
				    , "Specific of Entrez gene"
				    , "Specific of a discovery gene"
				    , "Ambiguous"
				    , "Genome or predicted gene"
				    , "Antisense"
				    , "Unmapped"
				    , 0 } ;
  const char **cqq, *qq[] = { "confirmed_gene && (NM_hit || NM_exact_hit)"
			     , "confirmed_gene && !(NM_hit || NM_exact_hit) && Confirmed_geneid"
			     , "confirmed_gene && !(NM_hit || NM_exact_hit) && !Confirmed_geneid"
			     , "probe_grade3"
			     , "probe_grade4"
			     , "probe_grade5"
			     , "probe_grade6"
			     , 0 } ;
  const char *titles2 [] = {  "Gene touched by a gene specific probe"
			      , "Gene with geneid touched by a gene specific probe"
			      , "Discovery gene touched by a gene specific probe"
			      , "Gene validated titrating A > B"
			      , "Gene validated titrating A < B"
			      , "Gene with geneid validated titrating A > B"
			      , "Gene with geneid validated titrating A < B"
			      , "Discovery gene validated titrating A > B"
			      , "Discovery gene validated titrating A < B"
			      , 0 } ;
  const char *qq2[] = { "Find gene maqc_probe"
			, "Find gene Geneid && maqc_probe"
			, "Find gene ! Geneid && maqc_probe"
			, "Find gene MAQC_acdb_plus" 
			, "Find gene MAQC_acdb_minus"
			, "Find gene Geneid && MAQC_acdb_plus" 
			, "Find gene Geneid && MAQC_acdb_minus"
			, "Find gene ! Geneid && MAQC_acdb_plus" 
			, "Find gene ! Geneid && MAQC_acdb_minus"
			, 0 } ;


  if (! maqc->db)
    return FALSE ;
  /*
   * 1a : inside NM
   * 1b : in geneid gene, no nm
   * 1c : in novel gene
   * 2 : may Xhyb
   * 3 : more than one gene
   * 5 : antisense only
   * 4 : genome
   * 6 : unmapped
   */
  printf ("Platform") ;
  for (title = titles ; *title ; title++)
    printf ("\t%s", *title) ; 
  
  for (lab = myLabs ; *lab ; lab++)
    { 
      printf ("\n%s", *lab) ;
      for (cqq = qq ; *cqq ; cqq++)
	{
	  ks = ac_dbquery_keyset (maqc->db
				  , messprintf ("Find probe IS %s_* AND %s", *lab, *cqq)
				  , 0) ; 
	  nn = ac_keyset_count (ks) ;
	  if (!strcmp (*lab, "AFX")) nn/=11 ;
	  printf ("\t%d",  nn) ;
	  ac_free (ks) ;     
	}
      fflush (stdout) ;
    }
  printf ("\n") ;

  for (cqq = qq2, title = titles2 ; *cqq && *title; cqq++, title++)
    {
      ks = ac_dbquery_keyset (maqc->db
			      , messprintf ("%s", *cqq)
			      , 0) ; 
      nn = ac_keyset_count (ks) ;
      printf ("%s\t%d\n", *title,  nn) ;
      ac_free (ks) ;     
    }
 
  

  return TRUE ;
} /* maqcProbeMappingPerGradePerLab_c6 */

/*************************************************************************************/
/*************************************************************************************/
/* maqcGlobalProbeMapping_c8 for the MS3 auxiliary material */
static BOOL maqcGlobalProbeMapping_c8 (MAQC *maqc)
{
  printf ("Use the code wacext/maggie -exportMaqcMapping\n") ;
  return TRUE ;
} /* maqcGlobalProbeMapping_c8 */

/*************************************************************************************/
/*************************************************************************************/

static void  maqcDoExportProbeClass2 (AC_OBJ probe, int pcl)
{
  int oldPcl = ac_tag_int (probe, "Probe_class", 0) ;
  int grade = (pcl%100)/10 ;

  freeOutf ("Probe %s\nProbe_class %d // %d\nProbe_grade%d\n\n"
	    , freeprotect (ac_name(probe)), pcl, oldPcl, grade
	    ) ;
} /* maqcDoExportProbeClass2 */

/*************************************************************************************/

static BOOL maqcExportProbeClass2 (MAQC *maqc)
{
  BOOL ok = FALSE ;
  float score, bestScore, bestScore2 ;
  int pcl, n, ne, napp, len, nErr, nErr2, bestnErr, ir, bestIr, bestIr2 ;
  AC_HANDLE h = ac_new_handle (), h1 = 0 ;
  AC_ITER iter ;
  AC_OBJ probe = 0, bestMrna, mrna, bestTg, tg, bestGene, gene, gene1, gene2 ;
  AC_TABLE exactGenes, appGenes, mHits = 0, sExactGene = 0 ;
  AC_KEYSET ks1, ks2 ;

  iter = ac_dbquery_iter (maqc->db
			  , messprintf ("Find probe IS %s"
					, maqc->template ? maqc->template : "*")
			  , h) ;
  while (ac_free (h1), ac_free (probe), probe = ac_iter_obj (iter))
    {
      h1 = ac_new_handle () ;
      mHits = 0 ;
      exactGenes = ac_tag_table (probe, "Exact_gene", h1) ;
      sExactGene = ac_tag_table (probe, "Single_exact_gene", h1) ;
      appGenes = ac_tag_table (probe, "Approximate_gene", h1) ;
      ne = exactGenes ? exactGenes->rows : 0 ;
      napp = appGenes ? appGenes->rows : 0 ;
      if (ne == 0 && napp == 0)
	continue ;
      if ((ne == 1 && napp == 0) ||
	  (ne ==2 && napp == 0 && sExactGene)) /* same map */
	{ maqcDoExportProbeClass2 (probe, 10) ; continue ; } 

      gene1 = gene2 = 0 ;
      if (ne == 2 && napp == 0) /* test for exact repeats */
	{
	  gene1 = ac_table_obj (exactGenes, 0, 0, h1) ;
	  gene2 = ac_table_obj (exactGenes, 1, 0, h1) ;
	}
      if (ne == 0 && napp == 2) /* test for exact repeats */
	{
	  gene1 = ac_table_obj (appGenes, 0, 0, h1) ;
	  gene2 = ac_table_obj (appGenes, 1, 0, h1) ;
	}
      if (gene2)
	{
	  ks1 = ac_objquery_keyset (gene1, ">transcribed_gene;>mrna;"
				    "{>Probe_hit} SETOR {>Probe_exact_hit}"
				    , h1) ;
	  
	  ks2 = ac_objquery_keyset (gene2, ">transcribed_gene;>mrna;"
				    "{>Probe_hit} SETOR {>Probe_exact_hit}"
				    , h1) ;
	  n = ac_keyset_xor (ks1, ks2) ;
	  if (n == 0)
	    { maqcDoExportProbeClass2 (probe, ne ? 110 : 112) ; continue ; } 
	}
      if (ne > 1)
	{ maqcDoExportProbeClass2 (probe, 30) ; continue ; } 
      if (napp == 0 && sExactGene && ac_has_tag (probe, "NM_exact_hit"))
	{ maqcDoExportProbeClass2 (probe, 11) ; continue ; } 

      len = ac_tag_int (probe, "Length", 1) ;
      /* in NM case, we may have the tag Single_approximate_gene and !Approximate_gene */
      if (ac_has_tag (probe, "Single_approximate_gene") &&
	  ! ac_has_tag (probe, "mRNA_hit") &&
	  ac_has_tag (probe, "NM_Hit"))
	{
	  mHits = ac_tag_table (probe, "NM_hit", h1) ;
	  nErr = ac_table_int (mHits, 0, 3, 9999) ;
	  nErr *= 100 ; nErr /= len ;
	  pcl = nErr < 5 ? 11 : 12 ;
	  maqcDoExportProbeClass2 (probe, pcl) ; 
	  continue ;
	}
      bestnErr = 999 ; bestIr = -1 ;
      mHits = ac_tag_table (probe, "mRNA_hit", h1) ;
      for (ir = 0 ; mHits && bestnErr && ir < mHits->rows ; ir++)
	{
	  n = ac_table_int (mHits, ir, 3, 9999) ;
	  if (n < bestnErr) { bestnErr = n ; bestIr = ir ; }
	} 
      bestnErr *= 100 ; bestnErr /= len ;
      if (bestnErr == 0) ;
      else if (bestnErr <= 5) bestnErr = 1 ;
      else bestnErr = 2;
      
       if (ne == 0 && napp == 1)
	{
	  pcl = bestnErr < 2 ? 11 : 12 ;
	  maqcDoExportProbeClass2 (probe, pcl) ; 
	  continue ;
	}
     
      if (bestnErr == 2)  /* best approximate is not accurate */
	{ maqcDoExportProbeClass2 (probe, 33) ;   continue ; }
      bestMrna = mHits ? ac_table_obj (mHits, bestIr, 0, h1) : 0 ;
      bestTg = bestMrna ? ac_tag_obj (bestMrna, "From_gene", h1) : 0 ;
      bestGene = bestTg ? ac_tag_obj (bestTg, "Gene", h1) : 0 ;
      bestScore = bestGene ? ac_tag_float (bestGene, "MAQC_avg_signal", -10) : -10 ;
      bestIr2 = -1 ; bestScore2 = -10 ;
      mrna = tg = gene = 0 ;
      for (ir = 0 ; mHits && ir < mHits->rows ; ir++, ac_free (mrna), ac_free (tg), ac_free (gene))
	{
	  if (ir == bestIr) continue ;
	  mrna = ac_table_obj (mHits, ir, 0, h1) ;
	  tg = ac_tag_obj (mrna, "From_gene", h1) ;
	  if (!tg || !ac_obj_cmp (tg, bestTg)) continue ;
	  gene = ac_tag_obj (tg, "Gene", h1) ;
	  if (!gene || !ac_obj_cmp (gene, bestGene)) continue ;
	  score = ac_tag_float (gene, "MAQC_avg_signal", -10) ;
	  if (score > bestScore2)
	    { bestIr2 = ir ; bestScore2 = score ; continue ; }
	}
      if (bestScore2 == -10)
	bestScore2 = bestScore ;
      if (mHits && bestIr2 >= 0)
	{ 
	  nErr2 = ac_table_int (mHits, bestIr2, 3, 9999) ;
	  nErr2 *= 100 ; nErr2 /= len ;
	  if (nErr2 == 0) ;
	  else if (nErr2 <= 5) nErr2 = 1 ;
	  else nErr2 = 2;
	}
      else 
	nErr2 = 2 ;

      pcl = 33 ;
      if (bestnErr == 0 && nErr2 == 0) /* first is perfect */
	pcl = 30 ;
      else if (bestnErr == 1 && nErr2 == 1) /* first is perfect */
	pcl = 31 ;
      else if (bestnErr == 0 && nErr2 == 1) /* first is perfect */
	{
	  if (bestScore > bestScore2 + 2) pcl = 21 ;
	  else if (bestScore < bestScore2 - 2) pcl = 34   ;
	  else pcl = 31 ;
	}
      else if (bestnErr == 0 && nErr2 == 2) /* first is perfect */
	{
	  if (bestScore > bestScore2 + 3) pcl = napp > 2 ? 22 : 12 ;
	  else if (bestScore > bestScore2 + 0) pcl = 22 ;
	  else pcl = 32 ;
	}
      else if (nErr2 == 2) /* first is ok second is bad */
	{
	  if (bestScore > bestScore2 + 2) pcl = 22 ;
	  else pcl = 32 ;
	} 
      maqcDoExportProbeClass2 (probe, pcl) ; 
    }

  ac_free (h1) ;
  ac_free (h) ;
  return ok ;
} /* maqcExportProbeClass2 */

/*************************************************************************************/
/*************************************************************************************/

static BOOL maqcExportProbeClass (MAQC *maqc)
{
  FILE *f, *fOut = 0 ;
  BOOL ok = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ;
  float score, exactScore, approxScore ;
  int level = 0, outLevel = 0, n, pcl, oldPcl ;
  int gene, probe, nExactGene, nApproxGene ;
  Array gene2signal = arrayHandleCreate (60000, float, h) ;
  KEYSET probe2pcl = keySetHandleCreate (h) ;
  KEYSET nProbe2exactGene = keySetHandleCreate (h) ;
  KEYSET nProbe2approxGene = keySetHandleCreate (h) ;
  Array probe2exactScore = arrayHandleCreate (600000, float, h) ;
  Array probe2approxScore = arrayHandleCreate (600000, float, h) ;
  DICT *dictProbe = dictHandleCreate (1000000, h) ;
  DICT *dictGene = dictHandleCreate (100000, h) ;

  dictAdd (dictGene, "toto", 0) ;
  dictAdd (dictProbe, "toto", 0) ;

  maqc->aa = arrayCreate (1000000, MM) ;

  f = filopen ("gene2signal","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      if (freefloat (&score))
	array (gene2signal, gene, float) = score ;
    }
  freeclose (level) ;

  f = filopen ("gene2signal","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      if (freefloat (&score))
	array (gene2signal, gene, float) = score ;
    }
  freeclose (level) ;

  f = filopen ("probe2exactGene","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictProbe, ccp, &probe) ;
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      score = array (gene2signal, gene, float) ;
      n = keySet (nProbe2exactGene, probe) ;
      keySet (nProbe2exactGene, probe) = n + 1 ;
      if (n)
	{
	  exactScore = array (probe2exactScore, probe, float) ;
	  if (score < exactScore)
	    score = exactScore ;
	}
      array (probe2exactScore, probe, float) = score ;	  
    }
  freeclose (level) ;

  f = filopen ("probe2approxGene","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictProbe, ccp, &probe) ;
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      score = array (gene2signal, gene, float) ;
      n = keySet (nProbe2approxGene, probe) ;
      keySet (nProbe2approxGene, probe) = n + 1 ;
      if (n)
	{
	  approxScore = array (probe2approxScore, probe, float) ;
	  if (score < approxScore)
	    score = approxScore ;
	}
      array (probe2approxScore, probe, float) = score ;	  
    }
  freeclose (level) ;

#ifdef JUNK
  f = filopen ("probe2len2hit","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      int len = 0, nerr ;
      float xnoise ;

      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictProbe, ccp, &probe) ;
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      if (!freeint (&len) || len < 1)
	continue ;
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictGene, ccp, &gene) ;
      ccp = freeword() ; /* the mrna */
      if (!ccp) continue ;
      if (!freeint (&nerr))
	continue ;
      nerr = 100 * nerr ; nerr /= len ;
      score = array (gene2signal, gene, float) ;
      n = keySet (nProbe2approxGene, probe) ;
      keySet (nProbe2approxGene, probe) = n + 1 ;
      if (n)
	{
	  approxScore = array (probe2approxScore, probe, float) ;
	  if (score < approxScore)
	    score = approxScore ;
	}
      array (probe2approxScore, probe, float) = score ;	  
    }
  freeclose (level) ;
#endif

  f = filopen ("probe2pclass","txt","r") ;
  if (!f)
    return FALSE ;
  level = freesetfile (f, 0) ;
  while (freecard (level))
    {
      ccp = freeword() ;
      if (!ccp) continue ;
      dictAdd (dictProbe, ccp, &probe) ;
      if (freeint (&pcl))
	keySet (probe2pcl, probe) = pcl ;
    }
  freeclose (level) ;

  fOut = filopen ("probe2pcl_new", "ace", "w") ;
  if (fOut)
    {
      outLevel = freeOutSetFile (fOut) ;
      for (probe = 1 ; probe <= dictMax (dictProbe) ; probe++) 
	{
	  nExactGene = keySet (nProbe2exactGene, probe) ;
	  nApproxGene = keySet (nProbe2approxGene, probe) ;
	  exactScore = array (probe2exactScore, probe, float) ;
	  approxScore = array (probe2approxScore, probe, float) ;
	  oldPcl = keySet (probe2pcl, probe) ;
	  pcl = 0 ;

	  if (nExactGene == 0)
	    {
	      if (nApproxGene == 1)
		pcl = 14 ;
	      else
		pcl = 35 ;
	    }
	  else if (nExactGene == 1)
	    {
	      if (nApproxGene == 0)
		pcl = 10 ;
	      else if (nApproxGene < 6 && approxScore < exactScore - 3)
		pcl = 11 ;
	      else if (nApproxGene < 3 && approxScore < exactScore - 1)
		pcl = 11 ;
	      else if (nApproxGene < 3 && approxScore > exactScore + 1)
		pcl = 32 ;
	      else if (nApproxGene < 3)
		pcl = 20 ;
	      else
		pcl = 31 ;
	    }
	  else /* nExactGene  > 1 */
	    {
	      pcl = 3 ;
	    }
	  if (pcl) freeOutf ("Probe %s\nProbe_class %d // old= %d\n\n"
			     , dictName (dictProbe, probe)
			     , pcl, oldPcl) ;
	}

      freeOutClose (outLevel) ;
      filclose (fOut) ;
      ok = TRUE ;
    }

  ac_free (h) ;
  
  return ok ;
} /* maqcExportProbeClass */

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: maqc [-title \"my title\"] [-db ACEDB] [-o outfile] [-t template] [-k keyset] [-c1] [-d1] [-m0 <float>] -nSigma <float>] [-allAgree]\n") ;

  fprintf (stderr, "// Example:  maqc -title \"this is a test\"\n") ; 
  fprintf (stderr, "//  -title \"text inside double quotes\" : the date [and title] are quoted at the top of output\n") ;
  fprintf (stderr, "// -o output_file : redirects the output") ;

  fprintf (stderr, 
	   "// Normalisation: on demand the program rescales the arrays (shifts the log-signal), \n"
	   "this is done globally for all the probes of the array. This renormalisation makes the  \n"
	   "arrays from diferent sites more comparable and has a strong effect on the global output \n"
	   "of the program, the various options are as follows\n") ;  
  fprintf (stderr, "// -mo target-scale : the default is 10\n") ;
  if (0) fprintf (stderr, "// -so normalizes the variance, it was worse when we tried, 0 means do not normalize\n") ;
  fprintf (stderr, "// [-mediane | -moyenne ] adjust the median or the average of the array to the target-scale\n") ;
  fprintf (stderr, "// [-basic | -nobasic] adjust the whole array but compute the shift so that the median/average of the probes specific of the genes listed in P2G/basicGenes.txt is shifted to the target-value\n") ;
  fprintf (stderr, "// -nocentre do not adjust the arrays.\n") ;
  fprintf (stderr, "// -flag     use flags (ignore even if default is not to)\n") ;
  fprintf (stderr, "// -noflag   do not  use flags (even if default is not to)\n") ;
  fprintf (stderr, "// -damir (use damir's files, -flag (use flags)\n") ;

  fprintf (stderr, "// -c1, -d1 , various analysis\n") ;
  fprintf (stderr, "// -c2, table of genes titrating in n platforms/tested in platforms\n") ;
  fprintf (stderr, "// -c3, export probe_class\n") ;
  fprintf (stderr, "// -c4, histo of signal per probe class\n") ;
  fprintf (stderr, "// -c5, GenesTouchedPerLab , data for MS-3, table 2\n") ;
  fprintf (stderr, "// -c6, ProbeMappingPerGradePerLab , data for MS-5, table 2\n") ;
  fprintf (stderr, "// -c7, GenesTouchedPerGrade , data for MS-5, table 3\n") ;
  fprintf (stderr, "// -c8, maqcGlobalProbeMapping_c8 for the MS3 auxiliary material\n") ;
  fprintf (stderr, "// -c9, per lab consistency of titrating probes\n") ;
  fprintf (stderr, "// -c10, not yet used\n") ;
  fprintf (stderr, "// -allAgree used in the sourire analysis to restrict to genes where all probes agree\n") ;
  fprintf (stderr, "//   -?? : do ??\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  float x ;
  const char *s = "ok" ;
  const char *ccp = 0, *outFileName = 0 ;
  const char *dbName = 0 ;
  int n, outlevel = 0 ;
  MAQC maqc ;
  AC_HANDLE h = 0 ;

  freeinit () ; 

  if (argc == 1)
    usage () ;

  getCmdLineOption (&argc, argv, "-title", &ccp) ;
  printf ("// start: %s\t%s\n", timeShowNow(), ccp && *ccp ? ccp : "") ;
  h = ac_new_handle () ;
  memset (&maqc, 0, sizeof (MAQC)) ;
  
  maqc.useNm = FALSE ;
  maqc.h = h ;
  maqc.s0 = maqc.m0 = 0 ;
  maqc.nSigma = 2 ;

  /* consume optional args */
  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "-o", &outFileName) ;
  getCmdLineOption (&argc, argv, "-t", &(maqc.template)) ;
  getCmdLineOption (&argc, argv, "-k", &(maqc.keysetName)) ;
  getCmdLineOption (&argc, argv, "-lab", &(maqc.labName)) ;
  maqc.c1 = getCmdLineOption (&argc, argv, "-c1", 0) ;
  maqc.d1 = getCmdLineOption (&argc, argv, "-d1", 0) ;
  maqc.c2 = getCmdLineOption (&argc, argv, "-c2", 0) ;
  maqc.c3 = getCmdLineOption (&argc, argv, "-c3", 0) ;
  maqc.c4 = getCmdLineOption (&argc, argv, "-c4", 0) ;
  maqc.c5 = getCmdLineOption (&argc, argv, "-c5", 0) ;
  maqc.c6 = getCmdLineOption (&argc, argv, "-c6", 0) ;
  maqc.c7 = getCmdLineOption (&argc, argv, "-c7", 0) ;
  maqc.c8 = getCmdLineOption (&argc, argv, "-c8", 0) ;
  maqc.c9 = getCmdLineOption (&argc, argv, "-c9", 0) ;
  maqc.c10 = getCmdLineOption (&argc, argv, "-c10", 0) ;

  if (getCmdLineOption (&argc, argv, "-damir", 0))
    maqc.damir = 1 ;
  if (getCmdLineOption (&argc, argv, "-damir2", 0))
    maqc.damir = 2 ;
  maqc.flag = FLAGBOF ;
  if (getCmdLineOption (&argc, argv, "-flag", 0))
    maqc.flag = FLAG ;
  if (getCmdLineOption (&argc, argv, "-noflag", 0))
    maqc.flag = NOFLAG ;

  maqc.centre = CENTREBOF ;
  if (getCmdLineOption (&argc, argv, "-moyenne", 0))
    maqc.centre = MOYENNE ;
  if (getCmdLineOption (&argc, argv, "-mediane", 0))
    maqc.centre = MEDIANE ;
  if (getCmdLineOption (&argc, argv, "-basic", 0))
    maqc.basic = BASIC ;
  if (getCmdLineOption (&argc, argv, "-nobasic", 0))
    maqc.basic = NOBASIC ;
  if (getCmdLineOption (&argc, argv, "-nocentre", 0))
    maqc.centre = NOCENTRE ;
  if (getCmdLineOption (&argc, argv, "-random", 0))
    maqc.random = TRUE ;

  maqc.allAgree = getCmdLineOption (&argc, argv, "-allAgree", 0) ;
  x = 0 ;
  x = 0 ;
  if (getCmdLineOption (&argc, argv, "-m0", &ccp) &&
      sscanf (ccp, "%f", &x) == 1)
    maqc.m0 = x ;
  x = 0 ;
  /* this option is disabled, adjusting the sigma of the platform was disastrous */
  if (0 &&
      getCmdLineOption (&argc, argv, "-s0", &ccp) &&
      sscanf (ccp, "%f", &x) == 1)
    maqc.s0 = x ;
  if (getCmdLineOption (&argc, argv, "-nSigma", &ccp) &&
      sscanf (ccp, "%f", &x) == 1)
    maqc.nSigma = x ;

  /* check the absolute args */
  if (argc != 1)
    usage () ;

  if (dbName)
    {
      maqc.db = ac_open_db (dbName, &s);
      if (!maqc.db)
	messcrash ("Failed to open db %s, error %s", dbName, s) ;
    }

  if (outFileName)
    {
      f = filopen (outFileName, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  if (maqc.c2)
    {
      maqcCountGenes_c2 (&maqc) ; 
      goto done ;
    }
  if (maqc.c3)
    {
      if (0) maqcExportProbeClass (&maqc) ; 
      if (1) maqcExportProbeClass2 (&maqc) ; 
      goto done ;
    }
  if (maqc.c4)
    {
      maqcSignalPerProbeClass_c4 (&maqc) ; 
      goto done ;
    }
  if (maqc.c5)
    {
      maqcGenesTouchedPerLab_c5 (&maqc) ; 
      goto done ;
    }
  if (maqc.c6)
    {
      maqcProbeMappingPerGradePerLab_c6 (&maqc) ; 
      goto done ;
    }
  if (maqc.c7)
    {
      maqcGenesTouchedPerGrade_c7 (&maqc) ; 
      goto done ;
    }
  if (maqc.c8)
    {
      maqcGlobalProbeMapping_c8 (&maqc) ; 
      goto done ;
    }
  if (maqc.c9)
    {
      maqcCountConfirmedGenesTitratingProbes_c9 (&maqc) ; 
      goto done ;
    }
  if (maqc.c10)
    {
      
      goto done ;
    }

  maqcGetData (&maqc) ;
  n = dictMax (maqc.dictGeneId) + 1 ;
  maqc.testedGeneId = keySetHandleCreate (h) ;
  maqc.signalGeneId = keySetHandleCreate (h) ;
  keySet (maqc.testedGeneId, n) = 0 ;
  keySet (maqc.signalGeneId, n) = 0 ;
  n = dictMax (maqc.dictGene) + 1 ;
  maqc.testedGene = keySetHandleCreate (h) ;
  maqc.signalGene = keySetHandleCreate (h) ;
  keySet (maqc.testedGene, n) = 0 ;
  keySet (maqc.signalGene, n) = 0 ;

  printf ("// data parsed, start statistics: %s\n", timeShowNow()) ;
  if (maqc.c1)
    { 
      LAB *lab ;
      KEYSET ks2 = keySetHandleCreate (h) ;
      KEYSET ks4 = keySetHandleCreate (h) ;
      KEYSET ksg2 = keySetHandleCreate (h) ;
      KEYSET ksg4 = keySetHandleCreate (h) ;
      KEYSET ksAllab = keySetHandleCreate (h) ;

      for (lab = &labos [0] ; *lab->lab ; lab++)
	{
	  if (maqc.labName && strcmp(maqc.labName,lab->lab) && strcmp (maqc.labName, lab->lab) &&
	      strcmp (maqc.labName, messprintf("%s_%s", lab->lab, lab->bench)))		
	    continue ;

	  printf ("\n##########################################################") ;
	  printf ("\n#### %s\n", lab->lab) ;
	  maqcC1 (&maqc, lab) ;  /* normalisation */
	  maqcC2 (&maqc, lab, ks2, ks4, ksg2, ksg4, ksAllab) ;
	  if (1)
	    {
	      maqcC3 (&maqc, lab) ; /* calcul des TM */
	      maqcC4 (&maqc, lab) ; 
	      maqcCorrelation (&maqc, lab) ;
	    }
	}
    }

  if (maqc.d1)
    {
      KEYSET ks2 = keySetHandleCreate (h) ;
      KEYSET ks4 = keySetHandleCreate (h) ;
      KEYSET ksg2 = keySetHandleCreate (h) ;
      KEYSET ksg4 = keySetHandleCreate (h) ;

      maqcD1 (&maqc) ; 
      printf ("\n##########################################################\n## NM analysis s=%.2f"
	      , maqc.s0) ;
      { int i ;
	for (i = 0 ; labos[i].id  ; i++)
	  if (big6 & labos[i].id) 
	    printf (" %s%s", labos[i].lab, labos[i].bench) ; 
      }
      printf ("\n") ;
      maqcD2 ("NM", &maqc, ks2, ks4, FALSE) ; 
      printf ("\n##########################################################\n## AceView analysis\n") ;
      maqcD2 ("AceView", &maqc, ksg2, ksg4, TRUE) ;
      maqcD3 ("NM A<>B", &maqc, ks2, ks4, FALSE, FALSE) ; 
      maqcD3 ("AceView A<>B", &maqc, ksg2, ksg4, FALSE, TRUE) ;
      maqcD3 (messprintf ("NM A<%.2f>B", maqc.nSigma), &maqc, ks2, ks4, TRUE, FALSE) ; 
      maqcD3 (messprintf ("AceView A<%.2f>B", maqc.nSigma), &maqc, ksg2, ksg4, TRUE, TRUE) ;
    }
  else if (!maqc.c1)
    usage () ;

 done:
  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  ac_db_close (maqc.db) ;
  if (maqc.aceSignalOut)
    filclose(maqc.aceSignalOut) ;
  printf ("// done: %s\n", timeShowNow()) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* melting temperature of a single stranded oligo acoording to Maniatis
 * LaDeana Hillier in OSp analysis.c uses the same formula
 *
 * input is a zero terminated string of atgcATGC 
 * non case sensitive, other char are just ignored
 */
/* the  oligoTm function is now in acedna.c which should be linked in */


/*************************************************************************************/
/*************************************************************************************/
