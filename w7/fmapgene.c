/*  File: fmapgene.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Dynamic programming gene finding
 * HISTORY:
 * Last edited: Dec 17 03:22 1998 (rd)
 * * Jul 16 14:27 1998 (edgrif): Introduce private header fmap_.h
 * * Jun 12 14:15 1994 (rd): change to select on BEGIN/END to
 		allow refinement of selections in regions
 * * Jun 12 14:15 1994 (rd): added exonScore, move to pg rules
 * * Jun  7 18:27 1994 (rd): split set Parms function into two
 * Created: Tue Aug  4 18:58:53 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmapgene.c,v 1.8 2017/02/15 20:36:39 mieg Exp $ */

#include "fmap_.h"

#include "freeout.h"
#include "dna.h"
#include "peptide.h"
#include "systags.h"
#include "session.h"
#include "peptide.h"
#include "sysclass.h"
#include "pick.h"		/* for pickMatch () */
#include "main.h"		/* for checkWriteAccess () */

/* #define DEBUG */

static void setTempGene (LOOK look, Array a, int k, int frame, 
			 BOOL startNotFound, BOOL endNotFound,
			 BOOL isByHand) ;

static KEY findPeptideMatch (KEY temp,int newchk,KEY key) ;
static int calculateAndSaveChecksum (KEY protein);

/*************************/

typedef enum {
  DUMMY=0, BEGIN, END, START, STOP, SPLICE_5, SPLICE_3, FRAME
  } ELTtype ;	/* STOP < SPLICE_5 required for correct sort */

#ifdef DEBUG
static char *eltTypeName[] = {
  "DUMMY",
  "BEGIN", 
  "END",
  "START",
  "STOP",
  "SPLICE_5",
  "SPLICE_3",
  "FRAME" } ;
#endif

typedef struct ELTstruct {
  ELTtype type ;
  int	x ;	/* base before splice5, after splice3 */
		/* first base of ATG, last of stop */
  int	frame ;
  float	localScore ;
  float score ; /* DP partial score */
  struct ELTstruct *ptr ;	/* DP back pointer */
  int	flag ;
} ELT ;

#define FLAG_REQUIRED		0x0001
#define FLAG_ANTI_REQUIRED	0x0002

/*********************************************************/

static KEY _Gene_assemble_method ;
static KEY _Intron_min, _Exon_min, _Intron_cost, _Inter_gene_cost ;
static KEY _GF_range, _GF_ATG_cutoff, _GF_5_cutoff, _GF_3_cutoff ;

static KEY M_GF_coding_seg, M_GF_ATG, M_GF_splice ;
static KEY M_hexExon, M_hexIntron, M_hexExon_span ;
		/* hexExon_span is for old style cumulative scoring */

static Array codeSegs = 0 ;

void initialise (void)
{ 
  static BOOL isDone = FALSE ;

  if (isDone) return ;
				/* tags */
  lexaddkey ("Intron_min", &_Intron_min, 0) ;
  lexaddkey ("Exon_min", &_Exon_min, 0) ;
  lexaddkey ("Intron_cost", &_Intron_cost, 0) ;
  lexaddkey ("Inter_gene_cost", &_Inter_gene_cost, 0) ;
  lexaddkey ("GF_range", &_GF_range, 0) ;
  lexaddkey ("GF_ATG_cutoff", &_GF_ATG_cutoff, 0) ;
  lexaddkey ("GF_5_cutoff", &_GF_5_cutoff, 0) ;
  lexaddkey ("GF_3_cutoff", &_GF_3_cutoff, 0) ;
  lexaddkey ("Gene_assemble_method", &_Gene_assemble_method, 0) ;

				/* methods */
  lexaddkey ("GF_ATG", &M_GF_ATG, _VMethod) ;
  lexaddkey ("GF_splice", &M_GF_splice, _VMethod) ;
  lexaddkey ("GF_coding_seg", &M_GF_coding_seg, _VMethod) ;
  lexaddkey ("hexExon_span", &M_hexExon_span, _VMethod) ;
  lexaddkey ("hexExon", &M_hexExon, _VMethod) ;
  lexaddkey ("hexIntron", &M_hexIntron, _VMethod) ;

				/* default display values */
  methodSet ("GF_ATG", ORANGE,
	     METHOD_FRAME_SENSITIVE | METHOD_STRAND_SENSITIVE |
	     METHOD_SCORE_BY_WIDTH,
	     4.2, 2, 'C', 0, 3) ;
  methodSet ("GF_splice", 0,
	     METHOD_STRAND_SENSITIVE,
	     5.6, 10, 0, -1, 4) ;
  methodSet ("GF_coding_seg", LIGHTGRAY, 
	     METHOD_FRAME_SENSITIVE | METHOD_STRAND_SENSITIVE |
	     METHOD_SCORE_BY_WIDTH,
	     4.3, 2, 'C', 2, 8) ;
  methodSet ("hexExon", ORANGE,
	     METHOD_FRAME_SENSITIVE | METHOD_STRAND_SENSITIVE |
	     METHOD_SCORE_BY_WIDTH,
	     4.6, 2, 'H', 10, 50) ;
  methodSet ("hexIntron", BROWN,
	     METHOD_STRAND_SENSITIVE | METHOD_SCORE_BY_WIDTH,
	     5.1, 2, 'I', 10, 40) ;

				/* default methods for gene assembly */
  method (0, M_GF_ATG)->gaFac = 1 ;
  method (0, M_GF_splice)->gaFac = 1 ;
  method (0, M_hexExon_span)->gaFac = 0.2 ;

  isDone = TRUE ;
}

/********************************************************/
/********** intron and exon scoring routines ************/

static int gf_range = 10000 ;
static int intronMin = 40 ;
static int intronRateMin = 55 ;
static float intronRate = -1.5 ; /* per log base bp beyond rate min */
static float intronBase = -3.0 ;
static int exonMin = 33 ;
static float interGeneCost = -5.0 ;
float intron3Cutoff = 0.0;
float intron5Cutoff = -0.5;
float atgCutoff = 0.0;

static float interGeneScore (ELT *elt1, ELT *elt2)
{
  return interGeneCost ;
}   

static float intronScore (ELT *elt1, ELT *elt2)
{
  int length = elt2->x - elt1->x - 2 ; /* -2 because ->x are exon endpoints */

  if (length < intronMin && 
      elt1->type != BEGIN && elt2->type != END)
    return -1000.0 ;
  else if (length <= intronRateMin)
    return intronBase ;
  else
    return intronBase + intronRate*log10 (length - intronRateMin) ;
}

static float exonScore (LOOK look, ELT *elt1, ELT *elt2, 
			int frame, GFINFO *gf)
{
  int i, x1 = elt1->x, x2 = elt2->x, s1, s2 ;
  float score = 0, bit = 0 ;
  KEY methodk = 0 ;
  SEG *seg ;
  METHOD *meth ;
  HOMOLINFO *hinf ;

  if (elt2->type == STOP) 
    x2 -= 2 ; /* added 2 so sorted correctly with splices */

  if (x2-x1+1 < exonMin &&
      elt1->type != BEGIN && elt1->type != START &&
      elt2->type != END && elt2->type != STOP)
    return -1000.0 ;

  while (x1 % 3 != frame) 
    ++x1 ;			/* so start in frame */
  while (x2 % 3 != frame) 
    --x2 ;			/* so end in frame */

  if (arrayExists (codeSegs)) for (i = 0 ; i < arrayMax (codeSegs) ; ++i)
    { seg = arr (codeSegs,i,SEG*) ;
      if (seg->x1 >= x2 || seg->x2 <= x1)
	continue ;
      if (seg->type == HOMOL || seg->type == HOMOL_UP)
	{ hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  methodk = hinf->method ;
	  bit = hinf->score ;
	}
      else if (seg->type == FEATURE || seg->type == FEATURE_UP)
	{ methodk = seg->key ;
	  bit = seg->data.f ;
	}
      meth = method (0, methodk) ;
      if (!meth->gaFac)
	continue ;
      if ((meth->flags & METHOD_FRAME_SENSITIVE) && seg->x1 %3 != frame)
	continue ;
      s1 = seg->x1 > x1 ? seg->x1 : x1 ;
      s2 = seg->x2 < x2 ? seg->x2 : x2 ;
      score += bit * meth->gaFac *
	 (s2 - s1 + 1) / (float) (seg->x2 - seg->x1 + 1) ;
    }

	/* hack for oldfashioned approach */
  meth = method (0, M_hexExon_span) ;
  if (meth && meth->gaFac) /* mhmp 08.09.98 */
    { x1 += 3 ;
      x2 -= 3 ;
      if (x1 < x2)
	score += meth->gaFac * 
	  (gf->cum[x2 - gf->min] - gf->cum[x1 - gf->min]) ;
    }

#ifdef DEBUG
  printf (" <%.2f>", score) ;
#endif 
  return score ;
}

#ifdef OLD_EXON_SCORE
  float min, max = 0, *cum = 0 ;

  cum = &gf->cum[x1 - gf->min] ;
  min = *cum ; cum += 3 ;
  for (i = x1 + 3 ; i <= x2 ; i += 3, cum += 3)
    { if (*cum < min) 
	min = *cum ;
      else if (*cum - min > max)
	max = *cum - min ;
    }

#ifdef COMPENSATE_FOR_LENGTH
  max -= log10 (x2-x1+1) ;
#endif

  if (max > 0)
    return max * codingCost ;
  else
    return 0 ;
}

static void flattenGfCum (LOOK look)
{
  int i ;
  char *cp = arrp (look->dna, look->gf.min, char) ;
  float max[3] ;
  char * tt = pepGetTranslationTable (look->seqKey, 0) ;

  for (i = 0 ; i < 3 ; ++i)
    max[i] = -1.0e20 ;

  for (i = look->gf.min ; i < look->gf.max ; ++i, ++cp)
    if (e_codon (cp, tt) == '*')
      max[i%3] = -1.0e20 ;
    else if (look->gf.cum[i] < max[i%3])
      look->gf.cum[i] = max[i%3] ;
    else
      max[i%3] = look->gf.cum[i] ;
}
#endif

/********************** set parameters ***********************/

static Graph parmsGraph = 0 ;
static char parmsName[64] = "" ;
static void parmsEditor (void) ;
static int parmsNameBox ;

static void parmsLoadName (char *s)
{ 
  KEY key ;
  OBJ obj ;
  BSMARK mark = 0 ;
  float x ;

  initialise () ;
  if (!lexword2key (s, &key, _VMethod))
    { messout ("Can't find Method object %s", parmsName) ;
      return ;
    }
  strcpy (parmsName, name (key)) ;
  if (! (obj = bsCreate (key)))
    { messout ("Failed to open Method object %s", name (key)) ;
      return ;
    }

  bsGetData (obj, _Intron_min, _Int, &intronMin) ;
  bsGetData (obj, _Exon_min, _Int, &exonMin) ;
  x = bsGetData (obj, _Intron_cost, _Float, &intronBase) &&
    bsGetData (obj, _bsRight, _Float, &intronRate) &&
    bsGetData (obj, _bsRight, _Int, &intronRateMin) ;
  bsGetData (obj, _Inter_gene_cost, _Float, &interGeneCost) ;

  bsGetData (obj, _GF_range, _Int, &gf_range) ;
  bsGetData (obj, _GF_ATG_cutoff, _Float, &atgCutoff) ;
  bsGetData (obj, _GF_5_cutoff, _Float, &intron5Cutoff) ;
  bsGetData (obj, _GF_3_cutoff, _Float, &intron3Cutoff) ;

  {				/* check displayable */
    float rp = 2.0 ;
    if (!bsGetData (obj, str2tag ("Right_priority"), _Float, 0))
      bsAddData (obj, str2tag ("Right_priority"), _Float, &rp) ;
    if (!bsGetKeyTags (obj, _Colour, 0))
      if (bsAddTag (obj, _Colour) && bsPushObj (obj))
	{ bsAddTag (obj, _BLUE) ;
	  bsGoto (obj, 0) ;
	}
  }

				/* clear {methodInfo}->gaFac */
  { int i ;
    METHOD **methp = arrp (methodInfo, 0, METHOD*) ;
    for (i = 0 ; i < arrayMax (methodInfo) ; ++i, ++methp)
      if (*methp) 
	 (*methp)->gaFac = 0 ;
  }
  if (bsGetKey (obj, _Gene_assemble_method, &key)) do
    { mark = bsMark (obj, mark) ;
      if (bsGetData (obj, _bsRight, _Float, &x))
	{ methodAdd (0, key) ;
	  method (0, key)->gaFac = x ;
	}
      bsGoto (obj, mark) ;
    } while (bsGetKey (obj, _bsDown, &key)) ;
  bsMarkFree (mark) ;

  bsDestroy (obj) ;

  if (parmsGraph)
    parmsEditor () ;
}

static void parmsLoad (void)
{ 
  if (graphPrompt ("Give name of method object to load from", parmsName, "wz"))
    parmsLoadName (freeword ()) ;
}

static void parmsSave (void)
{
  KEY key ;
  OBJ obj ;
  int i ;

  initialise () ;
  if (!graphPrompt ("Give name of method object to save into:",
		    parmsName, "wz"))
    return ;
  strncpy (parmsName, freeword (), 63) ;
  if (!*parmsName)
    { messout ("You must give a name") ;
      return ;
    }
  if (!graphCheckEditors (parmsGraph, TRUE))
    return ;
  if (!lexaddkey (parmsName, &key, _VMethod) &&
      !graphQuery (messprintf ("Overwrite method %s?", name (key))))
    return ;
  if (! (obj = bsUpdate (key)))
    { messout ("Failed to edit object Method:%s", name (key)) ;
      return ;
    }

  bsAddData (obj, _Intron_min, _Int, &intronMin) ;
  bsAddData (obj, _Exon_min, _Int, &exonMin) ;
  bsAddData (obj, _Intron_cost, _Float, &intronBase) ;
  bsAddData (obj, _bsRight, _Float, &intronRate) ;
  bsAddData (obj, _bsRight, _Int, &intronRateMin) ;
  bsAddData (obj, _Inter_gene_cost, _Float, &interGeneCost) ;

  bsAddData (obj, _GF_range, _Int, &gf_range) ;
  bsAddData (obj, _GF_ATG_cutoff, _Float, &atgCutoff) ;
  bsAddData (obj, _GF_5_cutoff, _Float, &intron5Cutoff) ;
  bsAddData (obj, _GF_3_cutoff, _Float, &intron3Cutoff) ;

  {				/* check displayable */
    float rp = 2.0 ;
    if (!bsGetData (obj, str2tag ("Right_priority"), _Float, 0))
      bsAddData (obj, str2tag ("Right_priority"), _Float, &rp) ;
    if (!bsGetKeyTags (obj, _Colour, 0))
      if (bsAddTag (obj, _Colour) && bsPushObj (obj))
	{ bsAddTag (obj, _BLUE) ;
	  bsGoto (obj, 0) ;
	}
  }

  if (bsFindTag (obj, _Gene_assemble_method))
    bsRemove (obj) ;

  { METHOD **methp = arrp (methodInfo, 0, METHOD*) ;
    for (i = 0 ; i < arrayMax (methodInfo) ; ++i, ++methp)
      if (*methp && (*methp)->gaFac)
	if (bsAddKey (obj, _Gene_assemble_method, KEYMAKE (_VMethod, i)))
	  bsAddData (obj, _bsRight, _Float, & (*methp)->gaFac) ;
  }

  bsSave (obj) ;

  if (graphActivate (parmsGraph))
    graphBoxDraw (parmsNameBox, -1, -1) ;
}

static MENUOPT parmsMenu[] =
{ { graphDestroy, "Quit" },
  { graphPrint, "Print" },
  { parmsLoad, "Load" },
  { parmsSave, "Save" },
  { 0, 0 }
} ;

void newMethod (char *text)
{ 
  KEY key ;
  float x ;

  if (!graphPrompt (messprintf ("Give score factor for method %s", text), 
		    "0", "fz"))
    return ;
  freefloat (&x) ;
  lexaddkey (text, &key, _VMethod) ; methodAdd (0, key) ;
  method (0, key)->gaFac = x ;
  parmsEditor () ;
}

static void parmsEditor (void)
{ 
  int line = 0 ;
  int i ;
  KEY key ;
  static char newMethodName[64] ;

  initialise () ;

  if (graphActivate (parmsGraph))
    { graphPop () ;
      graphClear () ;
    }
  else
    { parmsGraph = graphCreate (TEXT_SCROLL, "Gene assembly parameters", 
				0, 0, 0.4, 0.4) ;
      graphMenu (parmsMenu) ;
    }

  if (!*parmsName)
    if (lexword2key ("assembly-default", &key, _VMethod) && iskey (key) == 2)
      parmsLoadName ("assembly-default") ;

  parmsNameBox = graphBoxStart () ;
  graphText ("Name:", 1, line) ;
  graphTextPtr (parmsName, 7, line++, 32) ;
  graphBoxEnd () ;
  ++line ;
  graphText ("Assembly parameters", 1, line++) ;
  graphIntEditor ("Min intron length", &intronMin, 2, line++, 0) ;
  graphIntEditor ("Min exon length", &exonMin, 2, line++, 0) ;
  graphFloatEditor ("Intron base cost", &intronBase, 2, line++, 0) ;
  graphFloatEditor ("Intron rate", &intronRate, 2, line++, 0) ;
  graphIntEditor ("Intron rate min", &intronRateMin, 2, line++, 0) ;
  graphFloatEditor ("Inter-gene cost", &interGeneCost, 2, line++, 0) ;
  ++line ;
  graphText ("Assembly methods", 1, line++) ;
  { METHOD **methp = arrp (methodInfo, 0, METHOD*) ;
    for (i = 0 ; i < arrayMax (methodInfo) ; ++i, ++methp)
      if (*methp && (*methp)->gaFac)
	graphFloatEditor (name (KEYMAKE (_VMethod, i)), & (*methp)->gaFac,
			  2, line++, 0) ;
  }
  *newMethodName = 0 ;
  graphTextScrollEntry (newMethodName, 63, 20, 2, line++, newMethod) ;
  ++line ;
  graphText ("Genefinder parameters", 1, line++) ;
  graphIntEditor ("Features range (bp)", &gf_range, 2, line++, 0) ;
  graphFloatEditor ("3-splice cutoff", &intron3Cutoff, 2, line++, 0) ;
  graphFloatEditor ("5-splice cutoff", &intron5Cutoff, 2, line++, 0) ;
  graphFloatEditor ("ATG cutoff", &atgCutoff, 2, line++, 0) ;
  ++line ;
  graphText ("Save as \"assembly-default\" to set the default",
	     1, line++) ;

  graphRedraw () ;
}

/**************/

int eltOrder (const void *x, const void *y)
{
  int diff = ((const ELT*)x)->x - ((const ELT*)y)->x ;

  if (diff)
    return diff ;
  else
    return ((ELT*)x)->type - ((ELT*)y)->type ;
}

static void fMapDP (void)
{
  static Array zombie[3] ;	/* active elts in intron */
  static Array alive[3] ;	/* active elts in coding frame */
  static Array elts ;		/* complete set */
  static Array gene ;
  ELT *elt, *elt2, *dead, end, *best ; /* best needed for one-gene version */
  int i, j, jKeep, frame = 0, zframe ;
  int min, max, endFrame ;
  float score, scoreKeep, fac ;
  SEG *seg ;
  HOMOLINFO *hinf ;
  METHOD *meth ;
  char *cp ;
  char *tt = 0 ;
  FMAPLOOKGET ("fMapDP") ;

  if (!*parmsName)		/* load default parameters */
    { KEY key ;
      if (lexword2key ("assembly-default", &key, _VMethod) && iskey (key) == 2)
	parmsLoadName ("assembly-default") ;
    }
  else if (parmsGraph)
    graphCheckEditors (parmsGraph, TRUE) ;

  if (!look->gf.cum)
    { fMapAddGfSegs () ;
      if (!look->gf.cum)
	messout ("I cannot locate the genefinder tables, sorry.  "
		 "Either put them in $ACEDB/wgf, or set the environment "
		 "variable GF_TABLES.") ;
      return ;
    }

  	/* find working zone */
  tt = pepGetTranslationTable (look->seqKey, 0) ;
  fMapSetZone () ;
  min = look->zoneMin ;
  max = look->zoneMax + 1 ;
  if (look->gf.min > min)
    min = look->gf.min ;
  if (look->gf.max < max)
    max = look->gf.max ;
  if (max < min)
    { messout ("Available range is (%d,%d) is negative", min, max) ;
      return ;
    }

	/* make elts Array of features for dynamic programming */

  elts = arrayReCreate (elts, 1024, ELT) ;
  codeSegs = arrayReCreate (codeSegs, 1024, SEG*) ;
				/* elts from SEGs: SPLICES + STARTS */
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->x1 < min ||
	  seg->x2 >= max ||
	  assFind (look->antiChosen, SEG_HASH (seg), 0))
	continue ;
      switch (seg->type)
	{
	case ATG:
	  meth = method (0, seg->key) ;
	  if ((fac = meth->gaFac))
	    { elt = arrayp (elts, arrayMax (elts), ELT) ;
	      elt->type = START ;
	      elt->x = seg->x1 ;
	      elt->frame = seg->x1 % 3 ;
	      elt->localScore = seg->data.f * fac ;
	    }
	  break ;
	case SPLICE5:	/* need 3 - indexed by frame */
	  meth = method (0, seg->key) ;
	  if ((fac = meth->gaFac))
	    for (frame = 0 ; frame < 3 ; ++frame)
	      { elt = arrayp (elts, arrayMax (elts), ELT) ;
		elt->type = SPLICE_5 ;
		elt->x = seg->x1 ;
		elt->localScore = seg->data.f * fac ;
		elt->frame = frame ;
	      }
	  break ;
	case SPLICE3:	/* 3 indexed by zframe */
	  meth = method (0, seg->key) ;
	  if ((fac = meth->gaFac))
	    for (zframe = 0 ; zframe < 3 ; ++zframe)
	      { elt = arrayp (elts, arrayMax (elts), ELT) ;
		elt->type = SPLICE_3 ;
		elt->x = seg->x2 ;
		elt->localScore = seg->data.f * fac ;
		elt->frame = zframe ;
	      }
	  break ;
	case FEATURE:
	  meth = method (0, seg->key) ;
	  if (meth->gaFac && seg->data.f)
	    array (codeSegs, arrayMax (codeSegs), SEG*) = seg ;
	  break ;
	case HOMOL:
	  hinf = arrp (look->homolInfo, seg->data.i, HOMOLINFO) ;
	  meth = method (0, hinf->method) ;
	  if (meth->gaFac && hinf->score)
	    array (codeSegs, arrayMax (codeSegs), SEG*) = seg ;
	  break ;
	default:
	  break ;
	}
    }
				/* all STOPs */
  cp = arrp (look->dna, min, char) ;
  for (i = min ; i <= max-3 ; ++i)
    if (e_codon (cp++, tt) == '*')
      { elt = arrayp (elts, arrayMax (elts), ELT) ;
	elt->type = STOP ;
	elt->x = i+2 ;
	elt->frame = i % 3 ;
      }

  arraySort (elts, eltOrder) ;

	/* do the dynamic programming - first initialise to BEGIN */
  { ELT begin ;

    memset (&begin, 0 , sizeof(begin)) ;
    begin.type = BEGIN ; 
    begin.x = min ; 
    
    for (frame = 0 ; frame < 3 ; ++frame)
      { alive[frame] = arrayReCreate (alive[frame], 32, ELT*) ;
	array (alive[frame], 0, ELT*) = (ELT*) messalloc (sizeof (ELT)) ;
	*arr (alive[frame], 0, ELT*) = begin ;
	arr (alive[frame], 0, ELT*)->frame = frame ;
      }
    for (zframe = 0 ; zframe < 3 ; ++zframe)
      { ELT *myelt ; /* I need to allocate *myelt, otherwise some SGI compilo crash */
        zombie[zframe] = arrayReCreate (zombie[zframe], 32, ELT*) ;
	array (zombie[zframe], 0, ELT*) = myelt = 
	  (ELT*) messalloc (sizeof (ELT)) ;
	*myelt = begin ;  /* here myelt is absolutely needed */
	myelt->frame = zframe ;
      }
    dead = (ELT*) messalloc (sizeof (ELT)) ;
    *dead = begin ;
    dead->score = -interGeneCost ;
    best = 0 ;
  }

  for (i = 0 ; i < arrayMax (elts) ; ++i)
    { elt = arrp (elts, i, ELT) ;
#ifdef DEBUG
      printf ("%d_%d %4.1f %s ", COORD (look,elt->x), elt->frame, 
	      elt->localScore, eltTypeName[elt->type]) ;
#endif
      switch (elt->type)
	{
	case START:
	  { ELT *wasDead = dead ;
	    if (dead)
	      { elt->score = dead->score + elt->localScore + 
		  interGeneScore (dead, elt) ;
		if (dead->type != BEGIN)
		  elt->ptr = dead ;
	      }
	    if (assFind (look->chosen, HASH (ATG, elt->x + 2), 0))
	      { for (j = 0 ; j < 3 ; ++j)
		  { arrayMax (alive[j]) = 0 ;
		    arrayMax (zombie[j]) = 0 ;
		  }
		dead = 0 ;
	      }
	    if (wasDead)
	      array (alive[elt->frame], 
		     arrayMax (alive[elt->frame]), ELT*) = elt ;
	  }
	  break ;
	case STOP:
	  if (assFind (look->chosen, HASH (ORF, elt->x - 3), 0))
	    dead = best = 0 ;		/* to force choice of this frame */
	  if (!assFind (look->antiChosen, HASH (ORF, elt->x - 3), 0))
	    for (j = 0 ; j < arrayMax (alive[elt->frame]) ; ++j)
	      { elt2 = arr (alive[elt->frame], j, ELT*) ;
		score = elt2->score + 
		  exonScore (look, elt2, elt, elt->frame, &look->gf) ;
#ifdef DEBUG
	printf (" (%d_%d %4.1f)", COORD (look,elt2->x), elt2->frame, score) ;
#endif
		if (!dead || score > dead->score)
		  { elt->score = score ;
		    elt->ptr = elt2 ;
		    dead = elt ;
		  }
		if (!best || score > best->score)
		  { elt->score = score ;
		    elt->ptr = elt2 ;
		    best = elt ;
		  }
	      }
	  if (assFind (look->chosen, HASH (ORF, elt->x - 3), 0))
	    for (j = 0 ; j < 3 ; ++j) /* kill all except dead */
	      { arrayMax (alive[j]) = 0 ;
		arrayMax (zombie[j]) = 0 ;
	      }
	  else
	    arrayMax (alive[elt->frame]) = 0 ;
	  break ;
	case SPLICE_5:
	  zframe = (elt->x + 4 - elt->frame) % 3 ;
	  for (j = 0 ; j < arrayMax (alive[elt->frame]) ; ++j)
	    { elt2 = arr (alive[elt->frame], j, ELT*) ;
	      score = elt2->score + elt->localScore + 
		exonScore (look, elt2, elt, elt->frame, &look->gf) ;
#ifdef DEBUG
	printf (" (%d_%d %4.1f)", COORD (look,elt2->x), elt2->frame, score) ;
#endif
	      if (!elt->ptr || score > elt->score)
		{ elt->score = score ;
		  elt->ptr = elt2 ;
		}
	    }
	  if (assFind (look->chosen, HASH (SPLICE5, elt->x+1), 0))
	    { arrayMax (zombie[zframe]) = 0 ;
	      arrayMax (alive[elt->frame]) = 0 ;
	      dead = 0 ;
	    }
	  if (elt->ptr)		/* i.e. found something */
	    array (zombie[zframe], 
		   arrayMax (zombie[zframe]), ELT*) = elt ;
	  break ;
	case SPLICE_3:
	  frame = (3 - elt->frame + elt->x) % 3 ;
	  scoreKeep = -1000000.0 ;
	  jKeep = 0 ;

	  for (j = 0 ; j < arrayMax (zombie[elt->frame]) ; ++j)
	    { elt2 = arr (zombie[elt->frame], j, ELT*) ;
	      score = elt2->score + elt->localScore + 
		intronScore (elt2, elt) ;
#ifdef DEBUG
	printf (" (%d_%d %4.1f)", COORD (look,elt2->x), elt2->frame, score) ;
#endif
	      if (!elt->ptr || score > elt->score)
		{ elt->ptr = elt2 ;
		  elt->score = score ;
		}
				/* compress zombie[elt->frame] */
	      if (score > scoreKeep ||
		  elt->x - elt2->x < intronMin) /* keep */
		{ arr (zombie[elt->frame], jKeep++, ELT*) = elt2 ;
		  scoreKeep = score ;
		}
	    }
	  if (assFind (look->chosen, HASH (SPLICE5, elt->x), 0))
	    { arrayMax (alive[frame]) = 0 ;
	      arrayMax (zombie[elt->frame]) = 0 ;
	      dead = 0 ;
	    }
	  else
	    arrayMax (zombie[elt->frame]) = jKeep ;
	  if (elt->ptr)
	    { array (alive[frame], arrayMax (alive[frame]), ELT*) = elt ;
	      if (elt->ptr->type == BEGIN)
		elt->ptr = 0 ;
	    }
	  break ;
	case FRAME:		/* kill other frames */
	  for (frame = 0 ; frame < 3 ; ++frame)
	    { arrayMax (zombie[frame]) = 0 ;
	      if (frame != elt->frame)
		arrayMax (alive[frame]) = 0 ;
	      dead = 0 ;
	    }
	  break ;
	default:
	  break ;
	}
#ifdef DEBUG
      if (elt->ptr)
	printf (" | %d_%d %4.1f\n", 
		COORD (look,elt->ptr->x), elt->ptr->frame, elt->score) ;
      else
	printf ("\n") ;
#endif      
    }

			/* find the best end point */
			/* want an end elt for intronScore,exonScore */
  end.type = END ; end.x = max - 1 ; end.ptr = 0 ;
  end.score = -1000000 ; 
  endFrame = 0 ;
  for (frame = 0 ; frame < 3 ; ++frame)
    for (j = 0 ; j < arrayMax (alive[frame]) ; ++j)
      { elt2 = arr (alive[frame], j, ELT*) ;
	score = elt2->score + 
	  exonScore (look, elt2, &end, frame, &look->gf) ;
	if (score > end.score)
	  { end.score = score ;
	    end.ptr = elt2 ;
	    endFrame = frame ;
	  }
      }
  for (zframe = 0 ; zframe < 3 ; ++zframe)
    for (j = 0 ; j < arrayMax (zombie[zframe]) ; ++j)
      { elt2 = arr (zombie[zframe], j, ELT*) ;
	score = elt2->score + intronScore (elt2, &end) ;
	if (score > end.score)
	  { end.score = score ;
	    end.ptr = elt2 ;
	  }
      }
  if (dead && dead->type != BEGIN)
    { if (dead->score > end.score)
	{ end.score = dead->score ;
	  end.ptr = dead ;
	}
    }
  else if (best && best->score > end.score)
    { end.score = best->score ;
      end.ptr = best ;
    }

	/* now backtrack, reversing the order of features using a stack
	*/
  if (end.ptr)
    { int kgene = 0 ;
      BOOL startNotFound = FALSE;
      Stack stack = stackCreate (32) ;

      messout ("min, max are %d, %d\n"
	       "total score is %f\n", 
	       COORD (look,min), COORD (look,max), end.score) ;

      for (elt = end.ptr ; elt ; elt = elt->ptr)
	push (stack, elt, ELT*) ;

      gene = arrayReCreate (gene, 16, int) ;
      while (!stackEmpty (stack))
        { elt = pop (stack, ELT*) ;

#ifdef DEBUG
	  printf ("\t%s at %6d score %6.2f localScore %6.2f\n", 
		  eltTypeName[elt->type], COORD (look,elt->x), 
		  elt->score, elt->localScore) ;
#endif

	  if (!arrayMax (gene))
	    startNotFound = (elt->type != START) ;
	  if (arrayMax (gene) == 1) /* find frame */
	    frame = elt->frame ;
	  array (gene, arrayMax (gene), int) = elt->x ;
	  if (elt->type == STOP)
	    { setTempGene (look, gene, ++kgene, frame,
			   startNotFound, FALSE, FALSE) ;
	      gene = arrayReCreate (gene, 16, int) ;
	    }
	}

      if (arrayMax (gene))
	{ if (arrayMax (gene) % 2)
	    { if (arrayMax (gene) == 1)
		frame = endFrame ;
	      array (gene, arrayMax (gene), int) = max-1 ;
	    }
	  setTempGene (look, gene, ++kgene, frame,
		       startNotFound, TRUE, FALSE) ;
	}

      stackDestroy (stack) ;
    }
  else
    messout ("Couldn't find a gene through the chosen features") ;
}

static void fMapDP1 (void) 
{ float temp = interGeneCost ;
  interGeneCost = -10000 ;
  fMapDP () ;
  interGeneCost = temp ;
}

/********************************************/

FREEOPT fMapChooseMenu[] = {
  { 4, "Splice menu" },
  { 1, "Select" },
  { 2, "Antiselect" },
  { 3, "Unselect" },
  { 4, "Delete" }
} ;

void fMapChooseMenuFunc (KEY key, int box)
{ 
  int index ;
  SEG *seg ;
  FMAPLOOKGET ("spliceMenuFunc") ;

  if (! (index = arr (look->boxIndex, box, int)))
    return ;
  seg = arrp (look->segs, index, SEG) ;

  switch (key)
    {
    case 1:
      assRemove (look->antiChosen, SEG_HASH (seg)) ;
      assInsert (look->chosen, SEG_HASH (seg), assVoid (seg->x1)) ;
      graphBoxDraw (box, -1, GREEN) ;
      break ;
    case 2:
      assRemove (look->chosen, SEG_HASH (seg)) ;
      assInsert (look->antiChosen, SEG_HASH (seg), assVoid (seg->x1)) ;
      graphBoxDraw (box, -1, LIGHTGREEN) ;
      break ;
    case 3:
      assRemove (look->chosen, SEG_HASH (seg)) ;
      assRemove (look->antiChosen, SEG_HASH (seg)) ;
      graphBoxDraw (box, -1, WHITE) ;
      break ;
    case 4:
      assRemove (look->chosen, SEG_HASH (seg));
      assRemove (look->antiChosen, SEG_HASH (seg));
      arrayMax (look->segs) = look->lastTrueSeg;
      *seg = arr (look->segs, arrayMax (look->segs)-1, SEG);
      --arrayMax (look->segs);
      arraySort (look->segs, fMapOrder) ;
      look->lastTrueSeg = arrayMax (look->segs);
      fMapDraw (look, 0);
      break;
    }
}

/******************************/

static int compareChoseOrder (const void *x, const void *y)
{ 
  return UNHASH_X2 (* (const char**)x) - UNHASH_X2 (* (const char**)y) ;
}

static void fMapChooseToGene (void)
{
  int i ;
  int min, max ;
  int kgene = 0 ;
  char* v ;
  static Array chosen, gene ;
  FMAPLOOKGET ("fMapChooseToGene") ;

  min = look->zoneMin ;
  max = look->zoneMax + 1 ;

  chosen = arrayReCreate (chosen, 8, char*) ;
  v = 0 ;
  while (assNext (look->chosen, &v, 0))
    if (UNHASH_X2 (v) >= min && UNHASH_X2 (v) < max)
      array (chosen, arrayMax (chosen), char*) = v ;
  if (!arrayMax (chosen))
    { messout ("No chosen features in the active zone") ;
      return ;
    }
  arraySort (chosen, compareChoseOrder) ;

#define IGNORE  { messout ("Impossible chosen feature %s near %d - ignoring", \
			   fMapSegTypeName[UNHASH_TYPE (v)], \
			   COORD (look, UNHASH_X2 (v))) ; break ; }

  gene = arrayReCreate (gene, 8, int) ;
  if (UNHASH_TYPE (arr (chosen,0,char*)) != ATG)
    array (gene, 0, int) = min ;
  for (i = 0 ; i < arrayMax (chosen) ; ++i)
    { v = arr (chosen, i, char*) ;
      switch (UNHASH_TYPE (v))
	{
	case SPLICE3:
	  if (arrayMax (gene) % 2) IGNORE ;
	  array (gene, arrayMax (gene), int) = UNHASH_X2 (v) ; 
	  break ;
	case SPLICE5:
	  if (! (arrayMax (gene) % 2)) IGNORE ;
	  array (gene, arrayMax (gene), int) = UNHASH_X2 (v) - 1  ; 
	  break ;
	case ATG:
	  if (arrayMax (gene)) IGNORE ;
	  array (gene, arrayMax (gene), int) = UNHASH_X2 (v) - 2  ; 
	  break ;
	case ORF:
	  if (! (arrayMax (gene) % 2)) IGNORE ;
	  array (gene, arrayMax (gene), int) = UNHASH_X2 (v) + 3  ;
	  setTempGene (look, gene, ++kgene, 0, FALSE, FALSE, TRUE) ;
	  gene = arrayReCreate (gene, 8, int) ;
	  break ;
	default:
	  messout ("Unrecognized chosen feature type %s at %d",
		   fMapSegTypeName[UNHASH_TYPE (v)],
		   COORD (look, UNHASH_X2 (v))) ; break ;
	  return ;
	}
    }

  if (arrayMax (gene))
    { if (arrayMax (gene) % 2)
	array (gene, arrayMax (gene), int) = max ;
      setTempGene (look, gene, ++kgene, 0, FALSE, TRUE, TRUE) ;
    }
}

static void fMapChooseFromGene (LOOK look, KEY gene)
{ int i ;
  SEG *seg ;

  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->parent == gene)
	switch (seg->type)
	  {
	  case INTRON:
	    assInsert (look->chosen, HASH (SPLICE5, seg->x1),
		       assVoid (seg->x1 - 1)) ;
	    assInsert (look->chosen, HASH (SPLICE3, seg->x2 + 1),
		       assVoid (seg->x2)) ;
	    break ;
	  case CDS:
	    assInsert (look->chosen, HASH (ATG, seg->x1 + 2),
		       assVoid (seg->x1)) ;
	    assInsert (look->chosen, HASH (ORF, seg->x2 - 3),
		       assVoid (1)) ;
	    break ;
	  default:
	    break ;
	  }
    }

  fMapDraw (look, 0) ;
}

static void fMapSetChoose (void)
{
  FMAPLOOKGET ("fMapSetChoose") ;

  if (look->activeBox)
    fMapChooseFromGene (look, BOXSEG (look->activeBox)->parent) ;
}

/********************************/

static BOOL killGene (KEY key)
{
  KEY source ;
  OBJ obj, sourceObj ;

  if ((obj = bsCreate (key)))
    { if (bsGetKey (obj, _Source, &source))
	{ if (! (sourceObj = bsUpdate (source)))
	    { messout ("Couldn't update %s", name (source)) ;
	      bsDestroy (obj) ;
	      return FALSE ;
	    }
	  bsAddKey (sourceObj, _Subsequence, key) ;
	  bsRemove (sourceObj) ;	/* delete reference */
	  bsSave (sourceObj) ;
	}
      bsDestroy (obj) ;
    }

  if ((obj = bsUpdate (key)))
    { if (bsFindTag (obj, _From)) bsRemove (obj) ;
      if (bsFindTag (obj, _Start_not_found)) bsRemove (obj) ;
      if (bsFindTag (obj, _End_not_found)) bsRemove (obj) ;
      if (bsFindTag (obj, _CDS)) bsRemove (obj) ;
      bsSave (obj) ;
    }
  else
    { messout ("Couldn't update %s", name (key)) ;
      return FALSE ;
    }

  return TRUE ;
}

/***********************************/

static void killGeneMenu (void)
{
  KEY parent ;
  FMAPLOOKGET ("killGeneMenu") ;

  if (!look->activeBox || 
      ! (parent = BOXSEG (look->activeBox)->parent) ||
      class (parent) != _VSequence)
    { messout ("You must select a sequence") ;
      return ;
    }

  if (messQuery ("Do you want to destroy the sequence location "
	    "information attached to %s ?", name (parent)))
    { killGene (parent) ;
      fMapConvert (look, TRUE) ;
      fMapDraw (look, 0) ;
    }
}

/**************************/

static void setTempGene (LOOK look, Array a, int k, int frame,
			 BOOL startNotFound, BOOL endNotFound,
			 BOOL isByHand)
{
  KEY key, parentKey ;
  OBJ obj, parentObj = 0 ;
  int i, base = 0 ;
  int x, y ;

  if (!arrayMax (a))
    return ;

  if (arrayMax (a) % 2)
    { messerror ("setTempGene called with an odd "
		 "number of defining points") ;
      return ;
    }

  lexaddkey (messprintf ("temp_gene_%d", k), &key, _VSequence) ;

  if (!killGene (key))	/* can't clear previous entry */
    return ;

  if (! (obj = bsUpdate (key)))
    { messout ("Couldn't update temp gene object") ;
      return ;
    }
  for (i = base = 0 ; i < arrayMax (a) ; i += 2)
    { x = arr (a, i, int) ; y = arr (a, i+1, int) ;
      if (!i)
	base = x - 1 ;
      x -= base ; y -= base ;
      bsAddData (obj, _Source_Exons, _Int, &x) ;
      bsAddData (obj, _bsRight, _Int, &y) ;
    }
  bsAddTag (obj, _CDS) ;
  if (startNotFound)
    { bsAddTag (obj, _Start_not_found) ;
      frame = (3 + frame - ((base+1)%3)) % 3 + 1 ; /* CDS entry */
      if (frame != 1)
	bsAddData (obj, _CDS, _Int, &frame) ;
    }
  if (endNotFound)
    bsAddTag (obj, _End_not_found) ;

  { KEY mkey ;
    OBJ mobj ;
    float rp = 2.0 ;

    if (isByHand)
      lexaddkey ("hand_built", &mkey, _VMethod) ;
    else if (lexIsGoodName (parmsName))
      lexaddkey (parmsName, &mkey, _VMethod) ;
    else
      lexaddkey ("assembly-default", &mkey, _VMethod) ;
    mobj = bsUpdate (mkey) ;
    if (mobj)
      { if (!bsGetData (mobj, str2tag ("Right_priority"), _Float, 0))
	  bsAddData (mobj, str2tag ("Right_priority"), _Float, &rp) ;
	if (!bsGetKeyTags (mobj, _Colour, 0))
	  if (bsAddTag (mobj, _Colour) && bsPushObj (mobj))
	    bsAddTag (mobj, isByHand ? _DARKBLUE : _BLUE) ;
	bsSave (mobj) ;
      }
    bsAddKey (obj, str2tag ("Method"), mkey) ;
  }

  bsSave (obj) ;

  x = base + 2 ;
  y = base + y + 1 ;

  parentKey = key ;		/* pass this so don't get it back! */
  if (fMapFindSpan (look, &parentKey, &x, &y) &&
      (parentObj = bsUpdate (parentKey)))
    { bsAddKey (parentObj, _Subsequence, key) ;
      bsAddData (parentObj, _bsRight, _Int, &x) ;
      bsAddData (parentObj, _bsRight, _Int, &y) ;
      bsSave (parentObj) ;
      fMapConvert (look, TRUE) ;
      fMapChooseFromGene (look, key) ;
    }
  else
    messout ("Failed to update parent sequence object") ;
}

/****************************************/


void recalcChecksums (){ 
  KEY peptide = 0;

  messout ("Recalulating NEW checksums"); 
  while (lexNext (_VProtein,&peptide))
      calculateAndSaveChecksum (peptide);
}      

/****************************************/

static void fixTemp (void)
{
  char *newName;
  KEY old, source, temp ;
  OBJ oldObj = 0, tempObj = 0, sourceObj = 0, obj = 0 ;
  static KEY lastKey = 0 ;
  int x1, x2, i, cdsOffset ;
  BOOL isNoStart, isNoEnd ;
  static Array a = 0 ;
  BOOL exists = FALSE,thisgo=TRUE;
  Array newpep;
  int newchk;
  KEY samepep,key,wormkey;
  static BOOL pepManage, first = TRUE;
  static char classtype[40],format[15];
  static int start,finish;
  char *cp = "checksum type 2" ; /* can use 3 4 5 later*/
  KEY dummy= 0;
  FMAPLOOKGET ("fixTemp") ;

#define error(x)	{ messout (x) ; goto abort ;}

  if (!checkWriteAccess ())
    return;

  if (!look->activeBox || 
      ! (temp = BOXSEG (look->activeBox)->parent) ||
      class (temp) != _VSequence)
    { messout ("You must select a sequence to fix") ;
      return ;
    }
  if (! (tempObj = bsCreate (temp)))
    error ("Could not open sequence object to fix") ;
  if (!bsGetKey (tempObj, _Source, &source))
    error ("sequence to fix is not attached to anything") ;
  bsDestroy (tempObj) ;
  if (messPrompt (messprintf ("Fix %s to:", name (temp)),
		  lastKey ? name (lastKey) : "","wz")){
    newName = freeword () ;
    if (lexword2key (newName, &old, _VSequence))
      { if ((oldObj = bsCreate (old)))
	  { if (bsFindTag (oldObj, _From) &&
		 (!messQuery (messprintf ("Overwrite %s ?\n", 
					 name (old))) ||
		 !killGene (old)))
	      goto abort ;
	    bsDestroy (oldObj) ;
	  }
	
	if (! (sourceObj = bsUpdate (source)))
	  error ("Could not update source gene") ;
	if (bsAddKey (sourceObj, _Subsequence, temp))
	  { bsRemove (sourceObj) ;
	    error ("No cross-reference of temp in parent") ;
	  }
	if (!bsGetData (sourceObj, _bsRight, _Int, &x1) ||
	    !bsGetData (sourceObj, _bsRight, _Int, &x2))
	  error ("No offsets for temp in its parent") ;
	bsAddKey (sourceObj, _Subsequence, temp) ;
	bsRemove (sourceObj) ;
	bsAddKey (sourceObj, _Subsequence, old) ;
	bsAddData (sourceObj, _bsRight, _Int, &x1) ;
	bsAddData (sourceObj, _bsRight, _Int, &x2) ;
	bsSave (sourceObj) ;
	
	if (! (oldObj = bsUpdate (old)))
	  error ("Could not update old gene") ;
	
	if (! (tempObj = bsUpdate (temp)))
	  error ("Could not update temp gene") ;
	if (!bsFindTag (tempObj, _Source_Exons))
	  error ("No source exons in temp gene") ;
	a = arrayReCreate (a, 32, BSunit) ;
	bsFlatten (tempObj, 2, a) ;
	bsFindTag (tempObj, _From) ;
	bsRemove (tempObj) ;
	isNoStart = bsFindTag (tempObj, _Start_not_found) ;
	isNoEnd = bsFindTag (tempObj, _End_not_found) ;
	cdsOffset = 0 ; bsGetData (tempObj, _CDS, _Int, &cdsOffset) ;
	bsSave (tempObj) ;
	
	for (i = 0 ; i < arrayMax (a) ; i += 2)
	  { bsAddData (oldObj, _Source_Exons, _Int, 
		       & (arr (a,i,BSunit).i)) ;
	    bsAddData (oldObj, _bsRight, _Int, 
		       & (arr (a,i+1,BSunit).i)) ;
	  }
	if (isNoStart) 
	  bsAddTag (oldObj, _Start_not_found) ;
	else if (bsFindTag (oldObj, _Start_not_found))
	  bsRemove (oldObj) ;
	if (isNoEnd) 
	  bsAddTag (oldObj, _End_not_found) ;
	else if (bsFindTag (oldObj, _End_not_found))
	  bsRemove (oldObj) ;
	if (cdsOffset) 
	  bsAddData (oldObj, _CDS, _Int, &cdsOffset) ;
	else if (bsGetData (oldObj, _CDS, _Int, &cdsOffset))
	  bsRemove (oldObj) ;
	bsAddTag (oldObj, _CDS) ;
	bsSave (oldObj) ;
	lastKey = old ;
	fMapConvert (look, TRUE) ;
      }
    else
      { lexAlias (&temp, newName, TRUE, FALSE) ;
	lastKey = temp ;
      }
  }

  if (first)			/* establish whether to pepManage */
    { 
      char junk[255];
      char *cp = sessionFilName ("wspec/database","wrm","r") ;
      FILE *fpace = cp ? filopen (cp, 0, "r") : 0 ;

      pepManage = FALSE;
      first = FALSE;
      if (fpace)
	{ while (fgets (junk,255,fpace) != NULL)
	    if (strstr (junk,"PEPMANAGE"))
	      {
		if (sscanf (junk,"PEPMANAGE %s %s %d %d",
			   classtype, format, &start, &finish)
		    ==4)
		  pepManage = TRUE;
		else
		  messout ("Error: could only get PEPMANAGE %s,%s,%d,%d.\n",
			  classtype, format, start, finish);
	      }
	  fclose (fpace);
	}
    }
  
  if (pepManage)
    { if ((obj = bsUpdate (source)))
	{ if (bsFindTag (obj,str2tag ("Wormpep")))
	    { thisgo = TRUE;
	      messout ("Using new pepcode routines\n");
	    }
	  else
	    thisgo = FALSE;
	  bsDestroy (obj);
	}
      else
	thisgo = FALSE;

      for (key = 0 ; lexNext (_VClass, &key) ; )
	if (pickMatch (name (key), classtype)) 
	  break ;

      if (!key)
	{ messout ("Could not find the class %s. Therefore exiting FixTemp gene\n"
		  "Check that %s is in the subclasses.wrm", classtype) ;
	  return ;
	}

      if (thisgo)
	{ if (!lexword2key (cp, &dummy, 0))
	    { lexaddkey (cp, &dummy, 0) ; /* only once */
	      recalcChecksums ();
	    }

	  newpep = peptideGet (temp);
	  newchk = hashArray (newpep);
	  if ((samepep = findPeptideMatch (temp,newchk,key)))
	    exists = TRUE;
	  i=start;
	  while (lexword2key (messprintf (format,i),&key,_VProtein)) i++ ;
	  if (i < finish)
	    lexaddkey (messprintf (format,i),&key,_VProtein);
	  else
	    { messout ("ERROR: Limit of numbers for format has been reached see PEPMANAGE in database.wrm\n"
		      "exiting fix temp gene") ; 
	      return ;
	    }
	  if (exists)
	    { messout ("existing protein %s matches temp gene %s so setting Corresponding_protein accordingly",
		      name (samepep), name (temp)) ;
	      lastKey = temp ;      
	      obj = bsUpdate (temp);
	      if (!obj || !bsAddKey (obj,str2tag ("Corresponding_Protein"), samepep))
		error ("Error adding Corresponding Protein");
	      bsSave (obj);
	    }
	  else
	    { messout ("Does NOT match any other protein.  "
		      "Creating new protein object %s", messprintf (format,i)) ;
	      if (!lexword2key (messprintf (format,i), &wormkey, _VProtein))
		{ error ("lexwordtokey could not find protein") ;}
	      else if (peptideStore (wormkey,newpep))
		{ if (! (obj = bsUpdate (temp)))
		    error ("Could not update old gene") ;
		  if (!bsAddKey (obj,str2tag ("Corresponding_Protein"),wormkey))
		    error ("Error adding Corresponding Protein");
		  bsSave (obj);
		  if ((obj = bsUpdate (wormkey)))
		    { if (!bsAddTag (obj,str2tag (classtype)))
			error ("Error adding Wormpep tag");
		      bsSave (obj);
		    }
		}
	      else  
		error ("Unable to add Peptide");
	    }
	}
    }

  fMapDraw (look, lastKey) ;
  return ;
  
 abort:
  bsDestroy (oldObj) ;
  bsDestroy (tempObj) ;
  bsDestroy (sourceObj) ;
}

/****************************************/
/****************************************/

static void fMapHideGfColumns (void)
{ 
  FMAPLOOKGET ("fMapHideGfColumns") ;

  mapColSetByName ("ATG", FALSE) ;
  mapColSetByName ("ORF's", FALSE) ; 
  mapColSetByName ("Coding Frame", FALSE) ; 
  mapColSetByName ("GF_coding_seg", FALSE) ;
  mapColSetByName ("GF_splice", FALSE) ;
  mapColSetByName ("hexExon", FALSE) ;

  look->isGeneFinder = FALSE;
  fMapDraw (look, 0) ;
}

static void fMapShowGfColumns (void)
{ 
  FMAPLOOKGET ("fMapShowGfColumns") ;

  mapColSetByName ("ATG", TRUE) ;
  mapColSetByName ("ORF's", TRUE) ; 
  mapColSetByName ("Coding Frame", TRUE) ; 
  mapColSetByName ("GF_coding_seg", TRUE) ;
  mapColSetByName ("GF_splice", TRUE) ;
  mapColSetByName ("hexExon", TRUE) ;

  look->isGeneFinder = TRUE;
  fMapDraw (look, 0) ;
}

/****************************************/

static LOOK gfLook ;

void fMapAddGfSeg (int type, KEY key, int x1, int x2, float score)
{
  SEG *seg = arrayp (gfLook->segs, arrayMax (gfLook->segs), SEG) ;

  fMapInitialise () ;

  seg->parent = 0 ;
  seg->type = type ;
  seg->key = key ;
  seg->x1 = x1 + gfLook->gf.min ;
  seg->x2 = x2 + gfLook->gf.min ;
  seg->data.f = score ;
}

void fMapAddGfSite (int type, int pos, float score, BOOL comp)
{

  fMapInitialise () ;
  
  switch (type)
    {
    case '3':
      fMapAddGfSeg (comp ? SPLICE3_UP : SPLICE3, M_GF_splice,
		    pos-1, pos, score) ;
      break ;
    case '5':
      fMapAddGfSeg (comp ? SPLICE5_UP : SPLICE5, M_GF_splice,
		    comp ? pos-2 : pos,  comp ? pos-1 : pos+1, score) ;
      break ;
    case 'a':
      fMapAddGfSeg (comp ? ATG_UP : ATG, M_GF_ATG,
		    comp ? pos-3 : pos,  comp ? pos-1 : pos+2, score) ;
      break ;
    }
}

void fMapAddGfCodingSeg (int pos1, int pos2, float score, BOOL comp)
{

  fMapInitialise () ;
  
  fMapAddGfSeg (comp ? FEATURE_UP : FEATURE, M_GF_coding_seg, 
		pos1-1, pos2-1, score) ;
}

void fMapAddGfSegs (void)
{
  int i, len ;
  SEG *seg ;
  char safe=0 ;
  FMAPLOOKGET ("fMapAddGfSegs") ;

  initialise () ;

  method (0, M_GF_ATG)->flags |= METHOD_CALCULATED ;
  method (0, M_GF_coding_seg)->flags |= METHOD_CALCULATED ;
  method (0, M_GF_splice)->flags |= METHOD_CALCULATED ;
  method (0, M_hexIntron)->flags |= METHOD_CALCULATED ;
  method (0, M_hexExon)->flags |= METHOD_CALCULATED ;

	/* first delete temp segs */
  arrayMax (look->segs) = look->lastTrueSeg ;
	/* then current gf segs */
  for (i = 1 ; i < arrayMax (look->segs) ; ++i)
    { seg = arrp (look->segs, i, SEG) ;
      if (seg->key == M_GF_coding_seg ||
	  seg->key == M_GF_ATG ||
	  seg->key == M_GF_splice ||
	  seg->key == M_hexExon ||
	  seg->key == M_hexIntron)
	{ *seg = arr (look->segs, arrayMax (look->segs)-1, SEG) ;
	  --arrayMax (look->segs) ; /* RD must do after, for ARR_CHECK */
	  --i ;			/* so test the new item */
	}
    }

  gfLook = look ;

  look->gf.min = look->map->centre - gf_range ;
  if (look->gf.min < 0) 
    look->gf.min = 0 ;
  if (look->gf.min > look->zoneMin)
    look->gf.min = look->zoneMin ;

  look->gf.max = look->map->centre + gf_range ;
  if (look->gf.max > arrayMax (look->dna)) 
    look->gf.max = arrayMax (look->dna) ;
  if (look->gf.max < look->zoneMax)
    look->gf.max = look->zoneMax ;

  len = look->gf.max - look->gf.min ;

  if (look->gf.cum) messfree (look->gf.cum) ;
  look->gf.cum = (float*) halloc (sizeof (float)* (len+1), graphHandle ()) ;
  if (look->gf.revCum) messfree (look->gf.revCum) ;
  look->gf.revCum = (float*) halloc (sizeof (float)* (len+1), graphHandle ()) ;

				/* ensure 0 termination is safe */
  array (look->dna, arrayMax (look->dna), char) = 0 ;
  --arrayMax (look->dna) ;

  if (look->gf.max < arrayMax (look->dna))
    { safe = arr (look->dna, look->gf.max, char) ;
      arr (look->dna, look->gf.max, char) = 0 ;
    }
  else
    safe = 0 ;
  geneFinderAce (arrp (look->dna,look->gf.min,char), &look->gf) ;

  if (look->gf.max < arrayMax (look->dna))
    arr (look->dna, look->gf.max, char) = safe ;

  hexAddSegs ("wgf/intron", FEATURE, M_hexIntron, 1, 10, FALSE,
	      arrp (look->dna,look->gf.min,char), len, 0) ;
  hexAddSegs ("wgf/intron", FEATURE_UP, M_hexIntron, 1, 10, TRUE,
	      arrp (look->dna,look->gf.min,char), len, 0) ;

	/* exon must come second, to fill cum properly */
  hexAddSegs ("wgf/cds", FEATURE, M_hexExon, 3, 10, FALSE,
	      arrp (look->dna,look->gf.min,char), len, look->gf.cum) ;
  hexAddSegs ("wgf/cds", FEATURE_UP, M_hexExon, 3, 10, TRUE,
	      arrp (look->dna,look->gf.min,char), len, look->gf.revCum) ;

  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs) ;

  fMapProcessMethods (look) ;
  fMapShowGfColumns () ;	/* also does fMapDraw () */
}

/**********************************************************/
/************** code to add ATG, splice sites *************/

static void addATG (void)
{
  SEG *seg;
  int x;
  FMAPLOOKGET ("addATG");

  x = look->DNAcoord + 1;

  initialise ();

  if (!look->activeBox)
    return;
  seg = BOXSEG (look->activeBox);
  if (seg->type != DNA_SEQ)
    return;

  arrayMax (look->segs) = look->lastTrueSeg;
  seg = arrayp (look->segs, arrayMax (look->segs), SEG);
  seg->type = ATG;
  seg->key = M_GF_ATG;
  seg->x1 = x-1;
  seg->x2 = x+1;
  seg->data.f = 4.0;
  assInsert (look->chosen, SEG_HASH (seg), assVoid (seg->x1)) ;
  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs);
  fMapProcessMethods (look);
  fMapDraw (look, 0);
}

static void addSplice3 (void)
{
  SEG *seg;
  int x;
  FMAPLOOKGET ("addATG");

  x = look->DNAcoord + 1;

  initialise ();

  if (!look->activeBox)
    return;
  seg = BOXSEG (look->activeBox);
  if (seg->type != DNA_SEQ)
    return;

  arrayMax (look->segs) = look->lastTrueSeg;
  seg = arrayp (look->segs, arrayMax (look->segs), SEG);
  seg->type = SPLICE3;
  seg->key = M_GF_splice;
  seg->x1 = x-2;
  seg->x2 = x-1;
  seg->data.f = 4.0;
  assInsert (look->chosen, SEG_HASH (seg), assVoid (seg->x1)) ;
  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs);
  fMapProcessMethods (look);
  fMapDraw (look, 0);
}

static void addSplice5 (void)
{ SEG *seg;
  int x;
  FMAPLOOKGET ("addATG");

  x = look->DNAcoord +1;

  initialise ();

  if (!look->activeBox)
    return;
  seg = BOXSEG (look->activeBox);
  if (seg->type != DNA_SEQ)
    return;

  arrayMax (look->segs) = look->lastTrueSeg;
  seg = arrayp (look->segs, arrayMax (look->segs), SEG);
  seg->type = SPLICE5;
  seg->key = M_GF_splice;
  seg->x1 = x-1;
  seg->x2 = x;
  seg->data.f = 4.0;
  assInsert (look->chosen, SEG_HASH (seg), assVoid (seg->x1)) ;
  arraySort (look->segs, fMapOrder) ;
  look->lastTrueSeg = arrayMax (look->segs);
  fMapProcessMethods (look);
  fMapDraw (look, 0) ;
}

/****************************************/

static void showSelect (void) ;

MENUOPT fMapGeneOpts[] = {
  { fMapAddGfSegs, "Genefinder Features"},
  { fMapDP1, "Autofind one gene"},
  { fMapDP, "Autofind genes"},
  { parmsEditor, "Parameters"},
  { fMapSetChoose, "Gene -> Selected"},
  { fMapChooseToGene, "Selected -> temp_gene"},
  { fixTemp, "Fix temp_gene"},
  { killGeneMenu, "Remove Gene"},
  { showSelect, "Show Selected"},
  { addATG, "New Gene Start (ATG)"},
  { addSplice3, "New Exon Start"},
  { addSplice5, "New Exon End"},
  { fMapHideGfColumns, "Hide Gf columns"},
  { fMapShowGfColumns, "Show Gf columns"},
  {0, 0}} ;

/**********************************************************************/
/************** code to show the selected set of features *************/

static MENUOPT selectMenu[] = {
  { graphDestroy, "Quit" },
  { graphPrint, "Print" },
  { showSelect, "Update" },
  { 0, 0 }} ;

void addLine (LOOK look, 
	      int *pStart, int *p3, int *p5, int *pStop, 
	      int line)
{
  static int frame ;
  static int oldLine = 10000 ;
  static float score ;
  int i ;
  SEG *seg ;
  static ELT elt1, elt2 ;
  float tscore = 0 ;

  elt1.type = DUMMY ;

  if (*pStart)
    { for (i = 0 ; i < arrayMax (look->segs) ; ++i)
	{ seg = arrp (look->segs, i, SEG) ;
	  if (seg->type == ATG && seg->x1 == *pStart)
	    { graphText (messprintf ("%5.2f", seg->data.f), 
			 22, line) ;
	      tscore += seg->data.f ;
	      break ;
	    }
	}
      frame = *pStart % 3 ;
      score = 0 ;
      elt2.type = DUMMY ;	/* prevents INTRON score */
      elt1.type = START ;
      elt1.x = *pStart ;
    }
  if (*p3)
    { for (i = 0 ; i < arrayMax (look->segs) ; ++i)
	{ seg = arrp (look->segs, i, SEG) ;
	  if (seg->type == SPLICE3 && seg->x2 == *p3)
	    { graphText (messprintf ("%5.2f", seg->data.f), 
			 22, line) ;
	      tscore += seg->data.f ;
	      break ;
	    }
	}
      elt1.type = SPLICE_3 ;
      elt1.x = *p3 ;
    }

  if (oldLine < line && elt1.type)
    { graphText (messprintf ("%6.2f", intronScore (&elt2, &elt1)),
		 50, line - 0.5) ;
      score += intronScore (&elt2, &elt1) ;
      frame += (elt1.x - elt2.x - 1) ; frame %= 3 ;
    }

  elt2.type = DUMMY ;
  if (*p5)
    { for (i = 0 ; i < arrayMax (look->segs) ; ++i)
	{ seg = arrp (look->segs, i, SEG) ;
	  if (seg->type == SPLICE5 && seg->x1 == *p5)
	    { graphText (messprintf ("%5.2f", seg->data.f), 
			 35, line) ;
	      tscore += seg->data.f ;
	    }
	}
      elt2.type = SPLICE_5 ;
      elt2.x = *p5 ;
    }
  if (*pStop)
    { elt2.type = STOP ;
      elt2.x = *pStop - 3 ;
    }

  if (elt1.type && elt2.type)
    { graphText (messprintf ("%6.2f", 
			     exonScore (look, &elt1, &elt2, 
					frame, &look->gf)),
		 28, line) ;
      tscore += exonScore (look, &elt1, &elt2, frame, &look->gf) ;
      oldLine = line ;
    }

  if (tscore)
    graphText (messprintf ("%6.2f", tscore), 42, line) ;

  score += tscore ;
  graphText (messprintf ("%6.2f", score), 58, line) ;

  *pStart = *p3 = *p5 = *pStop = 0 ;
}

/******************************/

static void showSelect (void)
{
  Array chosen ;
  char *v, *vx ;
  int i, line = 2 ;
  int xStart, x3, x5, xStop ;
  Graph oldGraph ;
  FMAPLOOKGET ("showSelect") ;

#define START_OFF 1
#define S3_OFF    3
#define S5_OFF	  11
#define STOP_OFF  13

  if (! (oldGraph = graphActivate (look->selectGraph)))
    { look->selectGraph = graphCreate (TEXT_SCROLL, 
				       "Gene structure", 
				       0.0, 0.0, 0.7, 0.2) ;
      graphMenu (selectMenu) ;
    }
  else
    { graphPop () ;
      graphClear () ;
    }

  graphText ("ATG/3'   5'/STOP    feat1 coding feat2    "
	     "Exon  Intron   Total",
	     2, 0.4) ;
  graphLine (0, 1.7, 64, 1.7) ;

  chosen = arrayCreate (32, char*) ;
  v = 0 ; 
  while (assNext (look->chosen, &v, &vx))
    { 
#ifdef GF_SEG_DEFINED
      x = assInt (vx) ;
      if (UNHASH_TYPE (v) == GF_SEG)
	v = HASH (GF_SEG, x + 3* ((UNHASH_X2 (v) - x)/6)) ;
#endif
      array (chosen, arrayMax (chosen), char*) = v ;
    }

  arraySort (chosen, compareChoseOrder) ;

  xStart = x3 = x5 = xStop = 0 ;

  for (i = 0 ; i < arrayMax (chosen) ; ++i)
    { v = arr (chosen, i, char*) ;
      switch (UNHASH_TYPE (v))
	{
	case ATG:
	  if (xStart || x3 || x5 || xStop)
	    addLine (look, &xStart, &x3, &x5, &xStop, line++) ;
	  xStart = UNHASH_X2 (v) - 2 ;
	  graphText (messprintf ("%7d", COORD (look, xStart)),
		     START_OFF, line) ;
	  break ;
	case SPLICE3:
	  if (xStart || x3 || x5 || xStop)
	    addLine (look, &xStart, &x3, &x5, &xStop, line++) ;
	  x3 = UNHASH_X2 (v) ;
	  graphText (messprintf ("%7d", COORD (look, x3)),
		     S3_OFF, line) ;
	  break ;
	case SPLICE5:
	  if (x5 || xStop)
	    addLine (look, &xStart, &x3, &x5, &xStop, line++) ;
	  x5 = UNHASH_X2 (v) - 1 ;
	  graphText (messprintf ("%7d", COORD (look, x5)),
		     S5_OFF, line) ;
	  break ;
	case ORF:
	  if (x5 || xStop)
	    addLine (look, &xStart, &x3, &x5, &xStop, line++) ;
	  xStop = UNHASH_X2 (v) + 3 ;
	  graphText (messprintf ("%7d", COORD (look, xStop)),
		     STOP_OFF, line) ;
	  break ;
	}
    }
  if (xStart || x3 || x5 || xStop)
    addLine (look, &xStart, &x3, &x5, &xStop, line++) ;

  arrayDestroy (chosen) ;

  graphLine (21, 0, 21, line) ;
  graphLine (49, 0, 49, line) ;
  graphTextBounds (64, line) ;
  graphRedraw () ;
}
 
/******************** autocreate peptide code ***********************/

static int calculateAndSaveChecksum (KEY protein){
  Array pep;
  int x1,i = 0;
  OBJ obj;
  KEY key;

   if (! (pep = peptideGet (protein)) || arrayMax (pep) < 2){
      return 0;
  }
  else 
    i = hashArray (pep);
  if (pep)
    arrayDestroy (pep);
  if (i){
    if ((obj = bsUpdate (protein))){
      bsFindTag (obj,_Peptide);
      bsGetKey (obj,_Peptide,&key); /* peptide */
      bsGetData (obj,_bsRight, _Int, &x1); /* length */
      bsAddData (obj,_bsRight, _Int,&i);   /* add the checksum */
      bsSave (obj);
    }
  }
  return i;
}

static KEY findPeptideMatch (KEY temp,int newchk,KEY key2) 
{
  KEY peptide = 0,key;
  int i,x1;
  Array pep1,pep2;
  BOOL found=TRUE;
  OBJ obj;
  char *p1,*p2;

/*  printf ("Searching for subclass %s\n",name (key2));*/
  while (lexNext (_VProtein,&peptide)){
    if (lexIsInClass (peptide,key2)){
      if ((obj = bsCreate (peptide))){
	if (bsFindTag (obj,str2tag ("Peptide"))){
	  if (bsGetKey (obj,_Peptide,&key)){              /* peptide */
	    if (!bsGetData (obj,_bsRight, _Int, &x1))     /* length */
	      continue;
	    x1 = 0;
	    bsGetData (obj,_bsRight, _Int, &x1);         /* checksum */
	    bsDestroy (obj);
	    if (x1==0)
	      x1 = calculateAndSaveChecksum (peptide);
	    if (x1==newchk && peptide != temp){ /* if same checksum and not the one we are testing against */
	      freeOut ("Same check sum found\n");
	      pep1 = peptideGet (temp);
	      pep2 = peptideGet (peptide);
	      found = TRUE;
	      if (arrayMax (pep1) == arrayMax (pep2)){ /* If they are the same length */
		freeOut ("They have the same length\n");
		p1 = arrp (pep1, 0, char) ;
		p2 = arrp (pep2, 0, char) ;
		for (i=0;i<arrayMax (pep1);i++){
		  if (pepDecodeChar[ (int)*p1++] != pepDecodeChar[ (int)*p2++]){
		    found = FALSE;
		    freeOut ("Not the same sequence though\n");
		    break;
		  }
		}
	      }
	      arrayDestroy (pep1);
	      arrayDestroy (pep2);
	      if (found)
		freeOut ("Identical peptides\n");
	      if (found)
		return peptide;
	    }
	  }
	}
	else
	  bsDestroy (obj);
      }
    }
  }
  return 0;
}
 
/*********************** end of file ************************/
 
 
 
