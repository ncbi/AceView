/*  File: geneannot.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
   Display cDNA
 * Exported functions:
 * HISTORY:
 * Created: Thu Jan   9 00:01:50 1999 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#include "display.h"

#include "bs.h"
#include "lex.h"
#include "systags.h"
#include "classes.h"
#include "dna.h"
#include "dnaalign.h"
#include "query.h"
#include "bindex.h"
#include "cdna.h"
#include "acembly.h"
#include "annot.h"

#include "cdna.h"


/*********************************************************/
/*********************************************************/

typedef struct SegStruct *SEG ;
typedef struct LookStruct* LOOK ;


typedef BOOL (*GASAVEFUNC) (LOOK look, SEG seg) ;
typedef void (*GACHECKFUNC) (char *text) ;
typedef void (*GADRAWFUNC) (LOOK look, SEG seg, int x, int *yp) ;

struct SegStruct 
{ 
  SEG right ;
  SEG down ;
  BOOL isTriangle ;
  BOOL open, fixedFont ; 
  int dx ;  /* skip dx char to the right */
  KEY tag, tag2 ;
  char *title, *tagName, *buf ;
  char *title2, *tagName2, *buf2 ;
  int nn, nn2 ;  /* length of tag text entry */
  GADRAWFUNC draw ;
  GACHECKFUNC check ;
  GASAVEFUNC save ;
} ;

struct LookStruct 
{
  void *magic;        /* == &GENEANNOT__MAGIC */
  KEY   key, tg ;
  Array boxes ;
  Array box2segs ;
  int activeBox ;
  Graph graph ;
  AC_HANDLE handle ;
  int openAll ;
  SEG seg ;
}  ;


#define LOOKGET(name) LOOK look ; \
                      if (!graphAssFind (&GENEANNOT_MAGIC,&look)) \
		        messcrash ("geneannot-graph not found in %s",name) ; \
		      if (!look) \
                        messcrash ("%s received a null look",name) ; \
                      if (look->magic != &GENEANNOT_MAGIC) \
                        messcrash ("%s received a non-magic look",name)

#define HOOK(_s0,_s2,_s) {if (_s2) {(_s2)->down = _s; _s2=(_s2)->down;}  \
	                else {(_s0)->right=(_s2)=_s;} \
                        while ((_s2)->down) _s2 = (_s2)->down ;\
                        (_s)->open=TRUE;}

/*********************************************************/
/*********************************************************/

static void lookDestroy (void) ;
static void lookPick (int box) ;
static void lookRedraw (void) ;
static void lookKbd (int k) ;
static void lookDraw (LOOK look) ;
static BOOL gaConvert (LOOK look) ;

/********** end of file **********/
 

#ifdef JUNK

static ANNOTATION  notes[] = 
{
  { 1, "Ready for submission",  "Ready_for_submission",  0, 'k', 0, {0}},
  { 1, "Fully edited",  "Fully_edited",  0, 'k', 0, {0}},
  { 1, "Alias", "Alias", 0, _Text, 1000, {0}},
  { 1, "Genefinder", 0, 0, 'k', 0, {0}},
  { 4,    "Identical to predicted gene", "Identical_prediction", 0, 'k', 0, {0}},
  { 4,    "Similar to predicted gene", "Compatible_prediction", 0, 'k', 0, {0}},
  { 4,    "Different", "Different", 0, 'k', 0, {0}},
  { 4,    "cDNA covers several predicted genes", "cDNA_covers_several_predicted_genes", 0, 'k', 0, {0}},
  { 4,    "Predicted covers several genes", "Predicted_covers_several_genes", 0, 'k', 0, {0}},
  { 1,  "Approximate genetic position", 0, 0, 'k', 0, {0}}, 
  { 4,    "Chromosome", "Chromosome", 0, 'k', 0, {0}},
  { 4,    "Position:", "Genetic_position", 0, 'f', 0, {0}},
  { 1,  "Spread over genomic", 0, 0, 'k', 0, {0}}, 
  { 4,    "from sequence", "from", 0, 'k', 0, {0}},
  { 4,    "popsition", "from:2", 0, 'i', 0, {0}},
  { 4,    "to sequence", "to_sequence", 0, 'k', 0, {0}},
  { 4,    "to sequence", "to_sequence:2", 0, 'i', 0, {0}},
  { 1, "Alternative products",  0, 0, 'k', 0, {0}},
  { 4,   "Possibly more",  "Possibly_more", 0, 'k', 0, {0}},
  { 1, "Sequencing error in genomic DNA",  "Sequencing_error_in_genomic_dna",  0, _Text, 1000, {0}},
  { 4,    "No Gap", "No_Gap", 0, 'k', 0, {0}},
  { 4,    "Gap in single open exon", "Single_gap", 0, 'k', 0, {0}},
  { 4,    "Gap", "Gaps", 0, 'k', 0, {0}},
  { 4,    "Single clone", "Single_clone", 0, 'k', 0, {0}},
  { 4,    "5-prime", 0, 0, 'k', 0, {0}},
  { 8,          "SL1",  "SL1", 0, 'k', 0, {0}},
  { 8,          "SL2",  "SL2", 0, 'k', 0, {0}},
  { 8,          "gccgtgctc",  "gccgtgctc", 0, 'k', 0, {0}},
  { 8,          "5p motif",  "5p_motif", 0, _Text, 1000, {0}},
  { 8,          "Many clones start here",  "Many_clones_start_here", 0, 'i', 0, {0}},
  { 8,          "5p alternatif",  "Alternative_5prime", 0, 'k', 0, {0}},
  { 4,    "Splicing", 0, 0, 'k', 0, {0}},
  { 8,          "Alternative splicing",  "Alternative_splicing", 0, 'k', 0, {0}},
  { 12,             "Coupled to alt 3p",  "coupled_to_alt_3p", 0, 'i', 0, {0}},
  { 12,             "Alternative exon",  "Alternative_exon", 0, 'k', 0, {0}},
  { 12,             "Overlapping alternative exons",  "Overlapping_alternative_exons", 0, 'k', 0, {0}},
  { 12,             "3_6_9 splice variations",  "3_6_9_splice_variations", 0, 'k', 0, {0}},
  { 8,          "Confirmed non gt_ag",  "Confirmed_non_gt_ag", 0, 'k', 0, {0}},
  { 4,    "3-prime", 0, 0, 'k', 0, {0}},
  /*   { 8,        "Poly A seen in 5p read",  "Poly_A_seen_in_5p_read", 0, 'k', 0, {0}}, */
  { 8,        "Alternative poly A",  "Alternative_poly_A", 0, 'k', 0, {0}},
  { 8,        "No Alternative poly A",  "No_Alternative_poly_A", 0, 'k', 0, {0}},
  { 8,        "Fake_internal poly A",  "Fake_internal_poly_A", 0, 'k', 0, {0}},
  /*
 { 1, "Translation", 0, 0, 'k', 0, {0}},
  { 4,    "Open from first to last exon",  "Open_from_first_to_last_exon", 0, 'k', 0, {0}},
  { 4,    "Nb introns after stop", "Nb_introns_after_stop", 0, _Int, 0, {0}},
  */
  { 1, "Problem",  0, 0, 'k', 0, {0}},
  { 4,    "Cosmid boundary",  "Cosmid_boundary", 0, 'k', 0, {0}},
  { 4,    "Lost reads_Bug",  "Lost_reads", 0, 'k', 0, {0}},
  { 4,    "Jack Pot Bug",  "Jack_Pot_Bug", 0, 'k', 0, {0}},
  { 4,    "Tag Bug",  "Tag_Bug", 0, 'k', 0, {0}},
  { 4,    "Other Program Bug",  "Other_Bug", 0, _Text, 1000, {0}},
  /*
*/
  { 0, 0, 0, 0, 0, 0, {0}}
} ;

#endif

static void* GENEANNOT_MAGIC ;

static MENUOPT gaMenu[] = {
  {graphDestroy, "Quit"},
  {help, "Help"},
  {graphPrint, "Print"},
  {displayPreserve, "Preserve"},
  {0, 0} } ;

KEY 
  _VGeneAnnot ;

static void geneAnnotInit (void)
{
  int i ; KEY key ;

  if (_Hit)
    return ;

  i = 256 ;  while (i--) B2[i] = 0 ;
  B2[A_] = 0x0 ; B2[T_] = 0x3 ;
  B2[G_] = 0x1 ; B2[C_] = 0x2 ;

  lexaddkey ("Gene_annotation", &key, _VMainClasses) ; 
  _VGeneAnnot = KEYKEY(key) ;
  lexaddkey ("Transcribed_gene", &key, _VMainClasses) ; 
  _VTranscribed_gene = KEYKEY(key) ;
  lexaddkey ("cDNA_clone", &key, _VMainClasses) ; 
  _VcDNA_clone = KEYKEY(key) ;
  lexaddkey ("Clone_Group", &key, _VMainClasses) ; 
  _VClone_Group = KEYKEY(key) ;
  lexaddkey ("Annotation", &key, _VMainClasses) ; 
  _VAnnotation = KEYKEY(key) ;

  _Transcribed_gene = str2tag ("Transcribed_gene") ;  
  _From_gene = str2tag ("From_gene") ;  
  _cDNA_clone = str2tag ("cDNA_clone") ; 
  _Hit = str2tag ("Hit") ; 
  _Forward = str2tag ("Forward") ;
  _Reverse = str2tag ("Reverse") ;
  _Confirmed_intron = str2tag ("Confirmed_intron") ;
  _Confirmed_exon = str2tag ("Confirmed_exon") ;
  _First_exon = str2tag ("First_exon") ;
  _Last_exon = str2tag ("Last_exon") ;
  _Splicing = str2tag ("Splicing") ;
  _Read = str2tag ("Read") ;
  _Clone_Group = str2tag ("Clone_Group") ;
  _Nb_possible_exons = str2tag ("Nb_possible_exons") ; 
  _Nb_alternative_exons = str2tag ("Nb_alternative_exons") ; 
  _Nb_confirmed_introns = str2tag ("Nb_confirmed_introns") ; 
  _Nb_confirmed_alternative_introns = str2tag ("Nb_confirmed_alternative_introns") ; 
  _Intron = str2tag ("Intron") ; 
  _Exon = str2tag ("Exon") ; 
  _Gap = str2tag ("Gap") ; 
  _Alternative_exon = str2tag ("Alternative_exon") ; 
  _Alternative_intron = str2tag ("Alternative_intron") ; 
  _Partial_exon = str2tag ("Partial_exon") ; 
  _Alternative_partial_exon = str2tag ("Alternative_partial_exon") ; 
  _Begin_not_found = str2tag ("Begin_not_found") ; 
  /*   _End_not_found = str2tag ("End_not_found") ;  in tags.h */
  _Total_length = str2tag ("Total_length") ; 
  _Total_intron_length = str2tag ("Total_intron_length") ; 
  _Total_gap_length = str2tag ("Total_gap_length") ; 
  _Transpliced_to = str2tag ("Transpliced_to") ; 
  _Transcribed_from = str2tag ("Transcribed_from") ; 
  _Genomic_sequence = str2tag ("Genomic_sequence") ; 
  _PolyA_after_base = str2tag ("PolyA_after_base") ;  
  _Intron_boundaries = str2tag ("Intron_boundaries") ;  
  _Discarded_cDNA = str2tag ("Discarded_cDNA") ; 
  _Discarded_from = str2tag ("Discarded_from") ; 
  _CTF_File = str2tag ("CTF_File") ;  
  /*   _ = str2tag ("") ;  */
}


/*********************************************************/
/********************************************************************/

static SEG gaNewSeg (LOOK look)
{
  return  (SEG) halloc (sizeof (struct SegStruct),look->handle) ;
}

/********************************************************************/
/*********************************************************/


BOOL geneAnnotDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  LOOK look = 0 ;
  KEY newName = 0 ;

  
  geneAnnotInit () ;


  if (!key && (!from || class(from) != _VTranscribed_gene))
    goto abort ;

  if (!from && !
      (from = keyGetKey (key, _Transcribed_gene)))
    goto abort ;

  look= (LOOK) messalloc(sizeof(struct LookStruct));
  look->magic = &GENEANNOT_MAGIC ;
  look->key = key;

  look->tg = from ;
      
  if (isOldGraph)
    { 
      messfree (look->handle) ;
      graphClear () ;
      graphGoto (0,0) ;
    }
  else
    { 
      if (!displayCreate ("DtGene_Annot") &&
	  !displayCreate (TREE))
	goto abort ;
      graphRegister (DESTROY, lookDestroy) ;
      graphRegister (PICK, lookPick) ;
      graphRegister (KEYBOARD, lookKbd) ;
      graphRegister (RESIZE, lookRedraw) ;
      graphMenu (gaMenu) ;
    }
  
  look->handle = handleCreate () ;
  look->graph = graphActive() ;
  graphAssociate (&GENEANNOT_MAGIC, look) ;
  look->box2segs = arrayHandleCreate (256, SEG, look->handle) ;
  newName = keyGetKey(from,str2tag("NewName")) ;
  graphRetitle (messprintf("%s: %s %s", className(from),name (from), 
			   newName ? name(newName) : "")) ;
  
  if (gaConvert (look))
    {
      look->boxes = 0 ;
      look->activeBox = 0 ;
      look->openAll = 0 ;
      lookDraw (look) ;
      
      return TRUE ;
    }
  
 abort:
  if (look)
    lookDestroy () ;
  return FALSE ;
}

/************************************************************/

static void lookDestroy (void)
{
  LOOKGET("lookDestroy") ;


  if (look && look->magic == GENEANNOT_MAGIC)
    {
      messfree (look->handle) ;
      look->magic = 0 ;
      messfree (look) ;
    }

  graphAssRemove (&GENEANNOT_MAGIC) ;
}

/************************************************************/

static void gaDrawTriangle (SEG seg, int x, int yp, BOOL global)
{
  Array aa = 0 ;
  int i, box ;
  int dx = 0, dy;
  float d1x = 0, d2x = 0, e1y = 0, e2y = 0 ;
  float width = 0, height = 0, aspect = 1.625;
  
  LOOKGET("gaDrawTriangle") ;

  if(graphActive()) 
    { graphTextInfo(&dx, &dy, &width, &height); 
    if(dx) aspect =  ((float) dy)/dx; }
  /* triangle equilateral de cote 0.8 hauteur de caractere */
  d1x = (global ? 1.6 : 0.8) * aspect ;
  e1y = (1 - 0.866 * d1x / aspect) * 0.5 ;
  d2x = d1x * 0.866 ;
  e2y = (1 - d1x / aspect) * 0.5 ;

  array(look->box2segs, box = graphBoxStart(), SEG) = seg ;
  seg->isTriangle = TRUE ;
  aa = arrayCreate (8, float) ; i = 0 ;
  if (seg->open)
    {
      array (aa, i++, float) = x ;
      array (aa, i++, float) = yp + e1y ;
      array (aa, i++, float) = x + d1x * 0.5 ;
      array (aa, i++, float) = yp + 1 - e1y ;
      array (aa, i++, float) = x + d1x ;
      array (aa, i++, float) = yp + e1y ; 
      array (aa, i++, float) = x ;
      array (aa, i++, float) = yp + e1y ; 
    }
  else
    {
      array (aa, i++, float) = x ;
      array (aa, i++, float) = yp + e2y;
      array (aa, i++, float) = x ;
      array (aa, i++, float) = yp + 1 - e2y;
      array (aa, i++, float) = x + d2x ;
      array (aa, i++, float) = yp + 0.5 ;
      array (aa, i++, float) = x ;
      array (aa, i++, float) = yp + e2y ;

    }
  graphColor (RED) ;
  graphPolygon (aa) ;
  graphColor (BLACK) ;
  arrayDestroy (aa) ;
  graphBoxEnd();
}

/*******************/

static void gaDrawSeg (LOOK look, SEG seg, int x, int *yp)
{
  if (look->openAll == 1)
    seg->open = FALSE ;
  else if (look->openAll == 2)
    seg->open = TRUE ;
  if (seg == look->seg)
    gaDrawTriangle (seg, x - 2, *yp, TRUE) ;
  if (seg->right)
    gaDrawTriangle (seg, x - 2, *yp, FALSE) ;
  if (seg->draw)
    (*seg->draw)(look, seg, x, yp) ;
  if (seg->right && seg->open)
    gaDrawSeg (look, seg->right, x + (seg->dx ? seg->dx : 2), yp) ;
  if (seg->down)
    gaDrawSeg (look, seg->down, x, yp) ;
}

/*******************/

static void lookDraw (LOOK look)
{
  int x = 3, y = 1 ;
  if (!graphActivate (look->graph))
    return ;

  graphClear () ;
  
  gaDrawSeg (look, look->seg, x, &y) ;
  look->openAll = 0 ;
  graphTextBounds (80, y) ;
  graphRedraw () ;
}

/*******************/

static void lookRedraw (void)
{
  LOOKGET("lookRedraw") ;

  look->box2segs = arrayReCreate (look->box2segs, 256, SEG) ;
  lookDraw(look) ;
}

/*****************/

static void lookKbd (int k)
{ 
  LOOKGET("lookKbd") ;

  if (isDisplayBlocked())
    return ;

  if (!look->activeBox) return ;

  switch (k)
    {
    case LEFT_KEY:
      break ;
    case RIGHT_KEY:
      break ;
    case DOWN_KEY:
      break ;
    case UP_KEY:
      break ;
    default:
      break ;
    }
}


/*********************************************************/

static void lookPick (int box)
{
  SEG seg ;

  LOOKGET("lookPick") ;

 if (!box || (box >= arrayMax (look->box2segs)))
   return ;
 else if ((seg = array(look->box2segs, box, SEG)) && seg)
   {
     if (seg->isTriangle)
       {
	 seg->open = !seg->open ;
	 if (seg == look->seg)
	   look->openAll = seg->open ? 2 : 1 ;
	 lookRedraw() ;
       }
     else if (seg->buf && *seg->buf)
       {
	 int i = strlen(seg->buf) ;
	 
	 if (i < 50) 
	   graphPostBuffer (seg->buf) ;
	 else
	   {
	     Stack s = stackCreate((52*i + 50)/50) ;
	     char *cp = seg->buf, cc ;
	     
	     while (strlen(cp) > 50)
	       {
		 cc = *(cp+50) ;
		 *(cp+50) = 0 ;
		 catText(s,cp) ;
		 catText(s,"\n") ;
		 *(cp+50) = cc ;
		 cp += 50 ;
	       }
	     catText (s,cp) ;
	     catText(s,"\n") ;
	     graphPostBuffer (stackText(s,0)) ;
	     stackDestroy (s) ;
	   }
	     
	 if (look->activeBox)
	   graphBoxDraw(look->activeBox, BLACK, WHITE) ;
	 look->activeBox = 0 ;
	 if (seg->buf)
	   {
	     look->activeBox = box ;
	     graphBoxDraw(look->activeBox, BLACK, LIGHTBLUE) ;
	   }
       }
   }
}

/*********************************************************/
/*********************************************************/

static void gaDrawHeader (LOOK look, SEG seg, int x, int *yp)
{
  KEYSET ks = queryKey (look->tg, ">Transcript") ; 
  KEY newName = 0 ;
  
  newName = keyGetKey(look->tg,str2tag("NewName")) ;
  graphText (messprintf ("Transcript: %s   %s   %d mRNA", 
			 name(look->tg),  newName ? name(newName) : "",
			 keySetMax(ks)), x, *yp) ;
  keySetDestroy (ks) ;
  *yp += 2 ;
}

/*********************************************************/

static void gaDrawTag (LOOK look, SEG seg, int x, int *yp)
{
  int box ; /*mhmp*/
  if (seg->title)
    {
      graphText (seg->title, x, *yp) ;
      x += strlen(seg->title) + 3 ;
    }
  else
    x += 3 ;
  if (seg->fixedFont)
    graphTextFormat (FIXED_WIDTH) ;
  array(look->box2segs, box = graphBoxStart(), SEG) = seg ;
  if (seg->nn && seg->save)
    graphTextEntry (seg->buf, seg->nn - 1, x, *yp, seg->check) ;
  else
    graphText (seg->buf, x, *yp) ;
  if (!seg->dx) (*yp)++ ;
  graphBoxEnd() ;
  graphBoxDraw (box, BLACK, TRANSPARENT) ;
  graphTextFormat (PLAIN_FORMAT) ;  
}

/*********************************************************/

static void gaDrawArray (LOOK look, SEG seg, int x, int *yp)
{
  int box ; /*mhmp*/
  int i, max, len = 50 ;
  char cc = 0, *cp ;
  
  if (seg->title)
    {
      graphText (seg->title, x, *yp) ;
      x += strlen(seg->title) + 3 ;
    }
  else
    x += 3 ;
  if (seg->fixedFont)
    graphTextFormat (FIXED_WIDTH) ;
  array(look->box2segs, box = graphBoxStart(), SEG) = seg ;
  max = strlen(seg->buf) ;

  for (i = 0, cp = seg->buf ; i < max ; i += len, cp += len)
    {
      if (max - i > len)
	{ 
	  cc = *(cp + len) ;
	  *(cp + len) = 0 ;
	}
      graphText (cp, x, (*yp)++) ;
      if (max - i > len)
	*(cp + len) = cc ;
    }
  if (!seg->dx) (*yp)++ ;
  graphBoxEnd() ;
  graphBoxDraw (box, BLACK, TRANSPARENT) ;
  graphTextFormat (PLAIN_FORMAT) ;  
}

/*********************************************************/

static void gaDrawDouble (LOOK look, SEG seg, int x, int *yp)
{
  if (seg->title)
    {
      graphText (seg->title, x, *yp) ;
      x += strlen(seg->title) + 3 ;
    }
  else
    x += 3 ;
  if (seg->nn && seg->save)
    graphTextEntry (seg->buf, seg->nn - 1, x, *yp, seg->check) ; 
  else if (seg->nn)
    graphText (seg->buf, x, *yp) ;
  x += seg->nn + 4 ;

 if (seg->title2)
    {
      graphText (seg->title2, x, *yp) ;
      x += strlen(seg->title2) + 3 ;
    }
  else
    x += 3 ; 
   if (seg->nn2 && seg->save)
     graphTextEntry (seg->buf2, seg->nn2 - 1, x, *yp, seg->check) ;
   else if (seg->nn2)
     graphText (seg->buf2, x, *yp) ;
  (*yp)++ ;
}
/*********************************************************/
/*********************************************************/
#ifdef JUNK
static BOOL gaCheckTrivial (char *text)
{
  return TRUE ;
}

/*********************************************************/
/*********************************************************/

static BOOL gaSaveTrivial (LOOK look, SEG seg) 
{
  return TRUE ;
}
#endif   /* JUNK */
/*********************************************************/
/*********************************************************/

static SEG gaConvertGenomic (LOOK look)
{
  SEG seg, seg2 ;
  KEY tag, cosmid = 0 ;
  OBJ obj = 0 ;
  BOOL ok = FALSE ;
  int a1, a2, len ;

  seg = gaNewSeg (look) ;
  seg->buf = halloc (12, look->handle) ;
  seg->title = "Begins in genomic sequence" ;
  seg->title2 = "at or before base" ;
  seg->nn = 12;
  seg->buf2 = halloc (12, look->handle) ;
  seg->draw = gaDrawDouble ;
  seg->save = 0 ;
  seg->nn2 = 12 ;

  seg2 = gaNewSeg (look) ;
  seg2->buf = halloc (12, look->handle) ;
  seg2->title = "Ends   in genomic sequence" ;
  seg2->title2 = "around base" ;
  seg2->nn = 12;
  seg2->buf2 = halloc (12, look->handle) ;
  seg2->draw = gaDrawDouble ;
  seg2->save = 0 ;
  seg2->nn2 = 12 ;

  seg->down = seg2 ;

  tag = str2tag ("Genomic_sequence") ;
  if ((obj = bsCreate (look->tg)))
    {
      if (bsGetKey (obj, tag, &cosmid) &&
	  bsGetData (obj, str2tag("Covers"), _Int, &len) &&
	  bsGetData (obj, _bsRight, _Text, 0) &&  /* jump from */
	  bsGetData (obj, _bsRight, _Int, &a1)  &&
	  bsGetData (obj, _bsRight, _Int, &a2))
	ok = TRUE ;
      bsDestroy (obj) ;
    }

  if (ok)
    {
      seg->nn = strlen (name (cosmid)) ;
      seg->buf = strnew ( name(cosmid), look->handle) ;
      strcpy (seg->buf2, messprintf("%d", a1)) ;

      seg2->nn = strlen (name (cosmid)) + 2 ; 
      seg2->buf = strnew ( name(cosmid), look->handle) ;
      strcpy (seg2->buf2, messprintf("%d", a2)) ;
    }

  return seg ; /* return seg2 because this si vetical column */
}

/*********************************************************/

static SEG gaConvertGmap (LOOK look)
{
  SEG seg ;
  KEY tag, cosmid ;
  char *cp, *cq, *cr, cc ;
  OBJ obj = 0 ;

  seg = gaNewSeg (look) ;
  seg->buf = halloc (4, look->handle) ;
  seg->title = "Chromosome" ;
  seg->title2 = "Position" ;
  seg->nn = 4 ;
  seg->buf2 = halloc (12, look->handle) ;
  seg->draw = gaDrawDouble ;
  seg->save = 0 ;
  seg->nn2 = 12 ;
  tag = str2tag ("Genomic_sequence") ;
  cosmid = keyGetKey (look->tg, tag) ;
  cq = 0 ;
  if (cosmid && (obj = bsCreate (cosmid)))
    {
      if (bsGetData (obj, str2tag("Full_name"), _Text, &cp))
	cq = strnew(cp, 0) ;
      bsDestroy (obj) ;
    }
  if (cq)
    {
      cp = cq ;  
      while (*cp && *cp != '(') cp++ ;
      if (*cp == '(')
	{
	  cr = cp + 1 ;
	  while (*cp && (*cp != ':' && *cp != ';')) cp++ ;
	  cc = *cp ;
	  *cp = 0 ;
	  strcpy (seg->buf, cr) ;
	  if (cc == ':' || cc == ';')
	    {
	      cp++ ; cr = cp ;
	      while (*cp && *cp != ')') cp++ ;
	      cc = *cp ; *cp = 0 ;
	      if (cc == ')') 
		strcpy (seg->buf2, cr) ;
	    }
	}
    }
  
  return seg ;
}

/*********************************************************/

static SEG gaConvertPG (LOOK look)
{
  SEG seg0, segb, seg ;
  KEYSET ks = 0 ;
  KEY gf = keyGetKey (look->tg, str2tag("Matching_genefinder_gene")) ;

  seg0  = gaNewSeg (look) ;
  seg0->title = "Observed transcript" ;
  seg0->open = TRUE ;
  seg0->draw = gaDrawTag ;
  segb= 0 ;
  
  seg = gaConvertGenomic (look) ;
  HOOK (seg0, segb, seg) ;
  
  ks = queryKey (look->tg, ">Annotations Compatible_prediction") ;
  if (keySetMax(ks))
      {
	seg  = gaNewSeg (look) ;
	seg->tag = str2tag ("Identical_prediction") ;
	seg->title = strnew (messprintf("Compatible to predicted gene %s", 
					gf ? name(gf) : ""),
			     look->handle) ;
	seg->draw = gaDrawTag ;
	seg->save = 0 ;
	HOOK (seg0, segb, seg) ;
      }
  else if (keyGetKey (look->tg, str2tag ("Compatible_prediction")))
      {
	seg  = gaNewSeg (look) ;
	seg->tag = str2tag ("Compatible_prediction") ;
	seg->title = strnew (messprintf("Compatible to predicted gene %s", 
					gf ? name(gf) : ""),
			     look->handle) ;
	seg->nn = 20 ;
	seg->buf = halloc (seg->nn, look->handle) ;
	seg->draw = gaDrawTag ;
	seg->save = 0 ;
	HOOK (seg0, segb, seg) ;
      }
  keySetDestroy (ks) ;
  return seg0 ;
}

/*********************************************************/

static SEG gaConvertClones (LOOK look)
{
  SEG seg0, segb, seg, segc ;
  KEY clone = 0 ;
  OBJ obj = 0 ;
  KEYSET ks = queryKey (look->tg, ">cdNA_clone"), ks2 = 0, ks3 = 0, ks4 = 0 ;

  seg0 = gaNewSeg (look) ;
  seg0->title = strnew(messprintf("Constructed from %d cDNA clone%s", 
		  keySetMax(ks), (keySetMax(ks) ? "s" : "")), look->handle) ;
  seg0->open = TRUE ;
  seg0->draw = gaDrawTag ;
  seg0->save = 0 ;
  segb = 0 ;

  ks2 = query (ks, "IS yk*") ;
  seg = gaNewSeg (look) ;
  seg->open = TRUE ;
  seg->title =  strnew(messprintf("of which %d from Y.Kohara libraries,",
				  keySetMax(ks2)), look->handle) ;
  seg->draw = gaDrawTag ;
  HOOK (seg0, segb, seg) ;

  ks3 = query (ks2, "Library != E") ;
  seg = gaNewSeg (look) ;
  seg->open = TRUE ;
  seg->title =  strnew(messprintf("    %d from mixed stage libraries",
				  keySetMax(ks3)), look->handle) ;
  seg->draw = gaDrawTag ;
  HOOK (seg0, segb, seg) ;

  ks4 = query (ks2, "Library = E") ;
  seg = gaNewSeg (look) ;
  seg->open = TRUE ;
  seg->title =  strnew(messprintf("    %d from embryonic libraries",
				  keySetMax(ks4)), look->handle) ;
  seg->draw = gaDrawTag ;
  HOOK (seg0, segb, seg) ;

  obj = bsCreate (look->tg) ;
  if (!obj) return seg0 ;

  segc = seg ; segb = 0 ;
  if (bsGetKey (obj, _cDNA_clone, &clone))
    do
      {
	seg = gaNewSeg (look) ;
	seg->title = "" ;
	seg->nn = strlen(name(clone)) ;
	seg->buf = strnew (name(clone), look->handle) ;			 
	seg->draw = gaDrawTag ;
	seg->save = 0 ;
	HOOK (segc, segb, seg) ;
      } while (bsGetKey (obj, _bsDown, &clone)) ;
  bsDestroy (obj) ;
  segc->open = FALSE ;
  keySetDestroy (ks) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
  keySetDestroy (ks4) ;
  return seg0 ;
}

/*********************************************************/

#ifdef JUNK
static SEG gaConvertLongestClone (LOOK look)
{
  SEG seg ;
  int len ;
  KEY clone = 0 ;
  OBJ obj = 0 ;

  seg = gaNewSeg (look) ;
  seg->title = "Longest length covered on genomic sequence:" ;
  seg->nn2 = 6 ;
  seg->buf2 = halloc (7, look->handle) ;

  obj = bsCreate (look->tg) ;
  if (!obj) return seg ;
  if (bsGetKey (obj, str2tag("Longest_cDNA_clone"), &clone))
    {
      seg->nn = strlen(name(clone)) ;
      seg->buf = strnew(name(clone),look->handle) ;
      if (bsGetData (obj, _bsRight, _Int, &len))
	strcpy (seg->buf2, messprintf("%d bp", len)) ;
      seg->draw = gaDrawDouble ;
    }
  return seg ;
}
#endif

/*********************************************************/

static SEG gaConvertSplicing (LOOK look)
{
  SEG seg0, segb, seg ;
  OBJ obj = 0 ;
  int x1, x2 ;
  char *cp, *cq ;
  KEY tag ;
  BSMARK mark = 0 ;
  
  seg0 = gaNewSeg (look) ;
  seg0->title = "Splicing" ;
  seg0->draw = gaDrawTag ;
  seg0->save = 0 ;  
  seg0->open = FALSE ;
  
  segb = 0 ;
  if ((obj = bsCreate (look->tg)))
    if (bsGetData (obj, _Splicing, _Int, &x1))
      do
	{
	  tag = 0 ; cp = 0 ;
	  mark = bsHandleMark (obj, mark, 0) ;  
	  if (bsGetData (obj, _bsRight, _Int, &x2))
	    {
	      if (bsPushObj (obj) && bsGetKeyTags (obj, _bsRight, &tag))
		bsGetData (obj, _bsRight, _Text, &cp) ;
	      
	      seg = gaNewSeg (look) ;
	      seg->title = "" ;
	      seg->draw = gaDrawTag ;
	      seg->save = 0 ;
	      seg->dx = 0 ;
	      cq = messprintf("%5d %5d %15s %s ", 
			      x1, x2, (tag ? name(tag) : ""), (cp ? cp : "")) ;
	      seg->nn = strlen(cq) ;
	      seg->buf = strnew (cq, look->handle) ;
	      HOOK (seg0, segb, seg) ;
	    }
	  bsGoto (obj, mark) ;
	} while (bsGetData (obj, _bsDown, _Int, &x1)) ;
  bsMarkFree (mark) ;
  bsDestroy (obj) ;
  return seg0 ;
}

/*********************************************************/

static char *gaConvertExonUsage (LOOK look, KEY seq)
{
  char *cp = 0 ;
  Stack s = stackCreate (80) ;
  int ii ;
  OBJ obj = bsCreate (seq) ;

  pushText (s, "Using exons ") ;
  if (obj)
    {
      if (bsGetKey (obj, str2tag("Derived_from_gene"), 0) &&
	  bsGetData (obj, _bsRight, _Text, 0) &&
	  bsGetData (obj, _bsRight, _Int, &ii))
	{
	  do 
	    {
	      catText(s, messprintf("%d ",ii)) ;  
	    } while (bsGetData (obj, _bsRight, _Int, &ii)) ;
	}
      bsDestroy (obj) ;
    }
  cp = strnew (stackText(s,0), look->handle) ;
  stackDestroy (s) ; 
  return cp ;
}

/*********************************************************/

static SEG gaConvertDnaMessager (LOOK look, KEY seq)
{ 
  SEG seg0, segb, seg ;
  char *cp = 0 ;
  Array dna = 0 ;
  KEY seq2 = keyGetKey (seq, str2tag("Spliced_sequence")) ;
  seg0  = gaNewSeg (look) ;
  seg0->title = "RNA" ;
  seg0->open = FALSE ;
  seg0->draw = 0 ;

  segb = 0 ;

  dna = dnaGet (seq2) ;
  if (dna)
    {
      seg0->title =  strnew (messprintf("Messenger RNA %s %d base pairs", name(seq), arrayMax(dna)),
			     look->handle) ; 
      seg0->title2 = gaConvertExonUsage (look, seq) ;
      rnaDecodeArray (dna) ;

      cp = arrp (dna, 0, char) ;
      seg = gaNewSeg (look) ;
      seg->buf = halloc (arrayMax(dna) + 10, look->handle) ;
      seg->draw = gaDrawArray ;
      seg->nn =  arrayMax(dna) + 10 ;
      seg->fixedFont = TRUE ; 
      strcpy (seg->buf, cp) ;
      HOOK (seg0, segb, seg) ;
      
      arrayDestroy (dna) ;
      seg0->draw = gaDrawDouble ;
    } 
  else
    seg0->title =  "Messenger RNA not determined" ;

  return seg0 ;
}

/*********************************************************/

static SEG gaConvert5pUTR (LOOK look, KEY gene, KEY seq)
{ 
  SEG seg0, segb, seg ;
  char *cp = 0 ;
  Array dna = 0 ;
  int ii, x1, x2, z1, z2, nn ;
  KEY  tt ;
  KEY seq2 = keyGetKey (seq, str2tag("Spliced_sequence")) ;
  OBJ obj = bsCreate (seq) ;
  Array aa = arrayCreate (120, BSunit) ;
  BSunit *uu ;

  seg0  = gaNewSeg (look) ;
  seg0->title = "5pUTR" ;
  seg0->open = FALSE ;
  seg0->draw = gaDrawTag ;

  segb = 0 ;

  nn = 0 ;
  x1 = x2 = z1 = z2 = -1 ;
  if (obj &&
      bsGetArray (obj, _Splicing, aa, 3))
    {
      for (ii = 0 ; ii < arrayMax(aa) ; ii += 3)
	{
	  uu = arrp(aa, ii, BSunit) ;
	  x1 = uu[0].i ;
	  x2 = uu[1].i ;
	  tt = uu[2].k ;
	  if (tt == str2tag("Gap") || tt == str2tag("Intron"))
	    nn += x2 - x1 + 1 ;
	  if (tt == str2tag("5prime_UTR"))
	    {
	      if (z1 == -1) { z1 = x1 - nn ; z2 = x2 - nn ; }
	      else z2 = x2 - nn ;
	    }
	}
      if (z1 > -1)
	{ 
	  BOOL isSL1 = FALSE ; 
	  KEY trs = keyGetKey (gene, _Transpliced_to) ;
	  if (trs && (!strcmp(name(trs), "SL1") || !strcmp(name(trs), "SL2")))
	    isSL1 = TRUE ;
	  
	  dna = dnaGet (seq2) ;
	  if (dna && arrayMax(dna) >= z2)
	    { 
	      if (z1 >= 1 && (isSL1 || z2 > 3))
		{
		  rnaDecodeArray (dna) ;
		  array(dna, z2, char) = 0 ;
		  seg0->title =  strnew (messprintf("5 prime UTR region %d base pairs ", z2 - z1 + 1),
					 look->handle) ; 
		  cp = arrp (dna, z1 - 1, char) ;
		  seg = gaNewSeg (look) ;
		  seg->buf = halloc (z2 + 10, look->handle) ;
		  seg->draw = gaDrawArray ;
		  seg->nn =  z2 + 10 ;
		  seg->fixedFont = TRUE ; 
		  strcpy (seg->buf, cp) ;
		  HOOK (seg0, segb, seg) ; 
		  seg0->draw = gaDrawDouble ;
		}
	      else
		seg0->title =  "5 prime UTR seems < 3 bp" ;
	    }
	}
    }
  if (z1 == -1)
    seg0->title =  "5 prime end probably not reached" ;

  arrayDestroy (dna) ;
  arrayDestroy (aa) ;
  bsDestroy (obj) ;
  return seg0 ;
}

/*********************************************************/

static SEG gaConvert3pUTR (LOOK look, KEY seq)
{
  SEG seg0, segb, seg ;
  char *cp = 0 ;
  Array dna = 0 ;
  int ii, x1=0, x2=0, z1, z2, nn ;
  KEY tt=0 ;
  KEY seq2 = keyGetKey (seq, str2tag("Spliced_sequence")) ;
  OBJ obj = bsCreate (seq) ;
  Array aa = arrayCreate (120, BSunit) ;
  BSunit *uu ;

  seg0  = gaNewSeg (look) ;
  seg0->title = "3pUTR" ;
  seg0->open = FALSE ;
  seg0->draw = gaDrawTag ;

  segb = 0 ;

  nn = 0 ;
  x1 = x2 = z1 = z2 = -1 ; 
  if (obj &&
      bsGetArray (obj, _Splicing, aa, 3))
    {
      for (ii = 0 ; ii < arrayMax(aa) ; ii += 3)
	{
	  uu = arrp(aa, ii, BSunit) ;
	  x1 = uu[0].i ;
	  x2 = uu[1].i ;
	  tt = uu[2].k ;
	  if (tt == str2tag("Gap") || tt == str2tag("Intron"))
	    nn += x2 - x1 + 1 ;
	  
	  if (tt == str2tag("3prime_UTR"))
	    {
	      if (z1 == -1) { z1 = x1 - nn ; z2 = x2 - nn ; }
	      else z2 = x2 - nn ;
	    }
	}
      if (z1 > -1)
	{ 
	  dna = dnaGet (seq2) ;
	  if (dna && arrayMax(dna) >= z2)
	    { 
	      if (z2 - z1  > 1)
		{
		  rnaDecodeArray (dna) ;
		  array(dna, z2, char) = 0 ;
		  seg0->title =  strnew (messprintf("3 prime UTR region %d base pairs ", z2 - z1 + 1),
					 look->handle) ; 
		  cp = arrp (dna, z1 - 1, char) ;
		  seg = gaNewSeg (look) ;
		  seg->buf = halloc (z2 + 10, look->handle) ;
		  seg->draw = gaDrawArray ;
		  seg->nn =  z2 + 10 ;
		  seg->fixedFont = TRUE ; 
		  strcpy (seg->buf, cp) ;
		  HOOK (seg0, segb, seg) ; 
		  seg0->draw = gaDrawDouble ;
		}
	      else
		seg0->title =  "3 prime UTR seems < 3 bp" ;
	    }
	}
    }
  if (z1 == -1)
    seg0->title =  "3 prime end probably not reached" ;

  arrayDestroy (dna) ;
  arrayDestroy (aa) ;
  bsDestroy (obj) ;
  return seg0 ;
}

/*********************************************************/
static Array myaa = 0 ;
static int  mycol ;
static int gaCloneLengthOrder (const void *va, const void *vb)
{
  int ia = *(const int*) va, ib = *(const int*)vb ;
  return -arr(myaa,13*ia+mycol,BSunit).i + arr(myaa,13*ib+mycol,BSunit).i ;
}

static SEG gaConvertCloneLengths (LOOK look, KEY seq)
{
  SEG seg0, segb, seg ;
  OBJ obj = 0 ;
  int i = 0, nn, iaa ;
  Array aa = arrayCreate (700, BSunit), zz = arrayCreate (64, int) ;
  BSunit *up ;
  Stack s = stackCreate (60) ;
  char *blank = "    " ;
  float xf = 0 ;

  seg0 = gaNewSeg (look) ;
  seg0->title = "Clone Lengths" ; 
  seg0->open = TRUE ;
  seg0->draw = gaDrawTag ;
  seg0->save = 0 ;  

  segb = 0 ;
  /* seg0->right = segb = gaConvertLongestClone (look) ; */

  seg = gaNewSeg (look) ;
  if (!i)
    seg->title = "Calculated lengths by alignment to the mRNA" ;
  seg->draw = gaDrawTag ;
  HOOK (seg0, segb, seg) ;

  if (seq && (obj = bsCreate (seq)) &&
      bsGetArray (obj, str2tag("Assembled_from_cDNA_clone"), aa, 13))
    {
      i = arrayMax(aa)/13 ;
      while (i--) array(zz,i,int) = i ;
      myaa = aa ; mycol = 10 ;
      arraySort (zz, gaCloneLengthOrder) ;
      for (i = 0 ; i < arrayMax(zz) ; i++)
	{
	  iaa = 13 * arr(zz,i,int) ;
	  up = arrp (aa, iaa, BSunit) ;
	  
	  seg = gaNewSeg (look) ;
	  seg->draw = gaDrawTag ;
	  seg->save = 0 ;
	  seg->dx = 0 ;
	  stackClear (s) ;
	  pushText (s, name(up[0].k)) ;
	  nn = strlen (stackText (s, 0)) % 4 ;
	  if (nn) catText (s, blank + nn) ;

	  if (up[10].i > 0 && up[12].f > 0)
	    {
	      int im = 1000 *  up[12].f, ic = up[10].i , pp = im - ic ;
	      char *signe ;
	      if (pp > 0) { signe = "+" ; }
	      else  { pp = - pp ;  signe = "-" ; }
	      pp = pp * 100 ; pp /= ic ;
	      catText (s, messprintf (" calculated %d bp measured %g kb %s%d%%  ",
				      up[10].i, up[12].f, signe, pp)) ;
	    }
	  else if (up[6].i > 0)
	    {
	      xf = up[6].i ; xf = xf/1000.0 ;
	      catText (s, messprintf (" Not calculated,  measured %g kb   ",
				    xf)) ;
	    }
	  seg->nn = strlen (stackText (s, 0)) ;
	  seg->buf = strnew (stackText (s, 0), look->handle) ;
	  
	  HOOK (seg0, segb, seg) ;
	}
    }
  bsDestroy (obj) ;
  arrayDestroy (aa) ;
  arrayDestroy (zz) ;
  stackDestroy (s) ;
  return seg0 ;
}

/*********************************************************/

static SEG gaConvertDnaCodant (LOOK look, KEY seq)
{ 
  SEG seg0, seg, segb = 0 ;
  Array dna = 0 ;
  int ii, nn, z1, z2, x1, x2 ;
  KEY tt, seq2 = keyGetKey (seq, str2tag("Spliced_sequence")) ;
  char *cp ;
  OBJ obj = bsCreate (seq) ;
  Array aa = arrayCreate (120, BSunit) ;
  BSunit *uu ;

  seg0  = gaNewSeg (look) ;
  seg0->open = FALSE ;
  seg0->draw = gaDrawTag ;
  seg0->save = 0 ;  

  nn = 0 ;
  x1 = x2 = z1 = z2 = -1 ; 
  if (obj &&
      bsGetArray (obj, _Splicing, aa, 3))
    {
      for (ii = 0 ; ii < arrayMax(aa) ; ii += 3)
	{
	  uu = arrp(aa, ii, BSunit) ;
	  x1 = uu[0].i ;
	  x2 = uu[1].i ;
	  tt = uu[2].k ;
	  if (tt == str2tag("Gap") || tt == str2tag("Intron"))
	    nn += x2 - x1 + 1 ;
	  
	  if (tt == str2tag("Exon"))
	    {
	      if (z1 == -1) { z1 = x1 - nn ; z2 = x2 - nn ; }
	      else z2 = x2 - nn ;
	    }
	}
      if (z1 > -1)
	{ 
	  dna = dnaGet (seq2) ;
	  if (dna && arrayMax(dna) >= z2)
	    { 
	      if (z2 - z1  > 1)
		{
		  rnaDecodeArray (dna) ;
		  array(dna, z2, char) = 0 ;
		  seg0->title =  strnew (messprintf("Coding region %d base pairs ", z2 - z1 + 1),
					 look->handle) ; 
		  cp = arrp (dna, z1 - 1, char) ;
		  seg = gaNewSeg (look) ;
		  seg->buf = halloc (z2 + 10, look->handle) ;
		  seg->draw = gaDrawArray ;
		  seg->nn =  z2 + 10 ;
		  seg->fixedFont = TRUE ; 
		  strcpy (seg->buf, cp) ;
		  HOOK (seg0, segb, seg) ; 
		  seg0->draw = gaDrawDouble ;
		}
	      else
		seg0->title =  "3 prime UTR seems < 3 bp" ;
	    }
	}
    }
 
  if (z1 == -1)
    seg0->title =  "Coding Region not determined" ;

  
  arrayDestroy (dna) ;
  arrayDestroy (aa) ;
  bsDestroy (obj) ;
  return seg0 ;
}

/*********************************************************/

static SEG gaConvertTrans (LOOK look, KEY seq)
{ 
  SEG seg0, segb, seg ;
  int j = 0 ;
  char *cp, *cq, *cr ;
  Array dna = 0 ;
  KEY seq2 = keyGetKey (seq, str2tag("Spliced_sequence")) ;
  KEY subseq = seq2 ? keyGetKey (seq2, str2tag("Subsequence")) : 0 ;
  int x1, x2, z1, z2, ii ;
  Array aa = arrayCreate (120, BSunit) ;
  BSunit* uu ;
  OBJ obj = 0 ;

  seg0  = gaNewSeg (look) ;
  seg0->title = "Protein    length %d aa (including the stop)" ;
  seg0->open = FALSE ;
  seg0->draw = 0 ;

  segb = 0 ;
  dna = seq2 ? dnaGet (seq2) : 0 ;
  z1 = z2 = -1 ;
  if (dna)
    {
      if ((obj = bsCreate (subseq)) &&
	  bsGetArray (obj, str2tag("Source_Exons"), aa, 2))
	{
	  for (ii = 0 ; ii < arrayMax(aa) ; ii += 2)
	    {
	      uu = arrp(aa, ii, BSunit) ;
	      x1 = uu[0].i ;
	      x2 = uu[1].i ;
	      if (z1 == -1) { z1 = x1  ; z2 = x2  ; }
	      else z2 = x2 ;
	    }
	}
      bsDestroy (obj) ;

      cr = halloc(arrayMax(dna)/3 + 2, look->handle) ;
      cp = cr ;
      cq = arrp (dna, z1 - 1, char) ;
      for (j = z1 - 1; j < z2 ; j +=3, cq +=3)
	*cp++ = e_codon (cq, 0) ; 
      if (*(cp - 1) == '*')
	*(--cp) = 0 ; 
      seg0->title = strnew(messprintf ("Protein    length %d aa (%c to %c included)",
				cp - cr , *cr, *(cp-1)), 0);

      seg = gaNewSeg (look) ;
      seg->buf = halloc (arrayMax(dna)/3 + 10, look->handle) ;
      seg->draw = gaDrawArray ;
      seg->fixedFont = TRUE ; 
      seg->nn = arrayMax(dna)/3 + 10 ;
      strcpy (seg->buf, cr) ;
      HOOK (seg0, segb, seg) ;
      
      arrayDestroy (dna) ;
      seg0->draw = gaDrawTag ;
    }
  arrayDestroy (aa) ;
  return seg0 ;
}

/*********************************************************/
/*********************************************************/


static BOOL gaConvert (LOOK look)
{
  SEG  seg ;
  int i ;
  KEYSET ks = queryKey (look->tg, ">Transcript") ;

  seg = look->seg = gaNewSeg (look) ;
  seg->draw = gaDrawHeader ;
  seg->save = 0 ;

  seg->down = gaConvertGmap (look) ;
  seg = seg->down ;
 
  seg->down = gaConvertPG (look) ;
  seg = seg->down ;
 
  seg->down = gaConvertClones (look) ;
  seg = seg->down ;

  seg->down = gaConvertSplicing (look) ;
  seg = seg->down ;

  for (i = 0 ; i < keySetMax(ks) ; i++)
    {
      KEY messenger = keySet (ks, i) ;

      seg->down = gaConvertDnaMessager (look, messenger) ;
      seg = seg->down ;
      seg->down = gaConvert5pUTR (look, look->tg, messenger) ;
      seg = seg->down ;
      seg->down = gaConvert3pUTR (look, messenger) ;
      seg = seg->down ;
      
      seg->down = gaConvertCloneLengths (look, messenger) ;
      seg = seg->down ;
      
      seg->down = gaConvertDnaCodant (look, messenger) ;
      seg = seg->down ;
      
      seg->down = gaConvertTrans (look, messenger) ;
      seg = seg->down ;
    }
  keySetDestroy (ks) ;
  return TRUE ;
}

/*********************************************************/
/*********************************************************/
