/*  File: glocdisp.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the genetic map
 * Exported functions:
       glocDisplay
 * HISTORY:

 * Created: July 15 2002 (mieg)
 *-------------------------------------------------------------------
 */

/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "bump.h"
#include "wac/ac.h"
#include "display.h"
#include "pick.h"
#include "session.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "wfiche/biolog.h"
#include "lex.h"
#include "a.h"
#include "graph.h"
#include "hseqdisp.h"
#include "query.h"
#include "cdna.h"

static AC_DB ac_db = 0 ;

typedef struct GLMapStruct {
  magic_t *magic;        /* == GLMap_MAGIC */
  AC_HANDLE h ;
  Graph graph ;

  float mag, yCurrent, leftMargin, offset, yChrom ; 
  int graphWidth, graphHeight ;
  int a1, max ;

  BOOL hideHeader ;

  int   activeBox ;
  int	messageBox ;
  char  messageText[128] ;

  Array segs ;      	/* array of SEG's from the obj */
  Array boxIndex ;    /* box -> int (index into look->segs) */
  int lastTrueSeg ;	/* arrayMax(look->segs) after Convert */
  float dh ;  /* unit of height while bumping */
  char cytoBuffer[60] ;
  void (*mapDraw) (void) ;
} GLMap ;

typedef struct GLocStruct {
  magic_t *magic;        /* == GLoc_MAGIC */
  AC_HANDLE h ;
  AC_DB db ;
  GLMap *map ;
  BOOL isSmall ;
  int type ; /* 1: gene, 2:mrna, 3:product */
  KEY intMap, key, gene, tg, mrna, product ;
  int a1, a2 ; /* int map coords of the outside object */
  Array menu ;
} *GLoc;


typedef struct gLocSegStruct
  { 
    KEY gene, map ;
    int col, bcol, boxNeeded ;
    int a1, a2, iy1, iy2 ;
    int quality, nClones, height ;
    BOOL isDown, cloud ;
    int gid ;
    BOOL biblio, interactions, conserved, regulation, disease ;
} SEG ;

static magic_t GLoc_MAGIC = "GLoc";
static magic_t GLMap_MAGIC = "GLMap";
static BOOL gLocConvert (GLoc look) ;

/************************************************************/
/************************************************************/
#define HMAP2GRAPH(_gLocMap,_x) (((_x) - (_gLocMap)->offset) * (_gLocMap)->mag + (_gLocMap)->leftMargin)

/************************************************************/

static void gLocShowSegs (Array segs)
{
  SEG *seg ;
  int ii ;
  
  if (arrayExists (segs))
    for (ii = 0 ; ii < arrayMax (segs) ; ii++)
      {
	seg = arrp (segs, ii, SEG) ;
	printf ("%4d %s %12s%12d%12d cl=%d\n"
		, ii
		, seg->isDown ? "+" : "-"
		, ac_key_name(seg->gene), seg->a1, seg->a2
		, seg->cloud
		) ;
      }
  return ;
}

/************************************************************/

static GLMap *gLocMapCurrent (char *caller)
{
  GLMap *map = 0 ;
  
  if (!graphAssFind (&GLMap_MAGIC,&map))
    messcrash("%s() could not find GLMap on graph", caller);
  if (!map)
    messcrash("%s() received NULL GLoc pointer", caller);
  if (map->magic != &GLMap_MAGIC)
    messcrash("%s() received non-magic GLMap pointer", caller);
  
  return map ;
} /* gLocMapCurrent */

/************************************************************/

static void gLocMapResize (void)
{
  GLMap *map = gLocMapCurrent ("hmapPick") ;
  
  pickRememberDisplaySize (DtGLOC) ;
  
  graphFitBounds (&map->graphWidth,&map->graphHeight) ;
  map->leftMargin = 1 ;
  map->mag = ( (map->graphWidth - map->leftMargin))/map->max ;
  
  map->mapDraw () ;
} /* gLocMapResize */

/************************************************************/

static void gLocMapKbd (int k)
{ 
  GLMap *map = gLocMapCurrent ("hmapPick") ;
  
  if (map)
    switch (k)
      {
      case LEFT_KEY :
	map->offset -= map->graphWidth/(4 * map->mag) ;
	break ;
      case RIGHT_KEY :
	map->offset += map->graphWidth/(4 * map->mag) ;
	break ;
      case UP_KEY :
	map->mag *= 1.414 ;
	break ;
      case DOWN_KEY :
	map->mag /= 1.414 ;
	break ;
      case HOME_KEY :
	gLocMapResize () ;
	return ;
      default:
	return ;
      }
  map->mapDraw () ;
}

/************************************************************/

static void gLocMapFollow (GLMap *map, int box)
{
  SEG *seg ;
  int iSeg ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg >= 0 && 
      (seg = arrp(map->segs, iSeg, SEG)) &&
      seg->gene)
    display (seg->gene, 0, TREE) ;
 
} /* gLocMapFollow */

/************************************************************/

static void gLocMapSelect (GLMap *map, int box)
{
  SEG *seg = 0 ;
  int iSeg ;

  if (map->activeBox > 0 &&
      map->activeBox < arrayMax(map->boxIndex) && 
      (iSeg = arr (map->boxIndex, map->activeBox, int)) &&
      iSeg >= 0 &&
      (seg = arrp(map->segs, iSeg, SEG)))
    {
      graphBoxDraw (map->activeBox, seg->col, seg->bcol) ;
    }
  
  map->activeBox = 0 ; seg = 0 ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg >= 0 &&
      (seg = arrp(map->segs, iSeg , SEG)))
    { 
      map->activeBox = box ;
      graphBoxDraw (map->activeBox, BLACK, CYAN) ;
    }

  if (map->activeBox && map->messageBox)
    { 
      strncpy (map->messageText, name (seg->gene), 100) ;
      strcat (map->messageText, messprintf(" %d  %d", seg->a1 - map->a1 + 1, seg->a2 - map->a1 + 1)) ;
      
      graphBoxDraw (map->messageBox, BLACK, PALEBLUE) ;
    }
} /* gLocMapSelect */

/*************************************/

static void gLocMapPick (int box,  double x , double y) 
{
  GLMap *map = gLocMapCurrent ("gLocMapPick") ;

  if (!box)
    return ;
  if (box < arrayMax (map->boxIndex) && arr(map->boxIndex,box,int)) /* a SEG */
    {
      if (box == map->activeBox) /* a second hit - follow it */
	gLocMapFollow (map, box) ;
      else
	gLocMapSelect (map, box) ;
    }

  return;
} /* gLocMapPick */

/************************************************************/

static GLMap *gLocMapCreate (void *look, AC_HANDLE handle, void (*mapDraw) (void))
{
  GLMap *map = (GLMap *) halloc (sizeof (struct GLMapStruct), handle) ;
  map->magic = &GLMap_MAGIC ;
  map->h = handle ; /* do not destroy, belongs to calling routine */
  map->boxIndex = arrayHandleCreate (64,int, handle) ;
  map->segs = arrayHandleCreate (64, SEG, handle) ;
  map->activeBox = 0 ;
  map->lastTrueSeg = arrayMax(map->segs) ;
  map->mapDraw = mapDraw ;

  return map ;
} /* gLocMapCreate */

/************************************************************/

static void gLocMapGraphAssociate (GLMap *map)
{
  graphRegister (KEYBOARD, gLocMapKbd) ;
  graphRegister (PICK, gLocMapPick) ;
  graphAssRemove (&GLMap_MAGIC) ;
  graphAssociate (&GLMap_MAGIC, map) ; /* attach map to graph */
  graphRegister (RESIZE, gLocMapResize) ;
  map->graph = graphActive() ;

  gLocMapResize () ; /* redraws */
  return ;
} /* gLocMapCreate */

/************************************************************/
/************************************************************/

static GLoc currentGLoc (char *caller)
{
  GLoc gLoc;

  if (!graphAssFind (&GLoc_MAGIC,&gLoc))
    messcrash("%s() could not find GLoc on graph", caller);
  if (!gLoc)
    messcrash("%s() received NULL GLoc pointer", caller);
  if (gLoc->magic != &GLoc_MAGIC)
    messcrash("%s() received non-magic GLoc pointer", caller);
  
  return gLoc;
} /* currentGLoc */


/************************************************************/
/***************** Registered routines *********************/

static void gLocDestroy (void)
{
  GLoc look = currentGLoc("gLocDestroy") ;

  if (look)
    {
      ac_free (look->h) ;
      
      graphAssRemove (&GLoc_MAGIC) ; /* detach look from graph */
      
      look->magic = 0 ;
      ac_free (look) ;
    }

  return;
} /* gLocDestroy */

/*************************************************************/
#ifdef JUNK
static void gLocRecalculate(void)
{
  GLoc look = currentGLoc("gLocRecalculate") ;
  
  if (look)
    {
      look->map->segs = arrayReCreate (look->map->segs,64, SEG) ;
      gLocConvert (look) ;
      look->map->mapDraw () ;
    } 
}
#endif
/*************************************************************/
/************************************************************/

static int  gLocSegOrder(const void *a, const void *b)  /* for arraySort() call */ 
{
  const SEG *sa = (const SEG*)a, *sb = (const SEG*) b ;
  int n ;

  n = sa->isDown - sb->isDown ;
  if (n) return -n ;
  n = sa->a1 - sb->a1 ;
  if (n) return n ;
  if (sa->gene > sb->gene)
    return 1 ;
  if (sa->gene < sb->gene)
    return -1 ;
  return 0 ;
}

/************************************************************/

static AC_OBJ  gLoc_key_obj (GLoc look, KEY key)
{
  return ac_get_obj (look->db, className (key), name(key), look->h) ;
}

/************************************************************/

static BOOL gLocConvertOneGene (GMP *gmp, GLoc look, AC_OBJ Gene)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE imap ;
  AC_KEYSET ks ;
  SEG *seg ;
  const char *qq ;
  float u1, u2 ;
  int a1, a2, ii, col, nSeg = arrayMax (look->map->segs) ;
  KEY gene = ac_obj_key (Gene) ;
  BOOL ok = FALSE ;
  if (0) printf ("gene %s \n", ac_name(Gene)) ;
  imap = ac_tag_table (Gene, "IntMap", h) ;
  if (!imap)
    goto done ;
  a1 = ac_table_int (imap, 0, 1, -1) ;
  a2 = ac_table_int (imap, 0, 2, -1) ;
  u1 = HMAP2GRAPH (look->map, a1) ;
  u2 = HMAP2GRAPH (look->map, a2) ;

  if (0)
    {
      if (u1 < 1) u1 = 1 ;
      if (u2 < 1) u2 = 1 ;
      if (u1 >= look->map->graphWidth)
	u1 = look->map->graphWidth - 1 ;
      if (u2 >= look->map->graphWidth)
	u2 = look->map->graphWidth - 1 ;
      if (u1 == u2)
	goto done ;
    }

  seg = arrayp (look->map->segs, nSeg++, SEG) ;
  
  seg->gene = gene ;
  seg->map = ac_table_key (imap, 0, 0, 0) ;
  a1 = ac_table_int (imap, 0, 1, -1) ;
  a2 = ac_table_int (imap, 0, 2, -1) ;
  if (a1 < a2) { seg->isDown = TRUE ;  seg->a1 = a1 ; seg->a2 = a2 ; }
  else { seg->isDown = FALSE ;  seg->a1 = a2 ; seg->a2 = a1 ; }
  seg->cloud = ac_has_tag (Gene, "Cloud_gene") ;

  /* < 3 clones ;  < 10 clones ; < 100  ->epaisseur de la fleche */
  seg->nClones = ac_keyset_count (ac_objquery_keyset (Gene, ">transcribed_gene ; >cDNA_clone", h)) ;
  if (seg->nClones > 230)
    seg->height = 5 ;
  else if (seg->nClones > 50)
    seg->height = 4 ;
  else if (seg->nClones > 8)
    seg->height = 3 ;
  else if (seg->nClones > 3)
    seg->height = 2 ;
  else 
    seg->height = 1 ;


  /* pastilles */
  seg->gid = ac_has_tag (Gene, "GeneId") ;
  if (seg->gid && strncmp (ac_key_name (seg->gene), "LOC", 3))
    seg->gid = 2 ;
  if (seg->gene == look->gene)
    seg->gid = 3 ;
  if (seg->gid == 1) seg->boxNeeded++ ;

  seg->conserved = ac_has_tag (Gene, " pastille_conserved") ;
  if (seg->conserved) seg->boxNeeded++ ;

  seg->biblio = ac_has_tag (Gene, "Reference") ;
  if (seg->biblio) seg->boxNeeded++ ;

  ks = ac_objquery_keyset (Gene, "{Interacts} SETOR {>GeneId ; Interacts}", h) ;
   if (ks && ac_keyset_count (ks) > 0)
     seg->interactions =  TRUE ;
  if (seg->interactions) seg->boxNeeded++ ;

  seg->regulation = ac_has_tag (Gene, "Pastille_regulation") ;
  if (seg->regulation) seg->boxNeeded++ ;

  seg->disease = ac_has_tag (Gene, "Pastille_disease") ;
  if (seg->disease) seg->boxNeeded++ ;

  for (ii = 4 ; ! seg->quality && ii > 0 ; ii--)
    {
      switch (ii)
	{
	case 4:
	  qq = ">Product ; Good_product" ; 
	  col = MAGENTA ;
	  break ;
	case 3:
	  qq = "balise && nIntrons" ;
	  col = RED ;
	  break ;
	case 2:
	  qq = "putative OR balise" ;
	  col = DARKBLUE ;
	  break ;
	case 1:
	  qq = 0 ;
	  col = BLACK ;
	  break ;
	  
	}
      ks = qq ? ac_objquery_keyset (Gene, qq, h) : 0 ;
      if (!qq || (ks && ac_keyset_count (ks) > 0))
	{ 
	  seg->quality = ii ; 
	  seg->col = seg->bcol = col ;
	}
    }
  ok = TRUE ;

 done:
  ac_free (h) ;
  return ok ;
} /* gLocConvertOneGene */

/************************************************************/

static BOOL  gLocGeneConvert (GLoc look)
{
  AC_HANDLE h = ac_new_handle () ;
  GMP *gmp = 0 ;
  BOOL ok = FALSE ;
  KEY gene = look->gene ;
  AC_OBJ Gene = 0, Gene2 = 0 ;
  int a1 = 0, a2 = 0 ;
  int amax, xmax, x0, x1, x2 ;
  const char* newName, *ccp ;
  AC_ITER iter ;

  Gene = gLoc_key_obj (look, gene) ;
  if (!Gene)
    goto abort ;
  gmp = gmpCreate (look->db, Gene, 0, 0, 0, 0, 0, 0) ;

  {
    AC_TABLE imap = ac_tag_table (Gene, "IntMap", h) ;
    if (!imap)
      goto abort ;

    a1 = ac_table_int (imap, 0, 1, -1) ;
    a2 = ac_table_int (imap, 0, 2, -1) ;
  }

  if (0)
    {
      ccp = ac_tag_printable (Gene, "Cytogenetic", 0) ;
      if (!ccp)
	ccp = ac_tag_printable (Gene, "IntMap", "") ;
      strncpy (look->map->cytoBuffer, ccp, 59) ;
    }

  newName = ac_tag_printable (Gene, "Newname", 0) ;
  if (!newName || strlen (newName) < 3)
    goto abort ;

  look->map->segs = arrayHandleCreate (600, SEG, look->h) ;
  amax = xmax = 1 ;
  /* assume max length da = 50 characters */
  xmax = 50 ;
  if (gmp->Spc == WORM)
    {
      const char *ccp ;
      char chrom, zone1, zone2, buf[64] ;
      int nn ;

      amax =  look->isSmall ? 50000 : 200000 ;
      strcpy (buf, "000000000") ;
      ccp = newName + 2 ;
      nn = strlen (ccp) ;
      memcpy (buf + strlen(buf) - nn, ccp, nn) ;
      x0 = atoi (buf) ;
      chrom = newName [0] ;
      zone1 =  newName [1] ;
      x1 = x0 - 1000 ; if (x1 < 0) { x1 += 1000 ; zone1-- ; }
      zone2 =  newName [1] ;
      x2 = x0 + 1000 ; if (x2 > 999) { x2 -= 1000 ; zone2++ ; }
      iter = ac_dbquery_iter (gmp->db
			      , messprintf ("Find newname IS > \"%c%c%d\" && IS < \"%c%c%d\" ; > gene"
					    , chrom, zone1, x1, chrom, zone2, x2)
			      , h) ;
    }
  else if (gmp->Spc == MOUSE || gmp->Spc == HUMAN || TRUE)
    {
      char *cp, buf[64], prefix [8] ;
      int nn ;

      amax = look->isSmall ? 200000 : 1000000 ;
      if (gmp->Spc == ARA) 
	amax =  look->isSmall ? 25000 : 200000 ;
      strncpy (prefix, newName, 7) ;
      cp = strstr (prefix, "_") ;
      if (cp) *cp = 0 ;

      strcpy (buf, "000000000") ;
      cp = strstr (newName, "_") ;
      if (!cp)
	goto abort ;
      cp++ ;
      nn = strlen (cp) ;
      if (!nn)
	goto abort ;
      memcpy (buf + strlen(buf) - nn, cp, nn) ;
      x0 = atoi (buf) ;
      x1 = x0 - 2*amax ; if (x1 < 0) x1 = 1 ;
      x2 = x0 + 2*amax ;
      iter = ac_dbquery_iter (gmp->db
			      , messprintf ("{Find newname IS > \"%s_%d\" && IS < \"%s_%d\" ; > gene } SETOR { Find gene IS \"%s\"}"
					    , prefix, x1, prefix, x2, ac_name(gmp->gene))
			      , h) ;
    }
  else
    goto abort ;

  look->map->max = amax ;
  look->map->offset = (a1 + a2)/2 - amax/2 ;
  look->map->a1 = a1 < a2 ? a1 : a2 ;
  while (ac_free (Gene2), Gene2 = ac_next_obj (iter))
    gLocConvertOneGene (gmp, look, Gene2) ;

  ok = TRUE ;
 abort:
  gmpDestroy (gmp) ;
  ac_free (h) ;
  return ok ;
} /* gLocGeneConvert */

/************************************************************/

static BOOL  gLocTgConvert (GLoc look)
{
  look->tg = look->key ;
  look->key = look->gene = keyGetKey (look->tg, _Gene) ;
  return gLocGeneConvert (look) ;
}  /* gLocTgConvert */

/************************************************************/

static BOOL  gLocMrnaConvert (GLoc look)
{
  look->mrna = look->key ;
  look->key = look->tg = keyGetKey (look->mrna, _From_gene) ;
  return gLocTgConvert (look) ;
} /* gLocMrnaConvert */

/************************************************************/

static BOOL gLocProductConvert (GLoc look)
{
  look->product = look->key ;
  look->key = look->mrna = keyGetKey (look->product, _mRNA) ;
  return gLocMrnaConvert (look) ;		  
}  /* gLocProductConvert */

/************************************************************/

static BOOL gLocConvert (GLoc look)
{
  BOOL ok = FALSE ;

  if (class(look->key) == _VGene)
    { look->gene = look->key ; look->type = 1 ; }
  else if (class(look->key) == _VTranscribed_gene)
    { look->mrna = look->key ; look->type = 2 ;}
  else if (class(look->key) == _VmRNA)
    { look->mrna = look->key ; look->type = 3 ;}
  else if (class(look->key) == _VmProduct)
    { look->product = look->key ; look->type = 4 ; }
  else 
    goto abort ;
  
  switch (look->type)
    {
    case 1:
      ok = gLocGeneConvert (look) ;
      break ;
    case 2:
      ok = gLocTgConvert (look) ;
      break ;
    case 3:
      ok = gLocMrnaConvert (look) ;
      break ;
    case 4:
      ok = gLocProductConvert (look) ;
      break ;
    default:
      break ;
    }
  arraySort (look->map->segs,  gLocSegOrder) ;
  arrayCompress (look->map->segs) ;
 abort:
  return ok ;
} /* gLocConvert */

/************************************************************/

static void gLocMenuSelect (GLoc look)
{
  int ii = 0 ;
  MENUOPT *mm ;
  
  arrayDestroy (look->menu) ;
  look->menu = arrayHandleCreate (30, MENUOPT, look->h) ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphDestroy ; mm->text = "Quit" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = help ; mm->text = "Help";

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphPrint ; mm->text = "Print" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = 0 ; mm->text = 0;


} /* gLocMenuSelect */

/************************************************************/

static void gLocBumpGene (GLoc look, BUMP bumper, BOOL boxNeeded, BOOL isDown, int *maxp)
{
  int iSeg, n, ln ;
  SEG *seg ;
  float u0, u1, u2, du ;

  for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG) ; iSeg  < arrayMax(look->map->segs) ; iSeg++, seg++) 
    {
      if (
	  (seg->boxNeeded && !boxNeeded) || 
	  (!seg->boxNeeded && boxNeeded) || 
	  (seg->isDown != isDown)
	  )
	continue ;
      
      u1 = HMAP2GRAPH (look->map, seg->a1) ;
      u2 = HMAP2GRAPH (look->map, seg->a2) ;
      u0 = (u1 + u2)/2.0 ;
      if (u2 > look->map->graphWidth - 1)
	u2 = look->map->graphWidth - 1 ;
      if (u1 < 1)
	u1 = 1 ;
      if (u1 >= u2)
	continue ;
      
      /* reserve a variable thickness */ 
      seg->iy1 = 0 ;
      n = 3 * boxNeeded +  seg->height ;
      n = 1 ; if (boxNeeded) n++ ;
      if (look->isSmall) 
	du = u2 - u1 ;
      else
	{
	  if (seg->gid >= 2 || boxNeeded)
	    du = 1 ;
	  else
	    du = .01 ;
	}
      bumpItem (bumper, n, du, &(seg->iy1), &u1) ;

      if (seg->gid >= 2)
	{
	  n++ ;
	  ln = strlen (ac_key_name(seg->gene)) ;
	  u1 = u0 - ln/2 - 1 ;
	  u2 = u0 + ln/2 + 1 ;
	  if (look->isSmall)
	    du = u2 - u1 ;
	  else
	    du = 1 ;
	  bumpItem (bumper, 1, du, &(seg->iy1), &u1) ;
	}
      if (boxNeeded) 
	{
	  n++ ;
	  ln = seg->boxNeeded ;
	  u1 = u0 - ln/2 - 1 ;
	  u2 = u0 + ln/2 + 1 ;
	  if (look->isSmall)
	    du = u2 - u1 ;
	  else
	    du = 1 ;
	  bumpItem (bumper, 1, du, &(seg->iy1), &u1) ;
	}
      if (seg->isDown)
	seg->iy1 += n ; /* top line of this gene */
      if (*maxp < seg->iy1)
	*maxp = seg->iy1 ;
    }
  if (!isDown) (*maxp)++ ;
} /* gLocBumpGene */

/************************************************************/

static BOOL gLocDrawPastilles (GLoc look, SEG *seg, float u0, float y)
{
  int box, oldFormat ;
  float old, dx = 1.1 ;

  u0 -= seg->boxNeeded * dx/2.0 ;
  if (u0 < 1 || u0 > look->map->graphWidth - 3)
    return FALSE ;
  old = graphTextHeight (FALSE && look->isSmall ? 1.2 : 1.0) ;
  oldFormat = graphTextFormat(BOLD) ;
  if (seg->gid == 1)
    {
      box = graphBoxStart () ;
      graphText ("G", u0, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, DARKBLUE, PALEBLUE) ;
      graphBoxSetPick (box, FALSE) ;
  
      u0 += dx ;
    }
  if (seg->disease)
    {
      box = graphBoxStart () ;
      graphText ("D", u0, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, CERISE, PALERED) ;
      graphBoxSetPick (box, FALSE) ;
      u0 += dx ;
    }
  if (seg->conserved)
    {
      box = graphBoxStart () ;
      graphText ("C", u0, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, PURPLE, YELLOW) ;
      graphBoxSetPick (box, FALSE) ;
      u0 += dx ;
    }
  if (seg->interactions)
    {
      box = graphBoxStart () ;
      graphText ("I", u0, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, DARKGREEN, PALEGREEN) ;
      graphBoxSetPick (box, FALSE) ;
      u0 += dx ;
    }
  if (seg->regulation)
    {
      box = graphBoxStart () ;
      graphText ("R", u0, y) ;
      graphBoxEnd () ;
      if (0) graphBoxDraw (box, MAGENTA, PALEMAGENTA) ; 
      graphBoxDraw (box, BLUE, PALECYAN) ;
      graphBoxSetPick (box, FALSE) ;
      u0 += dx ;
    }
  if (seg->biblio)
    {
      box = graphBoxStart () ;
      graphText ("P", u0, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, DARKVIOLET, PALEVIOLET) ;
      graphBoxSetPick (box, FALSE) ;
      u0 += dx ;
    }
  graphTextHeight (old) ;
  graphTextFormat (oldFormat) ;

  return TRUE ;
} /* gLocDrawPastilles */

/************************************************************/

static BOOL gLocDrawArrow (GLoc look, SEG *seg, float y)
{
  int type = 0 ;
  float *fp, uu, u1, u2, ua, uamin, ub, dy2 = 0, dy = seg->height * look->map->dh/10.0 ;
  Array pp = 0 ;

  y += look->map->dh ; /* tip of arrow */

  u1 = HMAP2GRAPH (look->map, seg->a1) ;
  u2 = HMAP2GRAPH (look->map, seg->a2) ; 
  if (u2 > look->map->graphWidth - 2)
    u2 = look->map->graphWidth - 2 ;
  if (u1 < 1)
    u1 = 1 ;
  if (u1 >= u2)
    return FALSE ;

  pp = arrayCreate (16, float) ;
  if (u2 - u1 < 1)
    {
      uu = (u1 + u2)/2 ;
      graphColor (TRANSPARENT) ;
      graphLine (uu - .5, y - .5 , uu + .5, y + .5) ;
    }
  ua = (3.0 * u2 + u1)/4.0 ; /* inflexion point of the arrow */
  ub = u2 ; /* tip of arrow */
  uamin = u2 - 3  ; /* impose a minimal angle */
  if (ua < uamin)
    ua = uamin ;
  uamin = u2 -  3 * dy ; /* impose a maximal angle */
  if (ua > uamin)
    { ua = uamin ; } /* clip the tail */
  if (ua == u1)
    { type = 1 ; ub = u2 ; } /* no tail left */
  else if (ua < u1) /* truncate the tip */
    { type = 2 ; ua = u1; ub = u2 ; dy2 = dy -  (ub - ua) ; }
  else
    ub = u1 + (u2 - ua) ;
  if (! seg->isDown) /* flip the arrow */
    { 
      uu = u1 ; u1 = u2 ; u2 = uu ;
      uu = ua ; ua = ub ; ub = uu ;
    }
  
  switch (type)
    {
    case 0: /* full fledged arrow */ 
      fp = arrayp (pp, 15, float) ; /* make room */
      fp = arrp (pp, 0, float) ;
      *fp++ = u1 ; *fp++ = y - dy ;
      *fp++ = ua ; *fp++ = y - dy ;
      *fp++ = ua ; *fp++ = y - 3 * dy ;
      *fp++ = u2 ; *fp++ = y ;
      *fp++ = ua ; *fp++ = y + 3 * dy ;
      *fp++ = ua ; *fp++ = y + dy ;
      *fp++ = u1 ; *fp++ = y + dy ;
      *fp++ = u1 ; *fp++ = y - dy ;
      break ;
    case 1: /* no tail */
      fp = arrayp (pp, 7, float) ; /* make room */
      fp = arrp (pp, 0, float) ;
      *fp++ = ua ; *fp++ = y - 1 * dy ;
      *fp++ = ub ; *fp++ = y ;
      *fp++ = ua ; *fp++ = y + 1 * dy ;
      *fp++ = ua ; *fp++ = y - 1 * dy ;
      break ;
    case 2: /* truncated arrow */
      dy2 = dy = .66 * dy ;
      fp = arrayp (pp, 9, float) ; /* make room */
      fp = arrp (pp, 0, float) ;
      *fp++ = ua ; *fp++ = y - 2 * dy ;
      *fp++ = ub ; *fp++ = y - 2 * dy2 ;
      *fp++ = ub ; *fp++ = y + 2 * dy2 ;
      *fp++ = ua ; *fp++ = y + 2 * dy ;
      *fp++ = ua ; *fp++ = y - 2 * dy ;
    }

  graphColor (seg->bcol) ;
  if (1) graphLineSegs (pp) ;
  else graphPolygon (pp) ;
  if (seg->col != seg->bcol)
    {
      graphColor (seg->bcol) ;
      graphLineSegs (pp) ;
    }
  graphColor (BLACK) ;
  
  arrayDestroy (pp) ;
  return TRUE ;
} /* gLocDrawArrow */

/************************************************************/

static float gLocDrawGene (GLoc look, BOOL boxNeeded, BOOL isDown)
{
  AC_HANDLE h = ac_new_handle () ;
  int box, iSeg ;
  KEYSET mrnas ;
  SEG *seg ;
  float u0, u1, u2, y, old, ymax = 0 ;
  const char *ccp ;
  vTXT buf = vtxtHandleCreate (h) ;

  for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG) ; iSeg  < arrayMax(look->map->segs) ; iSeg++, seg++) 
    {
      if (
	  (seg->boxNeeded && !boxNeeded) || 
	  (!seg->boxNeeded && boxNeeded) || 
	   (seg->isDown != isDown)
	  )
	continue ;
      
      u1 = HMAP2GRAPH (look->map, seg->a1) ;
      u2 = HMAP2GRAPH (look->map, seg->a2) ; 
      box = graphBoxStart() ;
      if (isDown)
	{
	  y = look->map->yChrom - look->map->dh * seg->iy1 ;
	  if (seg->gid >= 2)
	    {
	      ccp = messprintf ("%s b=%d", ac_key_name (seg->gene), seg->iy1) ;
	      ccp = ac_key_name (seg->gene) ;
	      u0 = (u1 + u2)/2.0 - strlen (ccp)/2.0 ;
	      if (u0 > 1 && u0 < look->map->graphWidth - 3)
		{
		  if (seg->gid == 3) 
		    { graphColor (RED) ;
		    old = graphTextHeight (1.3) ;
		    }
		  graphText (ccp, u0, y) ;
		  if (seg->gid == 3) 
		    {
		      graphColor (BLACK) ;
		      graphTextHeight (old) ;
		    }
		  y += 1.0 ;
		}
	    }

	  if (seg->boxNeeded)
	    {
	      u0 = (u1 + u2)/2.0 ;
	      if (gLocDrawPastilles (look, seg, u0, y))
		y += .8 ;
	    }

	  gLocDrawArrow (look, seg, y) ;
	  y += look->map->dh ;
	  if (ymax < y) ymax = y ;
	}
      else
	{
	  y = look->map->yChrom + look->map->dh * seg->iy1 ;
	  if (gLocDrawArrow (look, seg, y))
	    {
	      y += look->map->dh ;
	      
	      if (seg->boxNeeded)
		{ 
		  y += .4 ;
		  u0 = (u1 + u2)/2.0 ;
		  if (gLocDrawPastilles (look, seg, u0, y))
		    y += .8 ;
		}
	      if (seg->gid >= 2)
		{
		  ccp = ac_key_name (seg->gene) ;
		  y += .2 ;
		  u0 = (u1 + u2)/2.0 - strlen (ccp)/2.0 ;
		  if (u0 > 1 && u0 < look->map->graphWidth - 3)
		    {
		      old = 0 ;
		      if (seg->gid == 3) graphColor (RED) ;
		      else old = graphTextHeight (1.0) ;
		      graphText (ccp, u0, y) ;
		      if (seg->gid == 3) graphColor (BLACK) ;
		      else graphTextHeight (old) ; 
		      y += .8 ;
		    }
		}
	    }
	  if (ymax < y) ymax = y ;
	}
      graphBoxEnd () ;
      array (look->map->boxIndex, box, int) = iSeg ;
      vtxtClear (buf) ;
      vtxtPrintf (buf, "%s, %d accession%s"
		  , ac_key_name(seg->gene)
		  , seg->nClones
		  , seg->nClones > 1 ? "s" : ""
		  ) ;
     
      mrnas = queryKey (seg->gene, ">transcribed_gene; > mrna") ;
      if (mrnas && keySetMax (mrnas) > 1)
	vtxtPrintf (buf, ", %d variants", keySetMax (mrnas)) ;
      keySetDestroy (mrnas) ;
      
#ifdef JUNK
      if (0)  /* trop long de calculer les tissus */
	{
	  AC_KEYSET clones = ac_objquery_keyset (Gene, ">transcribed_gene ; >cdna_clone", h) ;
	  AC_KEYSET reads = ac_ksquery_keyset (clones, ">read ; tissue", h) ;
	  AC_TABLE cTbl = clones ? ac_keyset_table (reads, 0, -1, 0, h) : 0  ;
	  
	  ficheNewGeneExpressionTissue (buf, 0, cTbl, 3, 3, 0, ac_keyset_count (clones)) ;
	}
#endif

      graphBubbleInfo (box, seg->gene, "Gene", vtxtPtr (buf)) ;
    }
  ac_free (h) ;
  return ymax ;
}

/************************************************************/

static void gLocDrawChromosome (GLoc look)
{
  float x, y ;

  x = look->map->graphWidth - 2 ;
  if (x > 12)
    {
      y = look->map->yChrom = look->map->yChrom + .51 ;
      graphLine (3, y, x, y) ;
      y = look->map->yChrom = look->map->yChrom + .11 ;
    }

} /* gLocDrawChromosome */

/************************************************************/

static void gLocDrawScale (GLoc look)
{
  float u1, u2, u3, du, y1, y2 ;
  int sBox, i, sc = utMainRoundPart (look->map->max/10) ;
  int sc1 =  utMainPart (sc/5), sc3 ;
  
  sBox = graphBoxStart () ;
  u1 = .7 ;
  du = HMAP2GRAPH (look->map, sc) - HMAP2GRAPH (look->map, 0) ;
  u2 = u1 + du ;
  y2 = look->map->graphHeight - .07 ;
  graphLine (u1, y2, u2, y2) ;
  y1 = y2 - .4 ;
  graphLine (u1, y1, u1, y2) ; graphLine (u2, y1, u2, y2) ;
  y1 = y2 - .2 ;
  du = HMAP2GRAPH (look->map, sc1) - HMAP2GRAPH (look->map, 0) ;
  for (i = 1 ; i * sc1 < sc ; i++)
    { 
      u3 = u1 + i * du ;
      if (2*i*sc1 == sc)
	{
	  graphLine (u3, y2 - .4, u3, y2) ;
	  sc3 = i*sc1 ;
	  while (sc3 >= 1000)
	    sc3 /= 1000 ;
	  graphText (messprintf ("%d", sc3), u3 - .3 , y1 - 1.2) ;
	}
      else
	graphLine (u3, y1, u3, y2) ;
    }
  if (sc >= 1000000)
    graphText (messprintf ("%dMb", sc/1000000), u2 - .5 , y1 - 1.2) ;
  else if (sc >= 1000)
    graphText (messprintf ("%dkb", sc/1000), u2 - .5 , y1 - 1.2) ;
  else
    graphText (messprintf ("%dbp", sc), u2 - .5 , y1 - 1.2) ;

  graphText ("0",  u1 - .4, y1 - 1.2) ;
  graphText (" ",  1, y2 - .6) ;
  graphBoxEnd () ;
  graphBoxDraw (sBox, BLACK, WHITE) ;
}

/************************************************************/

static void gLocDraw (void)
{
  GLoc look = currentGLoc("gLocDraw") ;
  BUMP bumper = 0 ;
  int n1, n2, n3, n4 ;
  float ymax = 0 ;
  look->map->dh = .7 ; /* unit of height */
  gLocMenuSelect (look) ;

  graphClear () ;

  graphMenu (arrp(look->menu, 0, MENUOPT)) ;
  look->map->messageBox = graphBoxStart () ;
  graphBoxEnd () ;

  look->map->yCurrent = 1 ;

  arrayMax(look->map->boxIndex) = 0 ;
  if (arrayMax (look->map->segs))
    {
      n1 = n2 = n3 = n4 = 0 ;
      /* bump the genes to their proper place */
      bumpDestroy (bumper) ;
      bumper = bumpCreate (12,  0) ;
      gLocBumpGene (look, bumper, FALSE, TRUE, &n1) ; /* cloud, down */
      bumpDestroy (bumper) ;
      bumper = bumpCreate (12,  0) ;
      gLocBumpGene (look, bumper, TRUE, TRUE, &n2) ; /* non-cloud down */
 
      bumpDestroy (bumper) ; 
      bumper = bumpCreate (12,  0) ;      
      gLocBumpGene (look, bumper, FALSE, FALSE, &n3) ; /* cloud, up */
      bumpDestroy (bumper) ; 
      bumper = bumpCreate (12,  0) ;
      gLocBumpGene (look, bumper, TRUE, FALSE, &n4) ; /* non-cloud up */
 
      /* actually draw */
      look->map->yChrom = 0 ;

      if (n2)
	{
	  look->map->yChrom += look->map->dh * n2 ;
	  ymax = gLocDrawGene (look, TRUE, TRUE) ; /* non-cloud down */ 
	  look->map->yChrom -= 2 * look->map->dh ;
	}
      
      if (n1)
	{
	  look->map->yChrom += look->map->dh * n1 ;
	  ymax = gLocDrawGene (look, FALSE, TRUE) ; /* cloud, down */
	}
      gLocDrawChromosome (look) ;

      if (n3)
	{
	  ymax = gLocDrawGene (look, FALSE, FALSE) ; /* cloud, up */
	  look->map->yChrom += look->map->dh * n3 ;
	}
      if (n4)
	{
	  look->map->yChrom -= 1 ;
	  ymax = gLocDrawGene (look, TRUE, FALSE) ; /* non-cloud up */ 
	  look->map->yChrom += look->map->dh * n4 ; 
	}
    }
  if (n3 + n4 > 0 && isGifDisplay) /* adapt the height of the gif/flash screen */
   {
     swfGraphResize (look->map->graphWidth, ymax+2) ;
     graphFitBounds (&look->map->graphWidth,&look->map->graphHeight) ;
   }

  gLocDrawScale (look) ;
  bumpDestroy (bumper) ;
  graphRedraw () ;
}

/************************************************************/
/*************************************************************/

static BOOL glocDoDisplay (KEY key, KEY from, BOOL isOldGraph, BOOL isSmall)
{
  GLoc  oldlook = 0, look ;
  KEY gene = 0 ;
  AC_HANDLE handle = 0 ;
		/* centre on parent object if poss */
  if (!key)
    return FALSE ;

  cDNAAlignInit () ;
  if (!ac_db)
    ac_db = ac_open_db (0,0) ;

  if (
      isOldGraph && 
      graphAssFind (&GLoc_MAGIC, &oldlook) &&
      oldlook->magic == &GLoc_MAGIC &&
      oldlook->key == key &&
      oldlook->map &&
      graphActivate (oldlook->map->graph)
      )
    { 
      gLocDraw () ;
      return TRUE ;
    }
  if (oldlook)
    gLocDestroy () ;

  handle = handleCreate () ;
  look = (GLoc) halloc (sizeof (struct GLocStruct), 0) ;  
  look->magic = &GLoc_MAGIC;

  look->db = ac_db ;          /* cosmetic, since we are inside xace */
  look->h = handle ;
  look->key = key ;
  look->gene = gene ;
  look->isSmall = isSmall ;

  look->map = gLocMapCreate (look, look->h, gLocDraw) ;
  if (!gLocConvert (look) || ! (isOldGraph || displayCreate (DtGLOC)))
    { 
      handleDestroy (handle) ;
      ac_free (look) ;
      display (key,from,TREE) ;
      return FALSE ;
    }
    
  graphRegister (DESTROY, gLocDestroy) ;

  graphRetitle (name(key)) ;

  graphAssociate (&GLoc_MAGIC, look) ; /* attach look to graph */
  gLocMapGraphAssociate (look->map) ; /* redraws */

  return TRUE ;
} /* glocDoDisplay */

/************/

BOOL glocBigDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  return glocDoDisplay (key, from, isOldGraph, FALSE) ;
} /* glocBigDisplay */

/************/

BOOL glocDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  return glocDoDisplay (key, from, isOldGraph, TRUE) ;
} /* glocDisplay */

/*************************************************************/
/*************************************************************/
/*************************************************************/

