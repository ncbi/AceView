/*  File: geneexpdisp.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the AceDB/AceView genome database package, written by
 *	Danielle and Jean Thierry-Mieg (NCBI) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the expression index of a gene accross many tissues and species
 * Exported functions:
       geneExpDisplay   (interface for the acedb display package 
 * HISTORY:

 * Created: Aug 24, 2014 (mieg)
 *-------------------------------------------------------------------
 */

/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "ac.h"
#include "graph.h"

#include "bump.h"
#include "wfiche/gtitle.h"
#include "display.h"
#include "pick.h"
#include "session.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "bs.h"
#include "lex.h"
#include "a.h"
#include "graph.h"
#include "query.h"
#include "cdna.h"
#include "bitset.h"

static AC_DB ac_db = 0 ;

typedef struct GXDMapStruct {
  magic_t *magic;        /* == GXD_MAGIC */
  AC_HANDLE h ;
  Graph graph ;

  float mag, yCurrent, leftMargin, offset ; 
  int graphWidth, graphHeight ;
  int max ;
  int a1, a2 ;

  BOOL hideHeader ;

  int   activeBox ;
  int	messageBox ;
  char  messageText[128] ;

  int starBox ;
  char starBuffer [1000] ; /* a place to write the * when i touch an mrna */

  Array segs ;      	/* array of SEG's from the obj */
  Array introns, intronSupport ;
  Array boxIndex ;    /* box -> int (index into look->segs) */
  int lastTrueSeg ;	/* arrayMax(look->segs) after Convert */

  void (*mapDraw) (void) ;
} GXDMap ;

typedef struct GXDStruct {
  magic_t *magic;        /* == GXD_MAGIC */
  AC_HANDLE h ;
  AC_DB db ;
  GXDMap *map ;
  int type ; /* 1: gene, 2:mrna, 3:product */
  BOOL onlyNm, onlyTouched, allProbes, showNoSignalProbes, showOnlyExactProbes, showAmbiguousProbes ;
  KEY intMap, key, gene, tg, mrna, product, from ;
  int a1, a2 ; /* int map coords of the outside object */
  unsigned int flag ;
  int maxIntronClone ;
  DICT *tissues ;
  Array menu ;

  float geneMedian, geneAverage, geneSigma ; /* applies to index of look->gene */
  Array geneIndexHisto ; /* distrib of index of look->gene */
} *GXD ;

typedef enum { Tzero, Tgene, Tmrna
	       , MrnaExon, VeryGoodCodingExon, GoodCodingExon, PoorCodingExon, uORFExon, Utr5, Utr3
	       , MrnaIntron, MrnaIntronFuzzy, GAP
	       , Valid5p, V5Capped, V5SL, V5Aggregate, V5stop
	       , Valid3p, V3aataaa, V3variant
	       , MrnaProbe, GENE
} HSTYPE  ;
typedef struct
  { 
     KEY gene, map ;
    int col, bcol, boxNeeded ;
    int a1, a2, iy1, iy2 ;
    int quality, nClones, height ;
    BOOL isDown, cloud ;
    int gid ;
    BOOL biblio, interactions, conserved, regulation, disease ;
  } SEG ;

typedef struct {
  AC_HANDLE h ;
  DICT *titleDict ;
  DICT *runDict ;
  int runMax, titleMax ;
  KEYSET run2title, title2run ;
  Array aa ; /* iz = unsigned char [k] ; index = iz/10.0 + 5 ; */
  float median, average, sigma ;
  Array indexHisto ;
  int histoStep ; /* resolution in tenths of indexpoints */
} GXDATA ;

#define GXINDEX(_gene,_run) ((gxData->aa && (KEYKEY(_gene))*(gxData->runMax) + (_run) < arrayMax (gxData->aa)) ? 5.0 + arr(gxData->aa, (KEYKEY(_gene))*(gxData->runMax) + (_run), unsigned char)/10.0 : 0)

static GXDATA *gxData = 0 ;

static magic_t GXD_MAGIC = "GXD" ;
static magic_t GXDMap_MAGIC = "GXMap" ;
static BOOL gxdConvert (GXD look) ;

/************************************************************/
/************************************************************/

static BOOL gxdParseBinaryData (const char *fNam)
{
  BOOL ok = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreateFromFile (fNam, "rb", 0, h) ;
 
  if (ai)
    {
      int nn ;

      messStatus ("Binary Parsing expression") ;
      aceInBinary (ai, (char *)&nn, sizeof(int)) ;
      if (nn == 123456789) /* check big/little enddians */
	{
	  aceInBinary (ai, &nn, sizeof(int)) ;
	  if (nn == gxData->runMax)
	    {
	      aceInBinary (ai, &nn, sizeof(int)) ;
	      gxData->aa = arrayHandleCreate (gxData->runMax * (lexMax (_VGene) + 1), unsigned char, gxData->h) ;
	      array (gxData->aa, nn - 1, unsigned char) = 0 ; /* make room */
	      ok = aceInBinary (ai, arrp (gxData->aa, 0, unsigned char), nn) ; 
	      if (! ok) arrayDestroy (gxData->aa) ;	      
	    }
	}
    }

  ac_free (h) ;
 
  if (! ok)
    {
      arrayDestroy (gxData->aa) ;
      if (ai)
	unlink (fNam) ;
    }
  return ok ;
} /* gxdParseBinaryData  */      

/************************************************************/

static BOOL gxdWriteBinaryData (const char *fNam)
{
  int nn, err = 999 ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreateToFile (fNam, "wb", h) ;
 
  if (ao)
    {
      nn = 123456789 ; /* check big/little enddians */
      aceOutBinary (ao, &nn, sizeof(int)) ;
      
      nn = gxData->runMax ;
      aceOutBinary (ao, &nn, sizeof(int)) ;
      
      nn = arrayMax (gxData->aa) ;
      aceOutBinary (ao, &nn, sizeof(int)) ;
      err = aceOutBinary (ao, arrp (gxData->aa, 0, unsigned char), nn) ; 
    }
  ac_free (h) ;
  if (err && ao)
    unlink (fNam) ;
 
  return err == 0 ? TRUE : FALSE ; ;
} /* gxdWriteBinaryData  */      

/************************************************************/

static BOOL gxdParseRunTitleData (const char *tNam)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tNam, FALSE, h) ;
  BOOL ok = FALSE ;
  const char *ccp ;
  int run, title ;

  if (ai)
    {
      messStatus ("Parsing Run2SortingTitle") ;
      while (aceInCard (ai))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#') continue ;
	  dictAdd (gxData->runDict, ccp, &run) ;
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#') continue ;
	  dictAdd (gxData->titleDict, ccp, &title) ;
	  keySet (gxData->run2title, run) = title ;
	  keySet (gxData->title2run, title) = run ;
	}
      gxData->runMax = dictMax (gxData->runDict) + 1 ; 
      gxData->titleMax = dictMax (gxData->runDict) + 1 ; 
      ok = TRUE ;
    }
  ac_free (h) ;
  return ok ;
} /* gxdParseRunTitleData */

/************************************************************/

static BOOL gxdDataParse (void)
{
  BOOL ok = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array histo ;
  const char *ccp ;
  int run = 0 ;
  KEY gene = 0 ; 
  int step, i, nn ;
  const char *tNamB = "database/gxpRun2SortingTitle.txt_binary" ;
  const char *tNam  = "/home/mieg/GeneExpress/gxpRun2SortingTitle.txt" ;
  const char *fNamB = "database/gxp.av.GENE.u.binary" ;
  const char *fNam  = "/home/mieg/GeneExpress/gxp.av.GENE.u.ace.gz" ;

  gxData = (GXDATA *)messalloc (sizeof (GXDATA)) ;
  gxData->h = ac_new_handle () ;
  gxData->runDict = dictHandleCreate (1000, gxData->h) ;
  gxData->titleDict = dictHandleCreate (1000, gxData->h) ;
  gxData->run2title = keySetHandleCreate (gxData->h) ;
  gxData->title2run = keySetHandleCreate (gxData->h) ;
  gxData->indexHisto = histo = arrayHandleCreate (1024, int, gxData->h) ;
  gxData->histoStep  = step = 2 ;  /* number of points per unit */

  if (filName (tNam, 0, "r"))
    gxdParseRunTitleData (tNam) ;
  
  if (filName (tNamB, 0, "r") &&  filName (fNamB, 0, "rb"))
    {
      gxdParseRunTitleData (tNamB) ;
      ok = gxdParseBinaryData (fNamB) ;
    }

  if (! ok && filName (tNam, 0, "r"))
    gxdParseRunTitleData (tNam) ;

  if (! ok && filName (fNam, 0, "r"))
    {
      ai = aceInCreate (fNam, FALSE, h) ;
      if (ai)
	{
	  Array aa = gxData->aa = arrayHandleCreate (gxData->runMax * (lexMax (_VGene) + 1), unsigned char, gxData->h) ;
	  messStatus ("Parsing expression") ;
	  while (aceInCard (ai))
	    {
	      ccp = aceInWord (ai) ;
	      if (! ccp) gene = 0 ;
	      else if (! strcasecmp (ccp, "Gene"))
		{
		  ccp = aceInWord (ai) ;
		  gene = 0 ;
		  if (ccp) lexword2key (ccp, &gene, _VGene) ;
		  gene = KEYKEY (gene) ;
		}
	      else if (gene && (! strcasecmp (ccp, "Group_U") || ! strcasecmp (ccp, "Run_U")))
		{
		  float z = 0 ;
		  int iz ;
		  
		  ccp = aceInWord (ai) ;
		  run = 0 ;
		  if (ccp && dictFind (gxData->runDict, ccp, &run) && aceInFloat (ai, &z))
		    {
		      float seqs = 0 ;
		      int k = gene * gxData->runMax + run ;
		      z -= 5 ; if (z < 0) z = 0 ; iz = 10 * z + .499 ;  if (iz > 255) iz = 255 ;  
		      if (iz < 2) iz = 2 ; /* reserve 0/1 for NE/ NA/ */
		      aceInFloat (ai, &seqs) ;
		      if (seqs == 0) iz = 0 ;
		      else
			{
			  /* NA/NE is in column 29, this is horrible */
			  int i = 24 ;
			  while (i--)  aceInWord (ai) ;
			  
			  if ((ccp = aceInWord (ai)) && ! strcmp (ccp, "NA")) iz = 1 ;
			}
		      array (aa, k, unsigned char) = iz & 0xff ;
		    }
		}
	    }
	  ok = TRUE ;
	  gxdWriteBinaryData (fNamB) ;
	  {
	    vTXT txt = vtxtHandleCreate (h) ;
	    vtxtPrintf (txt, "\\cp %s ", filName(tNam, 0,"r")) ;
	    vtxtPrintf (txt, " %s/gxpRun2SortingTitle.txt_binary\n", sessionFilName ("database", 0, "rd")) ;
	    system (vtxtPtr (txt)) ;
	  }
	}
    }  
  
  /* normalize the histo surface to 100 */
  {
    Array aa = gxData->aa ;
    unsigned char * cp = arrp (aa, 0, unsigned char) ;

    nn = 0 ; i = arrayMax (aa) ;
    while (cp++, i--)
      if (*cp > 1)
	{
	  int iz = utArrondi ((50.0 + *cp) * step/10.0) ;
	  array (histo, iz, int) += 8 ;
	  array (histo, iz + 1, int) += 1 ;
	  array (histo, iz - 1, int) += 1 ;
	  nn++ ;
	}
  }

  for (i = 0 ; i < arrayMax (histo) ; i++)
    {
      int iz = array (histo, i, int) ;
      array (histo, i, int) = utArrondi (10 * iz/nn) ;
    }

  ac_free (h) ;
  fprintf (stderr, "gxdDataParse located %d values in %d runs\n", gxData->aa ? arrayMax(gxData->aa) : 0, gxData->runMax) ;
  return ok ;
} /* gxdDataParse */

/************************************************************/

static BOOL gxdDataInit (void)
{
  static BOOL ok = FALSE ;
  if (gxData)
    ok = TRUE ;
  else
    ok = gxdDataParse () ;

  return ok ;
} /* gxdDataInit */

/************************************************************/
/************************************************************/

#define GXD2GRAPH(_gxd,_x) (((_x) - (_gxd)->offset) * (_gxd)->mag + (_gxd)->leftMargin)

/************************************************************/
/************************************************************/
/************************************************************/

static void gxdShowSegs (Array segs)
{
  SEG *seg ;
  int ii ;
  
  if (segs && arrayExists (segs))
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

static GXDMap *gxdMapCurrent (char *caller)
{
  GXDMap *map = 0 ;
  
  if (!graphAssFind (&GXDMap_MAGIC,&map))
    messcrash("%s() could not find GXDMap on graph", caller);
  if (!map)
    messcrash("%s() received NULL GLoc pointer", caller);
  if (map->magic != &GXDMap_MAGIC)
    messcrash("%s() received non-magic GXDMap pointer", caller);
  
  return map ;
} /* gxdMapCurrent */

/************************************************************/

static void gxdMapResize (void)
{
  GXDMap *map = gxdMapCurrent ("hmapPick") ;
  
  pickRememberDisplaySize (DtGLOC) ;
  
  graphFitBounds (&map->graphWidth,&map->graphHeight) ;
  map->leftMargin = 1 ;
  map->mag = ( (map->graphWidth - map->leftMargin))/map->max ;
  
  map->mapDraw () ;
} /* gxdMapResize */

/************************************************************/

static void gxdMapKbd (int k)
{ 
  GXDMap *map = gxdMapCurrent ("hmapPick") ;
  
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
	gxdMapResize () ;
	return ;
      default:
	return ;
      }
  map->mapDraw () ;
}

/************************************************************/

static void gxdMapFollow (GXDMap *map, int box)
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
 
} /* gxdMapFollow */

/************************************************************/

static void gxdMapSelect (GXDMap *map, int box)
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
} /* gxdMapSelect */

/*************************************/

static void gxdMapPick (int box,  double x , double y) 
{
  GXDMap *map = gxdMapCurrent ("gxdMapPick") ;

  if (!box)
    return ;
  if (box < arrayMax (map->boxIndex) && arr(map->boxIndex,box,int)) /* a SEG */
    {
      if (box == map->activeBox) /* a second hit - follow it */
	gxdMapFollow (map, box) ;
      else
	gxdMapSelect (map, box) ;
    }

  return;
} /* gxdMapPick */

/************************************************************/

static GXDMap *gxdMapCreate (void *look, AC_HANDLE handle, void (*mapDraw) (void))
{
  GXDMap *map = (GXDMap *) halloc (sizeof (struct GXDMapStruct), handle) ;
  map->magic = &GXDMap_MAGIC ;
  map->h = handle ; /* do not destroy, belongs to calling routine */
  map->boxIndex = arrayHandleCreate (64,int, handle) ;
  map->segs = arrayHandleCreate (64, SEG, handle) ;
  map->activeBox = 0 ;
  map->lastTrueSeg = arrayMax(map->segs) ;
  map->mapDraw = mapDraw ;

  return map ;
} /* gxdMapCreate */

/************************************************************/

static void gxdMapGraphAssociate (GXDMap *map)
{
  graphRegister (KEYBOARD, gxdMapKbd) ;
  graphRegister (PICK, gxdMapPick) ;
  graphAssRemove (&GXDMap_MAGIC) ;
  graphAssociate (&GXDMap_MAGIC, map) ; /* attach map to graph */
  graphRegister (RESIZE, gxdMapResize) ;
  map->graph = graphActive() ;

  gxdMapResize () ; /* redraws */
  return ;
} /* gxdMapCreate */

/************************************************************/
/************************************************************/

static GXD currentGxd (char *caller)
{
  GXD gxd ;

  if (!graphAssFind (&GXD_MAGIC,&gxd))
    messcrash("%s() could not find GXD on graph", caller);
  if (!gxd)
    messcrash("%s() received NULL GXD pointer", caller);
  if (gxd->magic != &GXD_MAGIC)
    messcrash("%s() received non-magic GXD pointer", caller);
  
  return gxd;
} /* currentGxd */


/************************************************************/
/***************** Registered routines *********************/

static void gxdDestroy (void)
{
  GXD look = currentGxd("gxdDestroy") ;

  if (look)
    {
      ac_free (look->h) ;
      
      graphAssRemove (&GXD_MAGIC) ; /* detach look from graph */
      
      look->magic = 0 ;
      ac_free (look) ;
    }
  if (0) gxdShowSegs (0) ; /* for compiler happiness */

  return;
} /* gxdDestroy */

/*************************************************************/
/*
static void gXdRecalculate(void)
{
  GXD look = currentGxd("gXdRecalculate") ;
  
  if (look)
    {
      look->map->segs = arrayReCreate (look->map->segs,64, SEG) ;
      gxdConvert (look) ;
      look->map->mapDraw () ;
    } 
}
*/
/*************************************************************/
/************************************************************/

static BOOL gxdConvert (GXD look)
{
  BOOL ok = FALSE ;

  if (class(look->key) == _VGene)
    { look->gene = look->key ; look->type = 1 ; ok = TRUE ; }
  else
    goto done ;
  
 done :
    return ok ;
} /* gxdConvert */

/************************************************************/
/* intensity color scale */
static void gxScaleBox (const char *t, int color, int color1, float x, float y) 
{ 
  int box = graphBoxStart () ; 
  graphText (t, x, y) ; 
  graphBoxEnd () ; 
  
  graphBoxColor (box, color1, color) ; 
  if (0) graphText (t, x + 7, y) ;
} /*  gxScaleBox */

/************************************************************/
/* color scale */
static void gxColorScale (float x, float y)
{  /* caption */
  graphText ("Expression", x, y) ; y++ ;
  graphText ("quantiles", x+.5, y) ; y++ ;
  gxScaleBox(" None ", LIGHTGRAY, BLACK, x, y) ; y++ ;
  gxScaleBox(" Weak ", PALEGRAY, BLACK, x, y) ; y++ ;

  gxScaleBox("   1  ", GREEN1, BLACK, x, y) ; y++ ;
  gxScaleBox("   2  ", GREEN2, BLACK, x, y) ; y++ ;
  gxScaleBox("   3  ", GREEN3, BLACK, x, y) ; y++ ;
  gxScaleBox("   4  ", GREEN4, BLACK, x, y) ; y++ ;
  gxScaleBox("   5  ", GREEN5, BLACK, x, y) ; y++ ;
  gxScaleBox("   6  ", GREEN6, WHITE, x, y) ; y++ ;
  gxScaleBox("   7  ", GREEN7, WHITE, x, y) ; y++ ;
  gxScaleBox("   8  ", GREEN8, WHITE, x, y) ; y++ ;
  if(0)  gxScaleBox(" high ", GREEN8, YELLOW, x, y) ; y++ ;

} /*  gxColorScale */

/************************************************************/
/* index <--> fpkm table */
static void gxIndexFpkmScale (float x, float y)
{ 
  int a = 8 ;
  double b = .25 ;

  /* table */
  graphText ("Index" ,x+.5, y) ;  graphText ("sFPKM" , x + 5, y) ;
  for (a = 8, b = .25, y += 2 ; a < 21 ; a++, b = 2*b, y++)
    graphText (messprintf (b < 1 ? "  %2d         %.2f" :  "  %2d         %.0f", a, b), x, y) ;
  graphText ("  ---         ---", x, y) ;
  
  /* cadre */
  y += 1.7 ;
  graphLine (x + 0 , 1.5, x +  0, y) ;      
  graphLine (x + 4 , 1.5, x +  4, y) ;
  graphLine (x + 11, 1.5, x + 11, y) ;
  
  graphLine (x, 1.5, x + 11, 1.5) ;
  graphLine (x, 3.2, x + 11, 3.2) ;
  graphLine (x, y, x + 11, y) ;
} /*  gxIndexFpkmScale */

/************************************************************/

static void gxGetGeneIndexDistrib (GXD look, AC_HANDLE h)
{
  int run, nn, i, iz ;
  float z ;
  int step = gxData->histoStep ;
  Array a, histo ;

  a = arrayCreate (1024, float) ;
  histo = look->geneIndexHisto = arrayHandleCreate (1024, int, h) ;
 
  for (nn = 0, run = 1 ; run < gxData->runMax ; run++)
    {
      z = GXINDEX (look->gene,run) ;
      if (z > 5.05)
	{
	  array (a, nn++, float) = z ;
	  iz = utArrondi (z * step) ;
	  array (histo, iz, int) += 8 ;
	  array (histo, iz + 1, int) += 1 ;
	  array (histo, iz - 1, int) += 1 ;
	}
    } 
  floatVariance (a, &look->geneMedian, &look->geneAverage, &look->geneSigma) ;
  if (nn) /* normalize the surface to 100 */
    for (i = 0 ; i < arrayMax (histo) ; i++)
      {
	iz = array (histo, i, int) ;
	array (histo, i, int) = utArrondi (10 * iz/nn) ;
      }
  arrayDestroy (a) ;

  return ;
} /* gxGetGeneIndexDistrib */

/************************************************************/

static void gxPlotDistrib (float x0, float y0, GXD look)
{  
  Array aa ;
  Array aa0 = gxData->indexHisto ;
  Array aa1 = look->geneIndexHisto ;
  int i, iMax, box, h1, h2, color, pass ;
  float width = 14, height = 6, old ; /* diagram size */
  float dx, dy, ddx, z ;
  float a1, a2, b1, b2, aMax, bMax ;
  int step = gxData->histoStep ; /* step the histo every 5 tenth */

  y0 += height ;  /* set x0, y0 to the origin of the coordinates */
  iMax = arrayMax (aa0) ;		   
  for (bMax = 0, i = 0 ; i < iMax ; i++)
    {
      b1 = arr (aa0, i, int) ;
      if (bMax < b1) bMax = b1 ;
      b1 = arr (aa1, i, int) ;
      if (bMax < b1) bMax = b1 ;
    }
  aMax = utMainRoundPart (1.0 + arrayMax(aa0)) ;
  aMax = 25 ;
  if (arrayMax (aa1) < 20 * step) aMax = 20 ;
  dx = width / aMax ;  /* dx = movemennt per point of index */
  dy = height / (1 + bMax) ;
  
  /* Draw the exis 
   * y coordinates are oriented towards the bottom of the screen 
   */ 
  if (0) 
    graphLine (x0, y0, x0 + width, y0) ;
  if (0) graphLine (x0, y0, x0, y0 - height) ; 
  if (0)
    for (i = 4  ; i <= aMax ; i++)
      {
	a2 = x0 + (i - 4) * dx ;
	graphLine (a2, y0, a2, y0 + .2) ;
      }
  old = graphTextHeight (0) ;
  graphTextHeight (0.4 *old) ;
  for (i = 10 ; i <= aMax ; i += 10)
    {
      a2 = x0 + (i - 4) * dx ;
      graphLine (a2, y0, a2, y0 + .5) ;
      graphText (messprintf("%d", i), a2-1, y0 + .3) ;
    }
  graphTextHeight (old) ;

  /* Draw the histogram */
  ddx = dx / step ;  /* ddx = x step per exported value */
  for (pass = 0 ; pass < 2 ; pass++)
    {
      aa = pass ? aa1 : aa0 ;

      a1 = x0 ;
      i = 4 * step ;
      h1 = arr (aa, i, float) * ( 2 - pass) ;
      b1 = y0 - dy * h1 ; 
      
      for (++i ; i < aMax * step ; i++)
	{
	  h2 =  arr (aa, i, int) * ( 2 - pass) ;
	  a2 = a1 + ddx  ; b2 =  y0 - dy * h2 ;
	  box = graphBoxStart () ;
	  graphLine (a1, b1, a2, b2) ;
	  graphBoxEnd () ;
	  
	  if (0)
	    {
	      z = i/10.0 ;
	      color = WHITE ;
	      
	      if (z > 18) { color = GREEN1 + 7 ;  }
	      else if (z > 14.4) { color = GREEN1 + 7 ; }
	      else if (z > 13.5) { color = GREEN1 + 6 ; }
	      else if (z > 12.6) color = GREEN1 + 5 ;
	      else if (z > 11.7) color = GREEN1 + 4 ;
	      else if (z > 10.8) color = GREEN1 + 3 ;
	      else if (z > 9.9) color = GREEN1 + 2 ;
	      else if (z > 8.9) color = GREEN1 + 1 ;
	      else if (z > 5.1) color = GREEN1 + 0 ;
	      
	      else if (z >= 5.05) color =  PALEGRAY ;
	      else if (z == 5) color = LIGHTGRAY ;
	    }

	  color = pass ? GREEN6 : BLUE4 ;
	  graphBoxColor (box, color, TRANSPARENT) ;
	  a1 = a2 ; b1 = b2 ; h1 = h2 ;
	}
    }
  {
    box = graphBoxStart () ;
    graphText ("This gene", x0, y0 + 1) ;
    graphBoxEnd () ;
    graphBoxColor (box, GREEN6, TRANSPARENT) ;

    box = graphBoxStart () ;
    graphText ("All genes", x0, y0 + 1.8) ;
    graphBoxEnd () ;
    graphBoxColor (box, BLUE5, TRANSPARENT) ;
    graphText ("log2 distributions", x0, y0 + 2.8) ;
  }
  return ;
} /* gxPlotDistrib */

/************************************************************/

static void gxdDraw (void)
{ 
  AC_HANDLE h = ac_new_handle () ;
  GXD look = currentGxd("gxDraw") ;
  int run, title, box, x, y, globalBox ;
  const char *ccp ;
  char species[12] ;
  BitSet bx, by, bxy ;
  int xWidth = 5, xMax = 1, yMax = 1, titleWidth = 2 ;
  int color ;

  bx = bitSetCreate (gxData->titleMax, h) ;
  by = bitSetCreate (gxData->titleMax, h) ;
  bxy = bitSetCreate (gxData->titleMax * gxData->titleMax, h) ;

  for (title = 1 ; title <= dictMax (gxData->titleDict) ; title++)
    {
      ccp = dictName (gxData->titleDict, title) ;
      sscanf (ccp, "%02d_%02d_%3s_",&x,&y,species) ;
      bitSet (bx, x) ; bitSet (by, y) ;
    }
  graphClear () ;
  graphTextFormat(BOLD) ;

  gxGetGeneIndexDistrib (look, h) ;
  globalBox = graphBoxStart () ;
  box = graphBoxStart () ;
  {
    /*     double z = exp(look->geneAverage * log(2.0))/1000 ; */

    float old = graphTextHeight (0) ;
    graphTextHeight (1.2 * old ) ;
    graphText (name(look->gene), 1, 3) ;
    graphText (messprintf(" Gene expression in %d primates, %d tissues, from the NHPRTR project in sFPKM",  bitSetCount(bx), bitSetCount (by)), 20,3) ;
    graphTextHeight (old ) ;
  }
  graphBoxEnd () ;
  graphBoxColor (box, BLACK, YELLOW) ; 

  bx = bitSetReCreate (bx, 256) ;
  by = bitSetReCreate (by, 256) ;

  for (title = 1 ; title <= dictMax (gxData->titleDict) ; title++)
    {
      ccp = dictName (gxData->titleDict, title) ;
      x = strlen (ccp) ;
      if (x > titleWidth) titleWidth = x ;
    }
  titleWidth += 4 - xWidth  - 10 ;

  for (title = 1 ; title <= dictMax (gxData->titleDict) ; title++)
    {
      ccp = dictName (gxData->titleDict, title) ;
      sscanf (ccp, "%02d_%02d_%3s_",&x,&y,species) ;
      if (! bit(bx, x)) 
	{
	  color = WHITE ; 
	  if (x < 4) color = BLUE ;
	  else if (x < 10) color = DARKRED ;
	  else if (x < 12) color = BROWN ;
	  else if (x < 15) color = DARKGREEN ;
	  else  color = ORANGE ;

	  bitSet (bx, x) ; 
 
	  if (x > 3) x-- ;

 	  box = graphBoxStart () ;
	  graphText (species , titleWidth + xWidth * x, 4) ;
	  graphBoxEnd () ;

	  graphBoxColor (box, color, WHITE) ;
	}
      if (! bit(by, y)) 
	{
	  if (y == 1) color = GREEN7 ;
	  else if (y <= 4 ) color = DARKGREEN ;
	  else if (y <= 5) color = DARKBLUE ;
	  else if (y <= 7) color = MAGENTA ;
	  else if (y <= 8) color = DARKRED ;
	  else if (y <= 9) color = BROWN ;
	  else if (y <= 10) color = RED7 ;
	  else if (y <= 11) color = RED6 ;
	  else if (y <= 12) color = PURPLE ;
	  else if (y <= 13) color = DARKGRAY ;
	  else if (y <= 14) color = LAVANDER ;
	  else if (y <= 15) color = BLUE8 ;
	  else if (y <= 16) color = DARKVIOLET ;
	  else if (y <= 17) color = BLUE6 ;
	  else if (y <= 18) color = RED8 ;

	  else  color = BLACK ;

	  bitSet (by,y) ;

	  if (y >= 4) y-= 2 ;

	  box = graphBoxStart () ; 
	  graphText (ccp + 10 , 3, 4 + y) ;
	  graphBoxEnd () ;  

	  graphBoxColor (box, color, WHITE) ;

	  if (y > yMax) yMax = y ;
	}
    }

  graphTextFormat(PLAIN_FORMAT) ;

  for (run = 1 ; run < gxData->runMax ; run++)
    {
      if (0) graphText (dictName (gxData->runDict, run), 4, 3 + run) ;
      title = keySet (gxData->run2title, run) ;
      ccp = dictName (gxData->titleDict, title) ;
      
      sscanf (ccp, "%02d_%02d_%3s_",&x,&y,species) ;

      if (! bit (bxy, x * gxData->titleMax + y))
	{
	  float z = GXINDEX (look->gene,run) ;
	  int color1 = BLACK ;

	  color = WHITE ;

	  if (z > 18) { color = GREEN1 + 7 ; color1 = YELLOW ; }
	  else if (z > 14.4) { color = GREEN1 + 7 ; color1 = WHITE ; }
	  else if (z > 13.5) { color = GREEN1 + 6 ; color1 = WHITE ; }
	  else if (z > 12.6) color = GREEN1 + 5 ;
	  else if (z > 11.7) color = GREEN1 + 4 ;
	  else if (z > 10.8) color = GREEN1 + 3 ;
	  else if (z > 9.9) color = GREEN1 + 2 ;
	  else if (z > 8.9) color = GREEN1 + 1 ;
	  else if (z > 5.1) color = GREEN1 + 0 ;

	  else if (z >= 5.05) color =  PALEGRAY ;
	  else if (z == 5) color = LIGHTGRAY ;

	  else if (z > 16) color = 1 ;
	  
	  bitSet (bxy, x * gxData->titleMax + y) ;
	  box = graphBoxStart () ;

	  if (x >3) x-- ;
	  if (y >= 4) y-= 2 ;

	  if (z > 5.11)
	    {
	      double z1 = exp(z*log(2.0))/1000.0 ;
	      char *ff = "%.0f" ;

	      if (z1 < 1) ff = "%0.2f";
	      else if (z1 < 10) ff = "%.2f";
	      else if (z1 < 100) ff = "%.1f";
	      else if (z1 < 100) ff = "%.1f";
	      else if (z1 < 1000) ff = " %.0f";
	      graphText (messprintf (ff,  z1), titleWidth + xWidth * x, 4 + y) ;
	    }
	  else
	    graphText ("    ", titleWidth + xWidth*x, 4 + y) ;
	  graphBoxEnd () ;
	  if (x > xMax) xMax = x ;
	  if (y > yMax) yMax = y ;

	  graphBoxColor (box, color1, color) ;
	}
    }

  gxColorScale (titleWidth + xWidth*xMax + 8, 5) ;
  if (0) gxIndexFpkmScale (titleWidth + xWidth*xMax + 22, 4) ;
  gxPlotDistrib (titleWidth + xWidth*xMax +18, 4, look) ;

  graphBoxEnd () ; /* globalBox */

  if (1) /* adapt the height of the gif/flash screen */
    {
      float x1, x2, y1, y2 ;
      graphBoxDim (globalBox, &x1, &y1, &x2, &y2) ;
       if (isGifDisplay) swfGraphResize (x2 - x1 + 2, y2 - y1 + 4) ;
       if (isGifDisplay) svgGraphResize (x2 - x1 + 2, y2 - y1 + 4) ;
      graphFitBounds (&look->map->graphWidth,&look->map->graphHeight) ;
    }

  graphRedraw () ;
  ac_free (h) ;
  return ;
} /* gxdDraw */

/************************************************************/
/************************************************************/

BOOL  geneExpDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  GXD  oldlook = 0, look ;
  AC_HANDLE handle = 0 ;
		/* centre on parent object if poss */
  if (!key)
    return FALSE ;

  if (! gxdDataInit ())
    return FALSE ;

  cDNAAlignInit () ;
  if (!ac_db)
    ac_db = ac_open_db (0,0) ;

  if (
      isOldGraph && 
      graphAssFind (&GXD_MAGIC, &oldlook) &&
      oldlook->magic == &GXD_MAGIC &&
      oldlook->key == key &&
      oldlook->map &&
      graphActivate (oldlook->map->graph)
      )
    { 
      oldlook->from = from ;
      gxdDraw () ;
      return TRUE ;
    }
  if (oldlook)
    gxdDestroy () ;

  handle = handleCreate () ;
  look = (GXD) halloc (sizeof (struct GXDStruct), 0) ;  
  look->magic = &GXD_MAGIC;

  look->db = ac_db ;          /* cosmetic, since we are inside xace */
  look->h = handle ;
  look->from = from ;
  look->gene = look->key = key ;

  look->map = gxdMapCreate (look, look->h, gxdDraw) ;
  if (!gxdConvert (look) || ! (isOldGraph || displayCreate (DtGeneExp)))
    { 
      handleDestroy (handle) ;
      ac_free (look) ;
      display (key,from,TREE) ;
      return FALSE ;
    }
    
  graphRegister (DESTROY, gxdDestroy) ;

  pickRememberDisplaySize (DtGeneExp) ;

  graphRetitle (name(key)) ;

  graphAssociate (&GXD_MAGIC, look) ; /* attach look to graph */
  gxdMapGraphAssociate (look->map) ; /* redraws */

  return TRUE ;
}
