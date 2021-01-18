/*
 * Authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * June 2012
 * pg2pg.c
 *   Input: ace file containing sequence->Souce_Eons and IntMap
 *   Output: regroup the transcripts into gene models
*/


#include <regular.h>
#include <ac.h>
#include <acedna.h>
#include <aceio.h>
#include <dict.h>
#include <vtxt.h>
#include <math.h>

typedef struct gxStruct { 
  AC_HANDLE h ; 
  DICT *chromDict ;
  DICT *geneDict ;
  DICT *locusDict ;
  DICT *mrnaDict ;
  DICT *gidDict ;
  Array models, modelsSorted, segs, pairs, gx2gbx ;
  const char *title ;
  const char *pg2pgFileName ;
  const char *pg2gbxFileName ;
  const char *outFileName ;

  BOOL gzi, gzo ;  Stack info ;
} GX ; 

typedef struct mmStruct { int m, chrom, a1, a2, isDown, iseg1, iseg2, ln, gid, locus, nmid, gbx ; } MM ;
typedef struct segStruct { int m, chrom, a1, a2, x1, x2, isDown ; } SEG ;
typedef struct pairStruct { int m1, m2 ; } PP ;
typedef struct gbxStruct { int gid, locus, chrom, a1, a2, isDown ; } GBX ;

/*************************************************************************************/
/*************************************************************************************/

static void gxParseInit (GX *gx)
{
  gx->chromDict = dictHandleCreate (500, gx->h) ;
  gx->geneDict = dictHandleCreate (50000, gx->h) ;
  gx->mrnaDict = dictHandleCreate (50000, gx->h) ;
  gx->gidDict = dictHandleCreate (50000, gx->h) ;
  gx->locusDict = dictHandleCreate (50000, gx->h) ;
  gx->models = arrayHandleCreate (100000, MM, gx->h) ;
  gx->segs = arrayHandleCreate (100000, SEG, gx->h) ;
  gx->info = stackHandleCreate (10000, gx->h) ;
  pushText (gx->info, "toto") ;
 
} /* gxParseInit */

/*************************************************************************************/
/* parse the ace file run->deep and run->Length */
static int gxAceParse (GX *gx, const char* fileName)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  int pg, chrom, nseg = 0, gid, locus, a1, a2, x1, x2, nLine = 0 ;
  const char *ccp ;
  SEG *seg ;
  MM *mm ;
  BOOL isPg = FALSE ;
    
  ai = aceInCreate (fileName, gx->gzi, h) ;
  if (! ai)
    messcrash ("cannot open input file -deep %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  while (aceInCard (ai)) 
    { /* parse the ace file */
      ccp = aceInWord (ai) ;
      if (! (++nLine % 1000000))
	fprintf (stderr, "%d\n", nLine) ;
      if (! ccp || ! *ccp)
	{
	  isPg = FALSE ; pg = 0 ;
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Sequence"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      dictAdd (gx->geneDict, ccp, &pg) ;
	      mm = arrayp (gx->models, pg, MM) ; /* make room */
	      mm->m = pg ;
	      isPg = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "IntMap"))
	{
	  if (isPg && (ccp = aceInWord (ai)))
	    {
	      dictAdd (gx->chromDict, ccp, &chrom) ;
	      if (aceInInt (ai, &a1) && aceInInt (ai, &a2))
		{
		  mm = arrayp (gx->models, pg, MM) ; /* make room */
		  mm->chrom = chrom ;
		  mm->a1 = a1 ; mm->a2 = a2 ; 
		}
	    }
	  continue ;
	}
       if (ccp && !strcasecmp (ccp, "NM_id"))
	{
	  if (isPg && (ccp = aceInWord (ai)))
	    {
	      dictAdd (gx->gidDict, ccp, &gid) ;
	      mm = arrayp (gx->models, pg, MM) ; /* make room */
	      mm->nmid = gid ;
	    }
	  continue ;
	}
       if (ccp && !strcasecmp (ccp, "GeneId_pg"))
	{
	  if (isPg && (ccp = aceInWord (ai)))
	    {
	      mm = arrayp (gx->models, pg, MM) ; /* make room */
	      if (mm->gid) 
		ccp = messprintf ("%s:%s",dictName(gx->gidDict, mm->gid), ccp) ;
	      dictAdd (gx->gidDict, ccp, &gid) ;
	      mm->gid = gid ;
	    }
	  continue ;
	}
       if (ccp && !strcasecmp (ccp, "Locus"))
	{
	  if (isPg && (ccp = aceInWord (ai)))
	    {
	      mm = arrayp (gx->models, pg, MM) ; /* make room */
	      if (mm->locus) 
		ccp = messprintf ("%s:%s",dictName(gx->locusDict, mm->locus), ccp) ;
	      dictAdd (gx->locusDict, ccp, &locus) ;
	      mm->locus = locus ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Source_Exons"))
	{
	  if (isPg && aceInInt (ai, &x1) && aceInInt (ai, &x2))
	    {
	      seg = arrayp (gx->segs, nseg++, SEG) ;
	      seg->m = pg ;
	      seg->x1 = x1 ; 
	      seg->x2 = x2 ;
	    }
	  continue ;
	}
    }

  ac_free (h) ;
  fprintf (stderr, "// parsed %d models in %d chrom with %d exons\n"
	   , arrayMax (gx->models)
	   , dictMax(gx->chromDict)
	   , arrayMax (gx->segs)
	   ) ;
  return nseg ;
} /* gxAceParse */

/*************************************************************************************/

static int segA1Order (const void *a, const void *b)
{
  const SEG *up = (const SEG *)a, *vp = (const SEG *)b ;
  int n ;

  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->isDown - vp->isDown ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->m - vp->m ; if (n) return n ;

  return 0 ;
} /* segAiOrder */

/*************************************************************************************/

static int mmA1Order (const void *a, const void *b)
{
  const MM *up = (const MM *)a, *vp = (const MM *)b ;
  int n ;

  n = up->iseg1 - vp->iseg1 ; if (n) return n ;
  return 0 ;
} /* mmAiOrder */

/*************************************************************************************/

static int ppOrder (const void *a, const void *b)
{
  const PP *up = (const PP *)a, *vp = (const PP *)b ;
  int n ;

  n = up->m1 - vp->m1 ; if (n) return n ;
  n = up->m2 - vp->m2 ; if (n) return n ;
  return 0 ;
} /* ppOrder */

/*************************************************************************************/

static int gxIntMap (GX *gx)
{
  int ii, a1, a2 ;
  int nsegs = arrayMax (gx->segs) ;
  int nmodels = arrayMax (gx->models) ;
  SEG *seg ;
  MM *mm ;

  for (ii = 0, seg = arrayp (gx->segs, ii, SEG) ; ii < nsegs ; ii++, seg++)
    {
      mm = arrp (gx->models, seg->m, MM) ;
      seg->chrom = mm->chrom ;
      a1 = mm->a1 ; a2 = mm->a2 ;
      if (a1 < a2)
	{
	  seg->isDown = 1 ;
	  seg->a1 = a1 + seg->x1 - 1 ; 
	  seg->a2 = a1 + seg->x2 - 1 ; 
	}
      else
	{
	  seg->isDown = 0 ;
	  seg->a1 = a1 - seg->x2 + 1 ; 
	  seg->a2 = a1 - seg->x1 + 1 ; 
	}      
    }

  for (ii = 0, mm = arrayp (gx->models, ii, MM) ; ii < nmodels ; ii++, mm++)
    {
      a1 = mm->a1 ; a2 = mm->a2 ; mm->iseg1 = nsegs + 1 ;
      if (a1 < a2)
	{
	  mm->isDown = 1 ; /* do not use TRUE since we want to minus in order */
	}
      else
	{
	  int a0 = mm->a1 ; mm->a1 = mm->a2 ; mm->a2 = a0 ; 
	  mm->isDown = 0 ;
	}
    }
  arraySort (gx->segs, segA1Order) ;
  for (ii = 0, seg = arrayp (gx->segs, ii, SEG) ; ii < nsegs ; ii++, seg++)
    {
      mm = arrp (gx->models, seg->m, MM) ;
      if (mm->iseg1 > ii) mm->iseg1 = ii ;
      mm->iseg2 = ii ; mm->ln += seg->a2 - seg->a1 + 1 ;
    }
  gx->modelsSorted = arrayHandleCopy (gx->models,gx->h) ;
  arraySort (gx->modelsSorted, mmA1Order) ;
  return 0 ;
} /* gxIntMap */

/*************************************************************************************/

static int gxCompare (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, jj, iim, jjm, u1, u2, du, ddu, np = 0 ;
  BOOL ok ;
  int nmodels = arrayMax (gx->models) ;
  SEG *seg1, *seg2 ;
  MM *mm1, *mm2 ;
  PP *pp ;
  Array pairs ;

  pairs = gx->pairs = arrayHandleCreate (10000, PP, gx->h) ;
  for (iim = 0, mm1 = arrayp (gx->modelsSorted, iim, MM) ; iim < nmodels ; iim++, mm1++)
    {
      for (jjm = iim+1, mm2 = arrayp (gx->modelsSorted, jjm, MM) ; jjm < nmodels ; jjm++, mm2++)
	{
	  if (mm2->chrom != mm1->chrom || mm2->isDown != mm1->isDown || mm2->a1 > mm1->a2)
	    break ;
	  du = 0 ; ok = FALSE ; 
	  for (ii = mm1->iseg1, seg1 = arrayp (gx->segs, ii, SEG) ; !ok && ii <= mm1->iseg2 ; ii++, seg1++)
	    if (seg1->m == mm1->m)
	      {
		for (jj = mm2->iseg1, seg2 = arrayp (gx->segs, jj, SEG) ; !ok && jj <= mm2->iseg2 ; jj++, seg2++)
		  if (seg2->m == mm2->m)
		    {
		      u1 = seg1->a1 > seg2->a1 ? seg1->a1 : seg2->a1 ;
		      u2 = seg1->a2 < seg2->a2 ? seg1->a2 : seg2->a2 ;
		      ddu = u2 - u1 + 1 ;
		      if (ddu > 0) 
			du += ddu ;
		      if (seg1->a1 == seg2->a1 &&
			  seg1->a1 > mm1->a1 && seg2->a2 > mm2->a1
			  )
			ok = TRUE ; /* same acceptor */
		      if (seg1->a2 == seg2->a2 &&
			  seg1->a2 < mm1->a2 && seg2->a2 < mm2->a2
			  )
			ok = TRUE ; /* same donor */
		    }
	      }
	  if (!ok && (3*du > mm1->ln || 3*du > mm2->ln))
	    ok = TRUE ; /* sufficient exon intersection */
	  if (ok)
	    {      /* enter the pair in the 2 orientations */
	      pp = arrayp (gx->pairs, np++, PP) ;
	      pp->m1 = mm1->m ; pp->m2 = mm2->m ;
	      pp = arrayp (gx->pairs, np++, PP) ;
	      pp->m1 = mm2->m ; pp->m2 = mm1->m ; 
	      if(0) fprintf (stderr, "pair %d %d:%d %s %s\n", np-2, pp->m1, pp->m2, dictName(gx->geneDict, pp->m1),  dictName(gx->geneDict, pp->m2)) ;
	    }
	}
    }
  arraySort (pairs, ppOrder) ;
  arrayCompress (pairs) ;
  fprintf (stderr, "gxCompare found %d matching pairs\n", arrayMax(pairs)) ;
  ac_free (h) ;
  return 0 ;
} /* gxCompare */

/*************************************************************************************/

static void gxExportGeneId (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT aoGid = aceOutCreate (gx->outFileName, ".gid.txt", gx->gzo, h) ;
  ACEOUT aoNM = aceOutCreate (gx->outFileName, ".nmid.txt", gx->gzo, h) ;
  ACEOUT aoHit = aceOutCreate (gx->outFileName, ".hit.txt", gx->gzo, h) ;
  ACEOUT aoLocus = aceOutCreate (gx->outFileName, ".locus.txt", gx->gzo, h) ;
  PP *pp ;
  Array pairs = gx->pairs ;
  int ii, iiMax = arrayMax (pairs) ;
  MM *mm1, *mm2 ;

  for (ii = 0, pp = arrp (pairs, ii, PP) ; ii < iiMax ; ii++, pp++)
    {
      mm1 = arrp (gx->models, pp->m1, MM) ;
      mm2 = arrp (gx->models, pp->m2, MM) ;
      if (mm1->gid && ! mm2->gid) 
	aceOutf (aoGid, "%s\t%s\t%s\n", dictName(gx->geneDict, pp->m2), dictName(gx->geneDict, mm1->m), dictName(gx->gidDict, mm1->gid)) ;
      if (mm1->locus && ! mm2->locus) 
	aceOutf (aoLocus, "%s\t%s\t%s\n", dictName(gx->geneDict, pp->m2), dictName(gx->geneDict, mm1->m), dictName(gx->locusDict, mm1->locus)) ;
      if (mm1->nmid && ! mm2->gid) 
	aceOutf (aoNM, "%s\t%s\n", dictName(gx->geneDict, pp->m2),  dictName(gx->gidDict, mm1->nmid)) ;
      if (mm1 < mm2) 
	aceOutf (aoHit, "%s\t%s\n", dictName(gx->geneDict, pp->m1), dictName(gx->geneDict, pp->m2)) ;
    }
  ac_free (h) ;
} /* gxExportGeneId */

/*************************************************************************************/
/* sort the pg per gid, attribute them to a genebox if they have the exact same gid set
 * get the corrdinates of these gbx
 * recursivelly atribute the otherg pg
 */
static int gx2gbx (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, iim, nnpg = 0, gene ;
  int iiMax = arrayMax (gx->pairs) ;
  int nmodels = arrayMax (gx->models) ;
  MM *mm, *mm1, *mm2 ;
  PP *pp ;
  GBX *gbx ;
  Array gx2gbx ;
  BOOL  hasChanged = TRUE ;
  /* fit the genes with geneids into boxes */

  gx2gbx = gx->gx2gbx = arrayHandleCreate (10000, GBX, gx->h) ;
  while (hasChanged)
    {
      hasChanged = FALSE ;
      
      for (iim = 0, mm = arrayp (gx->models, iim, MM) ; iim < nmodels ; iim++, mm++)
	{
	  if ((! mm->locus && ! mm->gbx) || ! mm->chrom)
	    continue ;
	  nnpg++ ;
	  if (! (gene = mm->gbx))
	    gene = mm->locus ;
	  gbx = arrayp (gx2gbx, gene, GBX) ;
	  if (gbx->chrom)
	    {
	      if (gbx->chrom != mm->chrom)
		{
		  fprintf (stderr, "ERROR locus:%s pg:%s  mm->chrom=%s gbx->chrom=%s\n"
			   , dictName (gx->locusDict, gbx->locus)
			   , dictName (gx->geneDict, mm->m)
			   , dictName (gx->chromDict, mm->chrom)
			   , dictName (gx->chromDict, gbx->chrom)
			 ) ;
		  continue ;
		}
	      if (gbx->a1 > mm->a1) { gbx->a1 = mm->a1 ; hasChanged = TRUE ; }
	      if (gbx->a2 < mm->a2) { gbx->a2 = mm->a2 ; hasChanged = TRUE ; }
	    }
	  else
	    {
	      gbx->chrom = mm->chrom ;
	      gbx->isDown = mm->isDown ;
	      gbx->a1 = mm->a1 ; gbx->a2 = mm->a2 ; 
	      hasChanged = TRUE ;
	    }
	  mm->gbx = gene ;
	  if (mm->locus && ! gbx->locus)
	    gbx->locus = mm->locus ;
	  if (mm->gid && ! gbx->gid)
	    gbx->gid = mm->gid ;
	}
      
      /* use the contact between transcripts to propagate the gene boxes */
      for (ii = 0, pp = arrp (gx->pairs, ii, PP) ; ii < iiMax ; ii++, pp++)
	{
	  mm1 = arrp (gx->models, pp->m1, MM) ;
	  mm2 = arrp (gx->models, pp->m2, MM) ;
	  if (mm1->gbx || ! mm2->gbx)
	    continue ;
	  nnpg++ ;
	  mm1->gbx = mm2->gbx ;
	  hasChanged = TRUE ; 
	  gbx = arrayp (gx2gbx, mm1->gbx, GBX) ;	
	  if (gbx->a1 > mm1->a1) { gbx->a1 = mm1->a1 ; hasChanged = TRUE ; }
	  if (gbx->a2 < mm1->a2) { gbx->a2 = mm1->a2 ; hasChanged = TRUE ; }
	}
    }

  fprintf (stderr, "gx2gbx constructed %d gene-boxes containing %d pg\n", arrayMax(gx2gbx),nnpg) ;
  ac_free (h) ;
  return 0 ;
} /* gx2gbx */

/*************************************************************************************/

static void gxExportGeneBox (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT aoGbx = aceOutCreate (gx->outFileName, ".gbx.txt", gx->gzo, h) ;
  Array gx2gbx = gx->gx2gbx ;
  int ii, iim ;
  int nmodels = arrayMax (gx->models) ;
  GBX *gbx ;
  MM *mm ;

  for (iim = 0, mm = arrayp (gx->models, iim, MM) ; iim < nmodels ; iim++, mm++)
    if (mm->m && mm->gbx)
      {
	gbx = arrp (gx->gx2gbx, mm->gbx, GBX) ;
	aceOutf (aoGbx, "Gene %s\tGenefinder %s\n"
		 , dictName(gx->locusDict, gbx->locus)
		 , dictName(gx->geneDict, mm->m)
		 ) ;
      }
  for (ii = 0, gbx = arrp (gx2gbx, ii, GBX) ; ii < arrayMax (gx2gbx) ; ii++, gbx++)
    {
      if (gbx->locus)
	aceOutf (aoGbx, "Gene %s\tGeneId %s\tIntMap %s %d %d\nSequence %s\tGenes %s %d %d\n"
		 , dictName(gx->locusDict, gbx->locus)
		 , gbx->gid ? dictName(gx->gidDict, gbx->gid) : ""
		 , dictName(gx->chromDict, gbx->chrom)
		 , gbx->isDown ? gbx->a1 : gbx->a2
		 , gbx->isDown ? gbx->a2 : gbx->a1
		 , dictName(gx->chromDict, gbx->chrom)
		 , dictName(gx->locusDict, gbx->locus)
		 , gbx->isDown ? gbx->a1 : gbx->a2
		 , gbx->isDown ? gbx->a2 : gbx->a1
		 ) ;
    }
  ac_free (h) ;
} /* gxExportGeneBox */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// pg2pg: cluster sets of gene models, transfer the geneid and NM_id\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, March 2012, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a .ace file of transcript-models mapped to the genome\n"
	    "//    cluster them and transfer the geneid, to help name the genes\n"
	    "//\n"
	    "// Syntax:\n"
	    "// The options are all optional and may be specified in any order\n"
	    "// pg2pg [options]\n"
	    "//   -o output_file_prefix\n"
	    "//   -pg2pg file_name\n"
	    "//      expect a .ace file  Sequence->Source_Exons, ->IntMap, GeneId_pg, NM_id\n"
	    "// pg2gbx [options]\n"
	    "//   -o output_file_prefix\n"
	    "//   -pg2gbx file_name\n"
	    "//      expect a .ace file  Sequence->Source_Exons, ->IntMap, GeneId_pg, NM_id\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  geneindex -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argx, const char **argv)
{
  GX gx ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&gx, 0, sizeof (GX)) ;
  gx.h = h ;

  /* optional arguments */
  if (argx == 1 ||
      getCmdLineOption (&argx, argv, "-h", 0) ||
      getCmdLineOption (&argx, argv, "-help", 0) ||
      getCmdLineOption (&argx, argv, "--help", 0)
      )
    usage (0) ;

  /* input output redirections */
  gx.gzi = getCmdLineOption (&argx, argv, "-gzi", 0);
  gx.gzo = getCmdLineOption (&argx, argv, "-gzo", 0);
  getCmdLineOption (&argx, argv, "-pg2pg", &gx.pg2pgFileName) ;
  getCmdLineOption (&argx, argv, "-pg2gbx", &gx.pg2gbxFileName) ;
  getCmdLineOption (&argx, argv, "-o", &gx.outFileName) ;

  if (argx > 1)
    usage (messprintf ("COMMAND LINE ERROR: Unknown argument %s,", argv[1])) ;

  fprintf (stderr, "// geneindex start: %s\n", timeShowNow()) ;
  /* Input is an ace file Run->Deep */

  gxParseInit (&gx) ;
  if (gx.pg2gbxFileName)
    {
      gxAceParse (&gx, gx.pg2gbxFileName) ;
      gxIntMap (&gx) ;
      gxCompare (&gx) ; /* create the sets of connected pairs */ 
      gx2gbx (&gx) ; /* create a genbox surrounding the pg sharing a set of geneid */ 
      gxExportGeneBox (&gx) ; /* export the pg<-->gene-box table */
    }
  else if (gx.pg2pgFileName)
    {
      gxAceParse (&gx, gx.pg2pgFileName) ;
      gxIntMap (&gx) ;
      gxCompare (&gx) ; /* create the sets of connected pairs */ 
      gxExportGeneId(&gx) ; /* export the pg<-->geneid table */
    }

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
