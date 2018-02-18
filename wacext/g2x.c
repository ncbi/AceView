/*  File: g2x.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2014
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
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
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */
#define ARRAY_CHECK
#include <ac.h>

/* a data sructure to manage run/gene/expression values for sets of runs */

typedef struct gxStruct {
  AC_HANDLE h ;
  Array g2rs, r2gs ;
  DICT *geneDict, *runDict ;
  const char *root, *inFile, *aceInFile, *runName, *getGene ;
  const char *outFileName, *inFileName ;
  BOOL gzi, gzo ;
  BOOL init, save, test, dirty ;
} GX ;

/*************************************************************************************/

static void gxAceParse (GX *gx, const char * fileName)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (fileName, FALSE, h) ;
  int gene, run, state = 0 ;
  float x = 0 ;
  char *cp ;
  Array r2g, g2r ;

  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp)
	{
	  state = 0 ;
	  continue ;
	}
      if ((*cp == '/' && cp[1] == '/'))
	continue ;
      switch (state)
	{
	case 0:
	  if (strcmp (cp, "Gene"))
	    continue ;
	  cp = aceInWord (ai) ;
	  if (! cp)
	    continue ;
	  dictAdd (gx->geneDict, cp, &gene) ;
	  state = 1 ;
	  continue ;
	case 1:
	  if (strcmp (cp, "Run_U"))
	    continue ;
	  cp = aceInWord (ai) ;
	  if (! cp) continue ;
	  dictAdd (gx->runDict, cp, &run) ;
	  if (aceInFloat (ai, &x))
	    {
	      int gMax = dictMax (gx->geneDict) ;
	      int rMax = dictMax (gx->runDict) ;
	      unsigned short y ;

	      r2g = array (gx->r2gs, run, Array) ;
	      if (! r2g)
		r2g =  array (gx->r2gs, run, Array) = arrayHandleCreate (gMax + 1, unsigned short, gx->h) ;
	      x *= 100 ;
	      if (x < 0) x = 0 ; if (x > 4000) x = 4000 ; /* confined to 12 bits */
	      y = x ;
	      array (r2g, gene, unsigned short) = y ;

	      g2r = array (gx->g2rs, gene, Array) ;
	      if (! g2r)
		g2r =  array (gx->g2rs, gene, Array) = arrayHandleCreate (rMax + 1, unsigned short, gx->h) ;
	      array (g2r, run, unsigned short) = y ;
	      
	      gx->dirty = TRUE ;
	    }
	  continue ;
	}  
    }
  ac_free (h) ;
} /* gxAceParse */

/*************************************************************************************/

static void gxParse (GX *gx, const char * fileName)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (fileName, FALSE, h) ;
  int gene, run ;
  float x = 0 ;
  char *cp, buf[1024] ;
  Array r2g, g2r ;
  int gMax = dictMax (gx->geneDict) ;

  dictAdd (gx->runDict, gx->runName, &run) ;
  r2g = array (gx->r2gs, run, Array) ;
  if (! r2g)
    r2g =  array (gx->r2gs, run, Array) = arrayHandleCreate (gMax + 1, unsigned short, gx->h) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;
      strncpy (buf, cp, 1000) ;
      aceInStep (ai, '\t') ;
      if (aceInFloat (ai, &x))
	{
	  unsigned short y ;
	  
	  dictAdd (gx->geneDict, buf, &gene) ;
	  
	  x = 100 * x  ;
	  if (x < 0) x = 0 ; if (x > 4000) x = 4000 ; 
	  y = x ;
	  array (r2g, gene, unsigned short) = y ;
	  
	  g2r = array (gx->g2rs, run, Array) ;
	  array (g2r, run, unsigned short) = y ;
	  
	  gx->dirty = TRUE ;
	}  
    }

  ac_free (h) ;
} /* gxParse */

/*************************************************************************************/

static int gxGetGene (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int gene, run ;
  int nn = 0, rMax = dictMax (gx->runDict) ;
  Array histo = arrayHandleCreate (256, int, h) ;
  Array g2r ;
  ACEOUT ao = aceOutCreate (gx->outFileName, "", gx->gzo, h) ;

  if (dictFind (gx->geneDict, gx->getGene, &gene) &&
      (g2r = array (gx->g2rs, gene, Array))
      )
    {
      if (arrayMax (g2r) < rMax)
	rMax = arrayMax (g2r) ;
      for (run = 1 ; run <= rMax ; run++)
	{
	  float x = arr (g2r, run, unsigned short) ;
	  int ii = (x+24)/50, *ip ;
	  ip = arrayp (histo, ii, int) ; (*ip)++ ;
	  x = x/100.0  ;
	  aceOutf (ao, "%s\t%s\t%.2f\n"
		   , dictName (gx->geneDict, gene)
		   , dictName (gx->runDict, run)
		   , x
		   ) ;
	  nn++ ;
	}
    }

  if (1)
    {
      int i, *ip, ok = 0 ;
      
      for (i=1, ip = arrp (histo, i, int) ; i <= arrayMax (histo) ; i++, ip++)
	{
	  if (*ip > 0) ok = 1 ;
	  if (ok)
	    aceOutf (ao, "%.1f\t%d\n", i/2.0, *ip) ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* gxGetGene */

/*************************************************************************************/

static void gxSave (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao ;
  int gMax, rMax, run, gene ;
  Array r2g, g2r ;

  ao = aceOutCreate (hprintf (h, "%s/geneList", gx->root), "", FALSE, h) ;
  gMax = dictMax (gx->geneDict) ;
  for (gene = 1 ; gene <= gMax ; gene++)
    aceOutf (ao, "%s\n", dictName (gx->geneDict, gene)) ;
  ac_free (ao) ;

  ao = aceOutCreate (hprintf (h, "%s/runList", gx->root), "", FALSE, h) ;
  rMax = dictMax (gx->runDict) ;
  for (run = 1 ; run <= rMax ; run++)
    aceOutf (ao, "%s\n", dictName (gx->runDict, run)) ;
  ac_free (ao) ;

  ao = aceOutCreate (hprintf (h, "%s/g2rs", gx->root), "", FALSE, h) ;
  for (gene = 1 ; gene <= gMax ; gene++)
    {
      g2r = array (gx->g2rs, gene, Array) ;
      if (! g2r)
	g2r =  array (gx->g2rs, gene, Array) = arrayHandleCreate (rMax + 1, unsigned short, gx->h) ;
      if (arrayMax (g2r) < rMax + 1)
	array (g2r, rMax, unsigned short) = 0 ;
      if (arrayMax (g2r) > rMax + 1)
	arrayMax (g2r) = rMax + 1 ;
      aceOutBinary (ao, (void *)arrp (g2r, 0, unsigned short), 2*rMax) ;
      if (0)
	fprintf (stderr, "...\t%s\t%s\t%.1f\n"
		 , dictName (gx->geneDict, gene)
		 , dictName (gx->runDict, 3)
		 , arr (g2r, 3, unsigned short)/10.0
		 ) ;	       
    }
  ac_free (ao) ;

  ao = aceOutCreate (hprintf (h, "%s/r2gs", gx->root),  "", FALSE, h) ;
  for (run = 1 ; run <= rMax ; run++)
    {
      r2g = array (gx->r2gs, run, Array) ;
      if (! r2g)
	r2g =  array (gx->r2gs, run, Array) = arrayHandleCreate (gMax + 1, unsigned short, gx->h) ;
      if (arrayMax (r2g) < gMax + 1)
	array (r2g, gMax, unsigned short) = 0 ;
      if (arrayMax (r2g) > rMax + 1)
	arrayMax (r2g) = gMax + 1 ;
      aceOutBinary (ao, (void *)arrp (r2g, 0, unsigned short), 2*gMax) ;
    }
  ac_free (ao) ;

  gx->dirty = FALSE ;
  ac_free (h) ;
  return ;
}

/*************************************************************************************/
static BOOL filCheck (const char *cp, const char *type)
{
  return ! access (cp, (*type == 'w' ? W_OK : R_OK)) ;
}

static void gxInitFromTables (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  int gene, run, gMax, rMax ;
  Array g2r, r2g ;
  char *cp ;

  cp = hprintf (h, "%s/runList", gx->root) ;
  if (filCheck (cp, "r"))
    ai = aceInCreate (cp, FALSE, h) ;
  if (ai)
    {
      while (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (cp && *cp) dictAdd (gx->runDict, cp, 0) ;
	}
      ac_free (ai) ;
    }

  cp = hprintf (h, "%s/geneList", gx->root) ;
  if (filCheck (cp, "r"))
    ai = aceInCreate (cp, FALSE, h) ;
  if (ai)
    {
      while (aceInCard (ai))
	{
	  cp = aceInWord (ai) ;
	  if (cp && *cp) dictAdd (gx->geneDict, cp, 0) ;
	}
      ac_free (ai) ;
    }

  rMax = dictMax (gx->runDict) ;  gMax = dictMax (gx->geneDict) ;

  cp = hprintf (h, "%s/g2rs", gx->root) ;
  if (filCheck (cp, "r"))
    ai = aceInCreate (cp, FALSE, h) ;
  if (ai)
    {
      for (gene = 1 ; gene <= gMax ; gene++)
	{
	  g2r = array (gx->g2rs, gene, Array) = arrayHandleCreate (rMax + 1, unsigned short, gx->h) ;
	  array (g2r, rMax, unsigned short) = 0 ; /* make room */
	  aceInBinary (ai, (void *)arrp (g2r,0,unsigned short), 2*rMax) ;
	}
      ac_free (ai) ;
    }

  cp = hprintf (h, "%s/r2gs", gx->root) ;
  if (filCheck (cp, "r"))
    ai = aceInCreate (cp, FALSE, h) ;
  if (ai)
    {
      for (run = 1 ; run <= rMax ; run++)
	{
	  r2g = array (gx->r2gs, run, Array) = arrayHandleCreate (gMax + 1, unsigned short, gx->h) ;
	  array (r2g, gMax, unsigned short) = 0 ; /* make room */
	  aceInBinary (ai, (void *) arrp (r2g,0,unsigned short), 2*gMax) ;
	}
      ac_free (ai) ;
    }
  return ;
}  /* gxInitFromTables */

/*************************************************************************************/

static void gxInit (GX *gx) 
{
  gx->runDict = dictHandleCreate (10000, gx->h) ;
  gx->geneDict = dictHandleCreate (100000, gx->h) ;
  gx->g2rs = arrayHandleCreate (1024, Array, gx->h) ;
  gx->r2gs = arrayHandleCreate (1024, Array, gx->h) ;

  if (gx->root)
    gxInitFromTables (gx) ;
    
  return ;
}  /* gxInit */

/*************************************************************************************/

static void gxTest (GX *gx) 
{
  gx->runDict = dictHandleCreate (10000, gx->h) ;
  gx->geneDict = dictHandleCreate (100000, gx->h) ;
  gx->g2rs = arrayHandleCreate (1024, Array, gx->h) ;
  gx->r2gs = arrayHandleCreate (1024, Array, gx->h) ;

  if (gx->init)
    gxInitFromTables (gx) ;
    
  return ;
}  /* gxTest */

/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: clipalign  -t target_fasta_file -i probe_fasta_file [-errMax] -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  clipalign -p tags.fasta -t chromX.fasta -errMax 2\n"
	   "// -silent : suppress title lines and status reports from the output, just report the hits\n"
	   
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char *argv[])
{
  GX gx ;
  char commandBuf[1024] ;

  freeinit () ; 

  memset (&gx, 0, sizeof (GX)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  gx.h = ac_new_handle () ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  gx.test = getCmdLineOption (&argc, argv, "-test", 0) ;
  gx.init = getCmdLineOption (&argc, argv, "-init", 0) ;
  gx.save = getCmdLineOption (&argc, argv, "-save", 0) ;

  gx.gzi = getCmdLineOption (&argc, argv, "-gzi", 0) ;
  gx.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;
  getCmdLineOption (&argc, argv, "-i", &gx.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &gx.outFileName) ;

  getCmdLineOption (&argc, argv, "-db", &gx.root) ;
  getCmdLineOption (&argc, argv, "-run", &gx.runName) ;
  getCmdLineOption (&argc, argv, "-aceIn", &gx.aceInFile) ;
  getCmdLineOption (&argc, argv, "-getGene", &gx.getGene) ;

  gxInit (&gx) ;
 
  if (gx.inFile)
    gxParse (&gx, gx.inFile) ;
  if (gx.aceInFile)
    gxAceParse (&gx, gx.aceInFile) ;
  if (gx.getGene)
    gxGetGene (&gx) ;
  if (gx.dirty && gx.save) 
    gxSave (&gx) ;

  if (gx.test)
    gxTest (&gx) ;

  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
