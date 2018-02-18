/*  File: alibaba.c
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
#include "vtxt.h"
#include "dict.h"
#include "freeout.h"
#include "keyset.h"
#include <errno.h>

static BOOL debug = FALSE ;
static void usage (void) ;
/*************************************************************************************/
/***************************************************************************************/

typedef struct chipStruct { int id, type, m1, m2, xx1, xx2, len, dna1, dna2, nClone, nCloneTotal, nMrna, nMrnaTotal ; } CHIP ;

static void showChip (Array chips)
{
  CHIP *chip ;
  int i ;

  if (!debug)
    return  ;

  if (chips && arrayMax(chips))
    {
      printf ("//id\ttype\tlen\txx1\txx2\tnClone\tnMrna\t\n") ;
      for (i = 0 ; i < arrayMax (chips) ; i++)
	{
	  chip = arrp (chips, i, CHIP) ;
	  if (!chip->len)
	    continue ;
	  printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
		  , chip->id, chip->type, chip->len, chip->xx1, chip->xx2, chip->nClone, chip->nMrna) ;
	}
    }
} /* showChip */

static int chipOrder (const void *a, const void *b)
{
  const CHIP *ca = (const CHIP *)a, *cb = (const CHIP *) b; 
  int nn ;

  if (ca->len && cb->len)
      {
	  nn = cb->nClone - ca->nClone ; /* most clones supporting this exon */
	  if (nn) return  nn ;
	  
	  nn = cb->nMrna - ca->nMrna ; /* most alternative forms use this exon */
	  if (nn) return  nn ;
      }
  nn = cb->len - ca->len ; /* longest exon */
  if (nn) return  nn ;

  nn = cb->id - ca->id ; /* longest exon */
  return nn ;
}

#ifdef JUNK
static int chipOrderXX1 (const void *a, const void *b)
{
  const CHIP *ca = (const CHIP *)a, *cb = (const CHIP *) b; 
  int nn ;

  nn = ca->type - cb->type ; /* natural exon order */
  if (nn) return  nn ;
  nn = ca->xx1 - cb->xx1 ; /* natural exon order */
  return  nn ;
}
#endif
/*************************************************************************************/

static BOOL chipExport (AC_OBJ gene, AC_OBJ oCosmid, CHIP *chip, int zV)
{
#ifdef JUNK
  static int first = 1 ;
  int a1, a2, u1, u2, ir, iMrna, iChip ;
  char *txt ;
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ oMrna, tg = ac_tag_obj (gene, "Transcribed_gene", h) ;
  AC_TABLE gmSplicing, gMrnas = ac_tag_table (tg, "mRNA", h) ;

  if (gMrnas)
    for (iMrna = 0 ; iMrna < gMrnas->rows ; iMrna++)
      {
	oMrna = ac_table_obj (gMrnas, iMrna, 0, h) ;
	a1 = ac_table_int (gMrnas, iMrna, 1, 0) ; /* mrna coord in the gene */
	a2 = ac_table_int (gMrnas, iMrna, 2, 0) ;
	gmSplicing = ac_tag_table (oMrna, "Splicing", h) ;
	for (ir = 0 ; gmSplicing && ir < gmSplicing->rows ; ir++)
	  {
	    txt = ac_table_printable (gmSplicing, ir, 4, 0) ; 
	    if (!txt || ! strstr (txt, "xon"))
	      continue ;
	    u1 = a1 - 1 + ac_table_int (gmSplicing, ir, 0, 0) ; 
	    u2 = a1 - 1 + ac_table_int (gmSplicing, ir, 1, 0) ;
	    for (chip = arrayp (chips, 0, CHIP), iChip = 0 ; 
		 iChip <  arrayMax(chips) ; iChip++, chip++)
	      {
		chip->nMrnaTotal = gMrnas->rows ;
		if (u1 <= chip->xx1 && u2 >= chip->xx2)
		  chip->nMrna++ ;
	      }
	  }
      }
  ac_free (h) ;
#else
  messcrash ("chipExport is ill defined, look in the CVS for earlier versions") ;
#endif
  return TRUE ;
} /* chipExport */

/*************************************************************************************/

static int chipGeneZone (AC_OBJ gene, int z1, int z2, int zV)
{
  CHIP *chip, *oldChip ;
  Array chips ;
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ tg = ac_tag_obj (gene, "Transcribed_gene", h) ;
  AC_OBJ oCosmid = ac_tag_obj (tg, "Genomic_sequence", h) ;
  AC_TABLE gMap = ac_tag_table (tg, "IntMap", h) ;
  AC_TABLE gCovers = ac_tag_table (tg, "Covers", h) ;
  AC_TABLE gClone, gSplicing = ac_tag_table (tg, "Splicing", h) ;
  AC_TABLE gAssembled_from = ac_tag_table (tg, "Assembled_from", h) ; 
  int ir ;
  int oldN, iChip, jChip, nMrna ;
  int nChip = 0, oldnChip, nCloneTotal  ;
  const char *txt ;
  int ga1, ga2, gm1, gm2, xx1, xx2, len ;
	

  ga1 = ac_table_int (gCovers, 0, 2, 0) ;
  ga2 = ac_table_int (gCovers, 0, 3, 0) ;

  gm1 = ac_table_int (gMap, 0, 1, 0) ;
  gm2 = ac_table_int (gMap, 0, 2, 0) ;

  if (gSplicing && gAssembled_from)
    {} /*     maxRows = (gSplicing->rows) * (gAssembled_from->rows) ;  */
  else
    goto abort ;

  chips = arrayHandleCreate (256, CHIP, h) ;
  nChip = 0 ;
  gClone = ac_tag_table (tg, "cDNA_clone", h) ;
  nCloneTotal = gClone->rows ;
  /* find all exons */
  for (ir=0, oldnChip = 0 ; ir < gSplicing->rows ; ir++)
    {
      txt = ac_table_tag (gSplicing, ir, 2, "") ; 
      if (! strstr (txt, "xon"))
	continue ;
      xx1 = ac_table_int (gSplicing, ir, 0, 0) ;
      xx2 = ac_table_int (gSplicing, ir, 1, 0) ;
      if ((z1 || z2) && (xx1 > z2 || xx2 < z1))
	continue ; /* stay in zone */
      chip = arrayp (chips, nChip++, CHIP) ;
      chip->id = nChip ;
      chip->nCloneTotal = nCloneTotal ;
	/* exon */
      /* length & coordinate on Gene */
      {  /* get dna */
	chip->len = xx2 - xx1 + 1 ;
	chip->xx1 = xx1 ; chip->xx2 = xx2 ;
	if (ga1 < ga2)
	  {
	    chip->dna1 = ga1 + xx1 - 1 ; chip->dna2 = ga1 + xx2 - 1 ; 
	  }
	else
	  { chip->dna1 = ga1 - xx1 + 1 ; chip->dna2 = ga1 - xx2 + 1 ; }

	if (gm1 < gm2)
	  { chip->m1 = gm1 + xx1 - 1 ; chip->m2 = gm1 + xx2 - 1 ; }
	else
	  { chip->m1 = gm1 - xx1 + 1 ; chip->m2 = gm1 - xx2 + 1 ; }
	chip->type = oldnChip ; /* whole exon */
      }
      oldChip = chip ;
      if (oldChip->len > 200) /* prepare smaller bits */
	{
	  chip = arrayp (chips, nChip++, CHIP) ;
	  *chip = *oldChip ;
	  chip->id = nChip ;
	  chip->len = 200 ; chip->xx2 = chip->xx1 + 199 ;
	  if (chip->dna1 < chip->dna2)
	    chip->dna2 = chip->dna1 + 199 ;
	  else
	    chip->dna2 = chip->dna1 - 199 ;
	  if (chip->m1 < chip->m2)
	    chip->m2 = chip->m1 + 199 ;
	  else
	    chip->m2 = chip->m1 - 199 ;
	  chip->type = oldnChip ; /*exon top */
	}
      if (oldChip->len > 200) /* prepare smaller bits */
	{
	  chip = arrayp (chips, nChip++, CHIP) ;
	  *chip = *oldChip ;
	  chip->id = nChip ;
	  chip->len = 200 ; chip->xx1 = chip->xx2 - 199 ;
	  if (chip->dna1 < chip->dna2)
	    chip->dna1 = chip->dna2 - 199 ;
	  else
	    chip->dna1 = chip->dna2 + 199 ;
	  if (chip->m1 < chip->m2)
	    chip->m1 = chip->m2 - 199 ;
	  else
	    chip->m1 = chip->m2 + 199 ;
	  chip->type = oldnChip ; /* exon end */
	}
      if (oldChip->len > 300) /* prepare smaller bits */
	{
	  chip = arrayp (chips, nChip++, CHIP) ;
	  *chip = *oldChip ;
	  chip->id = nChip ;
	  chip->len = 300 ; chip->xx2 = chip->xx1 + 299 ;
	  if (chip->dna1 < chip->dna2)
	    chip->dna2 = chip->dna1 + 299 ;
	  else
	    chip->dna2 = chip->dna1 - 299 ;
	  if (chip->m1 < chip->m2)
	    chip->m2 = chip->m1 + 299 ;
	  else
	    chip->m2 = chip->m1 - 299 ;
	  chip->type = oldnChip ; /*exon top */
	}
      if (oldChip->len > 300) /* prepare smaller bits */
	{
	  chip = arrayp (chips, nChip++, CHIP) ;
	  *chip = *oldChip ;
	  chip->id = nChip ;
	  chip->len = 300 ; chip->xx1 = chip->xx2 - 299 ;
	  if (chip->dna1 < chip->dna2)
	    chip->dna1 = chip->dna2 - 299 ;
	  else
	    chip->dna1 = chip->dna2 + 299 ;
	  if (chip->m1 < chip->m2)
	    chip->m1 = chip->m2 - 299 ;
	  else
	    chip->m1 = chip->m2 + 299 ;
	  chip->type = oldnChip ; /* exon end */
	}
      
      if (oldChip->len > 3000) /* prepare smaller bits */
	{
	  chip = oldChip ;
	  chip->len = 3000 ; chip->xx2 = chip->xx1 + 3000 ;
	  if (chip->dna1 < chip->dna2)
	    chip->dna2 = chip->dna1 + 3000 ;
	  else
	    chip->dna2 = chip->dna1 - 3000 ;
	  if (chip->m1 < chip->m2)
	    chip->m2 = chip->m1 + 3000 ;
	  else
	    chip->m2 = chip->m1 - 3000 ;
	  chip->type = oldnChip ; /* exon end */
	}
      
      
      /* count the supporting Clone */
      for (iChip = oldnChip ; iChip < nChip ; iChip++)
	{
	  int iClone ;
	  int u1, u2 ;
	  chip = arrayp (chips, iChip, CHIP) ;
	  for (iClone = 0 ; iClone < gAssembled_from->rows ; iClone++)
	    {
	      u1 = ac_table_int (gAssembled_from, iClone, 0, 0) ;
	      u2 = ac_table_int (gAssembled_from, iClone, 1, 0) ;
	      
	      if (u1 <= chip->xx1  && u2 >= chip->xx2)
		chip->nClone++ ;
	    }
	}
      oldnChip = nChip ;
    }
  /* count the transcripts at once for all the exons */
  if (!nChip)
    goto abort ;
  /*   chipCountMrna (tg, chips) ; missing */
  showChip (chips) ;
  arraySort (chips, chipOrder) ;
  showChip (chips) ;
  
  /* select the relevant exons */
  for (oldN = iChip = 0 ;  iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (chip->nClone < 1)
	break ;
      if (oldN > 12 && chip->nClone < oldN/2)
	break ;
      if (!oldN) oldN = chip->nClone ;
    }
  nChip = arrayMax (chips) = iChip ;
  
  /* remove very short exons */
  for (iChip = 0 ; iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (chip->len < 80)
	chip->len = 0 ;
    }
  arraySort (chips, chipOrder) ;
  showChip (chips) ;
  for (iChip = 0 ; iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (!chip->len)
	break ;
    }
  nChip = arrayMax (chips) = iChip ;
  
  /* try to see if they are next to each other  and go over introns
   * not done
   *  arraySort (chips, chipOrderXX1) ;
   */

  /* keep a single piece of each exon, the longest among the well supported */
  
  /* find max number of supported mrna */
  for (iChip = nMrna = 0 ; iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (chip->len && nMrna < chip->nMrna)
	nMrna = chip->nMrna ;
    }
  
  for (oldN = len = iChip = 0 ; iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (! chip->len)
	continue ;
      for (jChip = iChip -1 ; chip->len && jChip >= 0 ; jChip--)
	{
	  oldChip = arrayp (chips, jChip, CHIP) ;
	  if (! oldChip->len)
	    continue ;
	  xx1 = chip->xx1 > oldChip->xx1 ? chip->xx1 : oldChip->xx1 ;
	  xx2 = chip->xx2 < oldChip->xx2 ? chip->xx2 : oldChip->xx2 ;
	  if (xx1 < xx2 - 80) /* intersects keep a single one */
	    {    
	      if (xx1 < xx2 - 120 ||
		  oldChip->nClone > 2 * chip->nClone)
		chip->len = 0 ;
	    }
	}
    }
  
  arraySort (chips, chipOrder) ;
  showChip (chips) ;
  
  /* find max number of supported mrna */
  for (iChip = nMrna = 0 ; iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (chip->len && nMrna < chip->nMrna)
	nMrna = chip->nMrna ;
    }
  
  /* select exons common to many forms */
  for (oldN = len = iChip = 0 ; len < 500 && iChip < nChip ; iChip++)
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (! chip->len)
	continue ;
      if (nMrna >= 6 && chip->nMrna < nMrna/2) continue ;
      if (chipExport (gene, oCosmid, chip, zV)) /* export if goot GC/AT ratio */
	{ len += chip->len ; }
      chip->len = 0 ;
    }
  /* fall back on less representative exons */
  for (iChip = 0 ; len < 500 && iChip < nChip ; iChip++) 
    {
      chip = arrayp (chips, iChip, CHIP) ;
      if (! chip->len)
	continue ;
      if (chipExport (gene, oCosmid, chip, zV)) /* export if goot GC/AT ratio */
	{ len += chip->len ; } 
      chip->len = 0 ;
    }
  
  abort :
    ac_free (h) ;
  return nChip ;
} /* chipGeneZone */

/*************************************************************************************/
/* find the different zones with 300 aa if disjoint, or several locus id */
typedef struct zoneStruct { int z1, z2, len ; } ZONE ;

static int zoneOrder (const void *a, const void *b)
{
  const ZONE *za = (const ZONE*) a, *zb = (const ZONE*) b;
  return zb->len - za->len ; /* decreasing order */
}

static int chipGeneMultiZone (AC_OBJ gene, AC_HANDLE h)
{
  Array zones = 0 ;
  AC_OBJ mRna ;
  AC_OBJ tg = ac_tag_obj (gene, "Transcribed_gene", h) ;
  AC_TABLE gProducts, gMrnas ;
  AC_TABLE gCoding ;
  int nChip = 0, ir, ic, iz = 0, jz, m1, jr, p1, p2, b1, b2 ;
  ZONE *zz, *zz2 ;
  
  if (!tg || !(gMrnas = ac_tag_table (tg, "mRNA", h)))
    return 0 ;
  zones = arrayCreate (32, ZONE) ;
  for (ir = 0 ; ir < gMrnas->rows ; ir++)
    {
      mRna = ac_table_obj (gMrnas, ir, 0, h) ;
      m1 = ac_table_int (gMrnas, ir, 1, 0) ;
      gProducts = ac_tag_table (mRna, "Product", h) ;
      gCoding = ac_tag_table (mRna, "Coding", h) ;
      for (jr = 0 ; jr < gProducts->rows ; jr++)
	{
	  p1 = ac_table_int (gProducts, jr, 1, 0) ;
	  p2 = ac_table_int (gProducts, jr, 2, 0) ;
	  b1 = b2 = 0 ;
	  for (ic = 0 ; ic < gCoding->rows ; ic++)
	    { /* unsplice p1 */
	      if (ac_table_int (gCoding, ic, 2, 0) >= p1)
		{ b1 = ac_table_int (gCoding, ic, 0, 0) ;break ; }
	    }
	  for ( ; ic < gCoding->rows ; ic++)
	    { /* unsplice p2 */
	      if (ac_table_int (gCoding, ic, 3, 0) >= p2)
		{ b2 = ac_table_int (gCoding, ic, 1, 0) ;break ; }
	    }
	  zz = arrayp (zones, iz++, ZONE) ;
	  zz->z1 = m1 + b1 - 1 - 100 ; /* allow a piece of UTR */
	  zz->z2 = m1 + b2 - 1 + 200 ; /* allow a piece of UTR */
	  zz->len = p2 - p1 + 1 ;
	}
    }
  arraySort (zones, zoneOrder) ;
  for (iz = 0 ; iz < arrayMax (zones) ; iz++)
    {
      zz = arrayp (zones, iz, ZONE) ;
      p1 = zz->z1 ; p2 = zz->z2 ;
      for (jz = iz + 1 ; zz->len && jz < arrayMax (zones) ; jz++)
	{
	  zz2 = arrayp (zones, jz, ZONE) ;
	  if (zz2->len < 900 || (p1 < zz2->z2 && p2 > zz2->z1)) /* intersect */
	    zz2->len = 0 ;
	}
    }
  arraySort (zones, zoneOrder) ;
  for (iz = jz = 0 ; iz < arrayMax (zones) ; iz++) 
    {
      zz = arrayp (zones, iz, ZONE) ;
      if (zz->len) jz++ ;
    }
  
  if (jz > 1) /* export one chip for each zone */
    for (iz = jz = 0 ; iz < arrayMax (zones) ; iz++)
      {
	zz = arrayp (zones, iz, ZONE) ;
	if (zz->len)
	  nChip += chipGeneZone (gene, zz->z1, zz->z2, ++jz) ;
      }
  else /* export a single chip */
    nChip += chipGeneZone (gene, 0, 0, 0) ;

  arrayDestroy (zones) ;
  return nChip ;
} /* chipGeneMultiZone */

/*************************************************************************************/

static int chipGene (AC_OBJ gene)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gLL = ac_tag_table (gene, "LocusId", h) ;
  int nChip, nLL= gLL ? gLL->rows : 0 ;

  if (nLL < 0) /* i think always use multi is better, we may have 2 proteins and no locusId */
    nChip = chipGeneZone (gene, 0, 0, 0) ;
  else
    nChip = chipGeneMultiZone (gene, h) ;

  ac_free (h) ;
  return nChip ;
} /* chipGene */

/*************************************************************************************/

static void chipAnalyse (AC_DB db, const char *template, const char *keysetName) 
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ gene = 0 ;
  AC_KEYSET ks = 0 ;
  int nn = 0 ;
  
  if (keysetName)
    {
      ks = ac_read_keyset (db, keysetName, h) ;
      if (!ks || !ac_keyset_count (ks))
	{ 
	  fprintf(stderr, "cannot find keyset %s", keysetName) ;
	  ac_db_close (db) ;
	  usage () ;
	}
      iter = ac_keyset_iter (ks, TRUE, h) ;
    }
  else if (template)
    {
      iter = ac_dbquery_iter (db 
			      , messprintf ("Find gene %s", template)			    
			      , h) ;
    }
  else
    {
      iter = ac_dbquery_iter (db 
			      , "{Find gene ! Cloud_gene} $- {Find tg to_be_fused_with ; >to_be_fused_with ; > gene } "
			      , h) ;
    }
  while (ac_free (gene), gene = ac_next_obj (iter))
    {
      nn += chipGene (gene) ; 
    }
  
  printf ("// Exported %d genes\n\n", nn) ;
  fprintf (stderr, "// Exported %d genes\n\n", nn) ;
  ac_free (h) ;
  return ;
} /* chipAnalyse */

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: dnachip ACEDB [-o outfile] [-debug] [-t template] [-k keyset]\n") ;
  fprintf (stderr, "// Example:  dnachip locusid 384D8-2 \n") ;
  fprintf (stderr, "//   -o filename : redirects the output in that file, useful for batch jobs\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  int nn = 0 ;
  FILE *f = 0 ;
  const char *s = "ok" ;
  const char *outfilename = 0 ;
  const char *dbName ;
  const char *template = 0 ;
  const char *keysetName = 0 ;
  int outlevel = 0 ;
  AC_DB db ;
  
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-out", &outfilename ) ;
  /* BOOL xml =  getCmdLineOption (&argc, argv, "-x", 0) ; */
  getCmdLineOption (&argc, argv, "-t", &template) ;
  getCmdLineOption (&argc, argv, "-k", &keysetName) ;
  debug = getCmdLineOption (&argc, argv, "-debug", 0) ;
  /* read absolute args */
  dbName = argc>=2 ? argv[1] : 0 ;
  if (!dbName) 
    usage () ;

  db = ac_open_db (dbName, &s);
  if (!db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  if (outfilename)
    {
      f = filopen (outfilename, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  chipAnalyse (db, template, keysetName) ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  ac_db_close (db) ;

  fprintf (stderr, "Exported %d alignments\n", nn) ;
  sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

