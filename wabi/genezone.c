/*  File: genezone.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the genetic map
 * Exported functions:
       geneZoneExport: export the fasta file of the dna zones
 * HISTORY:

 * Created: July 15 2002 (mieg)
 *-------------------------------------------------------------------
 */

/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "wac/ac.h"
#include "freeout.h"
#include "cdna.h"

static AC_DB ac_db = 0 ;

/*************************************************************/
typedef enum { GZ_NULL, GZ_EXON, GZ_ZONE, GZ_INTRON } GZ_TYPE ;
static char *gzTypeName [] = { "null", "exon", "zone", "intron"} ;

typedef struct gzStruct { KEY gene, tg, mrna ; int a1, a2, x1, x2, m1, m2, tg1, tg2, zone ; GZ_TYPE type ; } GZ ;
/* x1, x2 coords of the exon/intron  in the tg
 * a1, a2 coords of the exon intron in the cosmid (minus on up strand 
 * tg1, tg2 are the coords of the tg on the cosmid
 * m1, m2 are the coords of the exon on the spliced mrna
 * we compress using a1, a2 to handle the shedded genes
 * we fasta dump using the tg dna
 */ 
static int geneZoneOrder (const void *va, const void *vb)
{
  const GZ *a = (const GZ *)va, *b = (const GZ *)vb ;
  int nn ;
  
  nn = a->gene - b->gene ;  if (nn) return nn ;
  nn = a->type - b->type ;  if (nn) return nn ; /* exon before intron */
  nn = a->a1 - b->a1 ;  if (nn) return nn ;
  nn = a->a2 - b->a2 ;  if (nn) return nn ; 
  nn = a->mrna - b->mrna ;  if (nn) return nn ; 

  return 0 ;
} /* geneZoneOrder */

/*************************************************************/

static int geneZoneGetTg (KEY geneKey, AC_OBJ tg, Array aa)
{ 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl, mrnas, spls ;
  int ii, jj, m1, m2 = 0, tg1, tg2 ;
  const char *typ ;
  GZ *gz ;
  AC_OBJ mrna ;

  tbl = ac_tag_table (tg, "Covers", h) ;
  tg1 = tbl ? ac_table_int (tbl, 0, 2, 0) : 0 ; 
  tg2 = tbl ? ac_table_int (tbl, 0, 3, 0) : 0 ; 
  mrnas = ac_tag_table (tg, "Mrna", h) ;
  for (ii = 0 ; mrnas && ii < mrnas->rows ; ii++)
    {
      mrna = ac_table_obj (mrnas, ii, 0, h) ;
      m1 = ac_table_int (mrnas, ii, 1, 0) ; 
      m2 = ac_table_int (mrnas, ii, 2, 0) ; 
    
      if (m2)
	{      
	  spls = ac_tag_table (mrna, "Splicing", h) ;
	  for (jj = 0 ; spls && jj < spls->rows ; jj++)
	    {
	      gz = arrayp (aa, arrayMax(aa), GZ) ;
	      gz->gene =geneKey ;
	      gz->tg = ac_obj_key (tg) ;
	      gz->mrna = ac_obj_key (mrna) ;
	      gz->x1 = ac_table_int (spls, jj, 0, 0) + m1 - 1 ; 
	      gz->x2 = ac_table_int (spls, jj, 1, 0) + m1 - 1 ; 
	      gz->m1 = ac_table_int (spls, jj, 2, 0) + m1 - 1 ; 
	      gz->m2 = ac_table_int (spls, jj, 3, 0) + m1 - 1 ; 
	      if (tg1 < tg2)
		{ 
		  gz->a1 = tg1 + gz->x1 - 1 ; 
		  gz->a2 = tg1 + gz->x2 - 1 ; 
		}
	      else
		{ /* minus, so they order correctly */
		  gz->a1 = -(tg1 - gz->x1 + 1) ; 
		  gz->a2 = -(tg1 - gz->x2 + 1) ; 
		}
	      gz->tg1 = tg1 ;
	      gz->tg2 = tg2 ;
	      typ = ac_table_printable (spls, jj, 4, "toto") ; 
	      if (strstr (typ, "xon")) gz->type = GZ_EXON ; 
	      else if (strstr (typ, "tron")) gz->type = GZ_INTRON ; 
	    }
	}
    }

  ac_free (h) ;
  return arrayMax (aa) ;
} /* geneZoneGetTg */

/*************************************************************/

static int geneZoneCompress (Array aa)
{
  GZ *gza, *gzb ;
  int ii, jj, nx, ni ;

  arraySort (aa, geneZoneOrder) ;
  ni = nx = 0 ;
  for (ii = jj = 0, gzb = 0 ; ii < arrayMax (aa) ; ii++)
    {
      gza = arrp (aa, ii, GZ) ;
      if (gzb && gza->gene == gzb->gene &&
	  gza->type == gzb->type &&
	  gza->a1 == gzb->a1 &&
	  gza->a2 == gzb->a2)
	continue ;
      gzb = arrp (aa, jj++, GZ) ;
      if (gza != gzb)
	*gzb = *gza ; 
      if (gzb->type == GZ_EXON)
	gzb->zone = ++nx ;
      if (gzb->type == GZ_INTRON)
	gzb->zone = ++ni ;
    }
  arrayMax (aa) = jj ;
  return jj ;
} /* geneZoneCompress */

/*************************************************************/

static BOOL geneZoneFastaDump (KEY tg, Array aa, char *dna)
{
  int ii, nn = 0 ;
  char *cp, *cq, cc ;
  GZ *gz ;

  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      gz = arrp (aa, ii, GZ) ;
      if (gz->tg != tg)
	continue ;
      nn++ ;
      freeOutf (">%s|%s_%d\n"
		, ac_key_name (gz->gene)
		, gzTypeName [gz->type]
		, gz->zone
		) ;
      if (gz->type == GZ_EXON)
	{
	  cp = dna + gz->x1 - 1 ;
	  cq = dna + gz->x2 ;
	  cc = *cq ;
	  *cq = 0 ;
	  freeOutf ("%s\n", cp) ;
	  *cq = cc ;
	}
      else if (gz->type == GZ_INTRON)
	{
	  cp = dna + gz->x1 - 1 - 17 ;
	  cq = dna + gz->x1 - 1 ; /* first base of intron */
	  cc = *cq ;
	  *cq = 0 ;
	  freeOutf ("%s", cp) ;
	  *cq = cc ;

	  cp = dna + gz->x2 ; /* first base of second exon */
	  cq = dna + gz->x2 + 17 ;
	  cc = *cq ;
	  *cq = 0 ;
	  freeOutf ("%s\n", cp) ;
	  *cq = cc ;
	}
    }
  return nn ;
} /* geneZoneFastaDump */

/*************************************************************/

static int geneZoneExportTg (AC_OBJ tg, Array aa)
{
  int nn = 0, a1, a2 ;
  char *dna ;
  AC_OBJ cosmid ;
  AC_TABLE tbl ;
  AC_HANDLE h = ac_new_handle () ;

  cosmid = ac_tag_obj (tg, "Genomic_sequence", h) ;
  tbl = ac_tag_table (tg, "Covers", h) ;
  a1 = tbl ? ac_table_int (tbl, 0, 2, 0) : 0 ; 
  a2 = tbl ? ac_table_int (tbl, 0, 3, 0) : 0 ; 
  if (a1 && a2 && cosmid)
    {
      dna = ac_zone_dna (cosmid, a1, a2, h) ;
      geneZoneFastaDump (ac_obj_key (tg), aa, dna) ;
    }
  ac_free (h) ;
  return nn ;
} /* geneZoneExportTg */

/*************************************************************/

static int geneZoneExportGene (AC_OBJ gene)
{
  int nn = 0, itg ;
  Array aa ;
  AC_OBJ tg ;
  AC_TABLE tgs ;
  AC_HANDLE h = ac_new_handle () ;

  aa = arrayHandleCreate (64, GZ, h) ;
  tgs = ac_tag_table (gene, "Transcribed_gene", h) ;
  for (itg = 0; tgs && itg < tgs->rows ; itg++)
    {
      tg = ac_table_obj (tgs, itg, 0, h) ;
      if (! ac_has_tag (tg, "Shedded_from"))
      geneZoneGetTg (ac_obj_key (gene), tg, aa) ;
    }

  geneZoneCompress (aa) ;
  for (itg = 0; tgs && itg < tgs->rows ; itg++)
    {
      tg = ac_table_obj (tgs, itg, 0, 0) ;
      geneZoneExportTg (tg, aa) ;
    }
  
  ac_free (h) ;
  return nn ;
} /* geneZoneExportGene */

/*************************************************************/
/* Exports on freeout the fasta file of gene zones
 * and the ace file of mrna signatures
 * returns the number of exported genes
 */
int geneZoneExport (KEYSET ks)
{
  int ii, nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter ;
  KEY key ;
  AC_OBJ gene = 0 ;
  char *cp ;
  FILE *f = 0 ;
  int levelOut = 0 ;

  while ((cp = freeword ()))
    {
      if (! strcmp (cp, "-o"))
	{
	  if ((cp = freeword ()))
	    {
	      f = filopen (cp, 0, "w") ;
	    }
	}
    }
  if (!ac_db)
    {
      cDNAAlignInit () ;
      ac_db = ac_open_db (0,0) ;
    }
  if (f)
    levelOut = freeOutSetFile (f) ;
  
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      iter = ac_dbquery_iter (ac_db, messprintf ("Find gene %s", freeprotect (name(key))), h) ;
      while (ac_free (gene), gene = ac_iter_obj (iter))
	{
	  nn += geneZoneExportGene (gene) ;
	}
      ac_free (iter) ;
    }
  if (levelOut)
    freeOutClose (levelOut) ;
  filclose(f) ;
  ac_free (h) ;
  freeOutf ("Exported %d genes as fasta-zones\n",  nn) ;
  return nn ;
}  /* geneZoneExport  */

/************************************************************/
/************************************************************/
/************************************************************/

