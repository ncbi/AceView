/*  File: maggie.c
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
 *  are available from http://www.aceview.org/
 */

#include "../wac/ac.h"
#include "freeout.h"
#include "aceio.h"
#include "vtxt.h"
#include <errno.h>
#include "../wfiche/gtitle.h"

/*
 Exportation des alignements pour genecard et pour pierre de la grange
*/

typedef struct mmStruct { 
  AC_DB db ;
  BOOL exportSummary, exportAliTissue, pierre, html ;
  const char *template, *keysetName ;
} MM ;

static int genecardExportOneAliTissue (MM *mm, AC_OBJ tg)
{
  int a1, a2, ii ;
  const char *locusid, *locuslink, *tissue ;
  AC_OBJ chrom, gene, est ;
  AC_TABLE im, af = 0 ;
  AC_HANDLE h = handleCreate () ;
  
  im = ac_tag_table(tg, "IntMap", h) ;
  if (im && im->cols >= 3)
    {
      gene = ac_tag_obj (tg, "Gene", h) ;
      locusid = gene ? ac_tag_printable (gene, "locusId", "NULL") : "NULL" ;
      
      locuslink = gene ? ac_tag_printable (gene, "LocusLink", "NULL") : "NULL" ;
      
      chrom = ac_table_obj (im, 0, 0, h) ;
      a1 = ac_table_int (im, 0, 1, 0) ;  
      a2 = ac_table_int (im, 0, 2, 0) ;  
      af = ac_tag_table (tg, "Read", h) ;
      
      for (ii = 0 ; af && ii < af->rows ; ii++)
	{
	  printf ("AceView\t%s\t%s\t", ac_name(tg), locusid) ;
	  printf ("%s\t", locuslink) ;
	  printf ("%s\t%d\t%d\t", ac_name(chrom), a1, a2) ;
	  est = ac_table_obj (af, ii, 0, h) ;
	  tissue = ac_tag_printable (est, "Tissue", "NULL") ;
	  printf ("%s\t%s\n", ac_name(est), tissue) ;
	}
    }
  ac_free (h) ;
  return 1 ;
} /* genecardExportOneAliTissue */

static int genecardExportAliTissue (MM *mm)
{
  AC_ITER iter = 0 ; 
  AC_OBJ est ;
  int nn = 0 ;
  
  iter = ac_query_iter (mm->db, TRUE, "find tg", 0, 0) ; /* all alignments */
  while ((est = ac_next_obj (iter)))
    {
      nn += genecardExportOneAliTissue (mm, est) ;
      ac_free (est) ;
    }
  ac_free (iter) ;
  
  return nn ;
} /* genecardExportAliTissue */

/*************************************************************************************/

static int genecardExportSummary (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = 0 ; 
  AC_OBJ gene = 0 ;
  int nn = 0 ;
  unsigned char *cp ;
  ACEOUT out = 0 ;

  iter = ac_query_iter (mm->db, TRUE
			, messprintf ("find gene %s"
				      , mm->template ? mm->template : "")
			, 0, h) ; /* all genes matching template */
  
  while (ac_free (gene), (gene = ac_next_obj (iter)))
    {
      cp = ac_command (mm->db
		       , messprintf("view -v summary -c gene -n \"%s\" %s"
				    , ac_name (gene)
				    , mm->html ? "-refseq" : "-s"
				    )
		       , 0, 0) ;
      if (cp)
	{
	  nn++ ; 
	  if (mm->html)
	    {
	      out = aceOutCreateToFile (messprintf ("%s.html", ac_name(gene)), "w", 0) ;
	      aceOutf (out, "<html><head> <title>AceView gene %s</title></head>\n<body>\n", ac_name (gene)) ;
	    }
	  else
	    { 
	      out = aceOutCreateToStdout (0) ;
	      aceOutf (out, "Gene \"%s\"\n", ac_name(gene)) ;
	    }
	  aceOutf (out, "%s", cp) ;
	  aceOutf (out, "\n\n") ;
	  if (mm->html)
	    aceOutf (out, "</body></html>") ;
	  aceOutDestroy (out) ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* genecardExportSummary */

/*************************************************************************************/
/*************************************************************************************/
/* export tag a a2
              b b2
              c c2

   as
          a;b;c  \t   a2;b2;c2
*/
static int genecardExportTag (AC_OBJ gene, char *tag, int cols, AC_HANDLE h)
{
  int ir=1, jr ;
  char sep = '\t' ;
  AC_TABLE tbl = ac_tag_table (gene, tag, h) ;
  
  if (! tbl || !tbl->rows)
    for (jr = 0 ; jr < cols ; jr++, sep=';')
      freeOutf ("%c%s", sep, "-") ;
  else
    for (jr = 0 ; jr < cols ; jr++)
      for (ir = 0, sep='\t' ; ir < tbl->rows ; ir++, sep=';')
	freeOutf ("%c%s", sep, ac_table_printable (tbl, ir, jr, "-")) ;
  return ir*cols ;
} /* genecardExportTag */

/******************/
typedef struct intrStruct { int n, a1, a2, nsup ; } INTR ;

static int genecardExportPierre (MM *mm)
{
  AC_HANDLE h = 0 ;
  AC_ITER iter = 0 ; 
  AC_OBJ gene, tg = 0, mrna, product ;
  AC_TABLE mrnas, products, spl, ibs ;
  Array introns = 0 ;
  INTR *up ;
  int nn = 0, imrna, ir ;
  int nii, i1, i2, j1, j2, jr, nsup, iup ;
  const char *cp ; 
  char sep ;
  vTXT txt = vtxtCreate () ;
  vTXT txt1 = vtxtCreate () ;
  iter = ac_query_iter (mm->db, TRUE
			, messprintf ("find tg %s ; gt_ag && mrna && !shedded_from"
				      , mm->template ? mm->template : "")
			, 0, 0) ; /* all genes matching template */
  
  while (ac_free (tg), (tg = ac_next_obj (iter)))
    {
      h = ac_new_handle () ;

      introns = arrayReCreate (introns, 256, INTR) ;
      gene = ac_tag_obj (tg, "Gene", h) ;
      ibs = ac_tag_table (tg, "Intron_boundaries", h) ;

      freeOutf ("\nGENE\t%s", ac_name(tg)) ;
      if (gene)
	{
	  genecardExportTag (gene, "GeneId", 1, h) ;
	  genecardExportTag (tg, "IntMap", 3, h) ;
	  genecardExportTag (tg, "Covers", 1, h) ;
	}
      else
	freeOutf ("\t\t\t\t") ;

      spl = ac_tag_table (tg, "Splicing", h) ;
      for (ir = nii = 0 ;spl && ir < spl->rows ; ir++)
	{ 
	  cp = ac_table_printable (spl, ir, 2, "") ;
	  if (! strstr (cp, "tron"))
	    continue ;
	  i1 = ac_table_int (spl, ir, 0, 0) ;
	  i2 = ac_table_int (spl, ir, 1, 0) ;
	  cp = ac_table_printable (spl, ir, 3, "-") ;
	  up = arrayp (introns, nii++, INTR) ;
	  up->a1 = i1 ; up->a2 = i2 ; up->n = nii ;
	  freeOutf ("\nINTRON\t%d\t%s\t%d;%d", nii, cp, i1, i2) ;

	  vtxtClear (txt) ;
	  vtxtClear (txt1) ;
	  for (sep = '\t', jr = 0 , nsup = 0 ; jr < ibs->rows ; jr++)
	    {
	      j1 = ac_table_int (ibs, jr, 2, 0) ;
	      j2 = ac_table_int (ibs, jr, 3, 0) ;
	      if (j1 == i1 && j2 == i2)
		{
		  cp = ac_table_printable (ibs, jr, 4, 0) ;
		  if (cp)
		    { vtxtPrintf (txt, "%c%s", sep, cp) ; sep = ';' ; nsup++ ;}
		}
	    }
	  up->nsup = nsup ;
	  freeOutf ("\t%d", nsup) ;
	  if (nsup)
	    freeOut (vtxtPtr (txt)) ;
	}
      
      mrnas = ac_tag_table (tg, "Mrna", h) ;
      for (imrna = 0 ; imrna < mrnas->rows ; imrna++)
	{
	  int p3 = 0, p4 = 0 ;  /* coords of best product in pre-mrna */
	  int m1 = 0, m2 = 0 ;  /* coords of the mrna in the tg */
	  int iproduct ;
	  int quality = 0 ;

	  mrna = ac_table_obj (mrnas, imrna, 0, h) ;
	  m1 = ac_table_int (mrnas, imrna, 1, 0) ;
	  m2 = ac_table_int (mrnas, imrna, 2, 0) ;

	  if (! ac_has_tag (mrna, "gt_ag"))
	    continue ;
	  nn++ ;
	  freeOutf ("\nmRNA\t%s", ac_name(mrna)) ;
	  genecardExportTag (gene, "GeneId", 1, h) ;
	  freeOutf ("\t%d\t%d", m1, m2) ; /* mrna on gene */
	  if (ac_has_tag (mrna, "Valid5p")) freeOutf ("\t%d", m1) ;
	  else freeOutf ("\t-") ;
	  genecardExportTag (mrna, "valid3p", 1, h) ;
	  if (ac_has_tag (mrna, "Mrna_covered_by"))
	    genecardExportTag (mrna, "Mrna_covered_by", 1, h) ;
	  else
	    genecardExportTag (mrna, "CDS_covered_by", 1, h) ;
	  spl = ac_tag_table (mrna, "Splicing", h) ;
	  if (spl && spl->rows > 1)
	    {
	      vtxtClear (txt) ;
	      vtxtClear (txt1) ;
	      for (ir = nsup = 0, sep = '\t' ; ir < spl->rows ; ir++)
		{
		  if (strstr (ac_table_printable (spl, ir, 4, ""), "tron"))
		    {
		      i1 = ac_table_int (spl, ir, 0, 0) + m1 - 1 ;
		      i2 = ac_table_int (spl, ir, 1, 0) + m1 - 1 ;
		      nsup++ ;
		      vtxtPrintf (txt, "\t%d;%d", i1, i2) ;
		      for (up = arrp (introns, 0, INTR), iup = 0 ; iup < arrayMax (introns) ; iup++, up++)
			if (up->a1 == i1 && up->a2 == i2)
			  {
			    vtxtPrintf (txt1, "%c%d", sep, up->n) ; sep = ';' ;
			    vtxtPrintf (txt, ";%d", up->nsup) ;
			  }
		    }
		}
	      if (nsup)
		{
		  freeOutf ("\nINTRONS\t%d%s", nsup, vtxtPtr (txt1) ? vtxtPtr (txt1) : "\t-") ;
		  freeOut (vtxtPtr (txt)) ;
		}
	    }


	  products = ac_tag_table (mrna, "Product", h) ;
	  p3 = p4 = 0 ; product = 0 ;
	  for (iproduct = 0 ; products && iproduct < products->rows ; iproduct++)
	    {
	      product = ac_table_obj (products, iproduct, 0, h) ;
	      if (! ac_has_tag (product, "Good_product") ||
		  ! ac_has_tag (product, "Best_product"))
		{ product = 0 ; continue ; }
	      quality = ac_tag_int (product, "Quality", 0) ;
	      /* coords of best product in mrna */
	      /*   p1 = ac_table_int (products, iproduct, 1, 0) ; */
	      p3 = ac_table_int (products, iproduct, 3, 0) ;
	      p4 = ac_table_int (products, iproduct, 4, 0) ;
	      break ;
	    }

	  if (product)
	    {
	      AC_TABLE kz = 0, cvb = ac_tag_table (product, "Covered_by", h) ;
	      int icvb ;

	      freeOutf ("\nPROTEIN\t%d", quality) ;
	      freeOutf ("\t%d\t%d", m1 + p3 - 1, m1 + p4 - 1) ; /* product on gene */
	      if (0)
		{
		  genecardExportTag (product, "Nb_introns_in_CDS", 1, h) ;
		  genecardExportTag (product, "Nb_Introns_outside_CDS", 1, h) ;
		}
	      
	      if (ac_has_tag (product, "NH2_complete"))
		{
		  if ((kz = ac_tag_table (product, "First_Kozak", h)) &&
			   kz->cols >= 2 &&
			   ac_table_int (kz, 0,0,-1) == 1 &&
			   (cp = ac_table_printable (kz, 0,1,0))
			   )
		    freeOutf ("\t%s", cp) ;
		  else if  ((kz = ac_tag_table (product, "First_atg", h)) &&
			   ac_table_int (kz, 0,0,-1) == 1
			   )
		    freeOut ("\tATG") ;
		  else
		    freeOut ("\t-") ;
		}
	      else
		freeOut ("\tNH2_partial") ;

	      if (ac_has_tag (product, "COOH_complete"))
		freeOut ("\tStop") ;
	      else
		freeOut ("\tCOOH_partial") ;

	      for (icvb = 0, sep='\t' ; cvb && icvb < cvb->rows ; icvb++)
		{
		  if (ac_table_int (cvb, icvb, 2, -1) == 0 &&
		      ac_table_int (cvb, icvb, 4, -1) == 0
		      )
		    {
		      freeOutf ("%c%s", sep, ac_table_printable (cvb, icvb, 0, "-")) ;
		      sep = ';' ;
		    }
		}
	    }
	}
      ac_free (h) ;
    }    
  freeOut ("\n") ;

  ac_free (txt) ;
  ac_free (txt1) ;
  arrayDestroy (introns) ;
  ac_free (iter) ;
  return nn ;
} /* genecardExportPierre */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: genecardexport -db ACEDB [-summary]\n") ; 
  fprintf (stderr, "//   -t template: a query on class tg\n") ;
  fprintf (stderr, "//   -k keyset: a keyset.list of class tg\n") ;
  fprintf (stderr, "//   -summary : export summaries\n") ;
  fprintf (stderr, "//   -exportAliTissue \n") ;
  fprintf (stderr, "//   -pierre : export for Pierre de la grange\n") ;
  fprintf (stderr, "//   -out filename : redirects the output in that file, useful for batch jobs\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  const char *s = "ok" ;
  const char *outfilename = 0 ;
  const char *dbName = "" ;
  int outlevel = 0 ;
  MM mm ;
  
  freeinit () ; 
  memset (&mm, 0, sizeof (MM)) ;
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-out", &outfilename ) ;
  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  mm.html = getCmdLineOption (&argc, argv, "-html", NULL) ;
  mm.pierre =  getCmdLineOption (&argc, argv, "-pierre", NULL) ;
  mm.exportSummary = getCmdLineOption (&argc, argv, "-summary", NULL) ;
  mm.exportAliTissue = getCmdLineOption (&argc, argv, "-exportAliTissue", NULL) ;
  getCmdLineOption (&argc, argv, "-t", &(mm.template)) ;
  getCmdLineOption (&argc, argv, "-k", &(mm.keysetName)) ;

  /* check absolute args */
  if (!dbName[0]) 
    usage () ;

  mm.db = ac_open_db (dbName, &s);
  if (!mm.db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  if (outfilename)
    {
      f = filopen (outfilename, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  if (mm.pierre)
    genecardExportPierre (&mm) ;
  else if (mm.exportSummary)
    genecardExportSummary (&mm) ;
  else if (mm.exportAliTissue)
    genecardExportAliTissue (&mm) ;
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

