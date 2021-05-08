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
 * This file is part of the GOLD genome database package, written by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This code works as a client against an acedb server
 *  and should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  The present code and the acedb schema for this project
 *  are available from http://www.aceview.org/GOLD
 */

#include "../wac/ac.h"
#include "vtxt.h"
#include "dict.h"
#include "freeout.h"
#include "keyset.h"
#include <errno.h>
#include "bitset.h"
#include "dna.h"
#include "peptide.h"

#define mat(_i,_j) (arrp(Cost,(_i - x1 + 1)*JMAX + (_j - a1 + 1), int))
static void gChromStats (AC_DB db, BOOL xml, char *template) ;
static void gLibraryStats (AC_DB db, BOOL xml) ;
static void gExportTable (AC_DB db, AC_KEYSET ks, char *dirTitle, int method, char *fileTitle, char *methodTag) ;
static int gExportTable2 (AC_DB db, AC_KEYSET ks0,  char *myQuery, char *fileTitle, char *tableNam) ;
static void gExportTable3 (AC_DB db, AC_KEYSET ks, char *dirTitle,  char *fileTitle, int method, int dummy, char *methodTag) ;
static void gExportTable4 (AC_DB db, AC_KEYSET ks, char *dirTitle, char *library, char *defect) ;
static void gExportTable5 (AC_DB db, AC_KEYSET ks, char *dirTitle, char *library, char *defect) ;
static void gExportTable6 (AC_DB db, AC_KEYSET ks, char *dirTitle, char *library, char *defect) ;
static void usage (void) ;

typedef struct { 
  /* data imported directly */
  AC_OBJ Est, Gene ;
  int iMethod, est, ln, pA, pfA, clipTop, clipEnd ;  /* from the est object */
  int start, ali, err, orf ; /* est->from_gene */
  int gene, method, chrom, a1, a2, x1, x2, covers, nexons ; /* from the tg */
  int c1, c2 ; /* eventual biggest cdna gap */

  /* data evaluated dynamically on one alignment */
  int score, best_score, bad_score, no_score ;
  int orf_bonus ; /* 3 pt for best length if this changes the choice of gold */
  int m_split, m_score, m_ali, m_err ;
  int g_split, g_score, g_ali, g_err ;
  int   /* all this  only if v_repeat == 0 */
    masked, /* this gene is part of a hyper variable region */
    excellent,  /* ali >= 99% && err <= .4% */
    good,  /* ali >= 99% && .4% < err <= 1.5% */
    partial, 
    partial1,      /* 99% > ali >= 60% && err <= .4% */
    partial2,      /* 99% > ali >= 60% && err > .4% && err <= 1.5% */
    dubious, 
    dubious1,  /* ali >= 99% && err > 1.5% && err <= 3.0% */
    dubious2,  /* 99% > ali >= 60% && err > 1.5% && err <= 3.0% */
    dubious3,  /* ali < 60% && err <= .4% */
    dubious4,  /* ali < 60% && err > .4% && err <= 1.5% */
    bad,    
    bad1,  /* ali >= 99% && err > 3.0% */
    bad2,  /* 99% > ali >= 60% && err > 3.0% */
    bad3,  /* ali < 60% && err > 1.5% && err <= 3.0% */
    bad4,  /* ali < 60% && err > 3.0% */
    bad_short,    /* ali < 60% && ali < 200bp  && err < 3.9 : WARNING a bizare subset of bad Union dubious */
    unali ;   /* */


  /* comparison of all alignments */
  int best, best5, best50, best500, best5000 ; /* compared to other ali */
  int local_best, m_best, alibaba, diamond ;   /* properties of the best */
  int kill_method ; /* to fuse methods A and B if needed */
  /* bugs in the individual alignments */
  int twin, dubious_strand, wrong_strand ;
  int micro_deletion, v_repeat, v_pileUp, mosaic, inversion, fuse_with, transposition, ambiguous ;
  int no_polya_signal, genomic_contamination ;
  int bad_intron, polyA_exon, polyA_first_exon, vector_exon, too_short_first_exon ;
  int g_major, g_minor, m_major, m_minor ;
  int wrong_map, dubious_map_b, dubious_map_c ; /* , amb ; */

  int repetition, repetition5, repetition50 ;
  int genome_deletion ;
  int identical_in_n, alternate ; 
  int too_long,  missed_exon, more_exon, wrong_exon ;
  int missed_5p, missed_central, missed_3p ;
  int missed_intron, more_intron ;
  int missed_gling, more_gling ;
  int more_5p, more_central, more_3p ;
  int identical_exon, exon_tested ;
  int Duplicate_of_read ;
  BOOL g_forward, g_reverse ;
} STAT ;

static BOOL debug = FALSE ;
static int gIntronsTg (AC_OBJ tg) ;

/*****************************************************************************/
/*************************************************************************************/

static void gOutSection (BOOL xml, char *cp)
{
  if (xml)
    freeOutf ("<h3>%s</h3><br>\n", cp) ;
  else
    freeOutf ("%s\n", cp) ;
}

/*************************************************************************************/

static void gOutTableStart (BOOL xml, char *title)
{
  if (xml)
    freeOutf ("\nTABLE_START %s\n<table border=1>\n", title) ;
  else
    freeOutf ("\n\n") ;
}

/*************************************************************************************/

static void gOutTableEnd (BOOL xml)
{
  if (xml)
    freeOutf ("\n</table><br>\nTABLE_END\n") ;
  else
    freeOutf ("\n\n") ;
}

/*************************************************************************************/

static void gOutText (BOOL xml,const  char *cp)
{
  if (xml)
    freeOutf ("<td>%s</td>", cp) ;
  else
    freeOutf ("\t%s", cp) ;
}

/*****************************************************************************/

static void gOutNewLine (BOOL xml)
{
  if (xml)
    freeOutf ("</tr>\n") ;
  else
    freeOutf ("\n") ;
}

/*****************************************************************************/

static void gOutBreak (BOOL xml)
{
  if (xml)
    freeOutf ("<br>\n") ;
  else
    freeOutf ("\n") ;
}

/*************************************************************************************/
static int tableLine = 0 ;
static void gOutTitle (BOOL xml, const char *cp)
{
  if (xml)
    {
      if ((tableLine++) % 2)
	freeOutf ("<tr bgcolor=#efefff>") ;
      else
	freeOutf ("<tr bgcolor=#ffffff>") ;
      freeOutf ("<td>%s</td>", cp) ;
    }
  else
    freeOutf ("%-15s", cp) ;
}

/*************************************************************************************/

static void gOutTitles (BOOL xml, const char *cp0)
{
  char buff[1000] ;
  char *cp = buff, *cq ;
  int n = 0 ;
  tableLine = 0 ;
  if (xml)
    freeOutf ("<tr bgcolor=#afafff>") ;

  strcpy (buff, cp0) ; /* on some machine i cannot alter a fixed string */
  while ((cq = strstr (cp, "\t")) && *(cq+1))
    {
      *cq = 0 ;
      if (!n++ && !xml)
	 freeOutf ("%-15s", cp) ;
      else
	gOutText (xml, cp) ;
      *cq = '\t' ;
      cp = cq + 1 ;
    }
  gOutText (xml, cp) ;
  gOutNewLine (xml) ;
}

/*****************************************************************************/

static void gOutFloat (BOOL xml, float x)
{
  if (xml)
    freeOutf ("<td>%.2f</td>", x) ;
  else
    freeOutf ("\t%.2f", x) ;
}

/*****************************************************************************/

static void gOutInt (BOOL xml, int nn, int percent, int total) 
{
  if (percent == 2)
    {
      if (xml)
	freeOutf ("<td>%.2f</td>", (100.0 * nn)/ total) ;
      else
	freeOutf ("\t%.2f", (100.0 * nn)/ total) ;
    }
  else if (percent == 1)
    {
      if (xml)
	freeOutf ("<td>%.2f%%</td>", (100.0 * nn)/ total) ;
      else
	freeOutf ("\t%.2f%%", (100.0 * nn)/ total) ;
    }
  else
    {
      if (xml)
	freeOutf ("<td>%5d</td>", nn) ;
      else
	freeOutf ("\t%5d", nn) ;
    }
}

/*****************************************************************************/

static void gOutIntLink2 (BOOL xml, int nn, BOOL percent, int total, char *dirTitle, char *fileTitle) 
{
  if (percent)
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%s.%s.%d.htm\">%.2f%%</a></td>", dirTitle, dirTitle, fileTitle, nn, (100.0 * nn)/ total) ;
	  else	
	    freeOutf ("\t0.00%%") ;
	}
      else
	freeOutf ("\t%.2f%%", (100.0 * nn)/ total) ;
    }
  else
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%s.%s.%d.htm\">%d</a></td>", dirTitle, dirTitle, fileTitle, nn, nn) ;
	  else
	    freeOutf ("<td>%d</td>", nn) ;
	}
      else
	freeOutf ("\t%d", nn) ;
    }
}

static void gOutIntLink3 (BOOL xml, int nn, BOOL percent, int total, const char *dirTitle, const char *library, const char *defect) 
{
  if (percent)
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%s/%s.%d.htm\">%.2f%%</a></td>", dirTitle, library, defect, nn, (100.0 * nn)/ total) ;
	  else	
	    freeOutf ("\t0.00%%") ;
	}
      else
	freeOutf ("\t%.2f%%", (100.0 * nn)/ total) ;
    }
  else
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%s/%s.%d.htm\">%d</a></td>", dirTitle, library, defect, nn, nn) ;
	  else
	    freeOutf ("<td>%d</td>", nn) ;
	}
      else
	freeOutf ("\t%d", nn) ;
    }
}

/*****************************************************************************/

static void gOutIntLink (BOOL xml, int nn, BOOL percent, int total, int method, char *dirTitle, char *fileTitle) 
{
  if (percent)
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%c.%s.%d.htm\">%.2f%%</a></td>", dirTitle, method, fileTitle, nn, (100.0 * nn)/ total) ;
	  else
	    freeOutf ("<td>%.2f%%</td>", 0.0) ;
	}
      else
	freeOutf ("\t%.2f%%", (100.0 * nn)/ total) ;
    }
  else
    {
      if (xml)
	{
	  if (nn)
	    freeOutf ("<td><a href=\"%s/%c.%s.%d.htm\">%d</a></td>", dirTitle, method, fileTitle, nn, nn) ;
	  else
	    freeOutf ("<td>%d</td>", nn) ;
	}
      else
	freeOutf ("\t%d", nn) ;
    }
}

/*************************************************************************************/
/******************************** FLAG ***********************************************/
/*************************************************************************************/

/* genebox touch on EST */
static BOOL gSameIntrons (STAT *sp, STAT *sq)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE s1, s2 ;
  int ir, x1, x2, y1, y2; 
  BOOL ok = TRUE ;
  int dp, d2 ;

  dp = sp->a2 - sq->a1 ;
  if (dp > 0) d2 = sq->a1 - sp->a1 ;
  else  d2 = sp->a1 - sq->a1 ;

  s1 = ac_tag_table (sp->Gene, "Splicing", h) ;
  s2 = ac_tag_table (sq->Gene, "Splicing", h) ;
  if (s1->rows != s2->rows)
    ok = FALSE ;
  for (ir = 0 ; ok && ir < s1->rows ; ir++)
    {
      if (!strstr (ac_table_printable (s1, ir, 2, ""), "ntron"))
	continue ;
      x1 = ac_table_int (s1, ir, 0, 0) ;
      x2 = ac_table_int (s1, ir, 1, 0) ;
      y1 = d2 + ac_table_int (s2, ir, 0, 0) ;
      y2 = d2 + ac_table_int (s2, ir, 1, 0) ;
      if (x1 != y1 || x2 != y2)
	ok = FALSE ;
    }
  ac_free (h) ;
  return ok ;
} /* gSameInttrons */

/*************************************************************************************/
/* s1=gold, get s2->missed_exon/wrong_exon/more_exon/missed_5p/missed_central/missed_3p */
static void gExonAnalysis (STAT *s1, STAT *s2)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE t1, t2 ;
  int icr, x1, x2, old1, old2, max1, max2, i1, i2, j1, j2, in_intron1, in_intron2 ;
  int dp, d2 ;
  BOOL gling1 = FALSE, gling2 = FALSE ;

  if (s1->chrom != s2->chrom)
    goto done ;
  dp = s2->a2 - s1->a1 ;
  if (dp > 0) d2 = s2->a1 - s1->a1 ;
  else  d2 = s1->a1 - s2->a1 ;
  
  t1 = ac_tag_table (s1->Gene, "Splicing", h) ;
  t2 = ac_tag_table (s2->Gene, "Splicing", h) ;
  if (!t1->rows || !t2->rows)
    goto done ;

  i1 = i2 = 0 ; /* line number in t tables */
  j1 = j2 = 0 ; /* exon number */
  in_intron1 = in_intron2 = 0 ; /* in_intron */
  old1 = old2 = d2 > 0 ? -99999 : d2 -9999 ; /* previous position */
  max1 = t1->rows - 1 ;
  max2 = t2->rows - 1 ; /* last possible i value */
  
  s1->exon_tested = s2->exon_tested = 1 ;
  while (1)
    {
      /* check for correct line and end of gene */
      while (i1 < max1 && !strstr (ac_table_printable (t1, i1, 2, ""), "xon"))
	i1++ ;
      while (i2 < max2 && !strstr (ac_table_printable (t2, i2, 2, ""), "xon"))
	i2++ ;
      if (i1 > max1) /* reached end of gene 1 */
	{
	  if (in_intron2)
	    { i2++ ; in_intron2 = 0 ; continue ; }
	  if (i2 <= max2) /* else both are done */
	    {
	      s2->more_3p++ ;
	      s2->more_exon++ ; /* extra 3p exon */ 
	      i2++ ;
	      continue ;
	    }
	  goto done ;
	}
      if (i2 > max2) /* reached end of gene 2 */
	{
	  if (in_intron1)
	    { i1++ ; in_intron1 = 0 ; continue ; }
	  if (i1 <= max1) /* else both are done */
	    {
	      s2->missed_exon++ ; s2->missed_3p++ ;
	      i1++ ;
	      continue ; 
	    }
	  goto done ;
	}
      /* jump the gling=gling, treat them as a single exon , and slide to the end of the sceond exon */
      while (in_intron1 && i1 + 2 <= max1 && 
	  ac_table_int (t1, i1 + 2, 0, 0) < ac_table_int (t1, i1, 1, 0) + 40 &&
	  ac_table_int (t1, i1 + 2, 0, 0) > ac_table_int (t1, i1, 1, 0) - 40)
	{ i1 += 2 ; gling1 = TRUE ; }
      while (in_intron2 && i2 + 2 <= max2 && 
	  ac_table_int (t2, i2 + 2, 0, 0) < ac_table_int (t2, i2, 1, 0) + 40 &&
	  ac_table_int (t2, i2 + 2, 0, 0) > ac_table_int (t2, i2, 1, 0) - 40)
	{ i2 += 2 ; gling2 = TRUE ; }
      /* get and compare coordinates */
      /* flag as we leave the exon */
      x1 = ac_table_int (t1, i1, in_intron1, 0) ;
      x2 = d2 + ac_table_int (t2, i2, in_intron2, 0) ;
      switch (16 * in_intron1 + in_intron2)
	{
	case 0x00: /* exon exon */
	  break ;
	case 0x11: /*intron intron */
	  if (x1 != x2)
	    {
	      if (
		  (i1 < max1 && i2 < max2) || 
		  (old1 != old2 && j1 && j2)
		  ) s2->wrong_exon++ ;
	    }
	  else  /* x1 == x2 */
	    {
	      if (old1 == old2)	s2->identical_exon++ ;
	      else if (j1 && j2) s2->wrong_exon++ ;
	    }
	  if (gling1 && !gling2)
	    s2->missed_gling++ ;
	  if (!gling1 && gling2)
	    s2->more_gling++ ;
	  break ;
	case 0x01: /* exon intron */
	  if (x1 > x2 && old1 < old2) 
	    {
	      if (!j2) s2->more_5p++ ;
	      else s2->more_central++ ;
	      s2->more_exon++ ; 
	    }
	  if (x1 < x2 && old1 > old2 && j1)  
	    {
	      s2->missed_intron++ ;
	    }
	  break ;
	case 0x10: /* intron exon */
	  if (x1 < x2 && old1 > old2)  
	    {
	      if (i1 == max1) s2->missed_3p++ ;
	      else if (!j2) s2->missed_5p++ ; 
	      else s2->missed_central++ ; 
	      s2->missed_exon++ ;
	    }
	  if (x1 > x2 && old1 < old2 && j2)  
	    {
	      s2->more_intron++ ;
	    }
	  break ;
	}
      /* increment */
      icr = 0 ; /* falsify first because we may increment both */
      if (x1 <= x2) icr |= 1 ; /* increment x1 */
      if (x1 >= x2) icr |= 2 ; /* increment x2 */
      
      if (icr & 0x1)
	{ 
	  old1 = x1 ; gling1 = FALSE ;
	  if (! in_intron1) in_intron1++ ;
	  else { in_intron1 = 0 ; i1++ ; j1++ ; }
	}
      if (icr & 0x2)
	{
	  old2 = x2 ; gling2 = FALSE ;
	  if (! in_intron2) in_intron2++ ;
	  else { in_intron2 = 0 ; i2++ ; j2++ ; }
	}
    } /* loop along both genes */
 done:
  ac_free (h) ;
  return ;
} /* gExonAnalysis */

/*************************************************************************************/
/* genebox touch on EST */
static BOOL gGeneBoxOverlap (STAT *sp, STAT *sq, BOOL strandSensitive)
{
  int dp, dq, a1, a2, b1, b2, u1, u2, du ;

  if (sp->chrom != sq->chrom)
    return FALSE ;
  dp = sp->a2 - sp->a1 ;
  dq = sq->a2 - sq->a1 ;
  
  if (strandSensitive && ((dp > 0 && dq < 0) || (dq > 0 && dp < 0)))
    return FALSE ;

  if (dp > 0) { a1 = sp->a1 ; a2 = sp->a2 ; }
  else { a1 = sp->a2 ; a2 = sp->a1 ; }
  if (dq > 0) { b1 = sq->a1 ; b2 = sq->a2 ; }
  else { b1 = sq->a2 ; b2 = sq->a1 ; }

  u1 = a1 > b1 ? a1 : b1 ;
  u2 = a2 < b2 ? a2 : b2 ;
  if (dp < 0) dp = -dp ;
  if (dq < 0) dq = -dq ;
  du = u2 - u1 ;
  if (du < 0)
    return FALSE ;
  return TRUE ;
} /* gGeneBoxOverlap */

/*************************************************************************************/

/* overlaps on EST, returns size of actual small overlap */
static BOOL gCloneOverlap (STAT *sp, STAT *sq, int *dxp)
{
  int a1, a2, b1, b2, u1, u2, du, v1, v2, dv = 0 ;

  /* if not same zone of cdna, they both survive */
  a1 = sp->start ; a2 = a1 + sp->ali + sp->c2 - sp->c1 ;
  b1 = sq->start ; b2 = b1 + sq->ali + sq->c2 - sq->c1;
  u1 = a1 > b1 ? a1 : b1 ;
  u2 = a2 < b2 ? a2 : b2 ;
  du = u2 - u1 + 1 ;
  if (dxp) *dxp = du ;
  if (du < 30 || (2*du < sp->ali && 2*du < sq->ali)) /* was 3 , changed mars 15-2004*/
    return FALSE ;
  /* check for a gap */
  if (4 * (sp->c2 - sp->c1) > du)
    {
      v1 = u1 > sp->c1  ? u1 : sp->c1 ;
      v2 = u2 < sp->c2  ? u2 : sp->c2 ;
      dv = v2 - v1 ;
    }
  else  if (4 * (sq->c2 - sq->c1) > du)
    {
      v1 = u1 > sq->c1  ? u1 : sq->c1 ;
      v2 = u2 < sq->c2  ? u2 : sq->c2 ;
      dv = v2 - v1 ;
    }
  if (4*dv > 3 *du) /* the intersect falls in the gap */
    { 
      du -= dv ;
      if (du < 0) du = 0 ;
      if (dxp) *dxp -= du ;
      return FALSE ;
    }
  return TRUE ;
} /* gCloneOverlap */

/*************************************************************************************/

/* actual overlap of exons */
static BOOL gExonOverlap (STAT *sp, STAT *sq, int *dxp)
{
  int du = 0 ;
  int a1, a2, b1, b2, db, u1, u2, ir, jr ;
  AC_TABLE gAf1, gAf2 ;
  
  gAf1 = ac_tag_table (sp->Gene, "Assembled_from", 0) ;
  gAf2 = ac_tag_table (sq->Gene, "Assembled_from", 0) ;
  /* if not same zone of cdna, they both survive */
  a1 = sp->a1 ; a2 = sp->a2 ;
  b1 = sq->a1 ; b2 = sq->a2 ;
  if (a1 < a2) db = b1 - a1 + 1 ;
  else db = a1 - b1 + 1 ;
  
  for (du = ir = 0 ; gAf1 && ir < gAf1->rows ; ir++)
    {
      a1 = ac_table_int (gAf1, ir, 0, 0) ;
      a2 = ac_table_int (gAf1, ir, 1, 0) ;
      for (jr = 0 ; gAf2 && jr < gAf2->rows ; jr++)
	{
	  b1 = db + ac_table_int (gAf2, jr, 0, 0) ;
	  b2 = db + ac_table_int (gAf2, jr, 1, 0) ;
	  u1 = a1 > b1 ? a1 : b1 ;
	  u2 = a2 < b2 ? a2 : b2 ;
	  if (u1 <= u2) du += u2 - u1 + 1 ;
	}
    }
    
  ac_free (gAf1) ;
  ac_free (gAf2) ;
  if (dxp) *dxp = du ;
  
  if (du < 30 || (3*du < sp->ali && 3*du < sq->ali))
    return FALSE ;
  return TRUE ;
} /* gExonOverlap */

/*************************************************************************************/

static int gGetOrf (char *dna, int *p1, int *p2, int *metp)
{
  int best = 0, max = dna ? strlen (dna) : 0 ;
  char *cp, *cq, *cr ;
  int frame, i, i0, len;
  
  if (max)
    { 
      for (frame = 0 ; frame < 3 ; frame++)
	{
	  for (i = i0 = frame, cp = dna + i, cq = cp+1, cr = cp+2 ; 
	       i < max - 2 + 3 ; cp += 3, cq += 3, cr += 3, i += 3)
	    {
	      if (i >= max - 2 ||
		  (
		   (*cp=='t' && *cq=='a'&& *cr=='a') || /* ocre */
		   (*cp=='t' && *cq=='a'&& *cr=='g') || /* ambre */
		   (*cp=='t' && *cq=='g'&& *cr=='a')  /* azur */
		   )
		  )
		{
		  len = (i - i0)/3 ; 
		  if (len > best)
		    { 
		      best = len ; 
		      if (p1) *p1 = i0 + 1 ; /* plato */ 
		      if (p2) { *p2 = i < max - 2 ? i + 3 : i ; }
		    }
		  i0 = i + 3 ;
		}
	    }
	}
    }
  
  if (metp) *metp = 0 ;
  if (p1 && p2 && metp && *p1 > 0)
    {
      int met, leu ;

      met = leu = 0 ;
      for (i = *p1 - 1, cp = dna + i, cq = cp+1, cr = cp+2 ; 
	   i < *p2 - 32 ; /* impose at least 10 AA */
	   cp += 3, cq += 3, cr += 3, i += 3)
	{
	  if (
	      *cp == 'a' && *cq == 't' && *cr == 'g'
	      )
	    {
	      met = i + 1 ; /* plato */
	      break ;
	    }
	}
      for (i = *p1 - 1, cp = dna + i, cq = cp+1, cr = cp+2 ; 
	   i < *p2 - 32 ; /* impose at least 10 AA */
	   cp += 3, cq += 3, cr += 3, i += 3)
	{
	  if (
	      *cq == 't' && *cr == 'g' /* && any cp = ATGC */
	      )
	    {
	      leu = i + 1 ; /* plato */
	      break ;
	    }
	}
      if (leu && leu < met - 30)
	met = leu ;
      *metp = met ;
    }

  return best ;
} /* gGetOrf */

/*************************************************************************************/

static char* gGetDNA (AC_OBJ Gene, AC_HANDLE h)
{
  int max, i, j, ir ;
  int a1, a2, x1, x2 ;
  char *dna = 0, *dna2 = 0, *cp, *cq ;
  const char *nam = ac_name (Gene) ;
  AC_OBJ Cosmid = ac_tag_obj (Gene, "Genomic_sequence", h) ;
  AC_TABLE spl, tgs ;

  if ((Cosmid = ac_tag_obj (Gene, "Genomic_sequence", h)) && 
      (spl = ac_tag_table (Gene, "Assembled_from", h)) &&
      spl->rows &&
      (tgs = ac_tag_table (Cosmid, "Transcribed_gene", h)) &&
      tgs->rows
      )
    {
      for (ir = 0 ; ir < tgs->rows ; ir++)
	{
	  if (!strcmp (nam, ac_table_printable (tgs, ir, 0, "")))
	    {
	      a1 = ac_table_int (tgs, ir, 1, 0) ;
	      a2 = ac_table_int (tgs, ir, 2, 0) ;
	      if (a1 && a2) 
		dna = ac_zone_dna (Cosmid, a1, a2, h) ;
	    }
	}
      /* dna is allocated on h, it is mine, i can splice it freely */
      max = dna ? strlen (dna) : 0 ;
      dna2 = halloc (2*max+10000, h) ;
      for (ir = j = 0, cp = dna2 ; dna && ir < spl->rows ; ir++)
	{
	  x1 = ac_table_int (spl, ir, 0, 0) ;
	  x2 = ac_table_int (spl, ir, 1, 0) ;
	  if (x1 < 1) x1 = 1 ; 
	  if (x2 > strlen (dna)) x2 = strlen (dna) ;
	  /*
	    y1 = ac_table_int (spl, ir, 3, 0) ;
	    y2 = ac_table_int (spl, ir, 4, 0) ;
	    if (y1 < clipTop && y2 >= clipTop)
	    x1 = x2 - (y2 - clipTop)
	    JUNK define correctly form where to were in genome coords we are in the insert
	  */

	  for (i = x1 - 1, cq = dna + i ; i < x2 && i < max ; i++, j++, cp++, cq++) 
	    *cp = *cq ;
	  if (i > max || j > 2*max+10000)
	    messcrash ("LenghtError in gGetDna gene %s max=%d i=%d j=%d", ac_name(Gene), max, i, j) ;
	}
      *cp = 0 ;
    }

  return dna2 ;
} /* gGetDNA */

/*************************************************************************************/

static int gGetOrfLen (AC_OBJ Gene)
{
  AC_HANDLE h = ac_new_handle () ;
  char *dna = 0 ;
  int best = 0 ;

  if (/*! (best = ac_tag_int (Gene, "ORF", 0)) && */
      (dna = gGetDNA (Gene, h)))
    best = gGetOrf (dna, 0, 0, 0) ;
  
  ac_free (h) ;
  return best ;
} /* gGetOrfLen */

/*************************************************************************************/

typedef struct {char gene[256], chrom[256] ; int a1, a2 ;} HYP ;

static int gFlagEst (int *methods, char **methodTitle, Array aa, DICT *dict, AC_OBJ est, int errCost)
{
  HYP hyp[] = {
    {"IGKC", "2", 89532338, 89058814 },
    {"IGKC_zone", "2", 89008603, 89532580 },

    {"HLA-DRB1", "c_28212469", 137073, 64642 },
    {"HLA-G", "6", 29902755, 29906891 },
    {"HLA-21", "6", 29999906, 30019495 },
    {"HLA-B&A", "6", 31429180, 31341270 },
    {"HLA-DRB1.1", "6", 32612273, 32534888 },
    {"HLA-DQA1", "6", 32655720, 32761934 },
    {"HLA-DQB1", "6", 32778325, 32673028 },
    {"HLA-DPA1", "6", 33095447, 33079235 },
    {"HLA-DPB1", "6", 33090601, 33144302 },



    {"IGHD", "14", 104985153, 104274775 },
    {"IGHD.1", "14", 105103157, 104274775 },
    {"IGHD_zone", "14", 103900000, 105000000 },

    {"IGLC1", "22", 20874931, 21586552 },
    {"IGLC1.1", "22", 20777701, 21589720 },
    {"IGLC1_zone", "22", 20700000, 21700000 },
    { "", "", 0, 0 }
  } ;
  STAT *sp, *sq, *sr ;
  const char *gene, *chrom ;
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Gene, Clone ;
  AC_TABLE gDup = 0, gAf = 0, imap = 0, iam = 0, fg = ac_tag_table (est, "From_gene", h) ;
  int ii, iii, i2,  ihyp, jj, j0, nn, jBest = -1, bestScore = -888888 ;
  int clipTop = 1, pA, pfA, nerr, penalty, nDeletion ; /* so no 999 becomes best */
  int 
    intronCost = 0, micro_exon_cost = 0, /* should be equal */
    badIntronCost = 3,  /* must be bigger than 1 error  */
    delCost = 1 ;       /* should be errCost - (aliCost == 1) */
  int bestUniqueScoreA = -999999, bestUniqueScoreB =  -999999, bestUniqueScoreC =  -999999, bestUniqueScoreD =  -999999, bestMScore = -99999 ;
  BOOL isDown = ac_has_tag (est, "Forward") ;

  for (ihyp = 0 ; *hyp[ihyp].gene ; ihyp++)
    if (hyp[ihyp].a1 >  hyp[ihyp].a2)
      { int aa = hyp[ihyp].a1 ; hyp[ihyp].a1 = hyp[ihyp].a2 ; hyp[ihyp].a2 = aa ; }
    
  Clone = ac_tag_obj (est, "cDNA_clone", h) ;
  /* loop on all alignments of the EST */
  for (ii = 0, j0 = jj = arrayMax (aa)  ; fg && ii < fg->rows ; ii++)
    {
      BOOL myDown = TRUE ;

      Gene = ac_table_obj (fg, ii, 0, h) ;
      gene = ac_name (Gene) ;
      if (!gene)
	continue ;
      for (iii = 0 ; methods[iii]; iii++)
	if (*gene == methods[iii])
	  break ;
      if (! methods [iii])
	continue ;

      sp = arrayp (aa, jj++, STAT) ;

      /* data imported directly */ 
      /* from the est object */
      dictAdd (dict, ac_name(est), &nn) ;
      sp->Est = est ;
      sp->est = nn ;
      if ((gDup = ac_tag_table (est, "Duplicate_of_read", h)))
	{
	  vTXT buf = vtxtCreate () ;
	  int iDup ;
	  for (iDup = 0 ; iDup < gDup->rows ; iDup++)
	    vtxtPrintf (buf, "%s%s", iDup ? ";" : "", ac_table_printable (gDup, iDup, 0, "")) ;
	  dictAdd (dict, vtxtPtr (buf), &nn) ;
	  sp->Duplicate_of_read = nn ;
	  vtxtDestroy (buf) ;
	}
      ac_free (gDup) ;
      sp->iMethod = iii ;

      /* collect the position of the vector and polyA */
      {
	AC_TABLE mvClip, mDna = ac_tag_table (est, "DNA", h) ;
	sp->ln = ac_table_int (mDna, 0, 1, 0) ;
	mvClip = ac_tag_table (est, "Vector_clipping", h) ;
	if (mvClip)
	  {
	    sp->clipTop = clipTop = ac_table_int (mvClip, 0, 0, 0) ;
	    sp->clipEnd = ac_table_int (mvClip, 0, 1, 0) ;
	    ac_free (mvClip) ;
	  }
	sp->pA = pA = ac_tag_int (est, "PolyA_after_base", 0) ;
	sp->pfA = pfA = 0 ;
	if ((iam = ac_tag_table (est, "Initial_PolyA", h)))
	  {
	    sp->pfA = pfA = ac_table_int (iam, 0, 0, 0) + ac_table_int (iam, 0, 1, 0) ;
	  }
      }
      /* est->from_gene */
      sp->start = ac_table_int (fg, ii, 2, 0) ;
      sp->ali = ac_table_int (fg, ii, 3, 0) ;
      sp->err = ac_table_int (fg, ii, 4, 999999) ;
      /* from the tg */
      sp->Gene = Gene ;
      sp->method = *gene ;
      dictAdd (dict, gene, &nn) ;
      sp->gene = nn ;

      sp->orf = gGetOrfLen (Gene) ;
      chrom = 0 ;
      imap = ac_tag_table (Gene, "IntMap", h) ;
      if (imap && imap->cols >= 3)
	{
	  chrom = ac_table_printable (imap, 0, 0, 0) ;
	  if (chrom)
	    {
	      dictAdd (dict, chrom, &nn) ;
	      sp->chrom = nn ;
	      sp->a1 = ac_table_int (imap, 0, 1, 0) ;
	      sp->a2 = ac_table_int (imap, 0, 2, 0) ;
	      sp->covers = sp->a2 > sp->a1 ? sp->a2 - sp->a1 + 1 : sp->a1 - sp->a2 + 1 ;
	    }
	}

     /* special hack to please weissenbach and remove hyper variable zones*/
      for (ihyp = 0 ; chrom && ! sp->masked && *hyp[ihyp].gene ; ihyp++)
	{
	  if (!strcmp (dictName(dict, sp->chrom), hyp[ihyp].chrom))
	    {
	      if (
		  ( sp->a1 > hyp[ihyp].a1 && sp->a1 < hyp[ihyp].a2 ) ||
		  ( sp->a2 > hyp[ihyp].a1 && sp->a2 < hyp[ihyp].a2 )
		  )
		sp->masked = 1 ;
	    }
	} 

      sp->nexons = ac_tag_int (Gene, "Nb_possible_exons", 0) ;
      sp->bad_intron = sp->nexons - 1 
	- ac_tag_int (Gene, "gt_ag", 0) - ac_tag_int (Gene, "gc_ag", 0) - ac_tag_int (Gene, "at_ac", 0) ;

      if (sp->nexons == 1 && ac_has_tag (Clone, "Possible_genomic_contamination"))
	sp->genomic_contamination = 1 ;
      /* data evaluated dynamically on one alignment */
      /* get the score */
      {
	int jr, x1, x2, a2Old, x2Old, a1, a2, da, dx, dx2 ;
	myDown = TRUE ;

	gAf = ac_tag_table (Gene, "Assembled_from", h) ;
	nerr = penalty = nDeletion = 0 ;
	x2Old = a2Old = 0 ;
	sp->ali = 0 ; sp->start = sp->ln ;
	for (jr = 0 ; gAf && jr < gAf->rows ; jr++)
	  {
	    pA = sp->pA ; pfA = sp->pfA ;
	    a1 = ac_table_int (gAf, jr, 0, 0) ;
	    a2 = ac_table_int (gAf, jr, 1, 0) ;
	    x1 = ac_table_int (gAf, jr, 3, 0) ;
	    x2 = ac_table_int (gAf, jr, 4, 0) ;
	    da = a2 - a1 + 1 ; 
	    dx = x2 - x1 ; if (dx < 0) dx = -dx ; dx++ ;
	    if (x1 < sp->start) sp->start = x1 ;
	    if (x2 < sp->start) sp->start = x2 ;
	    if (!jr)
	      {
		myDown = x1 < x2 ? TRUE : FALSE ;
		if (myDown != isDown) sp->pA = pA = sp->pfA = pfA = 0 ;
		sp->x1 = x1 ; sp->x2 = x2 ;
		if (
		    (myDown && pA < 200)  ||
		    (!myDown && pA > 200)
		    )
		  sp->pA = pA = 0 ;
		if (!myDown)
		  sp->pfA = pfA = 0 ;
		if (x1 < x2)
		  sp->g_forward = 1 ;
		else
		  sp->g_reverse = 1;
	      }
	    pfA = sp->pfA ;
	    clipTop = sp->clipTop ;
	    if (clipTop <= pfA) clipTop = 0 ;
	    if (pfA < clipTop) pfA = 0 ;
	    sp->x2 = x2 ;

	    /* penalty for short first exon */
	    if (jr == 1 && sp->ali - 2 * nerr < 15 && (sp->ali - 2 * nerr < 1 || (1 << 2*(sp->ali - 2 * nerr)) < 4 * a1))
	      { sp->too_short_first_exon++ ; penalty += 3 + (sp->ali > 0 ? sp->ali : 0) ; }
	    /* we have an intron, errors in short first exon are doubled */
	    if (jr == 1 &&  sp->ali < 50) /* we are starting exon 2 */
	      penalty += nerr * errCost ; /* we wish to double the errors of the first exon
				   but nerr must stay exactly smith-waterman in the from_gene line */
	    /* cdna topology */
	    /* cdna used twice, subtract 3 times the number of reused bp */
	    dx2 = jr ? (myDown ? x2Old - x1 + 1 : x1 - x2Old + 1 ) : 0 ;
	    if (dx2 > 0)
	      penalty += 3 * dx2 ;
	    /* reuse of the genome for the next piece of cdna */
	    dx2 = jr && a1 <= a2Old ? (myDown ? (a2Old - a1) - (x2Old - x1) : (a2Old- a1) - (x1 - x2Old)) : 0 ;
	    if (dx2 > 0)
	      sp->v_repeat = dx2 ;

	    /* flag micro introns, no penalty */
	    dx2 = jr ? a1 - a2Old : 0 ;
	    if (dx2 > 0 && dx2 < 65)
	      sp->micro_deletion++ ;
	    
	    /* genome deletion == unaligned loop in the cdna */
	    /* this goes in the error cost not the penalty so they do not count as excellent */
	    dx2 = jr ? (myDown ?  x1 - x2Old - 1 :  x2Old - x1 - 1) : 0 ;
	    if (dx2 > 0)
	      {
		nDeletion += dx2 > 6 ? 6 : dx2 ;
		if (dx2 >= 18) sp->genome_deletion += dx2 ;
	      }
	    if (dx2 > 18 && dx2 > sp->c2 - sp->c1)
	      { 
		if (myDown) { sp->c1 = x2Old + 1 ; sp->c2 = x1 - 1 ; } 
		else { sp->c1 = x1 - 1 ; sp->c2 = x2Old + 1 ;}
	      }
	    if (dx - da >= 18) /* albeit we have a single exon, it is longer on the cdna than on the genome */
	      sp->genome_deletion += dx - da ;
	    if (sp->a2 < sp->a1 + 12) /* micro_exon or gling-gling */
	      penalty += micro_exon_cost ;

	    x2Old = x2 ; a2Old = a2 ;
	    if (  /* big penalty for pure first polyA exon, and do  not count the polyA aligned bp */
		(pfA && 3 * (x2 - x1) < 4 * (pfA - x1))
		)
	      { penalty += 3 + x2 - x1 ; sp->polyA_first_exon = 1 ; }
	    else if (pfA &&  myDown && x1 < pfA) /* little penalty for overaligning the AAA in same  exon */
	      penalty += pfA - x1 ;
	    if (  /* big penalty for pure polyA exon */
		(pA &&  myDown && 3 * (x2 - x1) < 4 * (x2 - pA)) ||
		(pA &&  myDown && pA < x2 && pA - x1 < 4) ||
		(pA && !myDown && pA > x2 && x1 - pA < 4) ||
		(pA && !myDown && 3 * (x1 - x2) < 4 * (pA - x2))
		)
	      { penalty += 3 ; sp->polyA_exon = 1 ; } /* and we do not count the aligned bp */
	    else if (pA &&  myDown && x2 > pA) /* penalty: do not count bp belonging inside the polyA */
	      sp->ali += pA - x1 + 1 ;
	    else if (pA &&  !myDown && x1 < pA)
	      sp->ali += x2 - pA + 1 ;
	    else
	      sp->ali += x2 > x1 ? x2 - x1 + 1 : x1 - x2 + 1 ;
	    /* vector clipping penalty */
	    if ( clipTop > 1 &&  /* big penalty for pure vector exon */
		(clipTop &&  myDown && 3 * (x2 - x1) < 4 * (clipTop - x1))
		)
	      { penalty += 3 + x2 - x1 + 1 ; sp->vector_exon = 1 ; }
	    else if (clipTop > 1 &&  myDown && x1 < clipTop) /* penalty: do not count bp belonging inside the vector */
	      penalty += clipTop - x1 ;
	    nerr += ac_table_int (gAf, jr, 5, 999999) ;
	  }
	}

      sp->err = nerr ;
      if (nerr >= 999999 ||
	  100*nerr > 20*sp->ali || 
	  sp->start + sp->ali > sp->ln + 200) /* 200 because of possible polya clipping */
	sp->score = -999999 ;
      else
	sp->score = 
	  sp->ali                             /* nb of aligned bp */
	  - delCost * nDeletion               /* holes in the cdna */         
	  - errCost * nerr                    /* point errors or 1 bp indel, from Smith-Waterman */
	  - penalty                           /* micro terminal exons, vector exons, polyA exons */
	  - intronCost * (sp->nexons - 1)     /* each intron has a cost */
	  - badIntronCost *  sp->bad_intron ; /* non classic introns have an extra cost */
    
      if (sp->score < -20000)
	sp->no_score = 1 ;
      
      if (jBest < 0 || sp->score > bestScore)
	{ jBest = jj ; bestScore = sp->score ; }
 

      if (sp->clipEnd && sp->clipEnd < sp->ln)
	sp->ln = sp->clipEnd ;
      if (sp->pA && myDown && sp->pA < sp->ln)
	sp->ln = sp->pA ;
      if (pfA > 1 && myDown)
	sp->ln -= sp->pfA - 1 ;
      if (clipTop > 1)
	sp->ln -= clipTop - 1 ;
    }


  /* add the orf bonus */
  {
    int bestOrf = 0 ;
    int bonusNeeded = 0 ;

    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (bestOrf < sp->orf)
	  bestOrf = sp->orf ;
      }

    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (bestOrf > sp->orf + 30 &&  /* sp should win */
	    bestScore <= sp->score + 3) /* sp could win */
	  bonusNeeded = 1 ;
      }

    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ; 
	if (bonusNeeded)
	  {
	    if (bestOrf == sp->orf)
	      {
		sp->orf_bonus = 3 ;
		sp->score += 3 ;
		if (bestScore < sp->score)
		  bestScore = sp->score ;
	      }
	    else
	      sp->orf_bonus = sp->orf - bestOrf ;
	  }
	else
	  sp->orf_bonus = 1 ; /* do not export */
      }
  }
 
  /* choose between the various equivalent methods */
  for (iii = 0 ; methods [iii] ; iii++)
    {
      int bestUniqueScore = -77777 ;

      for (ii = j0 ; ii < jj ; ii++)
	{
	  sp = arrayp (aa, ii, STAT) ;
	  if (sp->method != methods[iii] || sp->no_score)
	    continue ;
	  if (sp->score > bestUniqueScore)
	    bestUniqueScore = sp->score ;
	}
      if (methods[iii] == 'A') bestUniqueScoreA = bestUniqueScore ;
      if (methods[iii] == 'B') bestUniqueScoreB = bestUniqueScore ;
      if (methods[iii] == 'C') bestUniqueScoreC = bestUniqueScore ;
      if (methods[iii] == 'D') bestUniqueScoreD = bestUniqueScore ;
    }

  if (bestUniqueScoreA && bestUniqueScoreA <= bestUniqueScoreB)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'A')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreB && bestUniqueScoreB < bestUniqueScoreA)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'B')
	  sp->kill_method = 1 ;
      }

  if (bestUniqueScoreB && bestUniqueScoreB <= bestUniqueScoreC)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'B')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreC && bestUniqueScoreC < bestUniqueScoreB)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'C')
	  sp->kill_method = 1 ;
      }

  if (bestUniqueScoreA && bestUniqueScoreA <= bestUniqueScoreC)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'A')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreC && bestUniqueScoreC < bestUniqueScoreA)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'C')
	  sp->kill_method = 1 ;
      }

  if (bestUniqueScoreA && bestUniqueScoreA <= bestUniqueScoreD)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'A')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreD && bestUniqueScoreD < bestUniqueScoreA)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'D')
	  sp->kill_method = 1 ;
      }

  if (bestUniqueScoreB && bestUniqueScoreB < bestUniqueScoreD)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'B')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreD && bestUniqueScoreD <= bestUniqueScoreB)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'D')
	  sp->kill_method = 1 ;
      }

  if (bestUniqueScoreC && bestUniqueScoreC < bestUniqueScoreD)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'C')
	  sp->kill_method = 1 ;
      }
  if (bestUniqueScoreD && bestUniqueScoreD <= bestUniqueScoreC)
    for (ii = j0 ; ii < jj ; ii++)
      {
	sp = arrayp (aa, ii, STAT) ;
	if (sp->method == 'D')
	  sp->kill_method = 1 ;
      }



  /* find the split alignments  + m_best */
  for (iii = 0 ; methods [iii] ; iii++)
    {
      int bestUniqueScore = -77777, bestSplitScore = -77777, iBest = -1, jBest = -1 ;
      double aliTotal ; /* cast to double to correctly get ratio */
      int overlap ;

      for (ii = j0 ; ii < jj ; ii++)
	{
	  sp = arrayp (aa, ii, STAT) ;
	  if (sp->method != methods[iii] || sp->no_score || sp->kill_method)
	    continue ;
	  if (sp->score > bestUniqueScore)
	    bestUniqueScore = sp->score ;

	  for (i2 = ii + 1 ; i2 < jj ; i2++)
	    {
	      sq = arrayp (aa, i2, STAT) ;
	      if (sp->method != sq->method || sq->no_score)
		continue ;
	      if (gCloneOverlap (sp, sq, &overlap)) /*  || gGeneBoxOverlap (sp, sq, TRUE)) */
		continue ;
	      aliTotal = sp->ali + sq->ali ;
	      if (overlap >= aliTotal)
		continue ;
	      if (overlap < 0) overlap = 0 ; /* do not favor holes */
	      if ((aliTotal - overlap) * (sp->score + sq->score) > aliTotal * bestSplitScore)
		{
		  bestSplitScore = (sp->score + sq->score) * (aliTotal - overlap)/ aliTotal ;
		  iBest = ii ; jBest = i2 ;
		}
	    }
	}
      bestSplitScore -= 6 ;
      if (bestSplitScore > bestUniqueScore)
	{
	  sp = arrayp (aa, iBest, STAT) ;
	  sq = arrayp (aa, jBest, STAT) ;
	  sp->m_split = jBest + 1 ;
	  sq->m_split = iBest + 1;
	  
	  sp->m_score = sq->m_score = bestSplitScore ;
	  i2 = sp->ali + sq->ali ;
	  sp->m_ali = sq->m_ali = i2 ;
	  i2 = sp->err + sq->err ;
	  sp->m_err = sq->m_err = i2 ;
	  sp->m_best = sq->m_best = 1 ;
	  if (sp->ali > sq->ali)
	    { sp->m_major = 1 ; sq->m_minor = 1 ; }
	  else
	    { sq->m_major = 1 ; sp->m_minor = 1 ; }

	  /* try now to add a third bits AF155662, forget about rationalizing */
	  for (ii = j0 ; ii < jj ; ii++)
	    {
	      if (ii == iBest || ii == jBest)
		continue ;
	      sr = arrayp (aa, ii, STAT) ;
	      if (sr->no_score || sr->score < 30)
		continue ; 
	      if (sr->method != methods[iii])
		continue ;
	      if (gCloneOverlap (sp, sr, 0) || /* gGeneBoxOverlap (sp, sr, TRUE) || */
		  gCloneOverlap (sq, sr, 0) ) /* gGeneBoxOverlap (sq, sr, TRUE)) */
		continue ;
	      i2 = sp->m_ali + sr->ali ;
	      sp->m_ali = sq->m_ali = sr->m_ali = i2 ;
	      i2 = sp->m_err + sr->err ;
	      sp->m_err = sq->m_err = sr->m_err = i2 ;
	      i2 = sp->m_score + sr->score - 6 ;
	      bestSplitScore = sp->m_score = sq->m_score = sr->m_score = i2 ;
	      sr->m_minor = 1 ; sr->m_split = iBest + 1;
	      break ; /* 3 is enough, else we need to rewrite the code */
	    }

	} 
      else
	{
	  for (ii = j0 ; ii < jj ; ii++)
	    {
	      sp = arrayp (aa, ii, STAT) ;
	      if (sp->method != methods[iii])
		continue ;
	      if (sp->score == bestUniqueScore)
		sp->m_best = 1 ;
	    }
	}
      if (bestSplitScore > bestMScore)  bestMScore = bestSplitScore ;
    }

  /* find the best any method split alignments */
  if (1)
    {
      int iBest = -1, jBest = -1, bestSplitScore = -77777 ;

      for (ii = j0 ; ii < jj ; ii++)
	{
	  double aliTotal ; /* cast to double to correctly get ratio */
	  int overlap ;
	  
	  sp = arrayp (aa, ii, STAT) ;
	  if (sp->no_score || sp->kill_method)
	    continue ;

	  for (i2 = ii + 1 ; i2 < jj ; i2++)
	    {
	      sq = arrayp (aa, i2, STAT) ;
	      if (sq->no_score || sq->kill_method)
		continue ;
	      if (gCloneOverlap (sp, sq, &overlap) ) /* || gGeneBoxOverlap (sp, sq, TRUE)) */
		continue ; 
	      aliTotal = sp->ali + sq->ali ;
	      if (overlap >= aliTotal)
		continue ;
	      
	      if ((aliTotal - overlap) * (sp->score + sq->score) > aliTotal * bestSplitScore)
		{
		  bestSplitScore = (sp->score + sq->score) * (aliTotal - overlap)/ aliTotal ;
		  iBest = ii ; jBest = i2 ;
		}
 
	    }
	}
      bestSplitScore -= 6 ; /* penalty for using 2 parts */
      if (iBest != -1 && bestSplitScore >= bestScore)
	{
	  bestScore = bestSplitScore ;
	  sp = arrayp (aa, iBest, STAT) ;
	  sq = arrayp (aa, jBest, STAT) ;
	  sp->g_split = jBest + 1 ;
	  sq->g_split = iBest + 1;
	  sp->g_score = sq->g_score = bestScore ;
	  i2 = sp->ali + sq->ali ;
	  sp->g_ali = sq->g_ali = i2 ;
	  i2 = sp->err + sq->err ;
	  sp->g_err = sq->g_err = i2 ;
	  if (sp->ali > sq->ali)
	    { sp->g_major = 1 ; sq->g_minor = 1 ; }
	  else
	    { sq->g_major = 1 ; sp->g_minor = 1 ; }

	  /* try now to add a third bits AF155662 */
	  for (ii = j0 ; ii < jj ; ii++)
	    {
	      if (ii == iBest || ii == jBest)
		continue ;
	      sr = arrayp (aa, ii, STAT) ;
	      if (sr->no_score || sr->kill_method)
		continue ;
	      if (gCloneOverlap (sp, sr, 0) || /* gGeneBoxOverlap (sp, sr, TRUE) || */
		  gCloneOverlap (sq, sr, 0) ) /* gGeneBoxOverlap (sq, sr, TRUE)) */
		continue ;
	      i2 = sp->g_ali + sr->ali ;
	      sp->g_ali = sq->g_ali = sr->g_ali = i2 ;
	      i2 = sp->g_err + sr->err ;
	      sp->g_err = sq->g_err = sr->g_err = i2 ;
	      bestScore = sp->g_score + sr->score ;
	      sp->g_score = sq->g_score = sr->g_score = bestScore ;
	      sr->g_minor = 1 ; sr->g_split = iBest + 1;
	      break ; /* 3 is enough, else we need to rewrite the code */
	    }
	} 
    }

  if (bestMScore > bestScore) /* may happen is K finds a plit and S overwalks the exon as in AK023340 */
    {
      bestScore = bestMScore ;
      for (ii = j0 ; ii < jj ; ii++)
	{
	  sp = arrayp (aa, ii, STAT) ;
	  if (sp->m_score == bestMScore)
	    {
	      sp->g_split = sp->m_split ;
	      sp->g_major = sp->m_major ;
	      sp->g_minor = sp->m_minor ;
	      sp->g_score = sp->m_score ;
	      sp->g_ali = sp->m_ali ;
	      sp->g_err = sp->m_err ;
	    }
	  else
	     {
	      sp->g_split = 0 ;
	      sp->g_major = 0 ;
	      sp->g_minor = 0 ;
	      sp->g_score = 0 ;
	      sp->g_ali = 0 ;
	      sp->g_err = 0 ;
	    }
	}
    }

  /* comparison of all alignments */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->no_score) /* it is ok to give repetition_5 flag to killed_method */
	continue ;
      sp->best_score = bestScore ;
      if (!sp->m_split)
	{
	  sp->m_score = sp->score ;
	  if (sp->score == bestScore)
	    sp->best = 1 ;
	  else if (100 * sp->err > 20 * sp->ali)
	    sp->bad_score = 1 ;
	  else if (sp->score >= bestScore - 5)
	    sp->best5 = 1 ;
	  else if (sp->score >= bestScore - 50)
	    sp->best50 = 1 ;
	  else if (sp->score >= bestScore - 500)
	    sp->best500 = 1 ;
	  else 
	    sp->best5000 = 1 ;
	}
      else if (sp->m_split && sp->m_major && !sp->kill_method) /* no split killed */
	{
	  if (sp->m_score == bestScore)
	    sp->best = 1 ;
	  else if (100 * sp->m_err > 20 * sp->m_ali)
	    sp->bad_score = 1 ;
	  else if (sp->m_score >= bestScore - 5)
	    sp->best5 = 1 ;
	  else if (sp->m_score >= bestScore - 50)
	    sp->best50 = 1 ;
	  else if (sp->m_score >= bestScore - 500)
	    sp->best500 = 1 ;
	  else 
	    sp->best5000 = 1 ;
	}  
    }


  /***** alibaba, diamond ***********************/
  /* we now start by saying all alignments are aliba */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (!sp->no_score && sp->score && !sp->kill_method)
	sp->alibaba = 2 ;
      else
	sp->alibaba = 0 ;
    }
  /* now we keep locally only a single local_best alibaba */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (!sp->alibaba)
	continue ;      
      for (i2 = ii + 1 ;  sp->alibaba && i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ;
	  if (sq->alibaba && gGeneBoxOverlap (sp, sq, FALSE)) /* exons touch */
	    {
	      if (sp->g_score)
		sq->alibaba = 0 ;
	      else if (sq->g_score)
		sp->alibaba = 0 ;
	      else if (sp->score < sq->score)
		sp->alibaba = 0 ; 
	      else if (sp->score > sq->score)
		sq->alibaba = 0 ; 
	      else if (sp->bad_intron > sq->bad_intron)
		sp->alibaba = 0 ;
	      else if (sp->bad_intron < sq->bad_intron)
		sq->alibaba = 0 ;
	      else if (sp->ali >  sq->ali)
		sq->alibaba = 0 ;
	      else if (sp->ali < sq->ali)
		sp->alibaba = 0 ;
	      else if (sp->covers < sq->covers)
		{ sq->alibaba = 0 ; sq->too_long = 1 ; }
	      else if (sp->covers > sq->covers)
		{ sp->alibaba = 0 ; sp->too_long = 1 ; }
	      else if (sp->method == 'A' && sq->method == 'B')
		sp->alibaba = 0 ;
	      else if (sp->method == 'B' && sq->method == 'A')
		sq->alibaba = 0 ;
	      else if (sp->method < sq->method)
		sq->alibaba = 0 ;
	      else 
		sp->alibaba = 0 ;
	    }
	}
    }

  /* remember the local best */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->alibaba)
	sp->local_best = 1 ;
    }

  /* now for every pair of cdna-ovelapping alibaba, we kill one of them */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (!sp->alibaba)
	continue ;

      for (i2 = ii + 1 ; sp->alibaba && i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ;
	  if (!sq->alibaba)  /* => no genebox overlap */
	    continue ;
	  if (gCloneOverlap (sp, sq, 0)) /* geneboxes touch */
	    {
	      if (sp->g_score)
		sq->alibaba = 0 ;
	      else if (sq->g_score)
		sp->alibaba = 0 ;
	      else if (sp->score < sq->score)
		sp->alibaba = 0 ; 
	      else if (sp->score > sq->score)
		sq->alibaba = 0 ; 
	      else if (sp->bad_intron > sq->bad_intron)
		sp->alibaba = 0 ;
	      else if (sp->bad_intron < sq->bad_intron)
		sq->alibaba = 0 ;
	      else if (sp->ali >  sq->ali)
		sq->alibaba = 0 ;
	      else if (sp->ali < sq->ali)
		sp->alibaba = 0 ;
	      else if (sp->covers < sq->covers)
		{ sq->alibaba = 0 ; sq->too_long = 1 ; }
	      else if (sp->covers > sq->covers)
		{ sp->alibaba = 0 ; sp->too_long = 1 ; }
	      else if (sp->method == 'A' && sq->method == 'B')
		sp->alibaba = 0 ;
	      else if (sp->method == 'B' && sq->method == 'A')
		sq->alibaba = 0 ;
	      else if (sp->method < sq->method)
		sq->alibaba = 0 ;
	      else 
		sp->alibaba = 0 ;
	    }	      
	}
    } 

  /* because g_major is a composite, it may not be the alibaba, but should be */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->m_minor)
	sp->alibaba = 0 ;
      if (sp->g_major)
	{
	  for (i2 = j0 ; i2 < jj ; i2++)
	    {
	      sq = arrayp (aa, i2, STAT) ;
	      sq->alibaba = 0 ;
	    }
	  sp->alibaba = 2 ;
	}
    }

  /* in rare cases we end up with 2 alibaba, one beeing very short inside the other big one */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->alibaba)
	{
	  for (i2 = j0 ; sp->alibaba && i2 < jj ; i2++)
	    {
	      sq = arrayp (aa, i2, STAT) ;
	      if (sq->alibaba)
		{
		  if (sp->ali >  sq->ali)
		    sq->alibaba = 0 ;
		  else if (sp->ali < sq->ali)
		    sp->alibaba = 0 ;
		}
	    }
	}
    }

  /* locate dubious strand */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->alibaba != 2 || sp->g_split)
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{
	  if (ii == i2)
	    continue ; 
	  sq = arrayp (aa, i2, STAT) ;
	  if (!gCloneOverlap (sp, sq, 0))
	    continue ;
	  if (sp->chrom == sq->chrom)
	    {
	      int a1, a2, b1, b2, u1, u2, da, db, du ;
 
	      a1 = sp->a1 ; a2 = sp->a2 ;
	      b1 = sq->a1 ; b2 = sq->a2 ;
	      if (a1 < a2 && b1 > b2)
		{
		  da = a2 - a1 ; db = b1 - b2 ;
		  u1 = a1 > b2 ? a1 : b2 ;
		  u2 = a2 < b1 ? a2 : b1 ;
		  du = u2 - u1 ;
		  if (5*du > 4*da && 5*du > 4 * db) 
		    sq->dubious_strand = TRUE ;
		}
	      else if (a1 > a2 && b1 < b2)
		{
		  da = a1 - a2 ; db = b2 - b1 ;
		  u1 = b1 > a2 ? b1 : a2 ;
		  u2 = b2 < a1 ? b2 : a1 ;
		  du = u2 - u1 ;
		  if (5*du > 4*da && 5*du > 4 * db) 
		    sq->dubious_strand = TRUE ;
		}
	    }
	}
    }

  /* repetitions, they also count as gold */
  for (ii = j0 ; ii < jj ; ii++)
    {
      int found = 0 ;
      
      sp = arrayp (aa, ii, STAT) ;
      if (sp->alibaba != 2 || sp->g_split)  /* here killed_method is ok */
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{
	  if (ii == i2)
	    continue ; 
	  sq = arrayp (aa, i2, STAT) ;
	  if (!sq->local_best || sq->m_split)
	    continue ;
	  if (!gCloneOverlap (sp, sq, 0))
	    continue ;
	  if (sp->chrom == sq->chrom)
	    {
	      int a1, a2, b1, b2, u1, u2, da, db, du ;
 
	      a1 = sp->a1 ; a2 = sp->a2 ;
	      b1 = sq->a1 ; b2 = sq->a2 ;
	      if (a1 < a2 && b1 > b2)
		{
		  da = a2 - a1 ; db = b1 - b2 ;
		  u1 = a1 > b2 ? a1 : b2 ;
		  u2 = a2 < b1 ? a2 : b1 ;
		  du = u2 - u1 ;
		  if (5*du > 4*da && 5*du > 4 * db) 
		    {
		      if (sp->score < sq->score) sp->wrong_strand = TRUE ;
		      else sq->wrong_strand = TRUE ;
		    }
		}
	      else if (a1 > a2 && b1 < b2)
		{
		  da = a1 - a2 ; db = b2 - b1 ;
		  u1 = b1 > a2 ? b1 : a2 ;
		  u2 = b2 < a1 ? b2 : a1 ;
		  du = u2 - u1 ;
		  if (5*du > 4*da && 5*du > 4 * db) 
		    {
		      if (sp->score < sq->score) sp->wrong_strand = TRUE ;
		      else sq->wrong_strand = TRUE ;
		    }
		}
	    }
	  if (sq->best && ! sq->kill_method) /* set alibaba = 1 so we do not count quadratically */
	    { found = 1 ; sp->repetition++ ; sq->repetition++ ; sq->alibaba = 1 ; }
	  else if (sq->best5) /* here killed_method is ok */
	    { sp->repetition5++ ; sq->repetition5++ ; /* sq->alibaba = 1 ; */ }
	  else if (sq->best50)
	    { sp->repetition50++ ; sq->repetition50++ ; }
	}
      sp->repetition += found ;
    }

  /* intrinsic quality of the gold alignment */
  /* attributed to the loners (including repeats) or to the g_major of the split */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      pA = sp->pA ;
      if (sp->no_score || !sp->alibaba)
	continue ;
      if (!sp->g_split)
	{
	  int ln = sp->ln ;
	  int flag = sp->alibaba ? 0x3 : 0x1 ;
	  if (100.0 * sp->ali >= 99 * ln)
	    {
	      if (1000.0 * sp->err <= 4 * sp->ali)
		sp->excellent = flag ; /* ali >= 99% && err <= .4% */
	      else  if (1000.0 * sp->err <= 15 * sp->ali)
		sp->good = flag ;/* ali >= 99% && .4% < err <= 1.5 */
	      else if (1000.0 * sp->err <= 30 * sp->ali)
		sp->dubious = sp->dubious1 = flag ;
	      else
		sp->bad = sp->bad1 = flag ;
	    }
	  else if (100.0 * sp->ali >= 60 * ln)
	    {
	      if (1000.0 * sp->err <= 4 * sp->ali)
		sp->partial = sp->partial1 = flag ;   /* 99% > ali >= 60% && err <= .4% */
	      else if (1000.0 * sp->err <= 15 * sp->ali)
		sp->partial = sp->partial2 = flag ;   /* 99% > ali >= 60% && err <= 1.5% */
	      else if (1000.0 * sp->err <= 30 * sp->ali)
		sp->dubious = sp->dubious2 = flag ; /* ( ali > 60 % && 1.5% < err <= 3%) */
	      else
		sp->bad = sp->bad2 = flag ;
	    }
	  else
	    {
	      if (1000.0 * sp->err <= 4* sp->ali)
		sp->dubious = sp->dubious3 = flag ; /*( ali < 60 % && err <= .4 %)  */
	      else if (1000.0 * sp->err <= 15 * sp->ali)
		sp->dubious = sp->dubious4 = flag ; /*( ali < 60 % && err <= 1.5 %)  */
	      else if (1000.0 * sp->err <= 30 * sp->ali)
		sp->bad = sp->bad3 = flag ;
	      else 
		sp->bad = sp->bad4 = flag ;
	    }
	  if (sp->ali < 200 && 100.0 * sp->ali < 60 * ln && 1000.0 * sp->err <= 30 * sp->ali)
	    sp->bad_short = flag ;
	}
      else if (sp->g_split && sp->g_major)
	{
	  int ln = sp->ln ;
	  if (100.0 * sp->g_ali >= 99 * ln)
	    {
	      if (1000.0 * sp->g_err <= 4 * sp->g_ali)
		sp->excellent = 1 ; /* ali >= 99% && err <= .4% */
	      else  if (1000.0 * sp->g_err <= 15 * sp->g_ali)
		sp->good = 1 ;/* ali >= 99% && .4% < err <= 1.5 */
	      else if (1000.0 * sp->g_err <= 30 * sp->g_ali)
		sp->dubious = sp->dubious1 = 1 ;
	      else
		sp->bad = sp->bad1 = 1 ; 
	    }
	  else if (100.0 * sp->g_ali >= 60 * ln)
	    {
	      if (1000.0 * sp->g_err <= 4 * sp->g_ali)
		sp->partial = sp->partial1 = 1 ;   /* 99% > ali >= 60% && err <= .4% */
	      else if (1000.0 * sp->g_err <= 15 * sp->g_ali)
		sp->partial = sp->partial2 = 1 ;   /* 99% > ali >= 60% && err <= 1.5% */
	      else if (1000.0 * sp->g_err <= 30 * sp->g_ali)
		sp->dubious = sp->dubious2 = 1 ; /* ( ali > 60 % && 1.5% < err < 3%) */
	      else
		sp->bad = sp->bad2 = 1 ;
	    }
	  else
	    {
	      if (1000.0 * sp->g_err <= 4 * sp->g_ali)
		sp->dubious = sp->dubious3 = 1 ; /*( ali < 60 % && ali > 200 && err < .4 %)  */
	      else if (1000.0 * sp->g_err <= 15 * sp->g_ali)
		sp->dubious = sp->dubious4 = 1 ; /*( ali < 60 % && ali > 200 && err < 1.5 %)  */
	      else if (1000.0 * sp->g_err <= 30 * sp->g_ali)
		sp->bad = sp->bad3 = 1 ;
	      else 
		sp->bad = sp->bad4 = 1 ;
	    }
	  if (sp->g_ali < 200 && 100.0 * sp->g_ali < 60 * ln && 1000.0 * sp->g_err <= 30 * sp->g_ali)
	    sp->bad_short = 1 ;
	}
    }
  
  /* now we detect the diamond genes */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (!sp->excellent)
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{ 
	  sq = arrayp (aa, i2, STAT) ;
	  if (sp != sq && (sq->best || sq->best5))
	    break ;
	}
      if (i2 >= jj)
	sp->diamond = 1 ;
    }
  
  /************* BUG search ********************/
  /* split, mosaic, inversion */
  for (ii = j0 ; ii < jj ; ii++)
    {
      int dx = 0 ;
      
      sp = arrayp (aa, ii, STAT) ;
      if (!sp->g_major)
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ;
	  if (!sq->g_minor)
	    continue ;
	  if (sp->chrom != sq->chrom ||
	      (sp->a1 - sq->a1 < -300000 || sp->a1 - sq->a1 > 300000)
	      )
	    sp->mosaic = sq->mosaic = 1 ;
	  else if ( /* 2 strands */
		   (sp->a1 > sp->a2 && sq->a1 < sq->a2) ||
		   (sp->a1 < sp->a2 && sq->a1 > sq->a2)
		   )
	    { sp->inversion += sq->ln ; sq->inversion += sq->ln ; }
	  /* else same chrom, same strand */
	  else if (gExonOverlap (sp, sq, &dx)) /* the 2 genes pile up on each other */
	    { sp->v_pileUp += dx ; sq->v_pileUp += dx ; }
	  else if ( /* they are in the correct order */
		   (sp->a1 < sp->a2 && sp->a2 < sq->a1 && sp->x1 < sp->x2 && sp->x2 < sq->x1 + 100) ||
		   (sp->a1 < sp->a2 && sq->a2 < sp->a1 && sp->x1 < sp->x2 && sq->x2 < sp->x1+ 100) ||
		   (sp->a1 < sp->a2 && sp->a2 < sq->a1 && sp->x1 > sp->x2 && sp->x2 > sq->x1+ 100) ||
		   (sp->a1 < sp->a2 && sp->a2 < sq->a1 && sp->x1 > sp->x2 && sp->x2 > sq->x1+ 100) ||
		   
		   (sp->a1 > sp->a2 && sp->a2 > sq->a1 && sp->x1 < sp->x2 && sp->x2 < sq->x1+ 100) ||
		   (sp->a1 > sp->a2 && sq->a2 > sp->a1 && sp->x1 < sp->x2 && sq->x2 < sp->x1+ 100) ||
		   (sp->a1 > sp->a2 && sp->a2 > sq->a1 && sp->x1 > sp->x2 && sp->x2 > sq->x1+ 100) ||
		   (sp->a1 > sp->a2 && sp->a2 > sq->a1 && sp->x1 > sp->x2 && sp->x2 > sq->x1+ 100) 
		   )
	    sp->fuse_with = sq->fuse_with = 1 ;
	  else
	    sp->transposition = sq->transposition = 1;
	}
    }
  
  /************* BUG search ********************/
  /* Find Twins */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->twin)
	continue ;
      for (i2 = ii + 1 ;  i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ;
	  if (sp->method == sq->method &&
	      sp->chrom == sq->chrom &&
	      sp->start == sq->start &&
	      sp->a1 == sq->a1 &&
	      sp->a2 == sq->a2)
	    sq->twin = 1 ;
	}
    }

  /* wrong_map, dubious_map */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ; /* the possible bad guy */
      if (sp->alibaba || sp->best || sp->best5 ||
	  sp->no_score || sp->bad_score || 
	  ! sp->m_best ||
	  sp->v_repeat || sp->v_pileUp || sp->m_split)
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ; /* the golden boy */
	  if (!sq->alibaba || sq->g_split || !sq->best)
	    continue ;
	  if (sp->method == sq->method || gGeneBoxOverlap (sp, sq, TRUE))
	    continue ;
	  {
	    /* check if another ali of method sp matches any other (repeated) alibaba */
	    STAT *sp2, *sq2 ;
	    int i3, i4 ;
	    BOOL ok = FALSE ;

	    for (i3 = j0 ; !ok && i3 < jj ; i3++)
	      {
		sp2 = arrayp (aa, i3, STAT) ;
		if (sp->method != sp2->method)
		  continue ;
		for (i4 = j0 ; !ok && i4 < jj ; i4++)
		  {
		    sq2 = arrayp (aa, i4, STAT) ;
		    if (!sq2->alibaba)
		      continue ;
		    if (gGeneBoxOverlap (sp2, sq2, TRUE))
		      ok = TRUE ;		
		  }
	      }
	    if (ok)
	      continue ;
	  }

	  if (sq->excellent)
	    sp->wrong_map = 1 ;
	  else if (sq->good)
	    sp->dubious_map_b = 1 ;
	  else /* if (sq->partial) */
	    sp->dubious_map_c = 1 ;
	}
    }
  
  /* ambiguous now we flag the ambiguous: same method aligns same cdna part in several places */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      for (i2 = ii + 1 ; i2 < jj ; i2++)
	{
	  sq = arrayp (aa, i2, STAT) ;
	  if (sq->method == sp->method && !sp->kill_method &&
	      gCloneOverlap (sp, sq, 0))  /* no ambiguous kill */
	    sp->ambiguous = sq->ambiguous = 1 ;
	}
    }

  /* alternate/repeated_in_n : meme genome, meme score mais introns different */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->no_score) continue ;
      for (i2 = ii + 1 ; i2 < jj ; i2++)
	{
	  if (ii == i2)
	    continue ; 
	  sq = arrayp (aa, i2, STAT) ; 
	  if (sq->no_score) continue ;
	  if (!sp->alibaba && !sq->alibaba)
	    continue ;
	  if (!gGeneBoxOverlap (sp, sq, TRUE))
	    continue ;
	  if (gSameIntrons (sp, sq))
	    {
	      if (1 || sp->alibaba)
		sp->identical_in_n++ ;
	      if (1 || sq->alibaba)
		sq->identical_in_n++ ;
	    }
	  else if (sp->score == sq->score)
	    { sp->alternate++ ; sq->alternate++ ; }
	}
    }

  /* now we analyse the defects */
  for (ii = j0 ; ii < jj ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      if (sp->no_score || sp->bad_score || sp->v_repeat || sp->v_pileUp || sp->m_split)
	continue ;
      for (i2 = j0 ; i2 < jj ; i2++)
	{
	  if (ii == i2)
	    continue ;
	  sq = arrayp (aa, i2, STAT) ;
	  
	  if (sq->alibaba && !sq->g_split && !sq->v_repeat && !sp->v_pileUp && !sq->m_split && gGeneBoxOverlap (sp, sq, TRUE))
	    {
	      /* get sp->missed_exon/wrong_exon/more_exon/missed_5p/missed_central/missed_3p */
	      gExonAnalysis (sq, sp) ; 
	    }
	}
    } 

  /* export */
  ac_free (h) ;
  return 1 ;
} /* gFlagEst */

/*************************************************************************************/

static void gFlagExport (AC_DB db, DICT *dict, Array aa, char **methodTitle)
{
  STAT *sp ;
  int ii ;
  /* export keyset */
  AC_HANDLE h = ac_new_handle () ;

  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      sp = arrp (aa, ii, STAT) ;
      freeOutf ("\nSequence %s\n", dictName (dict, sp->est)) ;
      freeOutf ("From_gene \"%s\" %d %d %d %d\n\n"
		, dictName (dict, sp->gene)
		, sp->ln, sp->start, sp->ali, sp->err) ;
      freeOutf ("\nTranscribed_gene \"%s\"\n", dictName (dict, sp->gene)) ;
      freeOutf ("-D Rank\n-D Bug\n-D Split\n") ;
      if (sp->alibaba)
	freeOutf ("Per_mille %d\n", (int)(1000.0 * (1.0 - ((float)sp->ln - sp->ali)/sp->ln))) ;
      if (sp->score)
	freeOutf ("Score %d %d\n", sp->score, sp->best_score - sp->score) ;
      if (sp->m_score)
	freeOutf ("M_score %d %d\n", sp->m_score, sp->best_score - sp->m_score) ;
      if (sp->bad_score)
	freeOutf ("Bad_score\n") ;
      if (sp->no_score)
	freeOutf ("No_score\n") ;
      if (sp->g_score && sp->alibaba)
	freeOutf ("Gold_score %d\n", sp->g_score) ;
      if (sp->best)
	freeOutf ("Best\n") ;
      if (sp->best5)
	freeOutf ("Best5\n") ;
      if (sp->best50)
	freeOutf ("Best50\n") ;
      if (sp->best500)
	freeOutf ("Best500\n") ;
      if (sp->best5000)
	freeOutf ("Best5000\n") ;

      if (sp->alibaba == 2 && ! sp->g_minor)
	freeOutf ("Alibaba \"%s\"\n", dictName (dict, sp->est)) ;
      if (sp->local_best)
	freeOutf ("Local_best\n") ;
      if (sp->m_best)
	freeOutf ("Best_in_method\n") ;
      if (sp->diamond)
	freeOutf ("Diamond\n") ;

      if (sp->kill_method)
	freeOutf ("-D method\n") ;
      else if (sp->method == 'A' || sp->method == 'B' || sp->method == 'C' || sp->method == 'D')
	freeOutf ("AceView\n") ;

      if (!sp->g_minor)
	{
	  if (sp->excellent)
	    freeOutf ("Excellent\n") ;
	  if (sp->good)
	    freeOutf ("Good\n") ;
	  if (sp->partial)
	    freeOutf ("Partial\n") ;
	  if (sp->partial1)
	    freeOutf ("Partial1\n") ;
	  if (sp->partial2)
	    freeOutf ("Partial2\n") ;
	  if (sp->dubious)
	    freeOutf ("Dubious\n") ;
	  if (sp->dubious1)
	    freeOutf ("Dubious1\n") ;
	  if (sp->dubious2)
	    freeOutf ("Dubious2\n") ;
	  if (sp->dubious3)
	    freeOutf ("Dubious3\n") ;
	  if (sp->dubious4)
	    freeOutf ("Dubious4\n") ;
	  if (sp->bad)
	    freeOutf ("Bad\n") ;
	  if (sp->bad1)
	    freeOutf ("Bad1\n") ;
	  if (sp->bad2)
	    freeOutf ("Bad2\n") ;
	  if (sp->bad3)
	    freeOutf ("Bad3\n") ;
	  if (sp->bad4)
	    freeOutf ("Bad4\n") ;
	  if (sp->bad_short)
	    freeOutf ("Bad_short\n") ;
	}

      if (sp->twin)
	freeOutf ("Twin\n") ;
      if (sp->wrong_strand)
	freeOutf ("Wrong_strand\n") ;
      if (sp->dubious_strand)
	freeOutf ("Dubious_strand\n") ;
      if (sp->masked)
	freeOutf ("Masked\n") ;
      if ((sp->v_repeat + sp->v_pileUp) && sp->alibaba)
	freeOutf ("Variable_repeat %d\n", sp->v_repeat + sp->v_pileUp) ;
      if (0 && sp->g_split)  /* no need to report this */
	freeOutf ("Split\n") ;
      if (sp->mosaic)
	freeOutf ("Mosaic\n") ;
      if (sp->micro_deletion)
	freeOutf ("Micro_intron %d\n", sp->micro_deletion) ;
      if (sp->inversion)
	freeOutf ("Inversion %d\n", sp->inversion) ;
      if (sp->transposition)
	freeOutf ("Transposition\n") ;
      if (sp->fuse_with)
	freeOutf ("Fuse_with\n") ;
      if (sp->v_pileUp)
	freeOutf ("Pile_up %d\n", sp->v_pileUp) ;
      if (!sp->g_split && !sp->m_split && sp->genomic_contamination)
	freeOutf ("Genomic_contamination\n") ;
      if (sp->ambiguous)
	freeOutf ("Ambiguous\n") ;

      if (sp->orf)
	freeOutf ("ORF %d AA\n", sp->orf) ;
      if (sp->orf_bonus < 0)
	freeOutf ("ORF_malus %d \n", - sp->orf_bonus) ;
      else if (sp->orf_bonus == 0)
	freeOutf ("ORF_bonus\n") ;

      if (sp->polyA_exon)
	freeOutf ("PolyA_exon\n") ;
      if (sp->polyA_first_exon)
	freeOutf ("PolyA_first_exon\n") ;
      if (sp->polyA_first_exon)
	freeOutf ("PolyA_first_exon\n") ;
      if (sp->vector_exon)
	freeOutf ("Vector_exon\n") ;
      if (sp->too_short_first_exon)
	freeOutf ("Too_short_first_exon %d\n", sp->too_short_first_exon) ;

      if (sp->g_major)
	freeOutf ("Major\n") ;
      if (sp->g_minor)
	freeOutf ("Minor\nGold\n") ;


      if (sp->g_forward)
	freeOutf ("G_forward\n") ;
      if (sp->g_reverse)
	freeOutf ("G_reverse\n") ;
    
      if (sp->Duplicate_of_read)
	freeOutf ("Duplicate_of_read \"%s\"\n", dictName (dict, sp->Duplicate_of_read)) ;
      
      if (sp->m_major)
	freeOutf ("M_major\n") ;
      if (sp->m_minor)
	freeOutf ("M_minor\n") ;

      if (sp->wrong_map)
	freeOutf ("wrong_map\n") ;
      if (sp->dubious_map_b)
	freeOutf ("dubious_map_B\n") ;
      if (sp->dubious_map_c)
	freeOutf ("dubious_map_C\n") ;

       if (sp->identical_in_n)
	freeOutf ("Identical_in_n %d\n", sp->identical_in_n) ;

       if (sp->alibaba || sp->g_minor) /* i.e. exact repeats and split pieces */
	 {
	   if (sp->excellent)
	     freeOutf ("AliName Gold_excellent\n") ;
	   else if (sp->good)
	     freeOutf ("AliName Gold_good\n") ;
	   else if (sp->partial)
	     freeOutf ("AliName Gold_partial\n") ;
	   else if (sp->dubious)
	     freeOutf ("AliName Gold_dubious\n") ;
	   else
	     freeOutf ("AliName Gold_bad\n") ;
	 }
       else
	 freeOutf ("AliName \"%s\"\n", methodTitle [sp->iMethod]) ;
       if (sp->alibaba == 2)
	 {
	   if (sp->alternate)
	     freeOutf ("Alternate %d\n", sp->alternate) ;
	   if (sp->repetition)
	     freeOutf ("repetition %d\n", sp->repetition) ;
	   if (sp->repetition5)
	     freeOutf ("Repetition_5 %d\n", sp->repetition5) ;
	   if (0 && sp->repetition50)
	     freeOutf ("Repetition_50 %d\n", sp->repetition50) ;
	 }
       else
	 {
	   if (sp->alternate)
	     freeOutf ("Alternate\n") ;
	   if (sp->repetition)
	     freeOutf ("repetition\n") ;
	   if (sp->repetition5)
	     freeOutf ("Repetition_5\n") ;
	 } 

       if (sp->genome_deletion)
	freeOutf ("Genome_deletion %d\n", sp->genome_deletion) ;

       if (sp->too_long)
	freeOutf ("Too_long\n") ;
      if (sp->missed_exon)
	freeOutf ("Missed_exon %d\n", sp->missed_exon) ;
      if (sp->more_exon)
	freeOutf ("More_exon %d\n", sp->more_exon) ;
      if (sp->wrong_exon)
	freeOutf ("Wrong_exon\n", sp->wrong_exon) ;

      if (sp->missed_5p)
	freeOutf ("Missed_5p\n", sp->missed_5p) ;
      if (sp->missed_central)
	freeOutf ("Missed_central\n", sp->missed_central) ;
      if (sp->missed_3p)
	freeOutf ("Missed_3p\n", sp->missed_3p) ; 

      if (sp->missed_intron)
	freeOutf ("Missed_intron %d\n", sp->missed_intron) ; 
      if (sp->more_intron)
	freeOutf ("More_intron %d\n", sp->more_intron) ; 

      if (sp->missed_gling)
	freeOutf ("Missed_micro_indel %d\n", sp->missed_gling) ;
      if (sp->more_gling)
	freeOutf ("More_micro_indel %d\n", sp->more_gling) ;

      if (sp->more_5p)
	freeOutf ("More_5p %d\n", sp->more_5p) ;
      if (sp->more_central)
	freeOutf ("More_central %d\n", sp->more_central) ;
      if (sp->more_3p)
	freeOutf ("More_3p %d\n", sp->more_3p) ; 
      if (sp->exon_tested)
	freeOutf ("Exon_tested\n") ;

      if (sp->identical_exon)
	freeOutf ("Identical_exon %d\n", sp->identical_exon) ; 
    }
  freeOutf ("\n// done \n\n") ;  

  ac_free (h) ;
}  /* gFlagExport */

/*************************************************************************************/

static int gFlagAll (int *methods, char **methodTitle, AC_DB db, char *template, int errCost)
{
  AC_ITER iter = 0 ; 
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (500000, h) ;
  AC_OBJ est ;
  int nn = 0, nn1 = 0 ;
  Array aa = arrayHandleCreate (500000, STAT, h) ;
  char *qq = messprintf ("find est IS \"%s\"", template) ;

  if (template && strstr (template, "ctf")) qq = "find est ctf_file" ;
  if (0) qq = "find tg  alibaba ; split  ; > read " ;
  iter = ac_query_iter (db, TRUE, qq, 0, h) ; /* all alignments */
  while ((est = ac_next_obj (iter)))
    {
      nn1 += gFlagEst (methods, methodTitle, aa, dict, est, errCost) ;
      ac_free (est) ;
      if (nn1 > 200)
	{
	  nn += nn1 ; nn1 = 0 ;
	  gFlagExport (db, dict, aa, methodTitle) ; /* export all acedb flags */
	  aa = arrayReCreate (aa, 500000, STAT) ;
	}
    }
  nn += nn1 ;
  gFlagExport (db, dict, aa, methodTitle) ; /* export all acedb flags */
  ac_free (h) ;

  return nn ;
} /* gFlagAll */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*
  sur les Gold, 
  export number of ali and of clones,  %, #split and repeated clones
  foreach qaulity
  excellent, good, partial, dubious,  Bad, Bad_short
	bad -----> show table
*/

static void gStatsAliBabaQuality (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  int iii, nr, ng, nsg, nrep ;
  int nNoScore = 0 ;
  int nrt, ngt, nrtot, nsgt, nrept ;
  AC_KEYSET ks, ks0, ksR ;
  char *quals[] = { "Excellent", "Good", "Partial"
		  , "Dubious", "Bad", 0 };
 
  
  gOutSection (xml, "Document 6: Quality of the gold alignments\n") ;
  if (xml)
    {
      freeOutf ("For convenience, the categories are directly linked as html tables. ") ;
      freeOutf ("But for a thourough analysis, we also provide a "
		"complete excel-compatible tab delimited "
		"<a href=\"DATA/gold_quality.txt.gz\"> document</a>.<br>") ;

    }
  gOutTableStart (xml,"GoldQuality") ;
  if (xml)
    gOutTitles (xml, "Quality\tGold<br>alignments<sup>1</sup>\tmRNA<br>accession\t%mRNA\tmRNA with<br> rearranged<br> GOLD\tmRNA with<br> exact<br>repetitions\n") ;
  else
    gOutTitles (xml, "Quality\tGold\tmRNA\t%mRNA\tSplit\tRepeat\n") ;
 
  ks = ac_dbquery_keyset (db, "find est", h) ;
  nrtot = ac_keyset_count (ks) ;
  nrt = ngt = nsgt = nrept = 0 ;
  if (1) ks0 = ac_dbquery_keyset (db, "find tg gold && !minor", h) ; 
  else ks0 = ac_dbquery_keyset (db, "find est  IS af155662 ; > from_gene ; Gold", h) ;

  for (iii = 0 ; quals[iii] ; iii++)
    {
      ks = ac_ksquery_keyset (ks0, quals[iii], h) ;
      ng = ac_keyset_count (ks) ;
      ksR = ac_ksquery_keyset (ks, ">read",  h) ; 
      if (xml) gExportTable (db, ksR, "Gold", 0, quals[iii], "Gold") ;
      nr = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, "repetition ; >read",  h) ;
      nrep = ac_keyset_count (ksR) ;
      if (xml) gExportTable (db, ksR, "Repetition", 0, quals[iii], "Gold") ;
      ksR = ac_ksquery_keyset (ks, "rear ; >read",  h) ;
      nsg = ac_keyset_count (ksR) ;
      if (xml) gExportTable (db, ksR, "Rearrangement", 0, quals[iii], "Gold") ;
     
      /* we have a few (5) no score genes, we add them in bad */
      if (!strcmp (quals[iii], "Bad"))
	{	  
	  ksR = ac_dbquery_keyset (db, "Find est from_gene && ! alibaba", h) ;
	  nNoScore = ac_keyset_count (ksR) ;
	  ng += nNoScore ; nr += nNoScore ;
	}
      gOutTitle (xml, quals[iii]) ;
      gOutInt (xml, ng, FALSE, 0) ;
      if (iii >= 0)
	{
	  gOutIntLink2 (xml, nr, FALSE, 1, "Gold", quals[iii]) ;
	  gOutInt (xml, nr, TRUE, nrtot) ;
	  gOutIntLink2 (xml, nsg, FALSE, 0, "Rearrangement", quals[iii]) ;
	  gOutIntLink2 (xml, nrep, FALSE, 0, "Repetition", quals[iii]) ;
	}
      else
	{
	  gOutInt (xml, nr, FALSE, 1) ;
	  gOutInt (xml, nr, TRUE, nrtot) ;
	  gOutInt (xml, nsg, FALSE, 0) ;
	  gOutInt (xml, nrep, FALSE, 0) ;
	}
      gOutNewLine (xml) ;
      ngt += ng ; nrt += nr ; nsgt += nsg ; nrept += nrep ;
    }
  
  ksR = ac_dbquery_keyset (db,"find est ! from_gene", h) ;
  nr = ac_keyset_count (ksR) ;
  if (xml && nr) 
    gExportTable2 (db, ksR, 0
		   , messprintf ("Gold/Gold.Unaligned.%d.htm", nr)
		   , "alibaba.est_list.def") ;
  nrt += nr ;
  
  gOutTitle (xml, "Unaligned") ;
      gOutInt (xml, 0, FALSE, 0) ;
      gOutIntLink2 (xml, nr, FALSE, 1, "Gold", "Unaligned") ;
      gOutInt (xml, nr, TRUE, nrtot) ;
      gOutInt (xml, 0, FALSE, 0) ;
      gOutInt (xml, 0, FALSE, 0) ;
      gOutNewLine (xml) ;

      gOutTitle (xml, "Total") ;
      gOutInt (xml, ngt, FALSE, 0) ;
      gOutInt (xml, nrt, FALSE, 0) ;
      gOutInt (xml, nrt, TRUE, nrtot) ;
      gOutInt (xml, nsgt, FALSE, 0) ;
      gOutInt (xml, nrept, FALSE, 0) ;
      gOutNewLine (xml) ;

      gOutTableEnd (xml) ;

  if (!xml && nNoScore)
    freeOutf("\n\nWe added %d no-score mRNAs in the category Bad\n\n", nNoScore) ;

  if (0) /* export bad table */
    {
      ac_command_keyset (db
     , messprintf ("table -active -x -href www.aceview.org/av.cgi?db=gold&a=gold -title -o COMPARE/AliBaba/%d.bad.%shtm -f TABLES/Est2TgCond.def %s"
		   , 0
		   , ""
		   , " \"IS *\" ")
			 , ksR, h) ;
    }

  freeOutf ("\n\n") ;
  ac_free (h) ;
} /* gStatsAliBabaQuality */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*
  sur les Gold, 
  export the number in each category of ali/error
  i.e. a subdivision of the tags bad/dubious
  do not consider repeats
*/

static void gStatsAliBabaQuality2 (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  int line, col, nr ;
  int nNoScore = 0 ;
  int nrt = 0, nrtot = 0, nline , ncol[4] ; 
  char *lineTitle[] = { "<=4 ", "4 to 15", "15 to 30" , "> 30", "Total"} ;

  AC_KEYSET ks, ks0, ksR ;
  char *qual ;
  char *quals[] = { "Excellent", "Partial1", "Dubious3",
		    "Good",      "Partial2", "Dubious4",
		    "Dubious1",  "Dubious2", "Bad3",
                    "Bad1",      "Bad2",     "Bad4" } ;
  
  gOutSection (xml, "Document 6: Quality of the gold alignments\n") ;
  if (xml)
    {
      freeOutf ("For convenience, the categories are directly linked as html tables. ") ;
      freeOutf ("But for a thourough analysis, we also provide a "
		"complete excel-compatible tab delimited "
		"<a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/DATA/gold_quality.txt.gz\"> document</a>.<br>") ;

    }
  gOutTableStart (xml,"GoldQuality2") ;
  gOutTitles (xml, "err\\ali\t>=99%\t>=60%\t<60%\t0%\tTotal") ;
 
  ks = ac_dbquery_keyset (db, "find est", h) ;
  nrtot = ac_keyset_count (ks) ;

  if (1) ks0 = ac_dbquery_keyset (db, "find tg alibaba", h) ; 
  else ks0 = ac_dbquery_keyset (db, "find est  IS af155662 ; > from_gene ; Gold", h) ;
  
  ncol[0] =  ncol[1] = ncol[2] = ncol[3] = 0 ;
  for (line = 0 ; line < 4 ; line++)
    {
      nline = 0 ;
      gOutTitle (xml, lineTitle [line]) ;
      for (col = 0 ; col < 3 ; col++) 
	{
	  qual = quals [3 * line + col] ;
	  if (3 * line + col == 11) /* we have a few (78) no score genes, we add them in bad4 */ 
	    {
	      ksR = ac_dbquery_keyset (db, "Find est from_gene && ! alibaba", h) ;
	      nNoScore = ac_keyset_count (ksR) ;/* so we can report this number as a Note */
	      /* as one set so we can table-export it */
	      ksR = ac_dbquery_keyset (db, "{Find tg alibaba && !minor && Bad4 ; > read } $| {Find est from_gene && !alibaba}", h) ;
	      nr = ac_keyset_count (ksR) ;
	    }
	  else
	    {
	      ks = ac_ksquery_keyset (ks0, qual, h) ;
	      nr = ac_keyset_count (ks) ;
	    }
	  nrt += nr ; nline += nr ; ncol[col] += nr ;
	  if (line || col)
	    {  
	      ksR = ac_ksquery_keyset (ks, ">read",  h) ; 
	      gExportTable (db, ksR, "Gold", 0, qual, "Gold") ;
	      gOutIntLink2 (xml, nr, FALSE, 1, "Gold", qual) ;
	    }
	  else
	    gOutInt (xml, nr, FALSE, 1) ;
	}
      if (line == 3) /* add the unaligned */
	{
	  ksR = ac_dbquery_keyset (db,"find est ! from_gene", h) ;
	  nr = ac_keyset_count (ksR) ;
	  nrt += nr ; nline += nr ; ncol[col] += nr ;
	  if (xml && nr) 
	    gExportTable2 (db, ksR, 0
			   , messprintf ("Gold/Gold.Unaligned.%d.htm", nr)
			   , "alibaba.est_list.def") ;
	  gOutIntLink2 (xml, nr, FALSE, 1, "Gold", "Unaligned") ; 
	}
      else
	gOutInt (xml, 0, FALSE, 1) ;
      gOutInt (xml, nline, FALSE, 1) ;
      gOutNewLine (xml) ;
    }

  gOutTitle (xml, lineTitle [line]) ;
  for (col = 0 ; col < 4 ; col++) 
    gOutInt (xml, ncol[col], FALSE, 1) ;
  gOutInt (xml, nrt, FALSE, 1) ;
  
  gOutTableEnd (xml) ;
  freeOutf("\n\n Total %d/%d\n", nrt, nrtot) ;
  if (!xml && nNoScore)
    freeOutf("\n\nWe added %d no-score mRNAs in the category Bad\n\n", nNoScore) ;
  

  freeOutf ("\n\n") ;
  ac_free (h) ;
} /* gStatsAliBabaQuality2 */

/*************************************************************************************/
/*
  sur les gold
  export number of repetitions
  foreach quality
  excellent, good, partial, dubious
	bad -----> show table

*/

static void gStatsAliBabaRepeats (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  int iii, nr, ng, nrtot, n2, nrand, n99 ;
  AC_KEYSET ks, ks0, ks1, ksR ;
  char *quals[] = { "Excellent", "Good", "Partial", "(Excellent || good || partial)",  0 };
  char *title[] = { "Excellent", "Good", "Partial", "Total",  0 };
 
  ks = ac_dbquery_keyset (db,"find tg alibaba && (Excellent || good || partial) && ! split", h) ;
  nrtot = ac_keyset_count (ks) ;
  gOutSection (xml
	       , messprintf ("Analysis of the exact repetitions,\n excluding split alignments\nout of %d mRNA in this class", nrtot)
	       ) ;
  gOutTableStart (xml, "GoldRepeats1") ;
  gOutTitles (xml, "Quality\tGold\tmRNA\t2 copies\t>2 copies\n") ;
 
  if (1) ks0 = ac_dbquery_keyset (db,"find tg gold", h) ;
  else ks0 = ac_dbquery_keyset (db,"find est bc029540 ; > from_gene ; Gold", h) ;

  for (iii = 0 ; quals[iii] ; iii++)
    {
      ks1 = ac_ksquery_keyset (ks0, quals[iii], h) ;
      ks = ac_ksquery_keyset (ks1, "repetition ; alibaba", h) ;
      ksR = ac_ksquery_keyset (ks, ">read",  h) ;
      nr = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, ">read; >from_gene ; repetition", h) ;
      ng = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, "COUNT {>read; >from_gene ; repetition} = 2",  h) ;
      n2 = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, "COUNT {>read; >from_gene ; repetition} > 2",  h) ;
      n99 = ac_keyset_count (ksR) ;

      gOutTitle (xml, title[iii]) ;
      gOutInt (xml, ng, FALSE, 1) ;
      if (iii < 3)
	gOutIntLink2 (xml, nr, FALSE, 0, "Repetition", title[iii]) ;
      else
	gOutInt (xml, nr, FALSE, 0) ;
      gOutInt (xml, n2, FALSE, 1) ;
      gOutInt (xml, n99, FALSE, 1) ;

      gOutNewLine (xml) ;
   }
  gOutTableEnd (xml) ;

  gOutSection (xml
	       , messprintf ("Analysis of the exact or nearly exact  repetitions,\n excluding split alignments\nOut of %d E/G/P mRNA", nrtot)
	       ) ;
  gOutTableStart (xml, "GoldRepeats2") ;
  gOutTitles (xml, "Quality\tAlignments\tIn random contig\tmRNA\t2 copies\t>2 copies\n") ;
 
  ks = ac_dbquery_keyset (db,"find est ", h) ;
  nrtot = ac_keyset_count (ks) ;
  if (1) ks0 = ac_dbquery_keyset (db,"find tg gold", h) ;
  else ks0 = ac_dbquery_keyset (db,"find est bc029540 ; > from_gene ; Gold", h) ;

  mkdir ("COMPARE/Repetition_5", 755) ;
  system ("chmod 755 COMPARE/Repetition_5") ;

  for (iii = 0 ; quals[iii] ; iii++)
    {
      ks1 = ac_ksquery_keyset (ks0, quals[iii], h) ;
 
      ks = ac_ksquery_keyset (ks1, "repetition_5  || repetition; alibaba", h) ;    
      ksR = ac_ksquery_keyset (ks, ">read",  h) ;
      nr = ac_keyset_count (ksR) ;   
      if (xml) 
	gExportTable2 (db, ksR, 0
		       , messprintf ("Repetition_5/Repetition_5.%s.%d.htm", title[iii], nr)
		       , "Est2TgCond.def \"(repetition_5  || repetition)\"") ;
      ksR = ac_ksquery_keyset (ks, ">read; >from_gene ; repetition_5 || repetition", h) ;
      ng = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, ">read; >from_gene ; repetition_5 || repetition ; !(IntMap == ? || IntMap = ?? )", h) ;
      nrand = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, "COUNT {>read; >from_gene ; repetition_5 || repetition} = 2",  h) ;
      n2 = ac_keyset_count (ksR) ;
      ksR = ac_ksquery_keyset (ks, "COUNT {>read; >from_gene ; repetition_5 || repetition} > 2",  h) ;
      n99 = ac_keyset_count (ksR) ;

      gOutTitle (xml, title[iii]) ;
      gOutInt (xml, ng, FALSE, 1) ;
      gOutInt (xml, nrand, FALSE, 1) ;
      if (!xml)
	gOutInt (xml, nr, FALSE, 1) ;
      else
	gOutIntLink2 (xml, nr, FALSE, 1, "Repetition_5", title[iii]) ;
      gOutInt (xml, n2, FALSE, 1) ;
      gOutInt (xml, n99, FALSE, 1) ; 

      gOutNewLine (xml) ;
    }
  gOutTableEnd (xml) ;

  freeOutf ("\n\n") ;
  ac_free (h) ;
} /* gStatsAliBabaRepeats */

/*************************************************************************************/
/*
  sur les gold split
*/

static void gStatsAliBabaSplit (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  int iii, nr, ng, ngtot, nrtot ;
  AC_KEYSET ks, ks0, ksR ;
  char *types[] = { "Mosaic", "Inversion", "Transposition", "Genome_deletion", "Variable_repeat", "Fuse_with", "rear", 0 } ;
  char *title[] = { "Mosaic", "Inversion", "Transposition", "Suspected genome deletion", "Variable tandem repeat number", "Pseudo-split or mRNA bridging two correctly ordered and oriented contigs", "Total rearrangements",  0 } ;
 
  gOutSection (xml, "Analysis of the split gold alignments, excellent good or partial,\n excluding repetitions") ;
  
  gOutTableStart (xml, "GoldRear") ;
  gOutTitles (xml, "Type\tmRNA\tParts\t%mRNA") ;
 
  ks = ac_dbquery_keyset (db,"find tg alibaba && (Excellent || good || partial) && !repetition ", h) ;
  nrtot = ac_keyset_count (ks) ;
  ks0 = ac_ksquery_keyset (ks,"Rear ; > read ; >from_gene ; gold", h) ;

  mkdir ("COMPARE/Split", 755) ;
  system ("chmod 755 COMPARE/Split") ;

  if (!nrtot) nrtot = 1 ;

   for (iii = 0 ; types[iii] ; iii++)
    {
      ks = ac_ksquery_keyset (ks0, types[iii], h) ;
      ng = ac_keyset_count (ks) ;
      ksR = ac_ksquery_keyset (ks, ">read",  h) ;
      nr = ac_keyset_count (ksR) ;
      if (xml) 
	gExportTable2 (db, ksR, 0
		       , messprintf ("Split/Split.%s.%d.htm", types[iii], nr)
		       , "Est2TgCond.def") ;      
      gOutTitle (xml, title[iii]) ;
      if (!xml) 
	gOutInt (xml, nr, FALSE, 1) ; 
      else
	gOutIntLink2 (xml, nr, FALSE, 1, "Split", types[iii]) ;
      gOutInt (xml, ng, FALSE, 1) ;
      gOutInt (xml, nr, TRUE, nrtot) ;

      gOutNewLine (xml) ;
    }

   ks0 = ac_dbquery_keyset (db,"find tg gold", h) ; /* Gold, !alibaba, since we count the parts */
   ngtot = ac_keyset_count (ks0) ;

   gOutTitle (xml,"Total alignments at this quality") ;
   gOutInt (xml, nrtot, FALSE, 1) ;
   gOutInt (xml, ngtot, FALSE, 1) ;
   gOutInt (xml, nrtot, TRUE, nrtot) ;
   gOutNewLine (xml) ;

   gOutTableEnd (xml) ;
  ac_free (h) ;
} /* gStatsAliBabaRepeats */

/*************************************************************************************/
/*************************************************************************************/
#ifdef JUNK
static void gDTag (BOOL xml, AC_DB db, char* wanted, char *g_qual, char *a_rank
		   , int *methods, char **methodTag
		   , char *fileTitle, char *tag, char *g_tag)
{
  AC_HANDLE h = ac_new_handle () ;
  int ng, iMethod ;
  AC_KEYSET ks0, ks1, ksR ;

  if (g_tag)
    ks0 = ac_dbquery_keyset (db
			     , messprintf ("find tg alibaba && %s; >read; >from_gene ; %s && %s && best_in_method && ! m_minor && COUNT {>read;>from_gene;alibaba  && %s && %s} > 0"
					   , g_tag, tag, a_rank, g_qual, g_tag)
			     , h) ; 
  else
    ks0 = ac_dbquery_keyset (db
			     , messprintf ("find tg %s && %s && best_in_method && ! m_minor && COUNT {>read;>from_gene;alibaba  && %s} > 0"
					   ,tag, a_rank, g_qual)
			     , h) ; 

  for (iMethod = 0 ; methodTag [iMethod] ; iMethod++)
    {
      ks1 = ac_ksquery_keyset (ks0
			       , methodTag [iMethod] /* messprintf ("IS %c_*", methods [iMethod]) */
			       , h) ;
      ng = ac_keyset_count (ks1) ; 
      
      if (xml)
	{ 
	  char htmName[800], bw[3] ;
	  ksR = ac_ksquery_keyset (ks1, ">Read", h) ;
      
	  strncpy (bw, wanted, 2) ;
	  bw[2] = 0 ;
	  if (bw[1] == ' ') bw[1] = 0 ;
	  if (bw[1] == '*') bw[1] = 0 ; 

	  sprintf (htmName, "COMPARE/%s", fileTitle) ;
	  mkdir (htmName, 755) ;
	  system (messprintf ("chmod 755 %s", htmName)) ;
	  sprintf (htmName,
		   "%s/%c.%s%s.%d.htm"
		   , fileTitle
		   , methods [iMethod]
		   , fileTitle
		   , bw
		   , ng
		   ) ;
      
      	  ac_command_keyset (db
			     , messprintf ("table -active -x -href www.aceview.org/av.cgi?db=gold&a=gold -title -o COMPARE/%s -f TABLES/Est2TgCond.def \"gold || %s\""
					   , htmName
					   , methodTag [iMethod])
			     , ksR, h) ;
	  
	  if (ng> 0)
	    freeOutf ("<td><a href=\"%s\">%5d</a></td>",htmName, ng ) ;
	  else
	    freeOutf ("<td>%5d</td>", ng ) ;
	} 
      else
	freeOutf ("\t%5d", ng) ;
    }
  ac_free (h) ;
} /* gDTag */

/*****************************************************************************/

static void gStatsAliBabaDoDefects (BOOL xml, AC_DB db
				    , int *methods, char **methodTag
				    , char **wantedList
				    , char *fileTitle, char *tag, char *g_tag)
{
  char g_qual_symbol [] = { 'A', 'B', 'C', 0 } ;
  char *g_qual [] = { "Excellent", "Good", "(Gold_rank && !Excellent && !good)", 0 } ;
  static int xColor = 0 ;
  char a_rank_symbol [] = { ' ', '1', '2', '3', '4', '5', 0 } ;
  char *a_rank [] = { "*", "u", "best5000", "best500", "best50", "best5", 0 } ;
  
  char wanted[3] ;
  int iQual, iRank, iWanted ; 
  
  xColor = (xColor + 1) % 2 ;
  for (iWanted = 0 ; wantedList[iWanted] ; iWanted++)
    for (iQual = 0 ; g_qual_symbol [iQual] ; iQual++)
      for (iRank = 0 ; a_rank_symbol [iRank] ; iRank++)
	{
	  if (a_rank_symbol [iRank] == '1') 
	    continue ;
	  wanted[0] = g_qual_symbol [iQual] ; 
	  wanted[1] = a_rank_symbol [iRank] ;
	  wanted[2] = 0 ;
	  if (strcmp ( wantedList[iWanted], wanted))
	    continue ;
	  if (xml)
	    {
	      if (xColor)
		freeOut ("<tr bgcolor=#efefff>") ;
	      else
		freeOut ("<tr>") ;
	      freeOutf ("<td>%s</td><td>%s</td>", fileTitle, wanted) ;
	    }
	  else freeOutf ("%12s\t%s", fileTitle, wanted) ;
	  gDTag (xml, db, wanted, g_qual [iQual], a_rank [iRank], methods, methodTag, fileTitle, tag, g_tag) ;

	  if (xml) freeOutf ("</tr>\n") ;
	  else freeOutf ("\n") ;
	}  
} /* gStatsAliBabaDoDefects */

/*****************************************************************************/
/*****************************************************************************/
/*       
        identical_in_n = 0, 1, 2, 3, 4, 5...  -----> show table if 0
        alternate     -----> show table of alibaba + alternate
*/

static void gStatsAliBabaDefects (AC_DB db, BOOL xml, int *methods, char **methodTag)
{ 
  gOutSection (xml, "Detailled tables of differences between each method and the gold alignments") ;
  gOutTableStart (xml, "PgmDefectsOld") ;
  
  { /* titles */
    vTXT blkp = vtxtHandleCreate (0) ;
    int iMethod ;
    
    vtxtPrintf (blkp, "Defect\tType") ;
    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
      vtxtPrintf (blkp, "\t%s", methodTag [iMethod]) ;
    gOutNewLine (xml) ;
  }
  

  if (!xml)
  { /* missed, should be A5 and B5 */
    /* done here because the query is very different : ! best-in_method  */
    char g_qual_symbol [] = { 'A', 'B', 'C', 0 } ;
    char *g_qual [] = { "Excellent", "Good", "(Gold_rank && !Excellent && !good)", 0 } ;
    int ng, iQual, iMethod ;
    AC_HANDLE h = ac_new_handle () ;
    AC_KEYSET ks1, ks2 ;

    for (iQual = 0 ; g_qual [iQual] ; iQual++)
      {
	freeOutf ("%12s\t%c", "Unaligned", g_qual_symbol [iQual]) ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  {
	    ks1 = ac_dbquery_keyset (db
				     , messprintf ("find tg alibaba && %s" 
						   , g_qual [iQual])
				     , h) ;
	    ks2 = ac_ksquery_keyset (ks1
				     , messprintf ("COUNT {>read; >from_gene ; %s } = 0"
						   , methodTag[iMethod])
				     , h) ;
	    ng = ac_keyset_count (ks2) ;
	    freeOutf ("\t%5d", ng) ;
	  }
	freeOutf ("\n") ;
      }
    ac_free (h) ;
  }

  if (0)
    { /* missed_repeat */
      char *wanted [] = { "A ", "B ", "C ", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "missed_repeat"
			      , "!ambiguous"
			      , "Repetition") ;
    }
  
  
  if (!xml)
    { /* redundant_losing */
      /* done here because the query is very different : ! best-in_method  */
      char g_qual_symbol [] = { 'A', 'B', 'C', 0 } ;
      char *g_qual [] = { "Excellent", "Good", "(Gold_rank && !Excellent && !good)", 0 } ;
      int ng, iQual, iMethod ;
      AC_HANDLE h = ac_new_handle () ;
      AC_KEYSET ks1, ks2 ;
      
      for (iQual = 0 ; g_qual [iQual] ; iQual++)
	{
	  freeOutf ("%12s\t%c", "Redundant", g_qual_symbol [iQual]) ;
	  for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	    {
	      ks1 = ac_dbquery_keyset (db
				       , messprintf ("find tg !best_in_method && !no_score && !bad_score && !m_split && %s" 
						     , methodTag[iMethod])
				       , h) ;
	      ks2 = ac_ksquery_keyset (ks1
				       , messprintf ("COUNT {>read; >from_gene ; alibaba ; !split && !repetition && %s } > 0"
						     ,  g_qual [iQual])
				       , h) ;
	      ng = ac_keyset_count (ks2) ;
	      freeOutf ("\t%5d", ng) ;
	    }
	  freeOutf ("\n") ;
	}
      ac_free (h) ;
    }
  
  { /* wrong_map */
    char *wanted [] = { "A ", "B ", "C ", 0 } ;
    gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			    , "Wrong_map"
			    , "(Wrong_map || dubious_map_b || dubious_map_c)"
			    , 0) ;
  }

  /* fin du tableau Hit or Miss */


  /* tableau exon choice defects */

  if (1)
    { /* missed_exon */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Missed_exon"
			      , "(missed_exon && !Missed_3p && !Missed_5p && !Missed_central)"
			      , 0) ;
    }

  if (1)
    { /* wrong_exon */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Inaccurate_exon", "Wrong_exon", 0) ;
    }
  
  if (1)
    { /* more_exon */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Too_many_exons"
			      , "(More_exon && !PolyA_exon && !PolyA_first_exon)"
			      , 0) ;
    }

  if (0)
    { /* losing */
      char *wanted [] = { "A2", "A3", "A4", "A5", "B2", "B3", "B4", "B5", "C2", "C3", "C4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Losing"
			      , "*"
			      , 0) ;
    }
  
  if (0)
    { /* missed_5p */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Missed_5p"
			      , "Missed_5p"
			      , 0) ;
    }
  
  if (0)
    { /* missed_central */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Missed_central", "Missed_central"
			      , 0) ;
    }

  if (0)
    { /* missed_3p */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Missed_3p", "Missed_3p"
			      , 0) ;
    }


  { /* vector_exon */
    char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
    gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			    , "vector_exon"
			    , "vector_exon"
			    , 0) ;
  }

  { /* polyA-exons */
    char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
    gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			    , "PolyA_exon"
			    , "( PolyA_exon || PolyA_first_exon)"
			    , 0) ;
  }

  if (0) /* contained in previous case */
  { /* PolyA_first_exon */
    char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
    gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			    , "PolyA_first_exon"
			    , "PolyA_first_exon"
			    , 0) ;
  }

  if (0)
    { /* missed_split */
      char *wanted [] = { "A ", "B ", "C ", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "missed_split"
			      , "! m_split"
			      , "Split") ;
    }

  if (0)
    { /* Too_short_first_exon */
      char *wanted [] = { "A2", "A3", "A4", "B2", "B3", "B4", 0 } ;
      gStatsAliBabaDoDefects (xml, db, methods, methodTag, wanted
			      , "Too_short_first_exon"
			      , "Too_short_first_exon"
			      , 0) ;
    }


  gOutTableEnd (xml) ;
} /* gStatsAliBabaDefects */
#endif

/*****************************************************************************/

static void gDoExportTable (AC_DB db, AC_KEYSET ks, char *filNam, char *tableNam)
{
  AC_HANDLE h = ac_new_handle () ;
  char buf [1024] ;

  sprintf (buf
	   , "table -active -x -href www.aceview.org/av.cgi?db=gold&a=gold -title -o COMPARE/%s -l -f TABLES/%s"
	   , filNam, tableNam
	   ) ;

  ac_command_keyset (db, buf, ks, h) ;
  ac_free (h) ;
} /* gDoExportTable */
	  
static void gExportTable (AC_DB db, AC_KEYSET ks, char *dirTitle, int method, char *fileTitle, char *methodTag)
{
  int ng = ac_keyset_count (ks) ;
  char buf[1000] ;
  char htmName [1024] ;
  
  sprintf (htmName, "COMPARE/%s", dirTitle) ;
  mkdir (htmName, 755) ;
  system (messprintf ("chmod 755 %s", htmName)) ;

  if (method)
    sprintf (htmName,
	     "%s/%c.%s.%d.htm"
	     , dirTitle
	     , method
	     , fileTitle ? fileTitle : dirTitle
	     , ng
	     ) ;
  else if (fileTitle)
    sprintf (htmName,
	     "%s/%s.%s.%d.htm"
	     , dirTitle, dirTitle
	     , fileTitle ? fileTitle : dirTitle
	     , ng
	     ) ;

  sprintf (buf, "Est2TgCond.def \"gold || %s\"" , methodTag) ;
  gDoExportTable (db, ks, htmName, buf) ;
  return ;
}

static int gExportTable2 (AC_DB db, AC_KEYSET ks0,  char *myQuery, char *fileTitle, char *tableNam)
{
  AC_KEYSET ks = ks0 ;
  int nn ;
  if (!ks) ks = ac_dbquery_keyset (db, myQuery, 0) ;
  nn = ac_keyset_count (ks) ;
  gDoExportTable (db, ks, fileTitle, tableNam) ;
  if (!ks0) ac_free (ks) ; 
  return nn ;
}

static void gExportTable3 (AC_DB db, AC_KEYSET ks, char *dirTitle,  char *fileTitle, int method, int dummy, char *methodTag)
{
  int ng = ac_keyset_count (ks) ;
  char buf[1000] ;
  char htmName [1024] ;
  
  sprintf (htmName, "COMPARE/%s", dirTitle) ;
  mkdir (htmName, 755) ;
  system (messprintf ("chmod 755 %s", htmName)) ;

  if (method)
    sprintf (htmName,
	     "%s/%c.%s.%d.htm"
	     , dirTitle
	     , method
	     , fileTitle ? fileTitle : dirTitle
	     , ng
	     ) ;
  else if (fileTitle)
    sprintf (htmName,
	     "%s/%s.%s.%d.htm"
	     , dirTitle, dirTitle
	     , fileTitle ? fileTitle : dirTitle
	     , ng
	     ) ;

  sprintf (buf, "Est2TgCond.def \"gold || %s\"" , methodTag) ;
  gDoExportTable (db, ks, htmName, buf) ;
  return ;
}

static void gExportTable4 (AC_DB db, AC_KEYSET ks, char *dirTitle,  char *library, char *defect)
{
  char buf[1000] ;
  char htmName [1024] ;
  
  if (ks && ac_keyset_count (ks))
    {
      sprintf (htmName, "COMPARE/%s", dirTitle) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName, "COMPARE/%s/%s", dirTitle, library) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName,"%s/%s/%s.%d.htm", dirTitle, library, defect, ac_keyset_count (ks)) ;
      
      sprintf (buf, "Est2TgCond.def gold") ;
      gDoExportTable (db, ks, htmName, buf) ;
    }
  return ;
}

static void gExportTable5 (AC_DB db, AC_KEYSET ks, char *dirTitle,  char *library, char *defect)
{
  char buf[1000] ;
  char htmName [1024] ;
  
  if (ks && ac_keyset_count (ks))
    {
      sprintf (htmName, "COMPARE/%s", dirTitle) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName, "COMPARE/%s/%s", dirTitle, library) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName,"%s/%s/%s.%d.htm", dirTitle, library, defect, ac_keyset_count (ks)) ;
      
      sprintf (buf, "alibaba.clone_list.def") ;
      gDoExportTable (db, ks, htmName, buf) ;
    }
  return ;
}

static void gExportTable6 (AC_DB db, AC_KEYSET ks, char *dirTitle,  char *library, char *defect)
{
  char buf[1000] ;
  char htmName [1024] ;
  
  if (ks && ac_keyset_count (ks))
    {
      sprintf (htmName, "COMPARE/%s", dirTitle) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName, "COMPARE/%s/%s", dirTitle, library) ;
      mkdir (htmName, 755) ;
      system (messprintf ("chmod 755 %s", htmName)) ;
      
      sprintf (htmName,"%s/%s/%s.%d.htm", dirTitle, library, defect, ac_keyset_count (ks)) ;
      
      sprintf (buf, "alibaba.read_list.def") ;
      gDoExportTable (db, ks, htmName, buf) ;
    }
  return ;
}

/*****************************************************************************/
/*****************************************************************************/

static void gStatsABND1 (AC_DB db, int *methods, char **methodTag, char **methodTitle)
{
  int nn, iMethod ;
  int n, nPgm = 7 ;
  AC_KEYSET ks1, ks2 ;

  ks1 = ac_dbquery_keyset (db, "find tg alibaba && excellent && !split", 0) ;
  nn = ac_keyset_count (ks1) ;
  
  freeOutf ( "<h3><span style=\"font-weight: 400\"><font size=\"6\" face=\"Verdana\">Lists of \n"
	     "alignments that could be improved...</font></span></h3><h3>\n"
	     "<span style=\"font-weight: 400\"><font face=\"Verdana\">We compared the alignments \n"
	     "to the GOLD alignments that were selected anonymously from any of the\n"
	     " %d programs who contributed. We \n"
	     "restricted the comparison to the %d \n"
	     "accessions aligning completely (&gt;99%%) and with very few errors (&lt;0.4%% \n"
	     "differences from the genome). Mosaic clones, or clones with inversions or \n"
	     "transpositions were also excluded: all the mRNAs here should align in a \n"
	     "straightforward fashion all the way with very few errors.</font></span></h3>\n"
	    , nPgm, nn
	     ) ;
   for (iMethod = n = 0 ; methods [iMethod] ; iMethod++)
     if (methods [iMethod] != 1)
       n++ ;
   for (iMethod = 0 ; n < 3 && methods [iMethod] ; iMethod++)
     if (methods [iMethod] != 1)
       {
	 ks2 = ac_ksquery_keyset (ks1, messprintf (">read ; tested_by_%s", methodTag[iMethod]), 0) ;
	 nn = ac_keyset_count (ks2) ;
	 ac_free (ks2) ;
	 freeOutf ( "<h3><span style=\"font-weight: 400\"><font face=\"Verdana\">%s sent us data for\n"
		   " %d of those clones. </font>\n"
		   "</span></h3>\n"
		   "<h3><font face=\"Verdana\"><span style=\"font-weight: 400\"><font size=\"5\">1:&nbsp; \n"
		    "Accuracy of  the exon choice</font><br>\n"
		   "&nbsp;</span></font></h3><br>\n"
		   , methodTitle [iMethod], nn
		   ) ;
       }
   ac_free (ks1) ;
} /* gStatsABND1 */

/*****************************************************************************/

static void gStatsABND2 (void)
{
  freeOutf ("%s",
	    "<p><font face=\"Verdana\"><span style=\"font-weight: 400\"><font size=\"4\">The line &quot;Total&quot; \n"
	    "gives the number of accessions that you mapped to the same place on the genome \n"
	    "as the GOLD, and that could thus be compared to the GOLD in terms of exon and \n"
	    "intron structure.&nbsp; An exon is &quot;missed&quot; if it is part of the GOLD but your \n"
	    "alignment has no overlapping exon. An exon is &quot;extra&quot; if it belongs to your \n"
	    "alignment but does not overlap any exon in the GOLD. The first exon is the \n"
	    "5'-most, the last the 3'-most, central is any exon in between the first and \n"
	    "last. An inaccurate exon overlaps a gold exon, but has at least one intron \n"
	    "boundary not optimally chosen according to our analysis. Sometimes, the gold \n"
	    "solution is very close... The second column gives the number of accessions with \n"
	    "the defect. The same accession may appear in more than one line in the table.</font></span></font></p>\n"
	    "<p><font face=\"Verdana\" size=\"4\">Clicking on any blue number brings the \n"
	    "corresponding list of accessions with some useful information. The table \n"
	    "compares, for each accession, the main features of your alignment to the GOLD. \n"
	    "The length, position and number of errors of the match, number of exons and \n"
	    "exact map position on the genome are reported. </font></p>\n"
	    "<p class=\"MsoNormal\" style=\"line-height: 150%; margin-left: 36pt\">\n"
	    "<span style=\"color: black; font-family: Arial\">&nbsp;</span></p>\n"
	    "<h3><font face=\"Verdana\"><span style=\"font-weight: 400\"><font size=\"5\">2:&nbsp; \n"
	    "Hit or miss: problems finding the map on the genome</font></span><br>\n"
	    "</font></h3>\n"
	    "<br>\n"
	    ) ;
} /* gStatsABND2 */

/*****************************************************************************/

static void gStatsABND3 (void)
{
  freeOutf ("%s",
	    "<p><font size=\"4\" face=\"Verdana\">This table lists all accessions where finding \n"
	    "the place in the genome where the cDNA aligns was problematic. If you tested an \n"
	    "accession, but failed to produce any alignment, the accession is &quot;Unaligned&quot;. An \n"
	    "accession that you aligned only in a position where the quality of the match was \n"
	    "far less good than the GOLD (50 points below) is qualified of &quot;Wrong map&quot;. An \n"
	    "accession for which you provided two or more independent alignments of the same \n"
	    "part of the clone, but one was far less good than the GOLD (hence its map \n"
	    "position is not correct) is called &quot;Redundant&quot;: it might be better to discard \n"
	    "it. </font></p>\n"
	    "<h3><font face=\"Verdana\"><span style=\"font-weight: 400\"><font size=\"5\">3: \n"
	    "Missing exact repetitions and not recognizing mosaics and other split alignments</font></span><br>"
	    "</font></h3>\n"
	    ) ;
} /* gStatsABND3 */

/*****************************************************************************/
/*       
        identical_in_n = 0, 1, 2, 3, 4, 5...  -----> show table if 0
        alternate     -----> show table of alibaba + alternate
*/

static void gStatsAliBabaNewDefects (AC_DB db, BOOL xml, int *methods, char **methodTag, char **methodTitle)
{
  int iMethod ;
  AC_KEYSET ks1, ks2, ks3, ks1EG ;
  AC_HANDLE h = ac_new_handle () ;
  vTXT blkp = vtxtHandleCreate (h) ;
  int percent ;
  int nOk[64], nTotal[64],
    nMissed[64], nMissedRepeats[64], nMissedSplit[64], 
    nMissed5p[64], nMissedCentral[64], nMissed3p[64], 
    nMore5p[64], nMoreCentral[64], nMore3p[64], 
    nWrongMap[64],  nRedundant[64], 
    nMissedExons[64], nWrongExons[64],
    /* nWrongExons5[64],  */
    nMoreExons[64],
    nMissedGling [64], nMoreGling [64]
    ; 
  
  vtxtPrint (blkp, "Defect\t") ;
  for (iMethod = 0 ; methods [iMethod] ; iMethod++)
    if (methods [iMethod] != 1)
      vtxtPrintf (blkp, "%s\t", methodTitle [iMethod]) ;
  

  ks1EG = ac_dbquery_keyset (db
			     , "find tg alibaba && !split && (Excellent) ; >read"
			     , h) ;
  /* table 1:
   * exon defects : against alibaba excellent good but !split (was before multiple)
   *  - missed exon: the gold contacts the tg but one exon of the gold does not contact the tg
   *  - too many exon:  the gold contacts the tg but one exon of the tg does not contact the gold
   *  - wrong nexon: the global scores differ and a gold exon contacts the tg exon but the intron boundaries differ
   */
 
  { /* total */
    ks1 = ks1EG ;    
    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
      {	
	if (methods [iMethod] == 1)
	  continue ;
	/* denominator */
	ks2 = ac_ksquery_keyset (ks1
				 , messprintf (">from_gene %s && exon_tested"
					       , methodTag[iMethod])
				 , h) ;
	nOk [iMethod] =  ac_keyset_count (ks2) ; 
	ac_free (ks2) ;
	
	/* missed_exon */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && missed_exon && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMissedExons [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (0 && xml) gExportTable3 (db, ks2, "Problem", "Missed_exon", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* missed_5p */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && missed_5p && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMissed5p [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Missed_5p", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* missed_central */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && missed_central && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMissedCentral [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Missed_central", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* missed_3p */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && missed_3p && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMissed3p [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Missed_3p", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	
	/* more_exon */ 
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && more_exon && %s && !best ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMoreExons [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (0 && xml) gExportTable3 (db, ks2, "Problem", "More_exon", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;

	/* more_5p */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && more_5p && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMore5p [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "More_5p", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* more_central */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && more_central && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMoreCentral [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "More_central", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* more_3p */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && more_3p && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMore3p [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "More_3p", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	
	/* wrong_exon */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && wrong_exon && %s && !best ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nWrongExons [iMethod] = ac_keyset_and (ks2, ks1) ;
	ac_free (ks2) ;

	/* wrong_exon 5
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && wrong_exon && !best && !best5 && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nWrongExons5 [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Inaccurate_exon", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
*/
	/* missed_gling */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && Missed_micro_indel && !best && !best5 && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMissedGling [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Missed_micro_indel", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;

	/* more_gling */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("Find tg best_in_method && More_micro_indel && !best && !best5 && %s ; >read"
					       , methodTag[iMethod])
				 , h) ;
	nMoreGling [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "More_micro_indel", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;

      }

    
    
  for (percent = 0 ; percent < (xml ? 1 : 2) ; percent++)
      {
	if (xml)
	  gStatsABND1 (db, methods, methodTag, methodTitle) ;
	else
	  {
	    gOutSection (xml, "Comparison of the alignment programs: 1: accuracy of of the exon choice") ;
	    freeOutf ("Only excellent non-split gold alignments are considered\n") ;
	    gOutBreak (xml) ;
	    freeOutf ("Total is the number of genes of this method overlapping the tested gold alignment\n") ;
	    gOutBreak (xml) ;
	  }
	
	gOutTableStart (xml, "PgmQuality") ; 
	gOutTitles (xml, vtxtPtr (blkp)) ;
	
	gOutTitle (xml, "Total") ; 
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutInt (xml, nOk[iMethod], FALSE, 1) ;
	gOutNewLine (xml) ;	
	
	if (!xml)
	  {
	    gOutTitle (xml, "Missed exon") ;
	    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	      if (methods [iMethod] != 1)
		gOutIntLink (xml, nMissedExons [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_exon") ;
	    gOutNewLine (xml) ;	
	  }
	gOutTitle (xml, "missed first exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissed5p [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_5p") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "missed central exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissedCentral [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_central") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "missed last exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissed3p [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_3p") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "inaccurate exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nWrongExons [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Inaccurate exon") ;
	gOutNewLine (xml) ;	

	if (! xml)
	  {
	    gOutTitle (xml, "Too many exon") ;
	    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	      if (methods [iMethod] != 1)
		gOutIntLink (xml, nMoreExons [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "More_exon") ;
	    gOutNewLine (xml) ;	
	  }
	gOutTitle (xml, "extra first exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMore5p [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "More_5p") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "extra central exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMoreCentral [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "More_central") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "extra last exon") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMore3p [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "More_3p") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "missed micro in/del") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissedGling [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_micro_indel") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "extra micro in/del") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMoreGling [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "More_micro_indel") ;
	gOutNewLine (xml) ;	
	
	
	
	gOutTableEnd (xml) ;
      }
  }


  /* table 2:
   * hit or miss table : relative to alibaba excellent || good and ! multiple)
   * - unaligned and 
   * - wrong map 
   * - redundant
   */

  { /* total */
    ks1 = ks1EG ;
    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
      {
	if (methods [iMethod] == 1)
	  continue ;

	ks2 = ac_ksquery_keyset (ks1
				 , messprintf ("COUNT {>from_gene ; %s } > 0"
					       , methodTag[iMethod])
				 , h) ;
	nOk [iMethod] = ac_keyset_count (ks2) ; 

	ks3 = ac_ksquery_keyset (ks1
				 , messprintf ("Tested_by_%s"
					       , methodTag[iMethod])
				 , h) ;
	nTotal [iMethod] =  ac_keyset_count (ks3) ; 
	/* missed */
	nMissed [iMethod] = ac_keyset_minus (ks3, ks2) ;
	if (xml) gExportTable3 (db, ks3, "Problem", "Unaligned", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	ac_free (ks3) ;

	if (nOk [iMethod] < 0) nOk [iMethod] = 1;

	/* wrong_map (alibaba = excellent) */
	/* dubious_map_b  (alibaba = good) */
	/* limited to lose by at least 50 */
	ks2 = ac_dbquery_keyset (db
				 , messprintf ("find tg  %s &&  ! best50 &&  (Wrong_map || dubious_map_b) ; > read"
					       , methodTag[iMethod])
				 , h) ;
	nWrongMap [iMethod] = ac_keyset_and (ks2, ks1) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Wrong_map", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;

	
	/* redundant */
 	ks2 = ac_ksquery_keyset (ks1
				 , messprintf ("COUNT {>from_gene !no_score && !bad_score && !m_split && %s && ! best_in_method && !best && !best5 && !best50 } >= 1 ; COUNT {>from_gene gold && !multiple } = 1"
					       , methodTag[iMethod])
				 , h) ;
	nRedundant [iMethod] = ac_keyset_count (ks2) ;
	if (xml) gExportTable3 (db, ks2, "Problem", "Redundant", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ; 
      }
  }
  
  for (percent = 0 ; percent < (xml ? 1 : 2) ; percent++)
    {
      if (xml)
	gStatsABND2 () ;
      else
	{
	  gOutSection (xml, "Comparison of the alignment programs 2: hit or miss") ;
	  freeOutf ("Only excellent non-split gold alignments are considered\n") ;
	  gOutBreak (xml) ;
	}
      gOutTableStart (xml, "PgmMaps") ; 
      gOutTitles (xml, vtxtPtr (blkp)) ;
      
      gOutTitle (xml, "Unaligned") ;
      for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nMissed [iMethod], percent, nTotal[iMethod], methods [iMethod], "Problem", "Unaligned") ;
      gOutNewLine (xml) ;	
      
      gOutTitle (xml, "Wrong map") ;
      for (iMethod = 0 ; methods [iMethod] ; iMethod++)	
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nWrongMap [iMethod], percent, nTotal[iMethod], methods [iMethod], "Problem", "Wrong_map") ;
      gOutNewLine (xml) ;
          
      gOutTitle (xml, "Redundant") ;
      for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nRedundant [iMethod], percent, nTotal[iMethod], methods [iMethod], "Problem", "Redundant") ;
      gOutNewLine (xml) ;	

      gOutTitle (xml, "mRNA tested") ; 
      for (iMethod = 0 ; methods [iMethod] ; iMethod++)	
	if (methods [iMethod] != 1)
	  gOutInt (xml, nTotal[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;	
      
      gOutTableEnd (xml) ;
    }

  /* table 3: against alibaba excellent or good
   * - Missed repetition (alibaba not split)
   * - missed split (alibaba is split but not the method)
   */
  
  { /* total */
    ks1 = ac_dbquery_keyset (db
			     , "find tg alibaba && !split && (Excellent) && repetition ; > read"
			     , h) ;
    
    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
      {
	if (methods [iMethod] == 1)
	  continue ;

	/* missed_repeat */
	ks2 = ac_ksquery_keyset (ks1
				 , messprintf ("Tested_by_%s"
					       , methodTag[iMethod])
				 , h) ;
	nTotal [iMethod] =  ac_keyset_count (ks2) ; 
	ks3 = ac_ksquery_keyset (ks2
				 , messprintf ("COUNT {>from_gene ; %s ; } < 2"
					       , methodTag[iMethod])
				 , h) ;
	nMissedRepeats [iMethod] = ac_keyset_count (ks3) ;
	if (xml) gExportTable3 (db, ks3, "Problem", "Missed_repeat", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	ac_free (ks3) ;
      }
    ac_free (ks1) ;


    ks1 = ac_dbquery_keyset (db
			     , "find tg alibaba && multiple && (Excellent) ; > read"
			     , h) ;
    for (iMethod = 0 ; methods [iMethod] ; iMethod++)
      {
	if (methods [iMethod] == 1)
	  continue ;
	/* missed_split */
	ks2 = ac_ksquery_keyset (ks1
				 , messprintf ("Tested_by_%s"
					       , methodTag[iMethod])
				 , h) ;
	nOk [iMethod] =  ac_keyset_count (ks2) ; 
	ks3 = ac_ksquery_keyset (ks2
				 , messprintf ("COUNT {>from_gene ; %s ; m_split} = 0"
					       , methodTag[iMethod])
				 , h) ;
	nMissedSplit [iMethod] = ac_keyset_count (ks3) ;
	if (xml) gExportTable3 (db, ks3, "Problem", "Missed_split", methods[iMethod], 0, methodTag[iMethod]) ;
	ac_free (ks2) ;
	ac_free (ks3) ;
      }
    
    ac_free (ks1) ;
  
    
  for (percent = 0 ; percent < (xml ? 1 : 2) ; percent++)
      {
	if (xml)
	  gStatsABND3 () ;
	else
	  {
	    gOutSection (xml, "Comparison of the alignment programs: 3: Repetions, Mosaics and other split alignments") ;
	    freeOutf ("Only excellent gold alignments are considered\n") ;
	    freeOutf (" Total is the number of genes of this method overlapping the tested gold alignment\n") ;
	    gOutBreak (xml) ;
	    freeOutf (" A few repetitions mapping on random contigs were considered untested by UCSC because we had a littleproblem with their coordinate system.") ;
	    gOutBreak (xml) ;
	  }
	
	
	gOutTableStart (xml, "PgmRepeats") ; 
	gOutTitles (xml, vtxtPtr (blkp)) ;
	
	gOutTitle (xml, "mRNA with missed repetition") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissedRepeats [iMethod], percent, nTotal[iMethod], methods [iMethod], "Problem", "Missed_repeat") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "mRNA with missed split") ;
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutIntLink (xml, nMissedSplit [iMethod], percent, nOk[iMethod], methods [iMethod], "Problem", "Missed_split") ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "#mRNA with repetition tested") ; 
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutInt (xml, nTotal[iMethod], FALSE, 1) ;
	gOutNewLine (xml) ;	
	
	gOutTitle (xml, "#mRNA which are split tested") ; 
	for (iMethod = 0 ; methods [iMethod] ; iMethod++)
	  if (methods [iMethod] != 1)
	    gOutInt (xml, nOk[iMethod], FALSE, 1) ;
	gOutNewLine (xml) ;	
	
	
	gOutTableEnd (xml) ;
      }
  }
  ac_free (h) ;
} /* gStatsAliBabaNewDefects */

/*****************************************************************************/
/*****************************************************************************/
/*************************************************************************************/
/*       
        #std intron 0 1 2 3 4...70
          (
            excellent
            good
            partial
	    dubious
	    bad
          ) 

        #non-std intron 0 1 2 3 4...70
          (
            excellent
            good
            partial
	    dubious
	    bad
          ) 

        # ali-baba with 0 std and > 0 non-std
        # ali-baba with 0 std and > 0 ct_ac
  
     
      )
    )
*/
static void gStatsAliBabaIntrons (AC_DB db, BOOL xml)
{ 
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = ac_query_iter (db, TRUE, "Find tg gold && (excellent)", 0, h) ; 
  AC_TABLE gOther ;
  AC_OBJ gene = 0 ;
  int ir, 
    nn, ni,
    gt, gc, at, other, exons, 
    ngt, ngc, nat, nother,
    nOneGood, nOneOther, nOneExon, nJustBad ;
  /*
    AC_ITER genes ;
    DICT suspect = dictHandleCreate (1000, h) ;
    
    genes = ac_query_iter (db, 1, "find cdna_clone Suspected_internal_deletion ; > read ; >from_gene", 0, h) ;
    while (ac_free(gene), (gene = ac_next_obj (genes)))
     dictAdd (suspect, ac_name (gene), 0) ;
    ac_free (genes) ;
  */
  nn = ni = ngt = ngc = nat = nother = nOneGood = nOneOther = nOneExon = nJustBad = 0 ;
  while (ac_free (gene), (gene = ac_next_obj (iter)))
    {
    /*
      if (dictFind (suspect, ac_name (gene), 0))
      continue ;
  */
      nn++ ;
      ngt += (gt = ac_tag_int (gene, "gt_ag", 0)) ;
      ngc += (gc = ac_tag_int (gene, "gc_ag", 0)) ;
      nat += (at = ac_tag_int (gene, "at_ac", 0)) ;
      other = 0 ; 
      if ((gOther = ac_tag_table (gene, "Other", h)))
	{ /* schema is :  Other Text Int, example: Other at_at 1 */
	  for (ir = 0 ; ir < gOther->rows ; ir++)
	    other += ac_table_int (gOther, ir, 1, 0) ;
	  ac_free (gOther) ;
	}
      nother += other + ac_tag_int (gene, "ct_ac", 0) ;

      ni += gt + gc + at + other ;

      exons = ac_tag_int (gene, "Nb_possible_exons", 0) ;

      if (gt + gc + at)
	nOneGood++ ;
      if (other)
	nOneOther++ ;
      if (exons == 1)
	nOneExon++ ;
      if (exons > 1 && !(gt + gc + at))
	nJustBad++ ;
    }

  gOutSection (xml, "Statistics on the presence or absence of introns in the GOLD alignments\n") ;
  freeOutf ("Out of the %d excellent GOLD alignments,", nn) ;
  gOutBreak (xml) ;
/*      freeOutf ("excluding suspected internal deletions") ;
	gOutBreak (xml) ;
*/
  freeOutf ("including exact duplicates (but not the almost exact duplicates)") ;
  gOutBreak (xml) ;
  gOutTableStart (xml, "GoldIntrons") ;
  gOutTitles (xml, "Alignment has\t#alignments\t%%alignments\n") ;

  gOutTitle (xml, "At least one standard intron") ; gOutInt (xml, nOneGood, FALSE, 0) ; gOutInt (xml, nOneGood, TRUE, nn) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "no intron and no alignment gap") ; gOutInt (xml, nOneExon, FALSE, 0) ; gOutInt (xml, nOneExon, TRUE, nn) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "only non-standard introns or gaps") ; gOutInt (xml, nJustBad, FALSE, 0) ; gOutInt (xml, nJustBad, TRUE, nn) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "at least one non-standard intron") ; gOutInt (xml, nOneOther, FALSE, 0) ; gOutInt (xml, nOneOther, TRUE, nn) ;
  gOutNewLine (xml) ;

  gOutTableEnd (xml) ;

  gOutSection (xml, "Statistics on the types of intron boundaries, on the same sample of excellent alignments.") ;
  freeOutf ("There are on average %2.1f intron per alignment which is not intronless.", ((float)ni) / (nn - nOneExon + .01)) ;
  gOutBreak (xml) ;
  freeOutf ("The partition by type of intron boundary is as follows:") ;
  gOutBreak (xml) ;
  gOutTableStart (xml, "GoldIntrons") ;
  gOutTitles (xml, "Intron boundary\tNumber of introns\t% of all introns\n") ;

  gOutTitle (xml, "Standard gt-ag") ; gOutInt (xml, ngt, FALSE, 0) ; gOutInt (xml, ngt, TRUE, ni) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "Standard gc-ag") ; gOutInt (xml, ngc, FALSE, 0) ; gOutInt (xml, ngc, TRUE, ni) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "Standard U12 at-ac") ; gOutInt (xml, nat, FALSE, 0) ; gOutInt (xml, nat, TRUE, ni) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "Any other non-standard type#") ; gOutInt (xml, nother, FALSE, 0) ; gOutInt (xml, nother, TRUE, ni) ;
  gOutNewLine (xml) ;
  gOutTitle (xml, "Total") ; gOutInt (xml, ni, FALSE, 0) ; gOutInt (xml, ni, TRUE, ni) ;
  gOutNewLine (xml) ;
	
  gOutTableEnd (xml) ;
  freeOutf ("# Interpretations for the high level of such anomalous pseudointrons are discussed in the article.") ;
  gOutBreak (xml) ;

  gOutBreak (xml) ;

  ac_free (h) ;
} /* gStatsAliBabaIntrons */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*       

    ========

    
  foreach method
  best5000: tuer les best5000 ressemblant a des ok de la meme methode
  no_score: tuer les no-score ressemblant a des ok de la meme methode

    ========

    
  foreach method A K S J N
  (
    best
    best + best5
    best + best5 + best50
    best + best5 + best50 + best500
    best5000 -----> show table
    
*/

static void gStatsCountBest (AC_DB db, BOOL xml, int *methods, char **methodTag, char **methodTitle, BOOL isInt)
{
  AC_HANDLE h = ac_new_handle () ;
  int 
    ngold, nall = 0, ngminor, iMethod, 
    nn0[64], nn5[64], nn50[64], nn500[64], nn5000[64], 
    nnoscore[64], nminor[64], nunaligned[64], nread[64], tnn[64] ; 
  vTXT blkp = vtxtHandleCreate (h) ;
  AC_KEYSET ks0, ks, ksB, ksR, ksR0 ;
  
  
  gOutSection (xml, "Performance of the methods: including redundant alignments") ;
  freeOutf ("Only considering mRNA with Excellent  non-split gold alignment\n") ;
  gOutBreak (xml) ;
  gOutTableStart (xml, "PgmQuality2") ;
  vtxtPrint (blkp, "Method\t") ;
  vtxtPrintf (blkp, "%s\t", "Gold") ;
  for (iMethod = 0 ; methods [iMethod] ; iMethod++)
    if ( methods [iMethod] != 1) 
      vtxtPrintf (blkp, "%s\t", methodTitle [iMethod]) ;
  gOutTitles (xml, vtxtPtr (blkp)) ;

  ksR0 = ac_dbquery_keyset (db, "find tg alibaba AND (excellent ) AND !split  ; > read ", h) ;
  nall = ac_keyset_count (ksR0) ;

  ks0 = ac_ksquery_keyset (ksR0, " >from_gene", h) ;
  for (iMethod = 0 ; methods [iMethod] ; iMethod++)
    {
      if (methods [iMethod] == 1)
	continue ;
      if (1) ks = ac_ksquery_keyset (ks0
				     , messprintf ("%s", methodTag [iMethod])
				     , h) ;
      else
	ks = ac_dbquery_keyset (db
				, messprintf ("find est bc029540 ; >from_gene  %s ", methodTag [iMethod])
				, h) ;
      ksB = ac_ksquery_keyset (ks, "best", h) ;
      nn0[iMethod] = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "best5", h) ;
      nn5[iMethod]  = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "best50 ; > read", h) ;
      if (xml) gExportTable (db, ksB, "Problem", methods[iMethod], "losing_6_50", methodTag [iMethod]) ;
      nn50[iMethod]  = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "best500; > read", h) ;
      if (xml) gExportTable (db, ksB, "Problem", methods[iMethod], "losing_51_500", methodTag [iMethod]) ;
      nn500[iMethod]  = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "best5000; > read", h) ;
      if (xml) gExportTable (db, ksB, "Problem", methods[iMethod], "losing_over_500", methodTag [iMethod]) ;
      nn5000[iMethod]  = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "no_score; > read", h) ;
      if (xml) gExportTable (db, ksB, "Problem", methods[iMethod], "no_score", methodTag [iMethod]) ;
      nnoscore[iMethod]  = ac_keyset_count (ksB) ;

      ksB = ac_ksquery_keyset (ks, "m_minor", h) ;
      nminor[iMethod]  = ac_keyset_count (ksB) ;
      	  
      ksR = ac_ksquery_keyset (ksR0, messprintf ("tested_by_%s", methodTag [iMethod]), h) ;
      nread[iMethod] = ac_keyset_count (ksR) ;
      
      ksB = ac_ksquery_keyset (ks, ">read ", h) ;	 
      nunaligned[iMethod]  = ac_keyset_minus (ksR, ksB) ;
      if (xml) gExportTable (db, ksR, "Problem", methods[iMethod], "unaligned", methodTag [iMethod]) ;
		
      tnn[iMethod] =
	nn0[iMethod] + nn5[iMethod] + nn50[iMethod] + nn500[iMethod] + nn5000[iMethod] 
	+ nnoscore[iMethod] + nunaligned[iMethod]   ;
    }

  if (isInt)
    {
      gOutTitle (xml, "Gold") ;
      ksB = ac_dbquery_keyset (db, "find tg gold AND (excellent )  && ! split ", h) ; /* was multiple */
      ngold = ac_keyset_count (ksB) ;
      gOutInt (xml, ngold, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn0[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses 1 to 5") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn5[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml,"loses 6 to 50") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 	
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nn50[iMethod], FALSE, 1, methods [iMethod], "Problem", "losing_6_50") ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses 51 to 500") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nn500[iMethod], FALSE, 1, methods [iMethod], "Problem", "losing_51_500") ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses >500") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 	
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nn5000[iMethod], FALSE, 1, methods [iMethod], "Problem", "losing_over_500") ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Unaligned") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutIntLink (xml, nunaligned[iMethod], FALSE, 1, methods [iMethod], "Problem", "unaligned") ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "No score") ;
      gOutInt (xml, 0, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nnoscore[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Total") ;
      gOutInt (xml, ngold, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 	
	if (methods [iMethod] != 1)
	  gOutInt (xml, tnn[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Split minor") ;
      ksB = ac_dbquery_keyset (db, "find tg Gold && Excellent && !split &&  minor", h) ;
      ngminor = ac_keyset_count (ksB) ;
      gOutInt (xml, ngminor, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nminor[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml,"Grand total") ;
      gOutInt (xml, ngold + ngminor, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, tnn[iMethod] + nminor[iMethod] - nunaligned[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;

      gOutTitle (xml,"mRNA tested total") ;
      gOutInt (xml, nall, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++)
 	if (methods [iMethod] != 1)
	  gOutInt (xml, nread[iMethod], FALSE, 1) ;
      gOutNewLine (xml) ;
    }
  else
    {
      gOutTitle (xml, "Gold") ;
      ksB = ac_dbquery_keyset (db, "find tg gold AND (excellent )  && ! split ", h) ; /* was multiple */
      ngold = ac_keyset_count (ksB) ;
      gOutInt (xml, ngold, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn0[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses 1 to 5") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn5[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml,"loses 6 to 50") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn50[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses 51 to 500") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 	
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn500[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "loses >500") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nn5000[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Unaligned") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nunaligned[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "No score") ;
      gOutInt (xml, 0, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 	
	if (methods [iMethod] != 1)
	  gOutInt (xml, nnoscore[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Total") ;
      gOutInt (xml, ngold, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, tnn[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml, "Split minor") ;
      ksB = ac_dbquery_keyset (db, "find tg Gold && Excellent && !split &&  minor", h) ;
      ngminor = ac_keyset_count (ksB) ;
      gOutInt (xml, ngminor, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nminor[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;
      
      gOutTitle (xml,"Grand total") ;
      gOutInt (xml, ngold + ngminor, TRUE, nall) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, tnn[iMethod] + nminor[iMethod] - nunaligned[iMethod], TRUE, nread[iMethod]) ;
      gOutNewLine (xml) ;

      gOutTitle (xml,"mRNA tested total") ;
      gOutInt (xml, nall, FALSE, 1) ;
      for (iMethod = 0 ; methods[iMethod] ; iMethod++) 
	if (methods [iMethod] != 1)
	  gOutInt (xml, nread[iMethod], FALSE, nread[iMethod]) ;
      gOutNewLine (xml) ;
    }

  gOutTableEnd (xml) ;
  ac_free (h) ;
} /* gStatsCountBest */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void  gDocsAccessionList (AC_DB db)
{
  AC_ITER iter ;
  AC_KEYSET ks ;
  AC_OBJ est = 0 ;
  int level, nn = 0 ;
  int nVectorClip, nPolyA, nInitialA, nAlternate ;
  AC_HANDLE h = ac_new_handle () ;
  FILE *f ;
  char *fNam = filName ("COMPARE/DATA/mrna_accession_list", "txt", "r") ;

  ks = ac_dbquery_keyset (db, "Find est", h) ;
  nn = ac_keyset_count (ks) ;

  if (fNam)
    goto done ;

  fNam = filName ("COMPARE/DATA/mrna_accession_list", "txt", "w") ;
  if (!fNam)
    messcrash ("cannot open COMPARE/DATA/mrna_accession_list") ;

  f = filopen (fNam, 0, "w") ;
  level = freeOutSetFile (f) ;
  

  freeOutf ("# List of %d mRNA Genbank accession used in this project\n", nn) ;

  iter = ac_keyset_iter (ks, 0, h) ;
  while (ac_free (est), (est = ac_next_obj (iter)))
    freeOutf ("%s\n", ac_name(est)) ;

  freeOutClose (level) ;
  filclose (f) ;

 done:
  freeOutf ("<h3> Document 1: <a href=\"DATA/mrna_accession_list.txt\">List</a>"
	    " of %d mRNA Genbank accession used in this project</h3><p>"
	    , nn ) ;

  freeOutf ("<h3> Document 2: Genome sequence"
	    " as used in this project</h3><p>") ;

  freeOutf ("We have used the genome data from NCBI build 34, august 2003, which also supports NCBI "
	    "annotations release 34_3 (feb 2004). The chromosome fragments as exported here were "
	    "reconstructed using the NCBI agp file, and are identical to the main chromosome fragments "
	    "available from UCSC, except that we have not remapped the floating NT contigs into the "
	    "random contigs used by UCSC.<br>\n"
	    "The chomosomes are exported here as a pair of large fasta.gz files (each expands a little below 2 Gb)"
	    "which should be downloaded by ftp."
	    "<ul><li><a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/GOLD/Genome/genome_34_3.fasta1.gz\">"
	    "genome_34_3.fasta1.gz</a>"
	    "<li><a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/GOLD/Genome/genome_34_3.fasta2.gz\">"
	    "genome_34_3.fasta2.gz</a>"
	    "</ul><p>\n") ;

  freeOutf ("<h3> Document 3:  Gold alignments resulting from this project</h3>") ;

  freeOutf ("<ul>There are 3 files,\n") ;
  freeOutf ("<li>The <a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/mRNA_alignments/gold_alignments.excellent.ali.gz\">excellent or good Gold alignments</a>\n") ;
  freeOutf ("<li>The <a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/mRNA_alignments/gold_alignments.other.ali.gz\">other Gold alignments</a> of lesser quality\n") ;
  freeOutf ("<li>The <a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/mRNA_alignments/aceview.ali.gz\">AceView alignments</a>\n") ;
  freeOutf ("<li>The <a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/mRNA_alignments/ucsc.ali.gz\">UCSC alignments</a>\n") ;
  freeOutf ("<li>The <a href=\"ftp://ftp.ncbi.nih.gov/repository/acedb/Align/mRNA_alignments/ncbi.ali.gz\">NCBI alignments</a>\n") ;
  freeOutf ("</ul><p>") ;

  mkdir ("COMPARE/DATA", 755) ;
  system ("chmod 755 COMPARE/DATA") ;

  mkdir ("COMPARE/Feature", 755) ;
  system ("chmod 755 COMPARE/Feature") ;

  ks = ac_dbquery_keyset (db, "Find est polyA_after_base", h) ;
  nPolyA = ac_keyset_count (ks) ;
  nInitialA = gExportTable2 (db, 0, "Find est Initial_polyA", "DATA/initial_polyA.htm", "alibaba.initial_polyA.def") ;
  freeOutf ("<h3> Document 4:  <a href=\"DATA/polyA.txt.gz\">%d</a> PolyA addition sites</h3>", nPolyA) ;
  gOutBreak (TRUE) ;
  freeOutf ("We also signal the existence of <a href=\"DATA/initial_polyA.htm\">%d</a> initial PolyA stretches, and clip them.", nInitialA) ;
  gOutBreak (TRUE) ;

  nVectorClip = gExportTable2 (db, 0, "Find est vector_clipping > 3", "DATA/vector_clipping.htm", "alibaba.vector_clipping.def") ;
  freeOutf ("<h3> Document 5:  <a href=\"DATA/vector_clipping.htm\">%d</a> Vector clipping </h3>", nVectorClip) ;

  nAlternate = gExportTable2 (db, 0, "Find tg alternate; > read", "Feature/alternate.htm", "Est2TgCond.def \"alternate || gold\"") ;
  freeOutf ("<h3> Document 6:  <a href=\"Feature/alternate.htm\">%d</a> Alternate </h3>", nAlternate) ;

  ac_free (h) ;
  return ;
} /* gDocsAccessionList  */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void gStats (int *methods, char **methodTag, char **methodTitle, AC_DB db, BOOL xml, char *template)
{
  mkdir ("COMPARE/DATA", 755) ;
  system ("chmod 755 COMPARE/DATA") ;

  if (xml)
    freeOut ("<document><body><font face=\"Verdana\">\n") ;
  freeOutf ("%s\n", timeShowNow()) ;

  if (1 && xml)
    {
      gDocsAccessionList (db) ;
    }
  if (1) gStatsAliBabaQuality (db, xml) ; /* figure 1 */
  if (1) gStatsAliBabaQuality2 (db, xml) ; /* figure 1 */

  if (1)
    {
      if (1) gStatsAliBabaRepeats (db, xml) ; /* figure 2a, 2b */
      if (1) gStatsAliBabaSplit (db, xml) ;   /* figure 3a */
      if (1) gStatsAliBabaIntrons (db, xml) ; /* figure 3b */
      if (1) gChromStats (db, xml, template) ; /* figure 4a, 4b */
      if (1) gLibraryStats (db, xml) ; /* figure 5a, 5b */

      if (1) gStatsCountBest (db, xml, methods, methodTag, methodTitle, TRUE) ; /* int */
      if (1 && !xml) gStatsCountBest (db, xml, methods, methodTag, methodTitle, FALSE) ; /* percents */
      if (1) gStatsAliBabaNewDefects (db, xml, methods, methodTag, methodTitle) ; /* figure 6b, 6c */
    }
  if (xml)
    freeOut ("\n\n</font></body></document>\n") ;
} /* gStats */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* These functions implement the Banded Smith Waterman algorithm
 * I do not know the original reference, sorry,
 * but just wrote this code off a blank page, so all
 * implementation errors are mine
 * It is not optimised, but seems fast enough for our purpose
 * since it will evaluate 10^4 alignments in about one hour
 *
 * If you have a better version of the same thing, i'd be grateful
 */

/*
 * Return the number of error + insertion + deletion along the
 * optimal path joining (a1,x1) to (a2,x2) but constrained to
 * stay at maximal distance w from the first diagonal (a1+i,x1+i)
 *
 * Assumes a1 < a2, x1 < x2
 */

static int smith_waterman_one_band (int w, char *longDna, int a1, int a2, char *shortDna, int x1, int x2)
{
  int z, i, j, snpCost = 1, indelCost = 1, cost, cost2, d /* direction */ ;
  char *cp, *cq ;
  int *ip, *jp ;
  int IMAX ;
  int JMAX ;
  static Array Cost = 0 ;
  
  /* do banded smith_waterman in remaining rectangle
  */

  IMAX = x2 - x1 + 2 ;
  JMAX = a2 - a1 + 2 ;
  Cost = arrayReCreate (Cost, IMAX*JMAX + 2, int) ;
  array (Cost, IMAX * JMAX + 1, int) = 0 ; /* make room */
  
  for (i = x1 - 1 ; i <= x2 ; i++)
    {
      ip = mat(i, a1 - 1) ; *ip = i - x1 + 1 ;
    }
  for (j = a1 - 1 ; j <= a2 ; j++)
    {
      ip = mat(x1 - 1, j) ; *ip = j - a1 + 1 ;
    }
  for (i = x1, cp = shortDna +  i-1 ; i <= x2 ; i++, cp++)
    for (j = a1, cq = longDna + j-1 ; j <= a2 ; j++, cq++)
      {
	if (j - a1 < i - x1 - w - 1)
	  continue ;
	if (j - a1 > i - x1 + w + 1)
	  break ;
	ip = mat(i, j) ;
	if (j - a1 < i - x1 - w || j - a1 > i - x1 + w)
	  { *ip = debug ? 99 : 9999999 ; continue ; }
	cost = (*cp == *cq) ? 0 : snpCost ;
	if (i == x1 && j == a1) { *ip = cost ; continue ; }
	*ip = -1 ;
	for (jp = 0, d = 0 ; d < 3 ; d++)
	  {
	    cost2 = (d == 0 ? cost : indelCost) ; 
	    switch (d)
	      {
	      case 0: jp = mat(i-1, j-1) ; break ;
	      case 1: jp = mat(i-1, j) ; break ;
	      case 2: jp = mat(i, j-1) ; break ;
	      }
	    if (jp)
	      {
		z = *jp + cost2 ;
		if (*ip < 0 || z < *ip) *ip = z ;
	      }
	  }
      }
  
  cost = *mat(x2, a2) ;
  if (debug)
    {
      printf("\nBANDED\n   :") ;
      for (j = a1 ; j <= a2 ; j++)
	printf ("  %c",  *(longDna+j-1)) ;
      printf ("\n") ;
      for (i = x1 ; i <= x2 ; i++)
	{
	  printf ("  %c:",  *(shortDna+i-1)) ;
	  for (j = a1 ; j <= a2 ; j++)
	    {
	      z = *mat(i, j) ;
	      printf ("%3d", z) ;
	    }
	  printf ("\n") ;
	}
      printf ("x2=%d a2=%d cost=%d\n", x2, a2, cost) ;
    }
  return cost ;
} /* smith_waterman_one_band */

/*********************************************************************/
/*
 * Return the number of error + insertion + deletion along the
 * optimal path joining (a1,x1) to (a2,x2)
 *
 * Assumes a1 < a2, x1 < x2
 * Dna coordinates start at 1 (not zero)
 */

static int smith_waterman_banded (char *longDna, int a1, int a2, char *shortDna, int x1, int x2)
{
  int w, w0, cost = debug ? 99 : 9999999 ;
  char *cp, *cq ;
  
  /* squeeze the error square starting from both corners */
  cq = longDna + a1 - 1 ; cp = shortDna + x1 - 1 ;
  if (!debug)
    {
      while (*cp == *cq && a1 < a2 && x1 < x2)
	{ cp++ ; cq++ ; a1++ ; x1++ ;}
      if (a1 == a2 && x1 == x2)
	return *cp == *cq ? 0 : 1 ;
      cq = longDna + a2 - 1 ; cp = shortDna + x2 - 1 ;
      while (*cp == *cq && a1 < a2 && x1 < x2)
	{ cp-- ; cq-- ; a2-- ; x2-- ;}
    }
  if (a1 == a2 && x1 == x2)
    return *cp == *cq ? 0 : 1 ;

  w0 = (a2 - a1) - (x2 - x1) ; /* the band must at least reach the other corner */
  if (w0 < 0) w0 = - w0 ;
  for (w = w0 ; w <= a2 - a1 || w <= x2 - x1 ; w = w ? 2 * w : 1)
    {
      cost = smith_waterman_one_band (w, longDna, a1, a2, shortDna, x1, x2) ;
      if (cost < 2 * w + 2 - w0)
	break ; 
    }
  return cost ;
} /* smith_waterman_banded */


/*********************************************************************/
/*
 * Return the number of error + insertion + deletion along the
 * optimal path joining (a1,x1) to (a2,x2)
 *
 * Assumes a1 < a2, x1 < x2
 * Dna coordinates start at 1 (not zero)
 *
 * This version is a bit simpler but quite slower than smith_waterman_banded
 */

static int smith_waterman (char *longDna, int a1, int a2, char *shortDna, int x1, int x2)
{
  int z, i, j, snpCost = 1, indelCost = 1, cost, cost2, d /* direction */ ;
  char *cp, *cq ;
  int *ip, *jp ;
  int IMAX ;
  int JMAX ;
  static Array Cost = 0 ;
  
  /* squeeze the error square starting from both corners */
  if (!debug)
    {
      cq = longDna + a1 - 1 ; cp = shortDna + x1 - 1 ;
      while (*cp == *cq && a1 < a2 && x1 < x2)
	{ cp++ ; cq++ ; a1++ ; x1++ ;}
      if (a1 == a2 && x1 == x2)
	return *cp == *cq ? 0 : 1 ;
      cq = longDna + a2 - 1 ; cp = shortDna + x2 - 1 ;
      while (*cp == *cq && a1 < a2 && x1 < x2)
	{ cp-- ; cq-- ; a2-- ; x2-- ;}
      if (a1 == a2 && x1 == x2)
	return *cp == *cq ? 0 : 1 ;
    }

  /* do full smith_waterman in remaining rectangle
     but it should not be hard to limit to a band
     at most delta off the diagonal where delta
     is the current cost on the diagonal
     or something like this
  */

  IMAX = x2 - x1 + 2 ;
  JMAX = a2 - a1 + 2 ;
  Cost = arrayReCreate (Cost, IMAX*JMAX + 2, int) ;
  array (Cost, IMAX * JMAX + 1, int) = 0 ; /* make room */
  
  for (i = x1 - 1 ; i <= x2 ; i++)
    {
      ip = mat(i, a1 - 1) ; *ip = i - x1 + 1 ;
    }
  for (j = a1 - 1 ; j <= a2 ; j++)
    {
      ip = mat(x1 - 1, j) ; *ip = j - a1 + 1 ;
    }

  for (i = x1, cp = shortDna +  i-1 ; i <= x2 ; i++, cp++)
    for (j = a1, cq = longDna + j-1 ; j <= a2 ; j++, cq++)
      {
	ip = mat(i, j) ;
	cost = (*cp == *cq) ? 0 : snpCost ;
	if (i == x1 && j == a1) { *ip = cost ; continue ; }
	*ip = -1 ;
	for (jp = 0, d = 0 ; d < 3 ; d++)
	  {
	    cost2 = (d == 0 ? cost : indelCost) ; 
	    switch (d)
	      {
	      case 0: jp = mat(i-1, j-1) ; break ;
	      case 1: jp = mat(i-1, j) ; break ;
	      case 2: jp = mat(i, j-1) ; break ;
	      }
	    if (jp)
	      {
		z = *jp + cost2 ;
		if (*ip < 0 || z < *ip) *ip = z ;
	      }
	  }
      }
  
  cost = *mat(x2, a2) ;
  if (debug)
    {
      printf("\nDIRECT\n   :") ;
      for (j = a1 ; j <= a2 ; j++)
	printf ("  %c",  *(longDna+j-1)) ;
      printf ("\n") ;
      for (i = x1 ; i <= x2 ; i++)
	{
	  printf ("  %c:",  *(shortDna+i-1)) ;
	  for (j = a1 ; j <= a2 ; j++)
	    {
	      z = *mat(i, j) ;
	      printf ("%3d", z) ;
	    }
	  printf ("\n") ;
	}
      printf ("x2=%d a2=%d cost=%d\n", x2, a2, cost) ;
    }
  return cost ;
} /* smith_waterman */


/*********************************************************************/
/*
 * Writes on the acedb stdout (freeOutf) an acefile
 * giving the Smith-Waterman cost of each exon alignment
 *
 * the schema for the tg object scanned and exported is
 * ?Transcribed_gene Genomic_sequence ?Sequence // the genome contig
 *                   Covers Int Text Int Int // coords in that contig
 *                     // eg   18324 bp 100000 118324
 *                   Smith_Waterman_done // to avoid recursions
 *                   Assembled_from Int Int ?Sequence Int Int Int
 *                               // a1  a2   mRNA     x1  x2  score
 */


static int gExportEstTgCost (AC_OBJ est, AC_OBJ tg, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  int dn, dnb, a1, a2, b1, b2, c1, c2, ii, jj, x1, x2, y1, y2, z1, z2, xmin, xmax, totalDx ;
  const char *cp, *cq ;
  char *longDna, *shortDna ;
  AC_OBJ cosmid ;
  AC_TABLE im, af = 0 ;
  DICT *dict = dictHandleCreate (2000, h) ;
  
  shortDna = ac_obj_dna (est, h) ;
  cosmid = ac_tag_obj (tg, "Genomic_sequence", h) ;
  im = ac_tag_table(tg, "Covers", h) ;
  if (im && im->cols >= 3)
    {
      a1 = ac_table_int (im, 0, 2, 0) ;  
      a2 = ac_table_int (im, 0, 3, 0) ;  
      freeOutf ("Transcribed_gene \"%s\"\n", ac_name(tg)) ;
      if (1)
	{
	  af = ac_tag_table (tg, "Assembled_from", h) ;
	  
	  for (ii = 0 ; af && ii < af->rows ; ii++)
	    {
	      est = ac_table_obj (af, ii, 2, h) ;
	      cp = ac_name (est) ;
	      if (!dictAdd (dict, cp, 0))
		continue ;
	      xmin = xmax = z1 = z2 = -999 ;
	      totalDx = 0 ;

	      for (jj = ii ; shortDna && jj < af->rows ; jj++)
		{
		  b1 = c1 = ac_table_int (af, jj, 0, 0) ;
		  b2 = c2 = ac_table_int (af, jj, 1, 0) ;
		  cq = ac_table_printable (af, jj, 2, "") ;
		  x1 = y1 = ac_table_int (af, jj, 3, 0) ;
		  x2 = y2 = ac_table_int (af, jj, 4, 0) ;
	
		  if (strcmp (cp, cq))
		    continue ;
		  if (a1 < a2)
		    { b1 = b1 + a1 - 1 ; b2 = b2 + a1 - 1 ; }
		  else
		    { b1 = a1 - b1 + 1 ; b2 = a1 - b2 + 1 ; }
		  if (x1 > x2)
		    { int x0 = x1 ; x1 = x2 ; x2 = x0 ; x0 = b1 ; b1 = b2 ; b2 = x0 ; }
		  if (xmin == -999) xmin = x1 ;
		  xmax = x2 ; longDna = 0 ;
		  if (x2 - x1 > c2 - c1 + 1000 || x2 - x1 < c2 - c1 - 1000)
		    dn = 99999 ; /* just some silly big number of errors */
		  else
		    {
		      totalDx += x2 - x1 + 1 ;
		      longDna = ac_zone_dna (cosmid, b1, b2, h) ;
		      switch (type)
			{
			case 1:
			  dn = longDna ? smith_waterman (longDna, 1, strlen ((char*)longDna), shortDna, x1, x2) : 99999 ;
			  break ;
			case 2:
			  dn = longDna ? smith_waterman_banded (longDna, 1, strlen (longDna), shortDna, x1, x2) : 99999 ;
			  break ;
			case 3:
			  dn = longDna ? smith_waterman (longDna, 1, strlen ((char*)longDna), shortDna, x1, x2) : 99999 ;
			  dnb = longDna ? smith_waterman_banded (longDna, 1, strlen (longDna), shortDna, x1, x2) : 99999 ;
			  if (dn != dnb)
			    messcrash ("smith_waterman_banded (%d) != smith_waterman (%d) in %s", dnb, dn, ac_name(tg)) ;
			  break ;
			default:
			  dn = 99999 ;
			  break ;
			}
		    }
		  freeOutf ("Assembled_from %d %d \"%s\" %d %d %d\n"
			  , c1, c2, cq, y1, y2, dn) ;
		  ac_free (longDna) ;
		}
	      freeOut ("-D Rank\nSmith_Waterman_done\n\n") ;
	    }
	}  
    }
  ac_free (h) ;
  return 1 ;
} /* gExportEstTgCost */

/*********************************************************************/
/*
 * Returns the number of alignments evaluated
 */

static int gExportTgCost (AC_OBJ tg, int type)
{
  int ii, nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE ests = ac_tag_table (tg, "Read", h) ;
  AC_OBJ est ;

  /* loop on all ESTs/mRNAs participating in this gene */
  for (ii = 0 ; ests && ii < ests->rows ; ii++)
    {
      est = ac_table_obj (ests, ii, 0, h) ;
      nn += gExportEstTgCost (est, tg, type) ;
    }
  ac_free (h) ;
  return nn ;
} /* gExportTgEstCost */


/*********************************************************************/
/*
 * Returns the number of alignments evaluated
 */

static int gExportCosts (AC_DB db, char *template)
{
  AC_ITER iter = 0 ; 
  AC_OBJ tg ;
  int nn = 0, orf ;
  char *qq ;
  int type = 2 ;/* 3: test and compare both systems, 
		   1: full
		   2: banded (3 times faster)
		*/
  /* testing phase is to verify that 'full' and 'banded' agree */
  if (type != 3) /* production */
    qq = messprintf ("find tg IS \"%s\" AND !Smith_Waterman_done", template) ;
  else  /* test, redo even if known */
    qq = messprintf ("find tg IS \"%s\" ", template) ;

  iter = ac_query_iter (db, TRUE, qq, 0, 0) ; /* all alignments */
  while ((tg = ac_next_obj (iter)))
    {
      if (1) gIntronsTg (tg) ;
      if (1)
	{
	  nn += gExportTgCost (tg, type) ;
	  orf = gGetOrfLen (tg) ;
	  freeOutf ("Transcribed_gene %s\nORF %d\n\n", ac_name(tg), orf) ;
	}
      ac_free (tg) ;
    }
  ac_free (iter) ;

  freeOut ("\n\n ") ;

  return nn ;
} /* gExportCosts */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

/*
 * Export on acedb stdout (freeOutf) the nature of all introns feet (e.g. gt_ag)
 */

static int gIntronsTg (AC_OBJ tg)
{
  int ii, max, a1, a2, b1, b2, g1, g2 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE im, splicing = ac_tag_table (tg, "Splicing", h) ;
  AC_OBJ cosmid ;
  const char *tag ;
  char *longDna, *cp, *cq,  buf[6]  ;
  int gt_ag=0, gc_ag=0, ct_ac=0, at_ac=0, nExon = 0 ;
  vTXT vOther = vtxtHandleCreate (h) ;
  
  cosmid = ac_tag_obj (tg, "Genomic_sequence", h) ;
  im = ac_tag_table(tg, "Covers", h) ;

  if (cosmid && splicing && im && im->cols >= 3 &&
      (g1 = ac_table_int (im, 0, 2, 0)) && 
      (g2 = ac_table_int (im, 0, 3, 0)) &&     
      (longDna = ac_zone_dna (cosmid, g1, g2, h)))
    {
      freeOutf ("\nTranscribed_gene \"%s\"\n", ac_name (tg)) ;
      freeOutf ("-D Intron_boundaries\n") ;
      freeOutf ("-D Nb_possible_exons\n") ;
      freeOutf ("-D Splicing\n") ;
      
      max = longDna ? strlen ((char*)longDna) : 0 ;
      for (ii = 0 ; longDna && ii < splicing->rows ; ii++)
	{
	  tag = ac_table_printable (splicing, ii, 2, "") ;
	  if (strstr (tag, "xon"))
	    nExon++ ;

	  a1 = b1 = ac_table_int (splicing, ii, 0, 0) ;
	  a2 = b2 = ac_table_int (splicing, ii, 1, 0) ;
	  /* fix a bug in import from kent */
	  if (ii == 1 && nExon == 2) /* 2 exons de suite */
	    {
	      a1 = ac_table_int (splicing, 0, 1, 0) + 1 ;
	      a2 = ac_table_int (splicing, 1, 0, 0) - 1 ;
	      tag = "Intron" ;
	    }
	  if (strstr (tag, "tron"))
	    {
	      if (a1 > 0 && a2 > 0 && a1 < max && a2 < max)
		{
		  cp = longDna + a1 - 1 ;
		  cq = longDna + a2 - 2 ;
		  memcpy (buf, cp, 2) ;	  buf[2] = '_' ; memcpy (buf+3, cq, 2) ; buf[5]=0;
		  if (a2 > a1 + 1)
		    {
		      freeOutf ("Splicing %d %d %s %s Length %d bp\n", a1, a2, tag, buf, a2 - a1 + 1) ;
		      if (!strcmp (buf, "gt_ag"))
			gt_ag++ ;
		      else if (!strcmp (buf, "gc_ag"))
			gc_ag++ ;
		      else if (!strcmp (buf, "ct_ac"))
			ct_ac++ ;
		      else if (!strcmp (buf, "at_ac"))
			at_ac++ ;
		      else
			{
			  int nOther = 0 ;
			  char *cOther ;
			  vtxtPrintf (vOther, "%s ", buf) ;
			  
			  cOther = vtxtPtr (vOther) ;
			  while ((cOther = strstr (cOther, buf)))
			    { nOther++ ; cOther += 6 ; }
			  freeOutf ("Other %s %d\n", buf, nOther) ;
			}
		    }
		  else
		    freeOutf ("Splicing %d %d %s %s Length %d bp\n", a1, a2, tag, "v-rep", a2 - a1 + 1) ;
		}
	    }
	  if (strstr (tag, "xon") || a1 != b1)
	    freeOutf ("Splicing %d %d Exon %d Length %d bp\n", b1, b2, nExon, b2 - b1 + 1) ;
	}
      if (gt_ag)  freeOutf ("gt_ag %d\n", gt_ag) ;
      if (gc_ag)  freeOutf ("gc_ag %d\n", gc_ag) ;
      if (ct_ac)  freeOutf ("ct_ac %d\n", ct_ac) ;
      if (at_ac)  freeOutf ("at_ac %d\n", at_ac) ;
      if (nExon)  freeOutf ("Nb_possible_exons %d\n", nExon) ;
      freeOutf ("\n") ;
    }
  ac_free (h) ;
  return 1 ;
}  /* gIntronsTg */

/*********************************************************************/
/* to be called manually from the  debugger */
static void showStats (DICT *dict, Array aa)
{
  STAT *sp ;
  int ii ;
 
  if (!arrayExists (aa)) return ;
  freeOutf ("# alignments %d\n", aa ? arrayMax (aa) : 0) ;
  freeOutf ("      %8s%8s%8s%8s%8s%8s%8s%6s%6s%6s\n"
	  ,"est","gene","-","-","-","-","-", "score","ali","err") ;
  for (ii = 0 ; aa && ii < arrayMax (aa) && ii < 24 ; ii++)
    {
      sp = arrayp (aa, ii, STAT) ;
      freeOutf ("%4d: %8s%8s%8s%8s%8s%8s%8s%6d%6d%6d %8d %8d\n", 
	      ii, dictName(dict, sp->est), dictName(dict, sp->gene),
	      sp->alibaba ? "Alibaba" : "-",
	      sp->g_split ? "g_split" : "-",
	      sp->best == 10 ? "Best" : "-",
	      sp->best >= 3 ? "Best1_5" : "-",
	      sp->best >= 2 ? "Best6_50" : "-",
	      sp->score, sp->ali, sp->err,
	      sp->a1, sp->a2) ;
    }
  freeOutf ("\n\n") ;
  showStats (0, 0) ; /* for compiler happiness */
} /* showStats */

/*************************************************************************************/
/********************* Chromosomes quality report ************************************/
/*************************************************************************************/

typedef struct chromStatStruct 
{
  /* ATTENTION, only ints are allowed because of gCumulate */
  int nam ; /* index in dictionary */
  int nAlibaba, nOk, nOkKb, nCumul, nExcellent, nGood, nPartial, nDubious, nBad, nMasked, nUnaligned, nReverse, nUnclipped, nGenomic, nNoIntron ;
  int cLen, cAli, cErr ;
  int n99, ali99, err99, asterix ; /* just those gold aligning at 99% */
  int nRepetition, nRepeatSameChrom, nRepeatOtherChrom, nRepeatInRandom, ngDel, nQuasiExcellent ;
  int nRepeatSameChromIntronBoth, nRepeatSameChromIntronOne, nRepeatSameChromNoIntron ;
  int nRepeatTandem, nRepeatPalindromeInward, nRepeatPalindromeOutward, nRepeatDistant ;
  int nRepeatOtherChromIntronBoth, nRepeatOtherChromIntronOne, nRepeatOtherChromNoIntron ;
  int nTransposition , nMosaic, nInversion, nMosaicOrInversion, nGenomeProblem, nV_repeat ;
  int nSuspectedInternalDeletion, nAtypicalIntron, nDuplicated ;
}  CHROM_STAT ;

/*********************************************************************/

static void gCumulate (DICT *dict, Array rr) 
{
  int i, j, ii = arrayMax (rr) ;
  CHROM_STAT *r, *s ;
  int *ip, *jp ;

  dictAdd (dict, "Total", &ii) ;
  s = arrayp (rr, ii , CHROM_STAT) ;
  
  for (i = 0 ; i < ii  ; i++)
    {
      r = arrayp (rr, i , CHROM_STAT) ;
      for (j = 0 ; j < sizeof (CHROM_STAT)/sizeof(int) ; j++)
	{
	  ip = (int*) s + j ;
	  jp = (int*) r + j ;
	  *ip += *jp ;
	}
    }
  s->nam = ii ;
  s->asterix = ' ' ;
}

/*********************************************************************/
/* DD = 7*8*9*10*11*12 is divisible by 2*(1,2,3,4,5,6,7,8,9,10,11,12) */
#define DD 665280
static void gGetChromStats (AC_DB db, DICT *dict, Array rr, char *template)
{
  AC_HANDLE h = ac_new_handle () ;
  int jj, ir, gp ;
  AC_KEYSET ks ;
  AC_ITER genes ;
  AC_OBJ est, gene = 0 ;
  AC_TABLE fg ;
  CHROM_STAT *r ;
  DICT *suspect = dictHandleCreate (1000, h) ;
  char *cp, cc, buff[128] ;

  genes = ac_query_iter (db, 1, "find cdna_clone Suspected_internal_deletion ; > read ; >from_gene", 0, h) ;
  while (ac_free(gene), (gene = ac_next_obj (genes)))
    dictAdd (suspect, ac_name (gene), 0) ;
  ac_free (genes) ;
  
  if (template && strcmp (template, "*"))
    genes = ac_query_iter (db, 1, messprintf ("find est %s ;>alibaba", template), 0, h) ;
  else
    genes = ac_query_iter (db, 1, "find tg alibaba", 0, h) ;
  
  while (ac_free(gene), (gene = ac_next_obj (genes)))
    {
      strncpy (buff, ac_tag_printable (gene, "IntMap", "toto"), 127) ;
      cp = buff ;
      while (*++cp && *cp != '|' && *cp != '_') ;
      *cp = 0 ;
      if (!dictFind (dict, buff, &jj))
	continue ;	  
      r = arrayp (rr, jj, CHROM_STAT) ;
      r->nam = jj ;
      if (ac_has_tag (gene, "Masked"))
	{
	  r->asterix = '*' ;
	  r->nMasked++ ;
	  dictAdd (dict, "Masked", &jj) ;
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  r->nam = jj ;
	  r->asterix = '*' ;
	}
      if (!r->asterix)
	r->asterix = ' ' ;
      r->nAlibaba++ ;
      gp = 0 ;
      if (ac_has_tag (gene, "excellent"))
	r->nExcellent++ ;
      if (ac_has_tag (gene, "Good"))
	r->nGood++ ;
      if (ac_has_tag (gene, "Partial"))
	{ gp = 2 ; r->nPartial++ ; }
      if (ac_has_tag (gene, "Dubious"))
	{ gp = 1 ; r->nDubious++ ; }
      if (ac_has_tag (gene, "Bad"))
	{ gp = 1 ; r->nBad++ ; }
      if (gp == 1)
	continue ;
      r->nOk++ ;

      if (ac_has_tag (gene, "Genome_deletion"))
	{ gp = 2 ; r->ngDel++ ; }
      if (ac_has_tag (gene, "Transposition"))
	{ gp = 2 ; r->nTransposition++ ; }

      gp = 0 ;
      if (ac_has_tag (gene, "Mosaic"))
	{ gp = 2 ; r->nMosaic++ ; }
      if (ac_has_tag (gene, "Inversion"))
	{ gp = 2 ; r->nInversion++ ; }

      if (ac_has_tag (gene, "Variable_repeat"))
	{ gp = 2 ; r->nV_repeat++ ; }

      if (ac_has_tag (gene, "Other")) /* &&  !dictFind (suspect, ac_name(gene), 0))  rm sept 16 2004*/
	/* if you add back keep in synch with the construction of the keyset */
	{ r->nAtypicalIntron++ ; }

      if (0)
	{ /* disabled mars 12, shows how many clones align 99% length whatever the error rate */
	  est = ac_tag_obj (gene, "Read", 0) ;
	  fg = ac_tag_table (est, "From_gene", h) ;
	  for (ir = 0 ; fg && ir < fg->rows ; ir++)
	    {
	      if (strcmp (ac_name(gene), ac_table_printable (fg, ir, 0,  "")))
		continue ;
	      r->cLen += ac_table_int (fg, ir, 1, 0) ;
	      r->cAli += ac_table_int (fg, ir, 3, 0) ;
	      r->cErr += ac_table_int (fg, ir, 4, 0) ;
	      if (100 * ac_table_int (fg, ir, 3, 0) >= 99 * ac_table_int (fg, ir, 1, 0))
		{
		  r->n99++ ;
		  r->ali99 += ac_table_int (fg, ir, 3, 0) ;
		  r->err99 += ac_table_int (fg, ir, 4, 0) ;
		}
	    }
	  ac_free (fg) ;
	  ac_free (est) ;
	}
      ac_free (gene) ;
    }


  /* repetitions are counted differently, using Gold, rather than alibaba */
  /* actually we count all pairs twice as seen from each end */
  /* therefore if a clone is aligned N times we count N(N-1) pairs
     by dividing by N-1 we get an exact result in cases
     N=2 or N=3
     or any case where there is a single category, 
     which happens to be true, as of may 2004 HINV3, for N >= 5
  */
  if (template ) 
    genes = ac_query_iter (db, 1, messprintf ("find est %s ; >from_gene; (Repetition || repetition_5)", template), 0, h) ;
  else
    genes = ac_query_iter (db, 1, "find tg (Repetition || repetition_5)", 0, h) ;

  while (ac_free(gene), (gene = ac_next_obj (genes)))
    {
      AC_ITER repeatIter ;
      AC_TABLE MapGene, MapRep ;
      AC_OBJ rep = 0 ;
      const char *repMap ;

      int a1, a2, b1, b2 ;
      int nRepeats ;
      int nRepeatSameChrom, nRepeatOtherChrom, nRepeatInRandom ;
      int nRepeatSameChromIntronBoth ;
      int nRepeatSameChromIntronOne ;
      int nRepeatSameChromNoIntron ;
      
      int nRepeatOtherChromIntronBoth ;
      int nRepeatOtherChromIntronOne ;
      int nRepeatOtherChromNoIntron ;

      int nRepeatTandem, nRepeatPalindromeInward, nRepeatPalindromeOutward, nRepeatDistant ;

      int dd ;
      BOOL hasIntron = ac_has_tag (gene, "gt_ag") || ac_has_tag (gene, "gc_ag") || ac_has_tag (gene, "at_ac") ;

      strncpy (buff, ac_tag_printable (gene, "IntMap", "toto"), 127) ;
      cp = buff ;
      while (*++cp && *cp != '|' && *cp != '_') ;
      cc = *cp ; *cp = 0 ;
      if (!dictFind (dict, buff, &jj))
	continue ;
      ks = ac_objquery_keyset (gene, ">read; >alibaba ; Excellent", 0) ;
      gp = ac_keyset_count (ks) ;
      ac_free (ks) ;
      if (!gp) 
	continue ;

      r = arrayp (rr, jj, CHROM_STAT) ;
      r->nam = jj ;
      if (!ac_has_tag (gene, "Excellent"))
	r->nQuasiExcellent++ ; /* denominator = excellent + quasi==repeat_5 of excellent */
      nRepeats  = 0 ;
      nRepeatSameChrom = nRepeatOtherChrom = nRepeatInRandom = 0 ;
      nRepeatSameChromIntronBoth = 0 ;
      nRepeatSameChromIntronOne = 0 ;
      nRepeatSameChromNoIntron = 0 ;
      
      nRepeatOtherChromIntronBoth = 0 ;
      nRepeatOtherChromIntronOne = 0 ;
      nRepeatOtherChromNoIntron = 0 ;

      nRepeatTandem = nRepeatPalindromeInward = nRepeatPalindromeOutward = nRepeatDistant = 0 ;

      MapGene = 0 ;
      repeatIter = ac_objquery_iter (gene, ">read;>from_gene ;  (Repetition || repetition_5)", h) ;
      while (ac_free(rep), (rep = ac_next_obj (repeatIter)))
	{
	  BOOL repHasIntron = ac_has_tag (rep, "gt_ag") || ac_has_tag (rep, "gc_ag")|| ac_has_tag (rep, "at_ac") ;
	  nRepeats++ ;

	  if (!strcmp (ac_name(rep), ac_name(gene))) continue ;  
	  /* notice that a repeat 5 has no tag excellent/good since it is NOT gold */
	   /* one or the other member of the pair is on a floating contig */
	  repMap = ac_tag_printable (rep, "IntMap", "") ;
	  if (cc || strlen (repMap) > 2)
	    nRepeatInRandom++ ;
	  else /* both are real */
	    { 
	      /* if the repeat on the same chromo */
	      if (!strcmp (buff, repMap)) /* if the repeat on the same chromo */
		{
		  nRepeatSameChrom++ ;
		  if (hasIntron && repHasIntron) nRepeatSameChromIntronBoth++ ;
		  else if (hasIntron != repHasIntron) nRepeatSameChromIntronOne++ ;
		  else nRepeatSameChromNoIntron++ ;
		  
		  MapRep = ac_tag_table (rep, "IntMap", 0) ;
		  if (!MapGene) MapGene = ac_tag_table (gene, "IntMap", 0) ;

		  a1 = ac_table_int (MapGene, 0, 1, -1) ;
		  a2 = ac_table_int (MapGene, 0, 2, -1) ;
		  b1 = ac_table_int (MapRep, 0, 1, -1) ;
		  b2 = ac_table_int (MapRep, 0, 2, -1) ;
		  if (a1 > 0 && a2 > 0 && b1 > 0 && b2 > 0)
		    {
		      if (a1 < b1 + 1000000 && a1 > b1 - 1000000)
			{
			  if ((a1 < a2 && b1 < b2) || (a1 > a2 && b1 > b2))
			    nRepeatTandem++ ;
			  else if (a1 < a2 && b1 > b2)
			    {
			      if (a1 < b1) nRepeatPalindromeInward++ ;
			      else nRepeatPalindromeOutward++ ;
			    }
			  else if (a1 > a2 && b1 < b2)
			    {
			      if (a1 > b1) nRepeatPalindromeInward++ ;
			      else nRepeatPalindromeOutward++ ;
			    }
			}
		      else
			nRepeatDistant++ ;
		    }
		  ac_free (MapRep) ;
		}
	      else
		{ 
		  nRepeatOtherChrom++ ;
		  if (hasIntron && repHasIntron) nRepeatOtherChromIntronBoth++ ;
		  else if (hasIntron != repHasIntron)  nRepeatOtherChromIntronOne++ ;
		  else nRepeatOtherChromNoIntron++ ;
		}
	    }
	}
      ac_free (repeatIter) ;
      ac_free (MapGene) ;
      if (nRepeats < 2) nRepeats = 2 ; /* protection against zerodivide */ 
      dd = DD/(nRepeats - 1) ; /* to always get integral numbers */
      r->nRepetition += DD ;

      r->nRepeatInRandom += dd * nRepeatInRandom ;
      r->nRepeatSameChrom += dd * nRepeatSameChrom ;
      r->nRepeatOtherChrom += dd * nRepeatOtherChrom ;

      r->nRepeatSameChromIntronBoth += dd * nRepeatSameChromIntronBoth ;
      r->nRepeatSameChromIntronOne += dd * nRepeatSameChromIntronOne ;
      r->nRepeatSameChromNoIntron += dd * nRepeatSameChromNoIntron ;

      r->nRepeatOtherChromIntronBoth += dd * nRepeatOtherChromIntronBoth ;
      r->nRepeatOtherChromIntronOne += dd * nRepeatOtherChromIntronOne ;
      r->nRepeatOtherChromNoIntron += dd * nRepeatOtherChromNoIntron ;

      r->nRepeatTandem += dd * nRepeatTandem  ;
      r->nRepeatPalindromeInward += dd * nRepeatPalindromeInward ;
      r->nRepeatPalindromeOutward += dd * nRepeatPalindromeOutward ;
      r->nRepeatDistant += dd * nRepeatDistant ;

      ac_free (gene) ;
    }

  ac_free (h) ;
  return ;
}  /* gGetChromStats */

/*********************************************************************/

static void gChromStats (AC_DB db, BOOL xml, char *template)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER maps = ac_query_iter (db, 0, "find map  (IS ? OR IS ?? ) && ! IS Un", 0, h) ;
  AC_OBJ map ;
  CHROM_STAT *r ;
  Array rr = arrayHandleCreate (32, CHROM_STAT, h) ;
  int jj = 0 ;
  DICT *dict = dictHandleCreate (500, h) ;

  while (( map = ac_next_obj (maps)))
    {
      dictAdd (dict, ac_name(map), &jj) ;
      r = arrayp (rr, jj++, CHROM_STAT) ;
      r->nam = jj ;
      ac_free (map) ;
    }
  dictAdd (dict, "Un", &jj) ;
  r = arrayp (rr, jj++, CHROM_STAT) ;

  gGetChromStats (db, dict, rr, template) ;
  gCumulate (dict, rr) ;

  /*
   * Table for the web, to hook to the cdna keysets 
   */
  if (1)
    {
      gOutSection (xml, "Number of gold alignments of different qualities, sorted by chromosome") ;
      gOutBreak (xml) ;
      gOutTableStart (xml, "ChromQuality") ;
      gOutTitles (xml, "Chrom\tcDNA\tExcellent\tGood\tPartial\tDubious\tBad\tMasked") ;
      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nAlibaba) continue ;
	  gOutTitle (xml, messprintf ("%s%c", dictName (dict, r->nam), r->asterix)) ;
	  gOutInt (xml, r->nAlibaba, FALSE, 1) ;
	  gOutInt (xml, r->nExcellent, FALSE, 1) ;
	  gOutInt (xml, r->nGood, FALSE, 1) ;
	  if (jj != arrayMax (rr) - 2 && jj != arrayMax (rr) - 3) /* exclude Un , mask */
	    {
	      gOutIntLink3 (xml, r->nPartial, FALSE, 1, "Chroms", dictName (dict, r->nam), "Partial") ;
	      gOutIntLink3 (xml, r->nDubious, FALSE, 1, "Chroms", dictName (dict, r->nam), "Dubious") ;
	      gOutIntLink3 (xml, r->nBad, FALSE, 1, "Chroms", dictName (dict, r->nam), "Bad") ;
	      gOutIntLink3 (xml, r->nMasked, FALSE, 1, "Chroms", dictName (dict, r->nam), "Masked") ;
	    }
	  else
	    {
	      gOutInt (xml, r->nPartial, FALSE, 1) ;
	      gOutInt (xml, r->nDubious, FALSE, 1) ;
	      gOutInt (xml, r->nBad, FALSE, 1) ;
	      gOutInt (xml, r->nMasked, FALSE, 1) ;
	    }
	    
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }

  /*
   * Table for the paper, drawing of plot: quality of the genome
   */

  if (!xml)
    {
      gOutSection (xml, "Percentage of gold alignments of different qualities, sorted by chromosome") ;
      freeOutf ("Masking the hyper variable regions") ;
      gOutBreak (xml) ;
      gOutTableStart (xml, "ChromQuality2") ;
      gOutTitles (xml, "Chrom\tGold\tExcellent\tGood\tPartial\tDubious\tBad\tMasked") ;
     
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  int total ;

	  r = arrayp (rr, jj, CHROM_STAT) ;
	  total = r->nAlibaba ;
	  if (total < 1)
	    continue ;
	  gOutTitle (xml, messprintf ("%s%c", dictName (dict, r->nam), r->asterix)) ;
	  gOutInt (xml, total, FALSE, 1) ;
	  gOutInt (xml, r->nExcellent, TRUE, total) ;
	  gOutInt (xml, r->nGood, TRUE, total) ;
	  gOutInt (xml, r->nPartial, TRUE, total) ;
	  gOutInt (xml, r->nDubious, TRUE, total) ;
	  gOutInt (xml, r->nBad, TRUE, total) ;
	  gOutInt (xml, r->nMasked, TRUE, total) ;
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }


  /*
   * Table for the web, to hook to the cdna keysets 
   */
  if (1)
    {
      gOutSection (xml, "Split alignments, sorted by chromosome") ; 
      freeOutf (" only excellent, good and partial alignments are counted") ;
      gOutBreak (xml) ;
      freeOutf (" hyper variable regions are excluded") ;
      gOutBreak (xml) ;
      gOutTableStart (xml, "ChromDefects") ;
      gOutTitles (xml, "Chromosome\tGold\tPartial alignment\tGenome deletion\tTransposition\tVariable repeat\tInversion\tMosaic\tAtypical intron") ;
      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nOk) continue ; 
	  
	  gOutTitle (xml, messprintf ("%s%c", dictName (dict, r->nam), r->asterix)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;
	  if (jj != arrayMax(rr) - 2 && jj != arrayMax(rr) - 3)
	    {
	      gOutIntLink3 (xml, r->nPartial, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Partial") ;
	      gOutIntLink3 (xml, r->ngDel, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Genome_deletion") ;
	      gOutIntLink3 (xml, r->nTransposition, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Transposition") ;
	      gOutIntLink3 (xml, r->nV_repeat, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Variable_repeat") ;
	      gOutIntLink3 (xml, r->nInversion, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Inversion") ;
	      gOutIntLink3 (xml, r->nMosaic, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Mosaic") ;
	      gOutIntLink3 (xml, r->nAtypicalIntron, FALSE, r->nOk, "Chroms", dictName (dict, r->nam), "Atypical_intron") ;
	    }
	  else
	    {
	      gOutInt (xml, r->nPartial, FALSE, r->nOk) ;
	      gOutInt (xml, r->ngDel, FALSE, r->nOk) ;
	      gOutInt (xml, r->nTransposition, FALSE, r->nOk) ;
	      gOutInt (xml, r->nV_repeat, FALSE, r->nOk) ;
	      gOutInt (xml, r->nInversion, FALSE, r->nOk) ;
	      gOutInt (xml, r->nMosaic, FALSE, r->nOk) ;
	      gOutInt (xml, r->nAtypicalIntron, FALSE, r->nOk) ;
	    }
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }
  
  if (xml) /* export the corresponding tables */
    {
      char *defect[] = {"Partial","Dubious","Bad","Masked","Genome_deletion","Transposition","Variable_repeat","Inversion","Mosaic", "Atypical_intron", 0} ;
      char *tag[] = {
	"Partial"
	, "Dubious"
	, "Bad"
	, "Masked"   /* MUST stay as iDefect == 3 */
	,"Genome_deletion && (excellent || good || partial)"
	,"Transposition && (excellent || good || partial)"
	,"Variable_repeat && (excellent || good || partial)"
	,"Inversion && (excellent || good || partial)"
	,"Mosaic && (excellent || good || partial)"
	, "Other && (excellent || good || partial)"
	, 0} ;
      char *chroms[] = {"Total", "1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18","19","20",
			"21","22","X","Y", 0} ;
      int idefect, ichrom ;
      AC_KEYSET ks1, ks2 ;

      for (idefect = 0 ; defect[idefect] ; idefect++)
	{
	  ks1 = ac_dbquery_keyset (db
				   , messprintf ( "find tg alibaba && (%s)", tag[idefect])
				   , h) ; 
	  for (ichrom = 0; chroms[ichrom];ichrom++)
	    {
	      if (ichrom && idefect != 3)
		ks2 = ac_ksquery_keyset (ks1, messprintf ("!masked && (IntMap = \"%s\" || IntMap = \"%s?NT*\") ; > read", chroms[ichrom], chroms[ichrom]), h) ;
	      else if (ichrom && idefect == 3)
		ks2 = ac_ksquery_keyset (ks1, messprintf (" (IntMap = \"%s\" || IntMap = \"%s?NT*\")  ; > read", chroms[ichrom], chroms[ichrom]), h) ;
	      else
		ks2 = ac_ksquery_keyset (ks1, "> read", h) ;
	      gExportTable4 (db, ks2, "Chroms", chroms[ichrom], defect[idefect]) ;
	      ac_free (ks2) ;
	    }
	  ac_free (ks1) ;
	}
    }

  /*
   * Table for the paper: problems per chromosome 
   */
  if (!xml)
    {
      gOutSection (xml, "Percentage of Split alignments, sorted by chromosome") ; 
      freeOutf (" only excellent, good and partial alignments are counted") ;
      gOutBreak (xml) ;
      freeOutf (" hyper variable regions are excluded") ;
      gOutBreak (xml) ;
      gOutTableStart (xml,"ChromDefects2") ;
      gOutTitles (xml, "Chromosome\tGold\tPartial alignment\tGenome deletion\tTransposition\tVariable repeat\tInversion\tMosaic\tAtypical intron") ;

      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nOk) continue ; 
	  
	  gOutTitle (xml, messprintf ("%s%c", dictName (dict, r->nam), r->asterix)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;
	  gOutInt (xml, r->nPartial, TRUE, r->nOk) ;
	  gOutInt (xml, r->ngDel, TRUE, r->nOk) ;
	  gOutInt (xml, r->nTransposition, TRUE, r->nOk) ;
	  gOutInt (xml, r->nV_repeat, TRUE, r->nOk) ;
	  gOutInt (xml, r->nInversion, TRUE, r->nOk) ;
	  gOutInt (xml, r->nMosaic, TRUE, r->nOk) ;
	  gOutInt (xml, r->nAtypicalIntron, TRUE, r->nOk) ;
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }

  /*
   * Table for the web, to hook to the cdna keysets 
   */
  if (1)
    {
      gOutSection (xml, "Exact and approximate repetitions, sorted by chromosome") ;
      freeOutf (" only excellent alignments are counted") ;
      gOutBreak (xml) ;
      gOutBreak (xml) ;
      gOutTableStart (xml,"ChromRepeats") ;
      gOutTitles (xml, "Chromosome\tcDNA\tRepetition\tRepetition on the same chromosome\tOn another chromosome\tRepetition involving isolated contig\tChromosome\tSame both intron\tSame one intron\tSame no intron\tOther both intron\tOther one intron\tOther no intron\tTandem\tInward palindrome\tOutward palindrome\tFurther than 1 Mb") ;
      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nRepetition) continue ;

	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nExcellent, FALSE, 1) ;
	  gOutInt (xml, r->nRepetition/DD, FALSE, 1) ;
	  gOutInt (xml, r->nRepeatSameChrom/DD, FALSE, r->nRepetition) ;
	  gOutInt (xml, r->nRepeatOtherChrom/DD, FALSE, r->nRepetition) ;
	  gOutInt (xml, r->nRepeatInRandom/DD, FALSE, r->nRepetition) ;
	  gOutText (xml, dictName (dict, r->nam)) ;
	  gOutFloat (xml, r->nRepeatSameChromIntronBoth/(float)DD) ;
	  gOutFloat (xml, r->nRepeatSameChromIntronOne/(float)DD) ;
	  gOutFloat (xml, r->nRepeatSameChromNoIntron/(float)DD) ;
	  gOutFloat (xml, r->nRepeatOtherChromIntronBoth/(float)DD) ;
	  gOutFloat (xml, r->nRepeatOtherChromIntronOne/(float)DD) ;
	  gOutFloat (xml, r->nRepeatOtherChromNoIntron/(float)DD) ;

	  gOutFloat (xml, r->nRepeatTandem/(float)DD) ;
	  gOutFloat (xml, r->nRepeatPalindromeInward/(float)DD) ;
	  gOutFloat (xml, r->nRepeatPalindromeOutward/(float)DD) ;
	  gOutFloat (xml, r->nRepeatDistant/(float)DD) ;
	  
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }
  /*
   * Table for the paper: repeats per chromosome 
   */
  if (!xml)
    {
      gOutSection (xml, "Exact and approximate repetitions, sorted by chromosome") ;
      freeOutf (" only excellent alignments are counted") ;
      gOutBreak (xml) ;
      gOutBreak (xml) ;
      gOutTableStart (xml, "ChromRepeats2") ;
      gOutTitles (xml, "Chromosome\tcDNA\tRepetition\tRepetition on the same chromosome\tOn another chromosome\tRepetition involving isolated contig\tChromosome\tSame both intron\tSame one intron\tSame no intron\tOther both intron\tOther one intron\tOther no intron\tTandem\tInward palindrome\tOutward palindrome\tFurther than 1 Mb") ;
           
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nOk) continue ; 
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nExcellent, FALSE, 1) ;
	  gOutInt (xml, r->nRepetition/DD, TRUE, r->nExcellent + r->nQuasiExcellent) ;
	  gOutInt (xml, r->nRepeatSameChrom/DD, TRUE, r->nExcellent + r->nQuasiExcellent) ;
	  gOutInt (xml, r->nRepeatOtherChrom/DD, TRUE, r->nExcellent + r->nQuasiExcellent) ;
	  gOutInt (xml, r->nRepeatInRandom/DD, TRUE, r->nExcellent + r->nQuasiExcellent) ;
	  gOutText (xml, dictName (dict, r->nam)) ;
	  /* the next guys may be half integers, so beware*/
	  gOutInt (xml, 2*r->nRepeatSameChromIntronBoth/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatSameChromIntronOne/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatSameChromNoIntron/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatOtherChromIntronBoth/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatOtherChromIntronOne/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatOtherChromNoIntron/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;

	  gOutInt (xml, 2*r->nRepeatTandem/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatPalindromeInward/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatPalindromeOutward/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  gOutInt (xml, 2*r->nRepeatDistant/DD, TRUE, 2*(r->nExcellent + r->nQuasiExcellent)) ;
	  
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }

  ac_free (h) ;
  return ;
} /* gChromStats */

/*********************************************************************/
/*********************************************************************/

static void gGetLibraryStats (AC_DB db, DICT *dict, Array rr, BOOL mgcRef, int fljRef)
{
  AC_HANDLE h = ac_new_handle () ;
  int jj, gp ;
  const char *lib ;
  AC_ITER ests ;
  AC_TABLE gDna ;
  AC_OBJ est = 0, gene = 0, clone = 0 ;
  CHROM_STAT *r ;

  ests = ac_query_iter (db, 1, "find est", 0, h) ;
  
  while (ac_free(gene), ac_free(clone), ac_free(est), (est = ac_next_obj (ests)))
    {
      /* can we get the lib */
      clone = ac_tag_obj (est, "cdna_clone", 0) ;
      if (!clone) continue ;
      lib = ac_tag_printable (clone, "Library", 0) ;
      if (!lib) continue ;
      if (mgcRef)
	{
	  if (strcmp (lib, "MGC"))
	    continue ;
	  if (!ac_has_tag (clone, "MGC_ref"))
	    continue ;
	  lib = "MGC.filtered" ;
	}
      if (0 && strcasecmp (lib, "FLJ"))
	continue ;
      switch (fljRef)
	{
	case 1:
	  if (!ac_has_tag (clone, "FLJ.1_1"))
	    continue ;
	  lib = "FLJ.1_1" ;
	  break ;
	case 2:
	  if (!ac_has_tag (clone, "FLJ.1_2"))
	    continue ;
	  lib = "FLJ.1_2" ;
	  break ;
	case 3:
	  if (!ac_has_tag (clone, "FLJ.2"))
	    continue ;
	  lib = "FLJ.2" ;
	  break ;
	case 4:
	  if (!ac_has_tag (clone, "FLJ.3"))
	    continue ;
	  lib = "FLJ.3" ;
	  break ;
	case 5:
	  if (!ac_has_tag (clone, "FLJ_HRIandIMS"))
	    continue ;
	  lib = "FLJ_HRIandIMS" ;
	  break ;
	case 6:
	  if (!ac_has_tag (clone, "FLJ_KDRI"))
	    continue ;
	  lib = "FLJ_KDRI" ;
	  break ;
	case 7:
	  if (!ac_has_tag (clone, "FLJ_HRIandIMS") && !ac_has_tag (clone, "FLJ_KDRI"))
	    continue ;
	  lib = "FLJ.filtered" ;
	  break ;

	default:
	  break ;
	}

      dictAdd (dict, lib, &jj) ;
     
      r = arrayp (rr, jj, CHROM_STAT) ;
      r->nam = jj ;
      r->nAlibaba++ ;

      /* is the clone aligned */
      gene = ac_tag_obj (est, "alibaba", 0) ;
      if (!gene)
	{
	  r->nUnaligned++ ;
	  continue ;
	}

      /*
	union of
	unali , genome del, transpo, inver, mosaic, partial, genome_conta, 
	 = any
      */

      /* special hack to please weissenbach and remove hyper variable zones*/
      gp = 0 ;
      if (ac_has_tag (gene, "Masked"))
	{ gp = 1 ; r->nMasked++ ; }
      else if (ac_has_tag (gene, "excellent"))
	r->nExcellent++ ;
      else if (ac_has_tag (gene, "Good"))
	r->nGood++ ;
      else if (ac_has_tag (gene, "Partial"))
	{  r->nPartial++ ; }
      else if (ac_has_tag (gene, "Dubious"))
	{ gp = 1 ; r->nDubious++ ; }
      else if (ac_has_tag (gene, "Bad"))
	{ gp = 1 ; r->nBad++ ; }
      if (gp)
	continue ;
      r->nOk++ ;  /* good clones */
      gDna = ac_tag_table (est, "DNA", 0) ;
      if (gDna)
	{ /* cumulated length of the good clones */
	  r->nOkKb += ac_table_int (gDna, 0, 1, 0) ;
	  ac_free (gDna) ;
	}

      /* count the problems */
      if (ac_has_tag (gene, "Variable_repeat"))
	{ r->nV_repeat++ ; }
      if (ac_tag_int (est, "Vector_clipping", 0) > 2)
	{ r->nUnclipped++ ; }
      if (ac_has_tag (est, "Duplicate_of_read"))
	{ r->nDuplicated++ ; }
      if (ac_has_tag (gene, "G_reverse") && (ac_has_tag (gene, "gt_ag") || ac_has_tag (gene, "gc_ag")) && !ac_has_tag (gene, "ct_ac") && ! ac_has_tag (gene, "rear"))
	{ r->nReverse++ ; }
      if (ac_has_tag (gene, "Genome_deletion"))
	{  r->ngDel++ ; }

      /* count the structural problems */
      gp = 0 ;
      if (ac_has_tag (gene, "Transposition"))
	{ gp = 2 ; r->nTransposition++ ; }

      if (ac_has_tag (gene, "Mosaic"))
	{ gp = 2 ; r->nMosaic++ ; }
      if (ac_has_tag (gene, "Inversion"))
	{ gp = 2 ; r->nInversion++ ; }
      /* note that we do not and should not count aas clone defestc the v_repeat and the v_pileUp pile_up */
      if (ac_has_tag (clone, "Suspected_genomic_contamination"))
	{ gp = 2 ; r->nGenomic++ ; }
      else if (! ac_has_tag (gene, "Intron_boundaries"))
	{  r->nNoIntron++ ; }
      if (ac_has_tag (clone, "Suspected_internal_deletion"))
	{ gp = 2 ; r->nSuspectedInternalDeletion++ ; }
      else if (ac_has_tag (gene, "Other"))
	{ r->nAtypicalIntron++ ; }
      
      if (gp)
	r->nCumul++ ;
    }

  ac_free (h) ;
  return ;
}  /* gGetLibraryStats */

/*********************************************************************/

static int orderByProject (const void *a,const  void *b)
{
  const CHROM_STAT *r = (const CHROM_STAT *) a ;
  const CHROM_STAT *s = (const CHROM_STAT *) b ;

  return  s->nAlibaba - r->nAlibaba ;
}

/*********************************************************************/

static void gLibraryStats (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  Array rr = arrayHandleCreate (32, CHROM_STAT, h) ;
  int jj = 0 ;
  CHROM_STAT *r ;
  DICT *dict = dictHandleCreate (500, h) ;

  dictAdd (dict, "KIAA", &jj) ;
  dictAdd (dict, "MGC.filtered", &jj) ;
  dictAdd (dict, "MGC", &jj) ;
  dictAdd (dict, "DKFZ", &jj) ;
  dictAdd (dict, "FLJ", &jj) ;
  dictAdd (dict, "FLJ.3", &jj) ;
  dictAdd (dict, "FLJ.1_2", &jj) ;
  dictAdd (dict, "FLJ.1_1", &jj) ;
  dictAdd (dict, "FLJ.2", &jj) ;
  dictAdd (dict, "FLJ_HRIandIMS", &jj) ;
  dictAdd (dict, "FLJ_KDRI", &jj) ;
  dictAdd (dict, "FLJ.filtered", &jj) ;
  dictAdd (dict, "CHGC", &jj) ;
  dictAdd (dict, "other", &jj) ;

  while (jj--)
    {
      r = arrayp (rr, jj, CHROM_STAT) ;
      r->nam = jj ;
    }

 
  gGetLibraryStats (db, dict, rr, FALSE, 0) ;
  if (0) arraySort (rr, orderByProject) ;
  gCumulate (dict, rr) ;
  gGetLibraryStats (db, dict, rr, TRUE, 0) ;
  if (0)
    {   /* flj details */
      gGetLibraryStats (db, dict, rr, FALSE,1) ;
      gGetLibraryStats (db, dict, rr, FALSE,2) ;
      gGetLibraryStats (db, dict, rr, FALSE,3) ;
      gGetLibraryStats (db, dict, rr, FALSE,4) ;
      gGetLibraryStats (db, dict, rr, FALSE,5) ;
      gGetLibraryStats (db, dict, rr, FALSE,6) ;
    }
  gGetLibraryStats (db, dict, rr, FALSE,7) ; /* flj filtered */

  /* nOkKb = cumulated dna length, we rationalize it in kb
     this means a clone of length 2kb count as 2, in the
     denominator, rather than one when counting defects
  */
  for (jj = 0 ; jj < arrayMax (rr) - 1 ; jj++)
    {
      r = arrayp (rr, jj, CHROM_STAT) ;
      r->nOkKb = r->nOkKb/ 1000 ;
      if (r->nOkKb < 1) r->nOkKb = 1 ;
    }
  /*
   * Table for the web, to hook to the cdna keysets 
   */
  if (1)
    {
      gOutSection (xml, "cDNA by quality of their best alignment, sorted by project\n\n") ;
      freeOutf ("considering absolutely all clones  of this library") ;
      gOutBreak (xml);

      gOutTableStart (xml,"ProjectQuality") ;
      gOutTitles (xml, "Project\tcDNA\tExcellent\tGood\tPartial\tDubious\tBad\tMasked\tUnaligned") ;
      
      for (jj = 0 ; jj < arrayMax (rr) - 2 ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nAlibaba) continue ;
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutIntLink3 (xml, r->nAlibaba, FALSE, 1,"Libraries", dictName (dict, r->nam), "All") ; /* whole library */
	  gOutIntLink3 (xml, r->nExcellent, FALSE, 1,"Libraries", dictName (dict, r->nam), "Excellent") ;
	  gOutIntLink3 (xml, r->nGood, FALSE, 1,"Libraries", dictName (dict, r->nam), "Good") ;
	  gOutIntLink3 (xml, r->nPartial, FALSE, 1,"Libraries", dictName (dict, r->nam), "Partial") ;
	  gOutIntLink3 (xml, r->nDubious, FALSE, 1,"Libraries", dictName (dict, r->nam), "Dubious") ;
	  gOutIntLink3 (xml, r->nBad, FALSE, 1,"Libraries", dictName (dict, r->nam), "Bad") ;
	  gOutIntLink3 (xml, r->nMasked, FALSE, 1,"Libraries", dictName (dict, r->nam), "Masked") ;
	  gOutIntLink3 (xml, r->nUnaligned, FALSE, 1,"Libraries", dictName (dict, r->nam), "Unaligned") ;
	  
	  gOutNewLine (xml) ;
	}
      for (; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (!r->nAlibaba) continue ;
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nAlibaba, FALSE, 1) ; /* all this library */
	  gOutInt (xml, r->nExcellent, FALSE, 1) ;
	  gOutInt (xml, r->nGood, FALSE, 1) ;
	  gOutInt (xml, r->nPartial, FALSE, 1) ;
	  gOutInt (xml, r->nDubious, FALSE, 1) ;
	  gOutInt (xml, r->nBad, FALSE, 1) ;
	  gOutInt (xml, r->nMasked, FALSE, 1) ;
	  gOutInt (xml, r->nUnaligned, FALSE, 1) ;
	  
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }
  
  if (xml)
    {
      char *defect[] = {"Excellent", "Good", "Partial", "Dubious", "Bad", "Masked", 0} ;
      char *tag[] = {"Excellent && !masked", "Good && !masked", "Partial && !masked", "Dubious && !masked", "Bad && !masked", "Masked", 0} ;
      char *libs[] = {"MGC && mgc_ref","MGC","FLJ", "FLJ && (FLJ_HRIandIMS || FLJ_KDRI)", "KIAA","DKFZ", 0} ; 
      char *libs2[] = {"MGC.filtered","MGC","FLJ", "FLJ.filtered", "KIAA","DKFZ", 0} ; 
      int idefect, ilib ;
      AC_KEYSET ks1, ks2 ;
      
      for (idefect = 0 ; defect[idefect] ; idefect++)
	{
	  ks1 = ac_dbquery_keyset (db
				   , messprintf ( "find tg alibaba && (%s) ; >read ;>cdna_clone", tag[idefect])
				   , h) ; 
	  for (ilib=0; libs[ilib];ilib++)
	    {
	      ks2 = ac_ksquery_keyset (ks1, messprintf ("library = %s ; >read", libs[ilib]), h) ;
	      gExportTable4 (db, ks2, "Libraries", libs2[ilib], defect[idefect]) ;
	      ac_free (ks2) ;
	    }
	  ac_free (ks1) ;
	}

      ks1 = ac_dbquery_keyset (db, "{ find est !from_gene} SETOR {FIND tg no_score; >read ; COUNT {>from_gene;!no_score} = 0 } ;> cdna_clone ", h) ;
      for (ilib=0; libs[ilib];ilib++)
	{
	  ks2 = ac_ksquery_keyset (ks1, messprintf ("library = %s ; >read", libs[ilib]), h) ;
	  gExportTable6 (db, ks2, "Libraries", libs2[ilib], "Unaligned") ;
	  ac_free (ks2) ;
	}
      ac_free (ks1) ;

     for (ilib=0; libs[ilib];ilib++)
	{
	  ks2 = ac_dbquery_keyset (db, messprintf ("find cdna_clone library = %s ", libs[ilib]), h) ;
	  gExportTable5 (db, ks2, "Libraries", libs2[ilib], "All") ;
	  ac_free (ks2) ;
	}
    }


  /*
   * Table for the paper, drawing of plot: quality of the projects
   */

  if (!xml)
    {
      gOutSection (xml, "cDNA by quality of their best alignment, sorted by project\n\n") ;
      freeOutf ("considering all clones  of this library, except masked") ;
      gOutBreak (xml);
      
      gOutTableStart (xml, "ProjectQuality2") ;
      gOutTitles (xml, "Project\tcDNA\tExcellent\tGood\tPartial\tDubious\tBad\tUnaligned") ;
      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  int cumul ;

	  r = arrayp (rr, jj, CHROM_STAT) ;
	  cumul = r->nAlibaba - r->nMasked ; 
	  if (cumul <= 0) 
	    continue ;

	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, cumul, FALSE, 1) ;  /* all members of this library, except masked */
	  gOutInt (xml, r->nExcellent, TRUE, cumul) ;
	  gOutInt (xml, r->nGood, TRUE, cumul) ;
	  gOutInt (xml, r->nPartial, TRUE, cumul) ;
	  gOutInt (xml, r->nDubious, TRUE, cumul) ;
	  gOutInt (xml, r->nBad, TRUE, cumul) ;
	  gOutInt (xml, r->nUnaligned, TRUE, cumul) ;
	  
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }

  /*
   * Table for the web: problems in the libraries
   */
  if (1)
    {
      gOutSection (xml, "Problematic clone alignments, sorted by project") ;
      freeOutf ("considering only excellent or good or partial alignments outside of the hyper variable regions") ;
      gOutBreak (xml);
      freeOutf ("The structural problems are:Partial/Transposition/Mosaic/Inversion/genomic_contamination/internal_deletion\n") ;
      gOutBreak (xml);
     
      gOutTableStart (xml, "ProjectDefects") ;
      gOutTitles (xml, "Project\tcDNA\tStructural problem\tUnaligned\tSuspected genome deletion\tTransposition\tInversion\tMosaic\tPartial\tSuspected genomic contamination\tNo intron\t") ;
      
      for (jj = 0 ; jj < arrayMax (rr) - 2; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;   /* not masked, not bad */
	  gOutInt (xml, r->nCumul, FALSE, 1) ; /* union of the following problems */
	  gOutIntLink3 (xml, r->nUnaligned, FALSE, 1,"Libraries", dictName (dict, r->nam), "Unaligned") ;
	  gOutIntLink3 (xml, r->ngDel, FALSE, 1,"Libraries", dictName (dict, r->nam), "Genome_deletion") ;
	  gOutIntLink3 (xml, r->nTransposition, FALSE, 1,"Libraries", dictName (dict, r->nam), "Transposition") ;
	  gOutIntLink3 (xml, r->nInversion, FALSE, 1,"Libraries", dictName (dict, r->nam), "Inversion") ;
	  gOutIntLink3 (xml, r->nMosaic, FALSE, 1,"Libraries", dictName (dict, r->nam), "Mosaic") ;
	  gOutIntLink3 (xml, r->nPartial, FALSE, 1,"Libraries", dictName (dict, r->nam), "Partial") ;
	  gOutIntLink3 (xml, r->nGenomic, FALSE, 1,"Libraries", dictName (dict, r->nam), "Genomic_contamination") ;
	  gOutInt (xml, r->nNoIntron, FALSE, 1) ;
	  
	  gOutNewLine (xml) ;
	}
      for (; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;   /* not masked, not bad */
	  gOutInt (xml, r->nCumul, FALSE, 1) ; /* union of the following problems */
	  gOutInt (xml, r->nUnaligned, FALSE, 1) ;
	  gOutInt (xml, r->ngDel, FALSE, 1) ;
	  gOutInt (xml, r->nTransposition, FALSE, 1) ;
	  gOutInt (xml, r->nInversion, FALSE, 1) ;
	  gOutInt (xml, r->nMosaic, FALSE, 1) ;
	  gOutInt (xml, r->nPartial, FALSE, 1) ;
	  gOutInt (xml, r->nGenomic, FALSE, 1) ;
	  gOutInt (xml, r->nNoIntron, FALSE, 1) ;
	  
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;

      if (xml)
	{
	  char *defect[] = {"Genome_deletion", "Transposition", "Inversion", "Mosaic",  "Variable_repeat", "Reverse", 0} ;
	  char *tag[] = {"Genome_deletion", "Transposition", "Inversion", "Mosaic",  "Variable_repeat", "G_reverse && !rear && (gt_ag || gc_ag) && !ct_ac", 0} ;
	  char *libs[] = {"MGC && mgc_ref","MGC","FLJ", "FLJ && (FLJ_HRIandIMS || FLJ_KDRI)", "KIAA","DKFZ", 0} ; 
	  char *libs2[] = {"MGC.filtered","MGC","FLJ", "FLJ.filtered", "KIAA","DKFZ", 0} ; 
	  int idefect, ilib ; 
	  AC_KEYSET ks1, ks2 ; 

	  for (idefect = 0 ; defect[idefect] ; idefect++)
	    {
	      ks1 = ac_dbquery_keyset (db
				       , messprintf ("find tg %s && alibaba && !masked && (excellent || good || partial) ; >read ; >cdna_clone", tag[idefect])
				       , h) ; 
	      for (ilib=0; libs[ilib];ilib++)
		{
		  ks2 = ac_ksquery_keyset (ks1, messprintf ("library = %s ; >read", libs[ilib]), h) ;
		  gExportTable4 (db, ks2, "Libraries", libs2[ilib], defect[idefect]) ;
		  ac_free (ks2) ;
		}
	      ac_free (ks1) ;
	    }
	}
    }

  /*
   * Table for the paper: problems in the libraries
   */
  if (! xml)
    {
      gOutSection (xml, "Percentage of problematic clone alignments, sorted by project") ;
      freeOutf ("considering only excellent or good or partial alignments outside of the hyper variable regions") ;
      gOutBreak (xml);
      freeOutf ("The structural problems are:Partial/Transposition/Mosaic/Inversion/genomic_contamination/internal_deletion\n") ;
      gOutBreak (xml);
     
      gOutTableStart (xml, "ProjectDefects2") ;
      gOutTitles (xml, "Project\tcDNA\tStructural problem\tUnaligned\tSuspected genome deletion\tTransposition\tInversion\tMosaic\tPartial\tSuspected genomic contamination\tNo intron\t") ;
      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;   /* not masked, not bad */
	  gOutInt (xml, r->nCumul, TRUE, r->nOk) ; /* union of the following problems */
	  gOutInt (xml, r->nUnaligned, TRUE, r->nOk) ;
	  gOutInt (xml, r->ngDel, 2, r->nOkKb) ;
	  gOutInt (xml, r->nTransposition, 2, r->nOkKb) ;
	  gOutInt (xml, r->nInversion, 2, r->nOkKb) ;
	  gOutInt (xml, r->nMosaic, 2, r->nOkKb) ;
	  gOutInt (xml, r->nPartial, 2, r->nOkKb) ;
	  gOutInt (xml, r->nGenomic, 2, r->nOkKb) ;
	  gOutInt (xml, r->nNoIntron, TRUE, r->nOk) ;
	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }
  /*
   * Table for the web: vector clippng
   */
  if (1)
    {
      gOutSection (xml, "Problematic clone alignments, sorted by project") ;
      freeOutf ("considering only excellent or good or partial alignments outside of the hyper variable regions") ;

      gOutBreak (xml);
     
      gOutTableStart (xml, "ProjectRepeats") ;
      gOutTitles (xml, "Project\tcDNA\tVariable repeat\tUnclipped vector\tSubmitted on wrong strand\tDuplicate clone\tSuspected internal deletion\tAtypical intron") ;
      
      for (jj = 0 ; jj < arrayMax (rr) - 2 ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;
	  gOutIntLink3 (xml, r->nV_repeat, FALSE, 1,"Libraries", dictName (dict, r->nam), "Variable_repeat") ;
	  gOutIntLink3 (xml, r->nUnclipped, FALSE, 1,"Libraries", dictName (dict, r->nam), "Unclipped_vector") ;
	  gOutIntLink3 (xml, r->nReverse, FALSE, 1,"Libraries", dictName (dict, r->nam), "Reverse") ;
	  gOutIntLink3 (xml, r->nDuplicated, FALSE, 1,"Libraries", dictName (dict, r->nam), "Duplicate_clone") ;
	  gOutIntLink3 (xml, r->nSuspectedInternalDeletion, FALSE, 1,"Libraries", dictName (dict, r->nam), "SuspectedInternalDeletion") ;
	  gOutIntLink3 (xml, r->nAtypicalIntron, FALSE, 1,"Libraries", dictName (dict, r->nam), "AtypicalIntron") ;

	  gOutNewLine (xml) ;
	}
      for (; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;
	  gOutInt (xml, r->nV_repeat, FALSE, 1) ;
	  gOutInt (xml, r->nUnclipped, FALSE, 1) ;
	  gOutInt (xml, r->nReverse, FALSE, 1) ;
	  gOutInt (xml, r->nDuplicated, FALSE, 1) ;
	  gOutInt (xml, r->nSuspectedInternalDeletion, FALSE, 1) ;
	  gOutInt (xml, r->nAtypicalIntron, FALSE, 1) ;

	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
      if (xml)
	{
	  char *defect[] = { "Genomic_contamination",  "Unclipped_vector", "Duplicate_clone", "SuspectedInternalDeletion", "AtypicalIntron", 0} ;
	  char *geneDefect[] = { "",  "", "", "", " && Other", 0} ;
	  char *readDefect[] = { "",  "vector_clipping > 2", "Duplicate_of_read", "", "", 0} ;
	  char *cloneDefect[] = { "Suspected_genomic_contamination",  "", "", "Suspected_internal_deletion", "! Suspected_internal_deletion", 0} ;
	  char *libs[] = {"MGC && mgc_ref","MGC","FLJ", "FLJ && (FLJ_HRIandIMS || FLJ_KDRI)", "KIAA","DKFZ", 0} ; 
	  char *libs2[] = {"MGC.filtered","MGC","FLJ", "FLJ.filtered", "KIAA","DKFZ", 0} ; 
	  int idefect, ilib ;
	  AC_KEYSET ks1, ks2 ;
	  
	  for (idefect = 0 ; defect[idefect] ; idefect++)
	    {
	      ks1 = ac_dbquery_keyset (db
				       , messprintf ("find tg alibaba && !masked && (excellent || good || partial) %s ; >read %s ; >cdna_clone %s", geneDefect[idefect], readDefect[idefect], cloneDefect[idefect])
				       , h) ; 
	      for (ilib=0; libs[ilib];ilib++)
		{
		  ks2 = ac_ksquery_keyset (ks1, messprintf ("library = %s ; >read", libs[ilib]), h) ;
		  gExportTable4 (db, ks2, "Libraries", libs2[ilib], defect[idefect]) ;
		  ac_free (ks2) ;
		}
	      ac_free (ks1) ;
	    }
	}
      
    }

  /*
   * Table for the paper: vector clipping
   */
  if (! xml)
    {
      gOutSection (xml, "Percentage of problematic clone alignments, sorted by project") ;
      freeOutf ("considering only excellent or good or partial alignments outside of the hyper variable regions") ;
      gOutBreak (xml);
     
      gOutTableStart (xml, "ProjectDefects2") ;
      gOutTitles (xml, "Project\tcDNA\tVariable repeat\tUnclipped vector\tSubmitted on wrong strand\tDuplicate clone\tSuspected internal deletion\tAtypical intron") ;      
      for (jj = 0 ; jj < arrayMax (rr) ; jj++)
	{
	  r = arrayp (rr, jj, CHROM_STAT) ;
	  if (r->nOk <= 0)
	    continue ;
	  
	  gOutTitle (xml, dictName (dict, r->nam)) ;
	  gOutInt (xml, r->nOk, FALSE, 1) ;
	  gOutInt (xml, r->nV_repeat, 2, r->nOkKb) ;
	  gOutInt (xml, r->nUnclipped, TRUE, r->nOk) ;
	  gOutInt (xml, r->nReverse, TRUE, r->nOk) ;
	  gOutInt (xml, r->nDuplicated, TRUE, r->nOk) ;
	  gOutInt (xml, r->nSuspectedInternalDeletion, 2, r->nOkKb) ;
	  gOutInt (xml, r->nAtypicalIntron, 2, r->nOkKb) ;

	  gOutNewLine (xml) ;
	}
      gOutTableEnd (xml) ;
    }
  ac_free (h) ;
  return ;
} /* gLibraryStats */

/*********************************************************************/

static void gChromQualityStats (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate (32, CHROM_STAT, h) ;
  int ir, x1, x2, err, imap ;
  long int cAli, cErr, nAlibaba ;
  double s = 0 ;
  CHROM_STAT *r ;
  DICT *dict = dictHandleCreate (500, h) ;
  AC_ITER iter = ac_dbquery_iter (db, "Find tg gold && (excellent || good) && !masked && ! multiple", h) ;
  AC_TABLE af ;
  AC_OBJ gene = 0 ;
  const char *map ;
  char *cp, *ch[] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y", 0} ;
  
  for (imap = 0 ; ch[imap] ; imap++)
    dictAdd (dict, ch[imap], 0) ;
  while (ac_free (gene), gene = ac_next_obj (iter))
    {
      map = ac_tag_printable (gene, "IntMap", "toto") ;
      cp = strstr (map, "|") ;
      if (cp) *cp = 0 ;
      if (strlen (map) > 2)
	continue ;
      dictAdd (dict, map, &imap) ;
      r = arrayp (aa, imap, CHROM_STAT) ;
      r->nAlibaba++ ;

      af = ac_tag_table (gene, "Assembled_from", h) ;
      for (ir = 0; af && ir < af->rows ; ir++)
	{
	  x1 = ac_table_int (af, ir, 0, 0) ;
	  x2 = ac_table_int (af, ir, 1, 0) ;
	  err = ac_table_int (af, ir, 5, 0) ;
	  if (x1 && x2)
	    {
	      r->cAli += x2 - x1 + 1 ;
	      r->cErr += err ;
	    }
	}
    }

  printf ("Global quality of the chromosome estimated from the Excellent or Good, not split and not masked GOLD alignments, per chromosome, including the NT random contis\n") ;  
  printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chrom", "Genes", ".    bp", "Err", "ratio/kb", "1.96*sigma", "err/100kb") ;
  cAli = cErr = nAlibaba = 0 ;
  for (imap = 0 ; ch[imap] ; imap++)
    {
      r = arrayp (aa, imap, CHROM_STAT) ;
      nAlibaba += r->nAlibaba ;
      cAli += r->cAli ;
      cErr += r->cErr ;
      s = r->cErr ; s = 1.96 * 1000.0 * sqrt (s) / (1 +r->cAli) ;
      printf ("%s\t%d\t%d\t%d\t%g\t%g\t%d\n", ch[imap], r->nAlibaba, r->cAli, r->cErr, 1000.0 * r->cErr/(1+r->cAli), s, (int) (100000.0 * r->cErr/(1+r->cAli))) ;
    }
  printf ("\n") ;
  s = cErr ; s = 1.96 * 1000.0 * sqrt (s) / (1 + cAli) ;
  printf ("%s\t%ld\t%ld\t%ld\t%g\t%g\t%d\n", "Total", nAlibaba, cAli, cErr, 1000.0 * cErr/(1+cAli), s, (int) (100000.0 * cErr/(1+cAli))) ;
  printf ("\n") ;

  ac_free (h) ;
} /* gChromQualityStats */

/*********************************************************************/

static void gLibQualityStats (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate (32, CHROM_STAT, h) ;
  int ir, x1, x2, err, ilib ;
  long int cAli, cErr, nAlibaba ;
  double s = 0 ;
  CHROM_STAT *r ;
  AC_ITER iter = 0 ;
  AC_TABLE af ;
  AC_OBJ gene = 0 ;
  /*
  char *libs[] = {"MGC && mgc_ref","MGC","FLJ", "FLJ && (FLJ_HRIandIMS || FLJ_KDRI)", "KIAA","DKFZ", 0} ; 
  char *libs2[] = {"MGC.filtered","MGC","FLJ", "FLJ.filtered", "KIAA","DKFZ", 0} ; 
  */
  char *libs[] = {"MGC","FLJ", "KIAA","DKFZ", "IS NM*", 0} ; 
  char *libs2[] = {"MGC","FLJ","KIAA","DKFZ", "NM", 0} ; 

  for (ilib = 0 ; libs[ilib] ; ilib++)
    {
      iter = ac_dbquery_iter (db
			      , messprintf ("Find cdna_clone %s; > read ; >alibaba ; gold && (excellent || good) && !masked && ! multiple", libs[ilib])
			      , h) ;
      while (ac_free (gene), gene = ac_next_obj (iter))
	{
	  r = arrayp (aa, ilib, CHROM_STAT) ;
	  r->nAlibaba++ ;
	  
	  af = ac_tag_table (gene, "Assembled_from", h) ;
	  for (ir = 0; af && ir < af->rows ; ir++)
	    {
	      x1 = ac_table_int (af, ir, 0, 0) ;
	      x2 = ac_table_int (af, ir, 1, 0) ;
	      err = ac_table_int (af, ir, 5, 0) ;
	      if (x1 && x2)
		{
		  r->cAli += x2 - x1 + 1 ;
		  r->cErr += err ;
		}
	    }
	}
      ac_free (iter) ;
    }
  printf ("Global quality of the libraries estimated from the Excellent or Good, not split and not masked GOLD alignments\n") ;  
  printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Library", "#clone", ".    bp", "Err", "ratio/kb", "1.96*sigma", "err/100kb") ;
  cAli = cErr = nAlibaba = 0 ;
  for (ilib = 0 ; libs[ilib] ; ilib++)
    {
      r = arrayp (aa, ilib, CHROM_STAT) ;
      nAlibaba += r->nAlibaba ;
      cAli += r->cAli ;
      cErr += r->cErr ;
      s = r->cErr ; s = 1.96 * 1000.0 * sqrt (s) / (1 +r->cAli) ;
      printf ("%s\t%d\t%d\t%d\t%g\t%g\t%d\n", libs2[ilib], r->nAlibaba, r->cAli, r->cErr, 1000.0 * r->cErr/(1+r->cAli), s, (int) (100000.0 * r->cErr/(1+r->cAli))) ;
    }
  printf ("\n") ;
  s = cErr ; s = 1.96 * 1000.0 * sqrt (s) / (1 + cAli) ;
  printf ("%s\t%ld\t%ld\t%ld\t%g\t%g\t%d\n", "Total", nAlibaba, cAli, cErr, 1000.0 * cErr/(1+cAli), s, (int) (100000.0 * cErr/(1+cAli))) ;
  printf ("\n") ;
  
  ac_free (h) ;
} /* gLibQualityStats */

/*********************************************************************/

static void gCountN (AC_DB db, char *template)
{
  AC_HANDLE h = 0 ; 
  AC_ITER iter = 0 ;
  AC_TABLE clips = 0 ;
  AC_OBJ est = 0 ;
  int nn, n2, x1, x2, x ;
  char *cp ;

  iter = ac_dbquery_iter (db, messprintf ("Find est %s", template ? template : ""), h) ;

  while (ac_free (est), est = ac_next_obj (iter))
    {
      h = ac_new_handle () ;
      clips = ac_tag_table (est, "vector_clipping", h) ;
      if (clips)
	{
	  x1 = ac_table_int (clips, 0, 0, 1) ;
	  x2 = ac_table_int (clips, 0, 1, 99999999) ;
	}
      else
	{ x1 = 1 ; x2 = 99999999 ; }
      nn = n2 = x = 0 ;
      cp = ac_obj_dna (est, h) ; 
      if (cp) while (x++, *cp)
	if (*cp++ == 'n')
	  {
	    nn++ ;
	    if (x >= x1 && x <= x2) n2++ ;
	  }
      if (nn)
	printf ("Sequence \"%s\"\nNumber_of_N %d %d\n\n", ac_name (est), nn, n2) ;
      ac_free (h) ;
    }
  ac_free (iter) ;
} /* gCountN */

/*********************************************************************/

static void gIntronDistrib (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER genes ;
  AC_TABLE Spl ;
  AC_OBJ gene = 0, read, clone ;
  Array Gg, Other, Suspect, Gdel, Vrep ;
  const char *cp ;
  int nGenes = 0 ;
  int state, ir, len, minLen = 0, maxLen = 0, c1, c2, c3, c4, c5, isClassic, isSuspect ;
  int nGtInFirst, nGtInLast, nGtInDoublet, nGtInAnyShort, nGtInAnyLong, nGcInAnyLong, nAtInAnyLong ;
  int nOtherInFirst, nOtherInLast, nOtherInDoublet, nOtherInAnyShort, nOtherInAnyLong ;
  int nOtherLongFlagged, nGcLongFlagged, nOtherShortFlagged, nGcShortFlagged ;
  int lGtLast, lGtNonLast = 0 ;
  int lOtherLast, lOtherNonLast = 0 ;

  lGtLast = lGtNonLast = 0 ;
  lOtherLast = lOtherNonLast = 0 ;
  nGtInFirst = nGtInLast = nGtInDoublet = nGtInAnyShort = nGtInAnyLong  = nGcInAnyLong  = nAtInAnyLong = 0 ;
  nOtherInFirst = nOtherInLast = nOtherInDoublet = nOtherInAnyShort = nOtherInAnyLong = 0 ;
  nOtherLongFlagged = nGcLongFlagged = nOtherShortFlagged = nGcShortFlagged = 0 ;
  
  Gg = arrayHandleCreate (2000, int, h) ;
  Other = arrayHandleCreate (2000, int, h) ;
  Suspect = arrayHandleCreate (2000, int, h) ;
  Gdel = arrayHandleCreate (2000, int, h) ;
  Vrep = arrayHandleCreate (2000, int, h) ;

  genes = ac_dbquery_iter (db, "Find tg gold && excellent && genome_deletion", h) ;
  while (ac_free (gene), gene = ac_next_obj (genes))
    {
      if ((c1 = ac_tag_int (gene, "Genome_deletion", 0))
	  && c1 > 0 && c1 <= 1000)
	array (Gdel, 1000 - c1, int)++ ; /* register as negative */
    }
  ac_free (genes) ;
  
  genes = ac_dbquery_iter (db, "Find tg  gold && excellent && (intron_boundaries || variable_repeat)", h) ;
  while (ac_free (gene), gene = ac_next_obj (genes))
    {
      Spl = ac_tag_table (gene, "Splicing", h) ;
      if (!Spl)
	continue ;
      read = ac_tag_obj (gene, "Read", h) ;
      clone = ac_tag_obj (read, "cDNA_clone", h) ;
      isSuspect = ac_has_tag (clone, "Suspected_internal_deletion") ;
      ac_free (clone) ; ac_free (read) ;
      nGenes++ ;
      for (state = 0, ir = Spl->rows - 1 ; ir >= 0 ; ir--)
	if (strstr (ac_table_printable (Spl, ir, 2, ""), "tron"))
	  {
	    len = ac_table_int (Spl, ir, 5, -1) ;
	    cp = ac_table_printable (Spl, ir, 3, "") ;
	    if (strstr (cp, "gt_ag") ||strstr (cp, "gc_ag") || strstr (cp, "at_ac"))
	      isClassic = 1 ; 
	    else if (1 || !strstr (cp, "ct_ac"))  /*  merge v_rep  and ct_ac into Other */
	      isClassic = 2 ;  
	    else
	      isClassic = 0 ;
	    if (len >= -100 && len <= 300)
	      {
		if (minLen > len) minLen = len ;
		if (maxLen < len) maxLen = len ;
		if (state == 0) /* state == 0 ->  last exon */
		  {
		    if (isClassic == 1)
		      {
			if (len > 0)
			  array (Gg, len + 1000, int)++ ;
			else if (len < 0)
			  array (Vrep, len + 1000, int)++ ;
		      }
		    else if (isClassic == 2)
		      {
			if (0 && isSuspect)
			  array (Suspect, len + 1000, int)++ ;
			else
			  {
			    if (len > 0)
			      array (Other, len + 1000, int)++ ;
			    else if (len < 0)
			      array (Vrep, len + 1000, int)++ ;
			  }
		      }
		  }
	      }
	    if (len < 65)
	      {
		if (isClassic == 1)
		  nGtInAnyShort++ ;
		if (isClassic == 2)
		  nOtherInAnyShort++ ;
		if (ir == 1 && state == 2)
		  {
		    if (isClassic == 1)
		      nGtInFirst++ ;
		    if (isClassic == 2)
		      nOtherInFirst++ ;
		  }
		else if (ir == 1 && Spl->rows == 3)
		  {
		    if (isClassic == 1)
		      nGtInDoublet++ ;
		    if (isClassic == 2)
		      nOtherInDoublet++ ;
		  }
		else if (state == 0)
		  {
		    if (isClassic == 1)
		      nGtInLast++ ;
		    if (isClassic == 2)
		      nOtherInLast++ ;
		  }
		if (isSuspect)
		  {
		    if (isClassic == 1)
		      nGcShortFlagged++  ;
		    if (isClassic == 2)
		      nOtherShortFlagged ++ ;
		  }
	      }
	    else
	      {
		if (isClassic == 1) 
		  {
		    if (strstr (cp, "gt_ag")) 
		      nGtInAnyLong++ ;
		    if (strstr (cp, "gc_ag")) 
		      nGcInAnyLong++ ;
		    if (strstr (cp, "at_ac")) 
		      nAtInAnyLong++ ;
		  }
		if (isClassic == 2)
		  nOtherInAnyLong++ ;
		if (state == 0)
		  {
		    if (isClassic == 1)
		      lGtLast++ ;
		    if (isClassic == 2)
		      lOtherLast++ ;
		  }
		else
		  {
		    if (isClassic == 1)
		      lGtNonLast++ ;
		    if (isClassic == 2)
		      lOtherNonLast++ ;
		  }
		if (isSuspect)
		  {
		    if (isClassic == 1)
		      nGcLongFlagged++  ;
		    if (isClassic == 2)
		      nOtherLongFlagged ++ ;
		  }
	      }
	    if (state == 0) state = 1 ;
	    if (state == 1 && len > 65) state = 2 ;
	  }	  
      ac_free (Spl) ;
    }
  freeOutf ("Found %d genes with introns  \n", nGenes) ;
  freeOutf ("The gene is cut on long >65bp introns\n"
	    "   first is above a long intron\n"
	    "   last is last so below last long intron\n"
	    "   doublet is the unique short intron of the gene\n"
	    ) ;

  freeOutf ("Type\tFirst\tLast\tDoublet\tOther\tAnySht\tlongLast\tlongNonLast\tLongGt\tLongGc\tLongAt\tLongFlagged\tShortFlagged\n") ;
  freeOutf ("Classic\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
	    , nGtInFirst, nGtInLast, nGtInDoublet
	    , nGtInAnyShort - nGtInFirst - nGtInLast - nGtInDoublet
	    , nGtInAnyShort, lGtLast, lGtNonLast
	    , nGtInAnyLong, nGcInAnyLong, nAtInAnyLong
	    , nGcLongFlagged, nGcShortFlagged
	    ) ;
  freeOutf ("Other\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n\n"
	    ,nOtherInFirst, nOtherInLast, nOtherInDoublet
	    , nOtherInAnyShort - nOtherInFirst - nOtherInLast - nOtherInDoublet
	    , nOtherInAnyShort, lOtherLast, lOtherNonLast
	    , nOtherInAnyLong , 0, 0
	    , nOtherLongFlagged, nOtherShortFlagged
	    ) ;

  /* equalize the lengths one above what i wish to export */
  if (maxLen < 301) maxLen = 301 ;
  array (Gg, maxLen + 1000, int) = 0 ;
  array (Other, maxLen + 1000, int) = 0 ;
  array (Suspect, maxLen + 1000, int) = 0 ;
  array (Gdel, maxLen + 1000, int) = 0 ;
  array (Vrep, maxLen + 1000, int) = 0 ;
  c1 = c2 = c3 = c4 = c5 = 0 ;
  freeOutf ("Distribution of the lengths of the last intron, distinguishing\n"
	    "classic gt-ag,gc-ag, at-ac intron, from non classic\n"
	    "#\tClassic\t!classic\tGdel\nVrepeat\n") ;
  if (1) for (ir = -100 ; ir <= 300 ; ir++)
    {
      c1 += array (Gg, ir + 1000, int) ;
      c2 += array (Other, ir + 1000, int) ;
      c3 += array (Suspect, ir + 1000, int) ;
      c4 += array (Gdel, ir + 1000, int) ;
      c5 += array (Vrep, ir + 1000, int) ;

      printf ("%d\t%d\t%d\t%d\t%d\n", ir
	      , array (Gg, ir + 1000, int)
	      , array (Other, ir + 1000, int)
	      , array (Gdel, ir + 1000, int)
	      , array (Vrep, ir + 1000, int)
) ;
    }
  printf ("cumul\t%d\t%d\t%d\t%d\n", c1, c2, c4, c5) ;
  printf ("Out of %d genes with a last intron in the same quality range\n",  nGenes) ;
  ac_free (h) ;
} /* gIntronDistrib */

/*********************************************************************/

static int gIntron6mersCount (char *dna)
{
  char cc, *cp, *cq ;
  int i, nn = 0, max = strlen (dna) ;
  DICT *dict = dictCreate (1024) ;

  for (i = 0, cp = dna, cq = dna + i + 6 ; i < max - 6 ; cp++, cq++, i++)
    {
      cc = *cq ;
      *cq = 0 ;
      if (!dictAdd (dict, cp, 0))
	nn++ ; /* already known */
      *cq = cc ;
    }
  
  dictDestroy (dict) ;
  return nn ;
} /* gIntron6mersCount */

/*********************************************************************/
/* look for repeated 6mers inside or near a blue intron */
static void gIntron6mersOne (AC_OBJ tg, Array histoLen,  Array histo6mers, int shift)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ cosmid = 0 ;
  AC_TABLE gSpl, gCovers ;
  int nn, ir, a1, a2, b1, b2, x1, x2 ;
  char *dna ;
  const char *cp ;
 
  cosmid = ac_tag_obj (tg, "Genomic_sequence", h) ;
  gCovers = ac_tag_table (tg, "Covers", h) ;
  a1 = ac_table_int (gCovers, 0, 2, 0) ;
  a2 = ac_table_int (gCovers, 0, 3, 0) ;
  gSpl = ac_tag_table (tg, "Splicing", h) ;
  for (ir = 0 ; gSpl && ir < gSpl->rows ; ir++)
    {
      if (!strstr (ac_table_printable (gSpl, ir, 2, ""), "xon"))
	continue ;
      cp = ac_table_printable (gSpl, ir, 3, "gt_ag") ;
      if (!strcmp (cp, "gt_ag") ||
	  !strcmp (cp, "gc_ag") ||
	  !strcmp (cp, "at_ac"))
	continue ;
      x1 = ac_table_int (gSpl, ir, 0, 1) ;
      x2 = ac_table_int (gSpl, ir, 1, 1) ;
      if (x1 < 6)
	continue ;
      if (x2 - x1 > 30)
	continue ;
      x1 -= 6 ; x2 += 6 ;
      x1 -= shift ; x2 -= shift ;
      if (a1 < a2) { b1 = a1 + x1 - 1 ; b2 = a1 + x2 - 1 ; }
      else { b1 = a1 - x1 + 1 ; b2 = a1 - x2 + 1 ; }
      if ((dna = ac_zone_dna (cosmid, b1, b2, h)))
	{
	  nn = gIntron6mersCount (dna) ;
	  array (histoLen, x2 - x1 + 1 - 12, int)++ ;
	  array (histo6mers, nn, int)++ ;
	}
    }

  ac_free (h) ;
  return ;
} /* gIntron6mers */

/*********************************************************************/
/* look for repeated 6mers inside or near a blue intron */
static void gIntron6mers (AC_DB db, BOOL xml)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = ac_dbquery_iter (db, "Find tg gold && other", h) ;
  AC_OBJ tg = 0 ;
  Array histoLen = arrayHandleCreate (1024, int, h) ;
  Array histo6mers = arrayHandleCreate (1024, int, h) ;
  Array histoLen2 = arrayHandleCreate (1024, int, h) ;
  Array histo6mers2 = arrayHandleCreate (1024, int, h) ;
  int i ;

  while (ac_free (tg), (tg = ac_next_obj (iter)))
    {
      gIntron6mersOne (tg, histoLen, histo6mers, 0) ;
      gIntron6mersOne (tg, histoLen2, histo6mers2, 300) ;
    }
  
  if (arrayMax (histoLen) < arrayMax (histo6mers))
    { 
      i = arrayMax (histo6mers) ;
      if (i)
	array (histoLen, i-1, int) = 0 ;
    }
  else
    { 
      i = arrayMax (histoLen) ;
      if (i)
	array (histo6mers, i-1, int) = 0 ;
    }

  printf ("    n  len  6mer\ttemoin\n") ;
  for (i = 0 ; i < arrayMax (histoLen) ; i++)
    printf ("%d\t%d\t%d\t%d\n", i, array (histoLen, i, int) , array (histo6mers, i, int), array (histo6mers2, i, int)) ;

  ac_free (h) ;
  return ;
} /* gIntron6mers */

/*********************************************************************/
/************ Utility to assign the gene to a contig *****************/
/*********************************************************************/

typedef struct { int map, tile, nextTile, a1, a2, a3, a4 ; } MAP ;
static int mapOrder (const void *a, const void *b)
{
  const MAP *ma = (const MAP *)a, *mb = (const MAP *)b ;
  if (ma->map != mb->map)
    return ma->map - mb->map ;
  if (ma->a1 != mb->a1)
    return ma->a1 - mb->a1 ;
  if (ma->a2 != mb->a2)
    return ma->a2 - mb->a2 ;
  return ma->tile - mb->tile ;
}

static int gAssignTg (AC_DB db, char *template)
{
  AC_ITER iter = 0 ; 
  AC_OBJ tile, nextTile, tg ;
  AC_TABLE im = 0, jm = 0 ;
  int nn = 0 ;
  const char *cp, *qq ;
  int ii, jj = 0, i, a1, a2, b1, b2 ;
  MAP* up ;
  AC_HANDLE h = ac_new_handle () ;
  Array hits = arrayHandleCreate (1000, MAP, h) ;
  DICT *dict = dictHandleCreate (100, h) ;
  KEYSET vTiles = keySetHandleCreate (h) ;

  /* first we collect the map position of each tile */
  iter = ac_query_iter (db, TRUE, "find sequence genomic", 0, 0) ; /* all alignments */
  while ((tile = ac_next_obj (iter)))
    {
      im = ac_tag_table (tile, "IntMap", 0) ;
      if (im && im->cols >= 3)
	{
	  cp = ac_table_printable (im, 0, 0, 0) ;
	  a1 = ac_table_int (im, 0, 1, 0) ;
	  a2 = ac_table_int (im, 0, 2, 0) ;
	  if (cp)
	    {
	      dictAdd (dict, cp, &ii) ;
	      dictAdd (dict, ac_name(tile), &i) ;
	      up = arrayp (hits, jj++, MAP) ;
	      up->map = ii ;
	      up->tile = i ;
	      up->a1 = a1 ; up->a2 = a2 ; 
	      nextTile = ac_tag_obj (tile, "Overlap_right", 0) ;
	      if (nextTile)
		{
		  dictAdd (dict, ac_name (nextTile), &i) ;
		  up->nextTile = i ;
		  jm = ac_tag_table (nextTile, "IntMap", 0) ;
		  if (jm && jm->cols >= 3)
		    {
		      up->a3 = ac_table_int (jm, 0, 1, 0) ;
		      up->a4 = ac_table_int (jm, 0, 2, 0) ;
		    }
		  ac_free (jm) ;
		}
	      ac_free (nextTile) ;
	    }
	}
      ac_free (im) ;
      ac_free (tile) ;
    }
  arraySort (hits, mapOrder) ;
  ac_free (iter) ;

  /* then we assign each tg to the correct tile */
  qq = messprintf ("find transcribed_gene  IS \"%s\" ; IntMap && !Genomic_sequence", template) ;
  iter = ac_query_iter (db, TRUE, qq, 0, 0) ; /* all tg */
  while ((tg = ac_next_obj (iter)))
    {
      im = ac_tag_table (tg, "IntMap", 0) ;
      if (im && im->cols >= 3)
	{
	  cp = ac_table_printable (im, 0, 0, 0) ;
	  if (cp && dictFind (dict, cp, &ii))
	    {
	      a1 = ac_table_int (im, 0, 1, 0) ;
	      a2 = ac_table_int (im, 0, 2, 0) ;
	      if (a1 < a2) { b1 = a1 ; b2 = a2 ; }
	      else { b1 = a2 ; b2 = a1 ; }
	      for (i = 0 ; i < arrayMax (hits) ; i += 300)
		{
		  up = arrp (hits, i, MAP) ;
		  if (up->map >= ii)
		    break ;
		}
	      i -= 300 ; if (i < 0) i = 0 ;
	      for (; i < arrayMax (hits) ; i++)
		{
		  up = arrp (hits, i, MAP) ;
		  
		  if (up->map > ii)
		    break ;
		  if (up->map == ii && up->a1 <= b1 && up->a2 > b1 && (!up->nextTile || up->a3 > b1))
		    {
		      nn++ ;
		      freeOutf ("Transcribed_gene  \"%s\"\nCovers %d \"bp from\" %d  %d // b1=%d  b2=%d\n\n"
			      , ac_name (tg), b2 - b1 + 1, a1 - up->a1 + 1, a2 - up->a1 + 1, b1, b2
			      ) ;
		      if (up->a1 <= b2 && up->a2 >= b2)
			{
			  freeOutf ("Sequence \"%s\"\nTranscribed_gene \"%s\" %d %d\n\n", 
				  dictName (dict, up->tile),
				  ac_name (tg), a1 - up->a1 + 1, a2 - up->a1 + 1) ;
			}
		      else if (up->a1 <= b2 && up->a4 >= b2)
			{
			  freeOutf ("Sequence \"%s_%s\"\nTranscribed_gene \"%s\" %d %d // b1=%d  b2=%d\n\n", 
				  dictName (dict, up->tile),
				  dictName (dict, up->nextTile),
				  ac_name (tg), a1 - up->a1 + 1, a2 - up->a1 + 1, b1, b2) ;
			}
		      else /* create a virtual genomic tile */
			{
			  AC_ITER iter ;
			  AC_TABLE im ;
			  int ir, v1, v2, v3, w3 ;
			  AC_OBJ source ;

			  freeOutf ("Sequence \"%s__v\"\nTranscribed_gene \"%s\" %d %d // b1=%d  b2=%d\n\n", 
				  dictName (dict, up->tile),
				  ac_name (tg), a1 - up->a1 + 1, a2 - up->a1 + 1, b1, b2) ;
			  iter = ac_query_iter (db, 1,
						messprintf ("Find sequence IS \"%s\"; > source", dictName (dict, up->tile))
						, 0, 0) ;
			  if ((source = ac_next_obj (iter))) 
			    {
			      char vTileNam [1024] ;
			      
			      sprintf (vTileNam, "%s__v", dictName (dict, up->tile)) ;
			      im = ac_tag_table (source, "Subsequence", 0) ;
			      for (ir = 0 ; im && ir < im->rows ; ir++)
				{
				  if (strcmp (dictName (dict, up->tile), ac_table_printable (im, ir, 0, "")))
				    continue ;
				  v1 = ac_table_int (im, ir, 1, 0) ;
				  v2 = ac_table_int (im, ir, 2, 0) ;
				  if (v1 > v2)
				    messcrash ("stupid tile coords") ;
				  v3 = v1 + b2 - up->a1 + 1 ;
				  w3 = keySet (vTiles, up->tile) ;
				  if (!w3) /* verify it does not allready exist in the database */
				    {
				      AC_ITER iter2 = ac_query_iter (db, 1,
							    messprintf ("Find sequence IS \"%s\" ; > source"
									, vTileNam
									), 0, 0) ;
				      AC_OBJ Source = iter2 ? ac_next_obj (iter2) : 0 ;
				      AC_TABLE km = 0 ;
				      int kr ;
				      if (Source && (km = ac_tag_table (Source, "Subsequence", 0)))
					{
					  for (kr = 0 ; kr < km->rows ; kr++)
					    {
					      if (!strcmp (vTileNam, ac_table_printable (km, kr, 0, "")))
						{
						  w3 = ac_table_int (km, kr, 2, 0) ;
						  break ;
						}						
					    }
					}
				      ac_free (km) ;
				      ac_free (Source) ;
				      ac_free (iter2) ;
				    }
				  if (v3 > w3) 
				    { w3 = keySet (vTiles, up->tile) = v3 ; }
				  freeOutf ("Sequence \"%s\"\nSubsequence \"%s\" %d %d\n\n"
					  , ac_name (source)
					  , vTileNam
					  , v1, w3) ;
				  break ;
				}
			      ac_free (im) ;
			      ac_free (source) ;
			    }
			  ac_free (iter) ;			  
			}
		      break ;
		    }
		}
	      
	    }
	}
      ac_free (im) ;
      ac_free (tg) ;
    }
  
  ac_free (h) ;
  
  fprintf (stderr, "// assigned %d genes \n", nn ) ;
  return nn ;
} /* gAssignTg */

/*************************************************************************************/
/****** Utility to change the starnd of some manually selected alignments ************/
/*************************************************************************************/

static int gDoReverseCtAc (AC_OBJ gene)
{
  AC_TABLE im = 0 ;
  const char *cp ;
  int ir, nn = 0, max, nexon, a1, a2, x1, x2 ;

  max = ac_tag_int (gene, "Covers", 0) ;

  
  /* flip the intmap */
  nn++ ;
  freeOutf ("Transcribed_gene \"%s\"\n", ac_name(gene)) ;
  im = ac_tag_table (gene, "IntMap", 0) ;
  if (im && im->cols >= 3)
    {
      cp = ac_table_printable (im, 0, 0, 0) ;
      a1 = ac_table_int (im, 0, 1, 0) ;
      a2 = ac_table_int (im, 0, 2, 0) ;
      freeOutf ("-D IntMap\n") ;
      freeOutf ("IntMap \"%s\" %d %d\n", cp, a2, a1) ;
      max = a2 - a1 + 1 ;
      if (max < 0) max = -max ;
    }
  ac_free (im) ;
  if (!max) return 0 ;
  /* flip Splicing */
  freeOutf ("-D Splicing\n") ;
  im = ac_tag_table (gene, "Splicing", 0) ;
  nexon = 0 ;
  for (ir = im->rows - 1 ; ir >= 0 ; ir--)
    {
      x1 = ac_table_int (im, ir, 0, 0) ;
      x2 = ac_table_int (im, ir, 1, 0) ;
      cp = ac_table_printable (im, ir, 2, 0) ;
      x1 = max - x1 + 1 ;
      x2 = max - x2 + 1 ;
      freeOutf ("Splicing %d %d %s ", x2, x1, cp) ;
      if (strstr (cp, "xon"))
	freeOutf ("%d Length %d bp", ++nexon, x1 - x2 + 1) ;
      freeOutf ("\n") ;
    }
  
  /* flip assembled_from */
  freeOutf ("-D Assembled_from\n") ;
  freeOutf ("-D Rank\n-D Bug\n-D Split\n") ;
  im = ac_tag_table (gene, "Assembled_from", 0) ;
  nexon = 0 ;
  for (ir = im->rows - 1 ; ir >= 0 ; ir--)
    {
      x1 = ac_table_int (im, ir, 0, 0) ;
      x2 = ac_table_int (im, ir, 1, 0) ;
      cp = ac_table_printable (im, ir, 2, 0) ;
      x1 = max - x1 + 1 ;
      x2 = max - x2 + 1 ;
      freeOutf ("Assembled_from %d %d \"%s\" ", x2, x1, cp) ;
      x1 = ac_table_int (im, ir, 3, 0) ;
      x2 = ac_table_int (im, ir, 4, 0) ;
      freeOutf ("%d %d\n", x2, x1) ;
    }
  
  
  /* ask for recomputation */
  freeOutf ("-D covers\n") ;
  freeOutf ("-D Wrong_strand\n") ;
  freeOutf ("-D G_strand\n") ;
  freeOutf ("-D Genomic_sequence\n") ;
  freeOutf ("-D Smith_Waterman_done\n") ;
  
  freeOutf ("\n") ;
  
  return 1 ;
} /* gDoReverseCtAc */

static int gDoReverseEBI (AC_OBJ gene)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE im = 0 ;
  AC_OBJ est ;
  const char *cp ;
  int ir, max = 0, x1, x2 ;

  est = ac_tag_obj (gene, "Read", h) ; 
  im = ac_tag_table (est, "DNA", h) ;
  if (im && im->cols >= 2)
    max = ac_table_int (im, 0, 1, 0) ;
  ac_free (im) ;
  
  freeOutf ("Transcribed_gene \"%s\"\n", ac_name(gene)) ;
  /* flip assembled_from */
  freeOutf ("-D Assembled_from\n") ;
  freeOutf ("-D Rank\n-D Bug\n") ;
  freeOutf ("Wrong_strand\n") ;
  im = ac_tag_table (gene, "Assembled_from", 0) ;

  for (ir = 0 ; ir < im->rows ; ir++)
    {
      x1 = ac_table_int (im, ir, 0, 0) ;
      x2 = ac_table_int (im, ir, 1, 0) ;
      cp = ac_table_printable (im, ir, 2, 0) ;
      freeOutf ("Assembled_from %d %d \"%s\" ", x1, x2, cp) ;
      x1 = ac_table_int (im, ir, 3, 0) ;
      x2 = ac_table_int (im, ir, 4, 0) ;
      if (ir == 0) max = x1 + ac_table_int (im, im->rows -1, 4, 0) - 1 ;
      x1 = max - x1 + 1 ;
      x2 = max - x2 + 1 ;
      freeOutf ("%d %d\n", x1, x2) ;
    }
  
  
  /* ask for recomputation */
  freeOutf ("-D Smith_Waterman_done\n") ;
  
  freeOutf ("\n") ;
  ac_free (h) ;
  return 1 ;
} /* gDoReverseEBI */

/*************************************************************************************/

static int gReverseCtAc (AC_DB db, char *template, BOOL ebi)
{
  AC_ITER iter = 0 ; 
  AC_OBJ gene ;
  int nn = 0, nntot = 0 ;

  if (ebi) iter = ac_dbquery_iter (db, "find tg IS F_E_*  ", 0) ;
  else iter = ac_dbquery_iter (db, "find tg wrong_strand", 0) ;

  while ((gene = ac_next_obj (iter)))
    {
      nntot++ ;
      if (ebi) nn += gDoReverseEBI (gene) ;
      else  nn += gDoReverseCtAc (gene) ;
      ac_free (gene) ;
    }
  ac_free (iter) ;
  
  fprintf (stderr, "// flipped %d/%d %sgenes\n", nn, nntot, ebi ? "(EBI)" : "(wrong_strand)") ;
  return nn ;
} /* gReverseCtAc */

/*************************************************************************************/

static int gDoCreatePg (AC_OBJ tg)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt ;
  const char *cp ;
  int ir ;

  cp = ac_tag_printable (tg, "Genomic_sequence", 0) ;
  if (!cp) return 0 ;

  printf ("\nSequence \"%s\"", cp) ;
  printf ("Subsequence \"pg_%s\"", ac_name (tg)) ;
  tt = ac_tag_table (tg, "Covers", h) ;
  if (tt)
    printf (" %d %d\n", ac_table_int (tt, 0, 2, 0), ac_table_int (tt, 0, 3, 0)) ;
 
  printf ("\nSequence \"pg_%s\"\nIs_predicted_gene\n", ac_name (tg)) ;
  printf ("Matching_transcribed_gene \"%s\"", ac_name (tg)) ;
  tt = ac_tag_table (tg, "IntMap", h) ;
  if (tt) 
    printf ("IntMap %s %d %d\n"
	    , ac_table_printable (tt, 0, 0, "toto")
	    ,  ac_table_int (tt, 0, 1, 0)
	    ,  ac_table_int (tt, 0, 2, 0)
	    ) ; 
  tt = ac_tag_table (tg, "Splicing", h) ;
  for (ir = 0 ; tt && ir < tt->rows ; ir++)
    if (strstr (ac_table_printable (tt, ir, 2, "toto"), "xon"))
      printf ("Source_exons %d %d\n"
	      ,  ac_table_int (tt, ir, 1, 0)
	      ,  ac_table_int (tt, ir, 2, 0)
	      ) ;
  printf ("\n") ;

  ac_free (h) ;
  return 1 ;
} /* gDoCreatePg */

/*************************************************************************************/
/* a hack to cluster the alignments into gene models */
static int gCreatePg (AC_DB db)
{
  AC_ITER iter = 0 ; 
  AC_OBJ gene ;
  int nn = 0 ;

  iter = ac_dbquery_iter (db, "find tg gold", 0) ;
  while ((gene = ac_next_obj (iter)))
    {
      nn += gDoCreatePg (gene) ;
      ac_free (gene) ;
    }
  ac_free (iter) ;
  
  fprintf (stderr, "// created %d pg\n", nn) ;
  return nn ;
} /* gCreatePg */

/*************************************************************************************/

static void gWorm (AC_DB db)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = ac_dbquery_iter (db, "Find Gene", h) ;
  AC_OBJ obj=0 ;
  int n, ng = 0, i, ii, a1, a2, b1, b2 ;
  const char *cp ;
  DICT *dict = dictHandleCreate (256, h) ;
  BitSet b, ba, bb[256], bb2[256]  ;
  AC_TABLE tt ;

  memset (bb, 0, 256 * sizeof (BitSet)) ;
  memset (bb2, 0, 256 * sizeof (BitSet)) ;
  while (ac_free (obj), (obj = ac_next_obj (iter)))
    {
      tt = ac_tag_table (obj, "IntMap", h) ;
      if (tt && 
	  (cp = ac_table_printable (tt, 0, 0, 0)) &&
	  (a1 = ac_table_int (tt, 0, 1, 0)) &&
	  (a2 = ac_table_int (tt, 0, 2, 0)) &&
	  a1 > 0 && a2 > 0)
	{
	  b1 = a1 < a2 ? a1 : a2 ;
	  b2 = a1 < a2 ? a2 : a1 ;
	  dictAdd (dict, cp, &ii) ;
	  if (!bb[ii]) { bb[ii] = bitSetCreate (20000000, h) ; bb2[ii] = bitSetCreate (20000000, h) ;}
	  if (ii < 256)
	    {
	      b = bb[ii] ;
	      for (i = b1 ; i <= b2 ; i++)
		bitSet (b, i) ;
	      b = bb2[ii] ;
	      if (a1 < a2)
		{
		  for (i = a1 - 1000 ; i <= a1 ; i++)
		    if (i > 0) bitSet (b, i) ;
		}
	      else
		{
		  for (i = a2 - 1000 ; i <= a2 ; i++)
		    if (i > 0) bitSet (b, i) ;
		}
	      ng++ ;
	    }
	}
    }
  
  for (n = ii = 0 ; ii < 256 ; ii++)
     if ((b = bb[ii]))
      n += bitSetCount (b) ;
  
  printf ("Found %d genes occupying %d bp\n", ng, n) ;

  for (n = ii = 0 ; ii < 256 ; ii++)
    if ((b = bb2[ii]))
      n += bitSetCount (b) ;
  printf ("Their 1kb promoter regions occupy %d bp\n", n) ;
 
  for (n = ii = 0 ; ii < 256 ; ii++)
    if ((b = bb[ii]) && (ba = bb2[ii]))
      n += bitSetMINUS (ba, b) ;
  printf ("%d bp of these promoter regions are not masked by another gene\n", n) ;

  ac_free (h) ;
} /* gWorm */

/*************************************************************************************/
/************************* export gold in .ali format ********************************/
/*************************************************************************************/

static int gExportAliExons (char *geneId, char *method, AC_OBJ est, AC_OBJ tg, int gold,
			    AC_OBJ chrom, int g1, int g2, 
			    AC_HANDLE h)
{
  const char *cp ;
  int ii, ii1, jj, jjTotal, nExon ; 
  int a1, a2, x1, x2 ; /* a single aligned exon in tg coords */
  int b1, b2 ; /* same in chrom coords */
  AC_TABLE tt = ac_tag_table (tg, "Assembled_from", h) ;
  /* count all exons to number them in the direction of the EST */
  for (ii = jjTotal = 0 ; ii < tt->rows ; ii++)
    {
      cp = ac_table_printable (tt, ii, 2, "") ;
      if (strcmp (cp, ac_name (est)))
	continue ;
      jjTotal++ ;
    }

  for (ii1 = jj = 0 ; ii1 < tt->rows ; ii1++)
    {
      /* ii1 -> ii: export in the direction of the genome */
      ii = ii1 ;
      cp = ac_table_printable (tt, ii, 2, "") ;
      if (strcmp (cp, ac_name (est)))
	continue ;
      a1 = ac_table_int (tt, ii, 0, -9999) ;
      a2 = ac_table_int (tt, ii, 1, -9999) ;
      x1 = ac_table_int (tt, ii, 3, -9999) ;
      x2 = ac_table_int (tt, ii, 4, -9999) ;
      if (a1 == -9999 ||
	  a1 == -9999 ||
	  a1 == -9999 ||
	  a1 == -9999)
	messcrash ("bad values in Assembled_from tg = %s, est = %s, ii = %d", 
		   ac_name(tg), ac_name (est), ii) ;

      if (g1 < g2)
	{
	  b1 = g1 + a1 - 1 ;
	  b2 = g1 + a2 - 1 ;
	}
      else
	{
	  b1 = g1 - a1 + 1 ;
	  b2 = g1 - a2 + 1 ;
	}
      nExon = ++jj ; 

      /*we expect, in the direction of the gene
	gene = $1 ;
	acc = $2 ;
	method = $3 ;
	nexon = $4 ; 
	x1 = $5 ; x2 = $6 ;    // coords on the cDna
	chrom = $7 ;
	b1 = $8 ; b2 = $9 ;    // coords on the chromosome
      */
      printf("%s\t", geneId) ;
      printf("%s\t", ac_name(est)) ;
      printf ("%s\t", method) ;
      printf("%d\t", nExon) ;
      printf("%d\t%d\t", x1,x2) ;
      printf ("%s\t", ac_name(chrom)) ;
      printf("%d\t%d\t", b1,b2) ;
      printf ("\n") ;
    }    
  return 1 ;
} /* gExportAliExons */

/*************************************************************************************/

static int gExportAliEstTg (AC_OBJ est, AC_OBJ tg, int ii, int gold, AC_HANDLE h)
{
  AC_OBJ chrom = 0 ;
  int g1, g2 ; /* coords of the gene on the chromosome */
  AC_TABLE tt = 0 ;
  char geneId[1024], *method = "Toto" ;

  /* get the intmap of the gene */
  tt = ac_tag_table (tg, "IntMap", h) ;
  if (!tt || ! tt->rows || tt->cols < 3)
    messcrash ("missing IntMap in transcribed_gene %s", ac_name (tg)) ;
  chrom = ac_table_obj (tt, 0, 0, h) ;
  g1 = ac_table_int (tt, 0, 1, 0) ;
  g2 = ac_table_int (tt, 0, 2, 0) ;
  
  /* create the gene_id */
  switch (gold)
    {
    case 1:
      if (ac_has_tag (tg, "Excellent"))
	{
	  method = "Gold_excellent" ;
	  if (ii)
	    sprintf (geneId, "%s_e_%d", ac_name(est), ii) ;
	  else
	    sprintf (geneId, "%s_e", ac_name(est)) ;
	}
      else if (ac_has_tag (tg, "Good"))
	{
	  method = "Gold_good" ;
	  if (ii)
	    sprintf (geneId, "%s_g_%d", ac_name(est), ii) ;
	  else
	    sprintf (geneId, "%s_g", ac_name(est)) ;
	}
      break ;
    case 2:
      {
	if (ac_has_tag (tg, "Partial"))
	  method = "Gold_partial" ;
	else if (ac_has_tag (tg, "Dubious"))
	  method = "Gold_dubious" ;
	else
	  method = "Gold_bad" ;
	if (ii)
	  sprintf (geneId, "%s_o_%d", ac_name(est), ii) ;
	else
	  sprintf (geneId, "%s_o", ac_name(est)) ;
      }
      break ;
    case 3:
      method = "AceView" ;
      memcpy (geneId, ac_name (tg), sizeof (geneId)) ;
      break ;
    case 4:
      method = "NCBI" ;
      memcpy (geneId, ac_name (tg), sizeof (geneId)) ;
      break ;
    case 5:
      method = "UCSC" ;
      memcpy (geneId, ac_name (tg), sizeof (geneId)) ;
      break ;
    }
    
  /* get the coding region in gene coordinates */
  gExportAliExons (geneId, method, est, tg, gold, chrom, g1, g2, h) ;
  
  return 1 ;
} /* gExportAliEstTg */

/*************************************************************************************/

static int gExportAliEst (AC_OBJ est, char *template, int gold)
{
  int ii, go, jj, j2, nn = 0 ;
  AC_HANDLE h = handleCreate () ;
  AC_TABLE tgs = ac_tag_table (est, "From_gene", h) ;
  AC_OBJ tg ;
  const char *cp ;

  for (go = jj = 0 ; go < 2 ; go++)
    for (ii = j2 = 0 ; ii < tgs->rows ; ii++)
      {
	tg = ac_table_obj (tgs, ii, 0, h) ;
	cp = ac_name(tg) ;
	if (template && ! pickMatch (cp, template))
	  continue ;
	switch (gold)
	  {
	  case 1:
	    if (! ac_has_tag (tg, "gold"))
	      continue ;
	    if (!ac_has_tag (tg, "excellent") && ! ac_has_tag (tg, "good"))
	      continue ;
	    break ;
	  case 2:
	    if (! ac_has_tag (tg, "gold"))
	      continue ;
	    if (ac_has_tag (tg, "excellent") || ac_has_tag (tg, "good"))
	      continue ;
	    break ;
	  case 3:
	    if (! ac_has_tag (tg, "AceView"))
	      continue ;
	    break ;
	  case 4:
	    if (! ac_has_tag (tg, "NCBH"))
	      continue ;
	    break ;
	  case 5:
	    if (! ac_has_tag (tg, "Kent"))
	      continue ;
	    break ;
	  case 6:
	    if (! ac_has_tag (tg, "GBIRC"))
	      continue ;
	    break ;
	  }
	if (!go)
	  jj++ ; /* how many j2  will i need */
	else
	  nn += gExportAliEstTg (est, tg, jj > 1 ? ++j2 : 0, gold, h) ;
      }

  
  ac_free (h) ;
  return nn ;
} /* gExportAliEst */

/*************************************************************************************/

static int gExportAli (AC_DB db, char *template, int gold)
{
  AC_ITER iter = 0 ; 
  AC_OBJ est ;
  int nn = 0, nest = 0 ;

  if (template) /* select some gene */
    iter = ac_query_iter (db, TRUE, messprintf ("find tg IS \"%s\"; >read ; from_gene", template), 0, 0) ;
  else if (gold == 0)
    iter = ac_query_iter (db, TRUE, "find est ABI && from_gene", 0, 0) ; /* all alignments */
  else if (gold == 1)
    iter = ac_query_iter (db, TRUE, "find est from_gene", 0, 0) ; /* all alignments */
  else if (gold == 2)
    iter = ac_query_iter (db, TRUE, "find tg !excellent && !good ; >read", 0, 0) ; /* all alignments */
  else if (gold >= 3)
    iter = ac_query_iter (db, TRUE, "find est from_gene", 0, 0) ; /* all alignments */
  else
    {
      fprintf (stderr, "// bizare value for the gold parameter %d\n", gold) ;
      usage () ;
    }
  while ((est = ac_next_obj (iter)))
    {
      nn += gExportAliEst (est, template, gold) ;
      ac_free (est) ;
      nest++ ;
    }
  ac_free (iter) ;
  fprintf (stderr, "// Considered %d ests\n", nest) ;

  return nn ;
} /* gExportAli */

/*************************************************************************************/
/*************************************************************************************/

static int gExportOneGoldTgDna (AC_DB db, AC_OBJ est, AC_OBJ tg, char *buf)
{
  AC_HANDLE h = ac_new_handle () ;
  char cc, *cp, *cq ;

  cp = gGetDNA (tg, h) ;
  
  printf ("> %s\n", buf) ;
  while (strlen (cp) > 50)
    {
      cq = cp + 50 ;
      cc = *cq ; *cq = 0 ;
      printf ("%s\n", cp) ;
      *cq = cc ;
      cp += 50 ;
    }
  printf ("%s\n", cp) ;

  ac_free (h) ;

  return 1 ; 
} /* gExportGoldTgDna */

/*************************************************************************************/

static int gExportOneGoldTgPep (AC_DB db, AC_OBJ est, AC_OBJ tg, char *buf)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, x1, x2, met ;
  char cc, *dna, *cp, *cq = 0, buf2[4] ;

  dna = gGetDNA (tg, h) ;
  gGetOrf (dna, &x1, &x2, &met) ;
  if (x1 <= x2)
    {
      cq = dna ;
      if (met)
	{
	  *cq++ = 'M' ; /* impose Met translation */
	  x1 = met + 3 ;
	}
      for (i = x1 - 1, cp = dna + i, cq = dna ; i < x2 ; i += 3, j++, cp += 3)
	{
	  buf2[0] = dnaEncodeChar [(int)*(cp+0)] ;
	  buf2[1] = dnaEncodeChar [(int)*(cp+1)] ;
	  buf2[2] = dnaEncodeChar [(int)*(cp+2)] ;
	  *cq++ = codon (buf2) ;
	}
      *cq = 0 ;
      
      cp = dna ;
      printf ("> %s\n", buf) ;
      while (strlen (cp) > 50)
	{
	  cq = cp + 50 ;
	  cc = *cq ; *cq = 0 ;
	  printf ("%s\n", cp) ;
	  *cq = cc ;
	  cp += 50 ;
	}
      printf ("%s\n", cp) ;
    }
  ac_free (h) ;

  return 1 ; 
} /* gExportGoldTgPep */

/*************************************************************************************/

static int gExportOneGoldTgAnnot (AC_DB db, AC_OBJ est, AC_OBJ tg, char *buf)
{
  AC_HANDLE h = ac_new_handle () ;
  unsigned const char *cp, *cq ;

  cp = ac_command (db, messprintf ("find tg %s", ac_name(tg)), 0, h) ;
  cp = ac_command (db, "show -a", 0, h) ;
  
  cq = (unsigned const char *) strstr ((const char *)cp, "Transcribed_gene") ; cq++ ;
  cq = (unsigned const char *) strstr ((const char *)cq, "\n") ; cq++ ;
  printf ("Transcribed_gene \"%s\"\n%s\n\n", buf, cq) ;

  ac_free (h) ;

  return 1 ; 
} /* gExportGoldTgAnnot */

/*************************************************************************************/

static int gExportOneGoldAnnot (AC_DB db, AC_OBJ est, int type)
{
  AC_TABLE tt = ac_tag_table (est, "From_gene", 0) ;
  AC_OBJ tg ;
  char buf[4096] ;
  int ir, nn, ntg ;

  for (ir = ntg = 0 ; tt && ir < tt->rows ; ir++)
    {
      tg = ac_table_obj (tt, ir, 0, 0) ;
      if (ac_has_tag (tg, "Gold")) ntg++ ;
    }

  for (ir = nn = 0 ; tt && ir < tt->rows ; ir++)
    {
      tg = ac_table_obj (tt, ir, 0, 0) ;
      if (!ac_has_tag (tg, "Gold")) 
	continue ;
      if (ntg == 1)
	sprintf (buf, "Gold_%s", ac_name(est)) ;
      else
	sprintf (buf, "Gold_%s_repetition_%d", ac_name(est), ++nn) ;
      switch (type)
	{
	case 1: /* annot */
	  gExportOneGoldTgAnnot (db, est, tg, buf) ;
	  break ;
	case 2: /* dna */
	  gExportOneGoldTgDna (db, est, tg, buf) ;
	  break ;
	case 3: /* pep */
	  gExportOneGoldTgPep (db, est, tg, buf) ;
	  break ;
	}
      ac_free (tg) ;
    }
  ac_free (tt) ;

  return ir > 0 ? 1 : 0 ;
} /* gExportGoldAnnot */

/*************************************************************************************/

static int gExportGoldAnnot (AC_DB db, char *template, int gold, int type)
{
  AC_ITER iter = 0 ; 
  AC_OBJ est ;
  int nn = 0, nest = 0 ;

  if (template) /* select some est */
    iter = ac_query_iter (db, TRUE, messprintf ("find est %s", template), 0, 0) ;
  else
    iter = ac_query_iter (db, TRUE, "find est", 0, 0) ; /* all alignments */

  while ((est = ac_next_obj (iter)))
    {
      nn += gExportOneGoldAnnot (db, est, type) ;
      ac_free (est) ;
      nest++ ;
    }
  ac_free (iter) ;
  fprintf (stderr, "// Considered %d ests\n", nest) ;

  return nn ;
} /* gExportGoldAnnot */

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: gintegracost ACEDB [-o outfile] [-x] [-assign | -reverse | -cost | -flag | -stats | -chroms ....] [template]\n") ;
  fprintf (stderr, "// Example:  gintegracost \'?_*42\'\n") ;
  fprintf (stderr, "//    assign: exports an ace file: given gene IntMap, this creates cosmid coords\n") ;
  fprintf (stderr, "//    reverse: exports an ace file which reverses genes flagged as wrong_strand\n") ;
  fprintf (stderr, "//    cost: compute by smith waterman the ali and er-number of each exon and each tg\n") ;
  fprintf (stderr, "//    flag: computes and exports all the flags in the Tg class\n") ;
  fprintf (stderr, "//    stats: exports global statistics on stdout + ancilary files in ./COMPARE\n") ;
  fprintf (stderr, "//    chroms: exports chromosome statistics\n") ;
  fprintf (stderr, "//    chrom_quality: exports error rate per 10kb for each chromosonme\n") ;
  fprintf (stderr, "//    lib_quality: exports error rate per 10kb for each library\n") ;
  fprintf (stderr, "//    count_n: report existence of N in EST\n") ;
  fprintf (stderr, "//    libs: exports library statistics\n") ;
  fprintf (stderr, "//    tg2pg: create a pg out of the tg to allow clustering\n") ;
  fprintf (stderr, "//    export: export all gold in .ali format\n") ;
  fprintf (stderr, "//       -gold [1 | 2] 1: export excellent_or_good gold, 2: all others, 3:AceView, 4:NCBI, 5:Kent\n") ;
  fprintf (stderr, "//    exportGoldAnnot: exports and ace file for all Gold annotations\n") ;
  fprintf (stderr, "//    exportGoldDna: The sequences of the GOLD mRNAs, as derived from the genome, in fasta format, limited to excellent, good and partial Gold alignments\n") ;
  fprintf (stderr, "//    exportGoldPep: The sequences, in aminoacids, of the GOLD longest ORF as seen on the genome\n") ;
  fprintf (stderr, "//   -x: export in ./COMPARE html tables linked to other tables, (valid only with -stats)\n") ;
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
  const char *action ;
  char *template ;
  const char *who = 0 ;
  int gold = 0, outlevel = 0, errCost = 2 ; /* points per error */
  BOOL xml ;

  int methods1[] = { (int)'A',  (int)'B',  (int)'C', (int)'D', (int)'K', (int)'E', (int)'F', (int)'H', (int)'S', (int)'G', (int)'J', (int)'N', 0 } ;
  char *methodTitle1[] = {"AceView", "AceView", "AceView", "AceView", "UCSC", "Exonerate", "Exonerate", "NCBI 34.3", "Splign", "JBIRC", "JBIRC old", "NEDO", 0 } ;

  int methods2[] = { (int)'A', (int)'K', (int)'S', (int)'H',  (int)'G', (int)'J', (int)'N', (int)'E', 0 } ;
  char *methodTag2[] = {"AceView", "Kent", "Splign", "NCBH",  "GBIRC", "JBIRC", "NEDO", "EBI", 0 } ;
  char *methodTitle2[] = {"AceView", "UCSC", "Splign", "NCBI 34.3", "JBIRC", "JBIRC old", "NEDO", "Exonerate", 0 } ;
  AC_DB db ;
  
  if (argc <2)
    usage () ;

  /* consume optional args */
  outfilename = getArgV (&argc, argv, "-out") ;
  xml =  getArg (&argc, argv, "-x") ;
  who = getArgV (&argc, argv, "-who") ;
  {
    const char *ccp = getArgV (&argc, argv, "-gold") ;
    gold =  ccp ? atoi(ccp) : 0 ;
  }
  /* read absolute args */
  dbName = argc>=2 ? argv[1] : "." ;
  action = argc>=3 ? argv[2] : "";
  template = argc>=4 ? strnew (argv[3], 0) : "*" ;

 
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

  if (*template == 'X')
    *template =  '?' ; /* too hard to transmit ? or * via bsub scripts */

  if (!strcmp (action, "-assign"))
    nn = gAssignTg (db, template) ;
  else if (!strcmp (action, "-reverse"))
    nn = gReverseCtAc (db, template, FALSE) ;
  else if (!strcmp (action, "-ebi_reverse"))
    nn = gReverseCtAc (db, template, TRUE) ;
  else if (!strncmp (action, "-cost", 5)) /* silently allow cost costs, do banded_smith_waterman */
    nn = gExportCosts (db, template) ;
  else if (!strcmp (action, "-flag"))
    gFlagAll (methods1, methodTitle1, db, template, errCost) ;
  else if (!strcmp (action, "-stats"))
    gStats (methods2, methodTag2,  methodTitle2, db, xml, template) ;
  else if (!strcmp (action, "-worm"))
    gWorm (db) ;
  else if (!strcmp (action, "-bests"))
    {
      /* methods2[] = { (int)'A', (int)'K', (int)'S', (int)'H',  (int)'G', (int)'J', (int)'N', (int)'E', 0 } ; */
      int iii ;
      for (iii = 0 ; methods2[iii] ; iii++)
	{
	  if (who && ! strchr (who, methods2[iii]))
	    methods2[iii]=1; 
	}
       if (1) gStatsCountBest (db, xml, methods2, methodTag2, methodTitle2, TRUE) ; /* int */
      if (1 && !xml) gStatsCountBest (db, xml, methods2, methodTag2, methodTitle2, FALSE) ; /* percents */
      if (1) gStatsAliBabaNewDefects (db, xml, methods2, methodTag2, methodTitle2) ; /* figure 6b, 6c */
    }
  else if (!strcmp (action, "-chroms"))
    gChromStats (db, xml, template) ;
  else if (!strcmp (action, "-libs"))
    gLibraryStats (db, xml) ;
  else if (!strcmp (action, "-chrom_quality"))
    gChromQualityStats (db, xml) ;
  else if (!strcmp (action, "-lib_quality"))
    gLibQualityStats (db, xml) ;
  else if (!strcmp (action, "-count_n"))
    gCountN (db, template) ;
  else if (!strcmp (action, "-introns"))
    {
      gStatsAliBabaIntrons (db, xml) ; /* figure 3b */
      if (1) gIntronDistrib (db, xml) ;
    }
  else if (!strcmp (action, "-6mers"))
    {
      gIntron6mers (db, xml) ;
    }
  else if (!strcmp (action, "-tg2pg"))
    gCreatePg (db) ;
  else if (!strcmp (action, "-export"))
    gExportAli (db, template, gold) ;
  else if (!strcmp (action, "-exportGoldAnnot"))
    gExportGoldAnnot (db, template, gold, 1) ;
  else if (!strcmp (action, "-exportGoldDna"))
    gExportGoldAnnot (db, template, gold, 2) ;
  else if (!strcmp (action, "-exportGoldPep"))
    gExportGoldAnnot (db, template, gold, 3) ;
  else usage () ;

  showStats (0, 0) ;
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


