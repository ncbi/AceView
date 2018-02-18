#include "../wac/ac.h"
#include "../wfiche/gtitle.h"
#include "../wfiche/biolog.h"
/*
#include "../wfiche/biolog.h"
#define vONLINEMAXURLSIZE 2048
#define vSTRMAXNAME 2048 

#define _00    "\0\0" 

#define TABLEBORDER "2"
*/
#define	TBB(_line,_col)		 array (Tbb,  (_col) + TbbDim * (_line), int)
#define	TPP(_vtxt, _line, _col)	 vtxtAt (_vtxt, TBB(_line,_col))
#define	ooo  "\1\1\1"

static void ficheGoldSquareTable (vTXT blkp, GMP *gmp, Array Tbb, vTXT bfr, int TbbDim, int nRows, int *cols )
{
  /*the output line up in columns anyway.
   */
  int w, x, y ;
  char *cp, white[256] ;
  int widths [1024] ;
  AC_HANDLE h = ac_new_handle () ;
  /* clean up the ooo separators */
  cp = vtxtPtr (bfr) ;
  while (*cp) { if (*cp == 1) *cp = 0 ; cp++ ; }
  cp = white ;
  for (cp = white, y = 0 ; y < 255 ; y++, cp++) *cp = ' ' ;
  /*
   * tables always begin on a new line
   */
  vtxtBreak (blkp) ;

  /*
   * count the max width of each.
   */
  for (x = 0; x < nRows ; x++ )
    {
      for ( y = 0; y < 1024 && cols[y] != -1; y++ )
	{
	  w = strlen(TPP(bfr, x, cols[y])) ;
	  if (widths[y] < w)
	    widths[y] = w;
	}
    }
  
  /*
int   * display the data
   */
  for ( x = 0; x < nRows; x++ )
    {
      if (gmp->style == 'x')
	{
	  if (x == 0) /* titles */
	    vtxtPrint (blkp, "<table border=1>\n<tr bgcolor=#afafff>\n") ;
	  else if (x % 2)
	    vtxtPrint (blkp, "<tr>\n") ;
	  else
	    vtxtPrint (blkp, "<tr bgcolor=#efefff>\n") ;
	}
      for ( y = 0; cols[y] != -1; y++ )
	{
	  cp = TPP (bfr, x, cols[y]) ;
	  if (gmp->style == 'x')
	    vtxtPrintf (blkp, "<td>%s</td>\n", cp) ;
	  else
	    {
	      w = strlen (cp);
	      if (y < 1024 && widths[y] < 25)
		vtxtPrint (blkp, white + 25 - (widths[y] - w)) ;
	      vtxtPrint (blkp, cp) ;
	      vtxtPrint (blkp, "\t") ; /* column separator */
	    }	  
	}
      if (gmp->style == 'x')
	vtxtPrint (blkp, "\n</tr>\n") ;
      else
	vtxtPrint (blkp, "\n" ) ;
    }
  if (gmp->style == 'x')
    vtxtPrint (blkp, "</table>\n") ;
  else
    vtxtPrint (blkp, "\n" ) ;

  messfree (h) ;
}

/***************************************************************************/

static void ficheGoldCloneTitleParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp, buf [64] ;
  const char *ccp ;
  AC_KEYSET ks ;
  AC_TABLE tt ;
  int ir, nn = 0, x1, x2, dx ;

  /* try to catch the vervion number */
  ccp = ac_name (gmp->est) ;
  tt = ac_tag_table (gmp->est, "Database", h) ;
  
  for (ir = 0 ; tt && ir < tt->rows ; ir++)
    if (! strcasecmp (ac_table_printable (tt, ir, 0, "toto"), "Genbank"))
      ccp = ac_table_printable (tt, ir, 1, ac_name (gmp->est)) ;

  vtxtPrintf (blkp, "GenBank mRNA %s", ccp) ;
  if (gmp->est)
    {
      ks = ac_objquery_keyset (gmp->est, ">from_gene ; gold && !minor", h) ;
      nn = ac_keyset_count (ks) ;
    }
  if (nn == 1)
    vtxtPrint (blkp, " has a single Gold alignment on the human genome") ;

  tt = gmp->tg ? ac_tag_table (gmp->tg, "IntMap", h) : 0 ;
  if (!tt)
    goto abort ;
  ccp = ac_table_printable (tt, 0, 0, 0) ;
  if (ccp)
    {
      vtxtDot (blkp) ;
      vtxtPrint (blkp, " It aligns on") ;
      strncpy (buf, ccp, 63) ;
      cp = strstr (buf, "|") ;
      if (cp)
	{
	  vtxtPrint (blkp, " a random contig of") ;
	  *cp = 0 ;
	}
      vtxtPrintf (blkp, " chromosome %s", buf) ;
      if ((x1 = ac_table_int (tt, 0, 1, 1)) && (x2 = ac_table_int (tt, 0, 2, 0)))
	{
	  dx = x1 < x2 ? x2 - x1 + 1 : x1 - x2 + 1 ;
	  vtxtPrintf (blkp, " over %d kb, from bp %d to %d", dx/1000, x1, x2) ;
	}
    }
     
  vtxtDot (blkp) ;
 
 abort:
  ac_free (h) ;
  return ;

} /* ficheGoldCloneTitleParagraphContent */

/***************************************************************************/

static void ficheGoldExonsParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn, ngt, ngc, nat, nex, ni ;

  nn = nex =  ac_tag_int (gmp->tg, "Nb_possible_exons", 0) ;
  vtxtPrintf (blkp, "The best alignment has %d exon%s", nn, nn > 1 ? "s" : "") ;
  nn = ngt = ac_tag_int (gmp->tg, "gt_ag", 0) ;
  if (nn)
    vtxtPrintf (blkp, ", %d intron%s %s gt_ag boundaries"
		, nn
		, nn > 1 ? "s" : ""
		, nn > 1 ? "have" : "has"
		) ;
  nn = ngc = ac_tag_int (gmp->tg, "gc_ag", 0) ;
  if (nn)
    vtxtPrintf (blkp, ", %d intron%s %s gt_ag boundaries"
		, nn
		, nn > 1 ? "s" : ""
		, nn > 1 ? "have" : "has"
		) ;

  nat = ngc = ac_tag_int (gmp->tg, "at_ac", 0) ;
  if (nn)
    vtxtPrintf (blkp, ", %d intron%s %s at_ac boundaries"
		, nn
		, nn > 1 ? "s" : ""
		, nn > 1 ? "have" : "has"
		) ;
  ni = nex - 1 ;
  nn = ni - ngt - ngc - nat ;
  if (nn)
    vtxtPrintf (blkp, ", %d intron%s %s non classic boundaries"
		, nn
		, nn > 1 ? "s" : ""
		, nn > 1 ? "have" : "has"
		) ;

  vtxtDot (blkp) ;

  ac_free (h) ;
  return ;
} /* ficheGoldExonsParagraphContent */

/***************************************************************************/

static void ficheGoldExonsTable (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  Array Tbb = arrayHandleCreate (1024, int, h) ;
  const char *ccp, *estNam ;
  AC_OBJ oCosmid = ac_tag_obj (gmp->tg, "Genomic_sequence", h) ;
  AC_TABLE gCovers = ac_tag_table (gmp->tg, "Covers", h) ;
  AC_TABLE Spl, Af ;
  int ir, jr,jt, x1, x2, a1, a2, b1, b2, err, ok, maxCols, TbbDim ;
  int jExon = 0, ga1, ga2 ;
  char *colNames[]={ " " 
		     ,"Length<br>& DNA"
		     ,"Coordinates<br>on mRNA"
		     ,"Coordinates<br>on genome"
		     ,"#Differences"
		     ,0} ; 
  enum {COL_INTREX = 0, COL_LENG, COL_X1X2, COL_A1A2, COL_DIFF} ;
  int cols[] = { COL_INTREX, COL_LENG, COL_X1X2, COL_A1A2, COL_DIFF, -1 } ;
 
  estNam = ac_name (gmp->est) ;

  Af = ac_tag_table (gmp->tg, "Assembled_from", h) ;
  Spl = ac_tag_table (gmp->tg, "Splicing", h) ;
  
  ga1 = ac_table_int (gCovers, 0, 2, 0) ;
  ga2 = ac_table_int (gCovers, 0, 3, 0) ;

   /* initialize TbbDim */
  for (maxCols=TbbDim=0 ; colNames[maxCols] ; TbbDim++, maxCols++) ;
   /* initialize the column titles */
  for (maxCols=0 ; colNames[maxCols] ; maxCols++)
    TBB(0, maxCols) = vtxtPrintf (bfr, "%s"ooo, colNames[maxCols]) ; 
  
  /* filling the table */
  /* loop on splicing */
  for (ir=jr=0, jt=1 ; ir < Spl->rows ; jt++, ir++)
    {
      b1 = ac_table_int (Spl, ir, 0, 0) ;
      b2 = ac_table_int (Spl, ir, 1, 0) ;
      a1 = a2 = x1 = x2 = err = 0 ;
      /* loop on assembled_from to find the corresponding entry */
      ok = 0 ;
      for (jr = 0 ; !ok && jr < Af->rows ; jr++)
	{
	  ccp = ac_table_printable (Af, jr, 2, "") ;
	  if (strcmp (estNam, ccp))
	    continue ;
	  a1 = ac_table_int (Af, jr, 0, 0) ;
	  if (a1 != b1)
	    continue ;
	  a2 = ac_table_int (Af, jr, 1, 0) ;
	  if (a2 != b2)
	    continue ;
	  x1 = ac_table_int (Af, jr, 3, 0) ;
	  x2 = ac_table_int (Af, jr, 4, 0) ;
	  err = ac_table_int (Af, jr, 5, 0) ;
	  ok = 1 ;
	}
      if (!ok)
	continue ;
      /* Exon/Intron type */
      {
	int iss ; 
	struct {char *lookFor, * whatToSay ; }xonTypes[]={
	  {"Exon", "Exon"}, 
	  {"Alternative_exon", "Alternative exon"}, 
	  {"Predicted_exon", "Predicted exon"}, 
	  {"Stolen_exon", "Inferred exon"}, 
	  {"Intron", "Intron"}, 
	  {"Predicted_intron", "Predicted intron"}, 
	  {"Alternative_intron", "Alternative intron"}, 
	  {"Stolen_intron", "Inferred intron"}, 
	  {"Gap", "Sequencing gap"}, 
	  {0, 0}} ; 
	const char *txt2, *txt = ac_table_tag (Spl, ir, 2, "") ; 

	for (iss=0 ; xonTypes[iss].lookFor ; iss++)
	  {
	    if (!strcasecmp (txt, xonTypes[iss].lookFor))break ; 
	  }
	/* exon exon-number*/
	if (strstr (txt, "xon"))
	  {
	    TBB (jt, COL_INTREX) = vtxtPrintf (bfr, "<b>%s</b>", xonTypes[iss].whatToSay) ; 
	    vtxtPrintf (bfr, " %s"ooo, ac_table_printable (Spl, ir, 3, "")) ; 
	  }
	/* Intron */
	else if (strstr (txt, "ntron"))
	  {
	    txt2 = ac_table_printable (Spl, ir, 3, "") ;
	    TBB (jt, COL_INTREX) = vtxtPrintf (bfr, "%s", xonTypes[iss].whatToSay) ; 
	    if (!strcasecmp (txt2, "Fuzzy"))vtxtPrintf (bfr, " Fuzzy"ooo) ; 
	    vtxtPrintf (bfr, " [%c%c-%c%c]"ooo, txt2[0], txt2[1], txt2[3], txt2[4]) ; 
	  }
	else
	  {
	    TBB (jt, COL_INTREX)=vtxtPrintf (bfr, "%s", xonTypes[iss].whatToSay) ; 
	    vtxtPrintf (bfr, ooo) ; 
	  }
      }

      /* length & coordinate on Gene */
      {  /* view -v fasta -c mrna -n mog-6 -xml -p 12 35 */
	char buf1[256], buf2[256] ;
	int u1, u2 ;
	int exonColor ;

	if (ga1 < ga2)
	  { u1 = ga1 + b1 - 1 ; u2 = ga1 + b2 - 1 ; }
	else
	  { u1 = ga1 - b1 + 1 ; u2 = ga1 - b2 + 1 ; }
	TBB (jt, COL_LENG) = vtxtPrintf (bfr, " ") ;

	if (strstr (ac_table_printable (Spl, ir, 2, ""), "xon"))
	  {
	    if (jExon++ % 2)
	      exonColor = 1 ;
	    else
	      exonColor = 2 ;
	  }
	else
	  exonColor = 0 ;
	sprintf (buf1, "DNA:%d:%d:%d", u1, u2, exonColor) ;
	sprintf (buf2, "%dbp", b2 - b1 + 1) ;
	gmpFakeObjLink (bfr, gmp, buf1, oCosmid, buf2) ;
      
	vtxtPrintf (bfr, ooo) ; 
      }
      
      TBB (jt, COL_A1A2) = vtxtPrintf (bfr, "%d to %d"ooo, a1, a2) ;
      TBB (jt, COL_X1X2) = vtxtPrintf (bfr, "%d to %d"ooo, x1, x2) ;
      TBB (jt, COL_DIFF) = vtxtPrintf (bfr, "%d"ooo, err) ;
      
  
    } /* Spl loop */

  vtxtBreak (blkp) ;
  vtxtPrintf (blkp, "Sequences of exons and introns are available by clicking on the lengths in column 3") ;
  vtxtBreak (blkp) ;
  ficheGoldSquareTable (blkp, gmp, Tbb, bfr, TbbDim, jt, cols) ;
  
  ac_free (h) ;
  return ;
} /* ficheGoldExonsTable */

/***************************************************************************/
/***************************************************************************/

void ficheGoldCloneParagraph (vTXT blkp, GMP *gmp) 
{
  vTXT bfr ; 
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheGoldCloneTitleParagraphContent (bfr, gmp) ;
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "Summary", "Summary for this clone") ; 
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtClear (bfr) ;
  
  ficheGoldExonsParagraphContent (bfr, gmp) ;
  ficheGoldExonsTable (bfr, gmp) ;
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "Splicing", "Splicing of the alibaba alignment") ; 
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheCdnaCloneGoldParagraph */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/


