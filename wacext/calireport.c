/*  File: calireport.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2011
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
 *  available from http://www.aceview.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

#include "ac.h"
#include "acedb.h"
/* #define ARRAY_CHECK */

typedef struct qStruct { 
  char *chapter ;
  char *title ;
  char RA ;     /* r: Run, a: Ali */
  char *tag ;    /* tag to be found in the Run or Ali */
  int col ;      /* column in behind tag */
  char type ; /* t: text, i: integer, f: float */
  char nonEmpty ;
} QL ;

typedef struct baStruct { 
  AC_HANDLE h ; 
  AC_DB db ;           /* ACEDB database containing the semantics of the experiment */
  ACEOUT ao ;
  BOOL gzo, html, RNA_seq ;
  const char *outFileName ;
  const char *run, *project, *title ;
  QL *ql ;
  Array cells ;
  Stack s ;
  KEYSET rowNames ;       /* offsetts in rowStack */
  Stack rowNameStack ;    /* names of the rows */
  DICT *target_classDict ;
  KEYSET runNonEmpty ;
  AC_TABLE runs ;
} BA ;

typedef struct cellStruct { 
  int irun, iql, /* run and query */
    txt, i, f ; /* txt as offset in ba->s, int, float, as needed for that cell */
} CL ;

/*************************************************************************************/
/*************************************************************************************/

static int cellOrder (const void *a, const void *b)
{
  const CL *up = (const CL *)a, *vp = (const CL *)b ;
  int n ;

  n = up->iql - vp->iql ; if (n) return n ;
  n = up->irun - vp->irun ; if (n) return n ;
  return 0 ;
} /* cellOrder */

/*************************************************************************************/

static void crInit (BA *ba)
{
   static QL ql[] = {
    { "MetaData", "Run or group descriptor (manually entered in the database)", 'r', "Title", 1, 'K', 0},
    { 0, "Molecule targetted", 'r', "Protocol", 1, 'T', 0}, 
    { 0, "Sequencing platform", 'r', "Platform", 1, 'T', 0},
    { 0, "Authors or laboratory", 'r', "Author", 1, 'K', 0},
    { 0, "Date of data file", 'r', "Date_received", 1, 'd', 0},
    { 0, "SRA or GEO accession numbers", 'r', "SRR", 1, 'K', 0}, 
    { 0, " list of runs", 'r', "Runs", 1, 'K', 0}, 
    {"Sample", "Sample id", 'r', "Sample", 1, 'k', 0}, 
    {0, "Tissue or cell type", 's', "Tissue", 1, 'K', 0},
    {0, "Details on sample", 's', "T", 1, 'T', 0},
    { "Sequence data	original reads", "Distinct sequences", 'a', "Raw_data", 1, 'f', 0} ,
    { 0, "Reads received", 'a', "Raw_data", 3, 'f', 0} ,
    { 0, "kb", 'a', "Raw_data", 5, 'f', 0} ,
    { 0, "Reads to align (minus low entropy/no insert)", 'a', "Accepted", 3, 'f', 0},
    { 0, "kb to align (minus low entropy/no insert)", 'a', "Accepted", 5, 'f', 0},
    {"Reads destiny", "Distinct sequences mapped", 'a', "h_Ali any", 2, 'f', 0},
    { 0, "kb mapped", 'a', "h_Ali any", 6, 'f', 0},
    { 0, 0, 0, 0, 0, 0, 0 }
   } ;
   ba->title = "Analysis summary tables, ENCODE compliant, draft" ;
   ba->ql = ql ; 
} /* crInit */

/*************************************************************************************/
/****************************** HTML utilities ***************************************/

static void crHtmlInit (vTXT txt, BA *ba)
{
  vtxtPrint (txt,
	     "<html>\n"
	     "<head>\n"
	     ) ;
  vtxtPrintf (txt, "  <title>AceView global alignment results for sequencing project %s</title>\n", ba->project) ;

  vtxtPrint (txt,
	     "<meta http-equiv=\"content-type\" content=\"text/html;charset=iso-8859-1\">\n"
	     "<meta name=\"author\" content=\"Danielle Thierry-Mieg and Jean Thierry-Mieg, NCBI/NLM/NIH, mieg@ncbi.nlm.nih.gov\">"
	     "<meta name=\"ROBOTS\" CONTENT=\"NOINDEX, NOFOLLOW\">\n"
	     "</head>\n\n"
	     "<body>\n "
	     ) ;
 } /* crHtmlOpen */


static void crHtmlClose (vTXT txt)
{
  vtxtPrintf (txt, "\n</body>\n</html>\n") ;
} /* crHtmlClose */

static void cellPrint (BA *ba, vTXT txt, const char *ccp)
{
  if (ba->html) vtxtPrint (txt, "\n  <td VALIGN=TOP>") ;
  else vtxtPrint (txt, "\t") ;
  vtxtPrint (txt, ccp) ;
  if (ba->html) vtxtPrint (txt, "</td>")  ;
} /* cellPrint */

/*************************************************************************************/
/****************************** Actual work ******************************************/
/*
select r,t from r in class "run" where r->project = "Test3", t in r->type
*/
static int crGetRuns (BA *ba)
{
  char *errors = 0 ;
  char *qq ;

  qq = messprintf ("select r,t from r in class \"run\" where r->project == \"%s\", t in r->type", ba->project) ;
  ba->runs = ac_bql_table (ba->db, qq, 0, "+2+1", &errors, ba->h) ;
  if (errors)
    messcrash ("cGetRuns failed: %s\n", errors) ;
  else
    fprintf (stderr, "Found %d runs in project %s\n", ba->runs->rows,  ba->project) ;
  return ba->runs->rows ;
} /* crGetRuns */

/*************************************************************************************/

static char *crGetCellData (BA *ba, AC_TABLE table, QL *ql)
{
  static vTXT txt = 0 ;
  int ir = 0, jr, n = 0, col =  ql->col > 0 ? ql->col - 1 : 0 ;
  const char *ccp ;
  char *cp, *ct ;
  float z ;

  ct = strchr (ql->tag, ' ') ;
  if (ct) ct++ ;
  if (!txt)
    txt = vtxtHandleCreate (ba->h) ;
  vtxtClear (txt) ;

  if (ct)
    for (ir = 0 ; ir < table->rows ; ir++)
      {
	ccp = ac_table_printable (table, ir, 0, 0) ;
	if (ccp && !strcasecmp (ccp, ct))
	  break ;
      }
  switch (ql->type)
    {
    case 'K':  /* loop on all identifiers */
      for (ir = 0 ; ir < table->rows ; ir++)
	{
	  ccp = ac_table_printable (table, ir, col, 0) ;
	  if (ccp)
	    vtxtPrintf (txt, "%s%s", n++ ? ", " : "", ccp) ;
	}
      break ;
    case 'T':  /* loop on all content */
      for (ir = 0 ; ir < table->rows ; ir++)
	for (jr = 0 ; jr < table->cols ; jr++)
	  {
	    ccp = ac_table_printable (table, ir, jr, 0) ;
	    if (ccp)
	      vtxtPrintf (txt, "%s%s", n++ ? ", " : "", ccp) ;
	  }
      break ;
    case 'f':
      z = ac_table_float (table, ir, col, 0) ;
      if (z > 100)
	vtxtPrintf (txt, "%.0f", z) ;
      else
	vtxtPrintf (txt, "%.2f", z) ;
      break ;
    case 'k':
    case 'i':
    case 'd':
    default:
      ccp = ac_table_printable (table, ir, col, 0) ;
      vtxtPrint (txt, ccp) ;
      break ;
    }
  for (cp = vtxtPtr(txt) ; cp && *cp ; cp++)
    if (*cp == '_') *cp = ' ' ;
  return vtxtPtr (txt) ;
}  /* crGetCellData */

/*************************************************************************************/

static void crGetRunData (BA *ba, int iRun)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Run = 0, Ali = 0, Sample = 0 ;
  AC_TABLE table = 0 ;
  QL *ql ;
  CL *cell ;
  int cellMax = arrayMax (ba->cells) ;
  char *oldTag = 0 ;
  int iql ;
  const char *ccp ;

  Run = ac_table_obj (ba->runs, iRun, 0, h) ;
  if (Run)
    {
      Ali = ac_tag_obj (Run, "Ali", h) ;
      Sample = ac_tag_obj (Run, "Sample", h) ;
    }
  if (Ali)
    for (iql = 0, ql = ba->ql ; ql->tag ; iql++, ql++)
      { 
	if (!oldTag || strcmp (oldTag, ql->tag))
	  { 
	    AC_OBJ obj = 0 ;
	    char buf[100], *cp ;

	    strncpy (buf, ql->tag, 99) ;
	    ac_free (table) ;
	    switch (ql->RA)
	      {
	      case 'r': obj = Run ; break ;
	      case 'a': obj = Ali ; break ;
	      case 's': obj = Sample ; break ;
	      }
	    cp = strchr (buf, ' ') ;
	    if (cp) *cp = 0 ;
	    table = obj ? ac_tag_table (obj, buf, h) : 0 ;
	  }
	oldTag = ql->tag ;
	cell = arrayp (ba->cells, cellMax++, CL) ;
	cell->irun = iRun ;
	cell->iql = iql ;
	ccp = table ? crGetCellData (ba, table, ql) : 0 ;
	if (ccp)
	  { 
	    cell->txt = stackMark (ba->s) ;
	    pushText (ba->s, ccp) ;
	    ql->nonEmpty = 1 ;
	    keySet (ba->runNonEmpty, iRun) = 1 ;
	  }
      }

  ac_free (h) ;
  return ;
} /* crGetData */

/*************************************************************************************/

static void crGetData (BA *ba)
{
  int irun ;
  
  for (irun = 0 ; irun < ba->runs->rows ; irun++)
    crGetRunData (ba, irun) ;

  arraySort (ba->cells, cellOrder) ;
  return ;
} /* crGetData */

/*************************************************************************************/
/*************************************************************************************/

static void crExportData (BA *ba, vTXT txt)
{
  int ii, jj, irun, oldirun, oldiql ;
  CL *cell ;
  int cellMax = arrayMax (ba->cells) ;
  const char *ccp = "Group or run" ;
  char *cp, buf[1000] ;

  if (ba->html) 
    vtxtPrintf (txt, "<table width=\"98%%\" border=\"TABLEBORDER\">\n<tr bgcolor=\"#d5d5ff\">\n  <td>Project %s</td><td>%s</td>\n", ba->project, ccp) ;
  else
    vtxtPrintf (txt, "\tProject %s\t%s", ba->project, ccp) ;
  for (irun = 0 ; irun < ba->runs->rows ; irun++)
    {
      if (! keySet (ba->runNonEmpty, irun))
	continue ;
      if (ba->html) vtxtPrint (txt, "\n  <td>") ;
      strncpy (buf, ac_table_printable (ba->runs, irun, 0, "-"), 999) ;
      for (cp = buf ; *cp ; cp++)
	if (*cp == '_') *cp = ' ' ;
      vtxtPrintf (txt, "\t%s", buf) ;
      if (ba->html) vtxtPrint (txt, "</td>") ;
    }
  if (ba->html) vtxtPrint (txt, "</tr>") ;

  for (ii = jj = 0, cell = arrayp (ba->cells, 0, CL), oldiql = oldirun = -1 ; ii < cellMax ; cell++, ii++)
    {
      if (! ba->ql[cell->iql].nonEmpty)
	continue ;
      if (! keySet (ba->runNonEmpty, cell->irun))
	continue ;
      if (cell->iql != oldiql)
	{
	  if (ba->html) 
	    {
	      if (jj++) vtxtPrint (txt, "</tr>\n") ;
	      vtxtPrintf (txt, "<tr bgcolor=\"%s\">\n", jj%2 ? "white" : "#efefff") ;
	    }
	  else 
	    vtxtPrint (txt, "\n") ;
	  cellPrint (ba, txt, ba->ql[cell->iql].chapter) ;
	  cellPrint (ba, txt, ba->ql[cell->iql].title) ;
	  oldirun = -1 ;
	}
      while (oldirun++ < cell->irun - 1)
	{
	  if (! keySet (ba->runNonEmpty, oldirun))
	    continue ;
	  if (ba->html) vtxtPrint (txt, "\n  <td></td>") ;
	  else vtxtPrint (txt, "\t") ;
	}
      cellPrint (ba, txt, cell->txt ? stackText (ba->s, cell->txt) : "-") ;
      oldirun = cell->irun ;
      oldiql = cell->iql ;
    }
  if (ba->html) vtxtPrint (txt, "\n</tr>\n") ;
  else 	vtxtPrint (txt, "\n") ;
  
  if (ba->html) 
    vtxtPrint (txt, "</table>\n") ;

  return ;
} /* crExportData */

/*************************************************************************************/
/*************************************************************************************/

static void crExport (BA *ba)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (ba->outFileName, ba->html ? ".htm" : ".txt", ba->gzo, h) ;
  vTXT txt = vtxtHandleCreate (h) ;

  if (ba->html)
    {
      vtxtMarkup (txt) ;    
      crHtmlInit (txt, ba) ;
    }

  vtxtPrintf (txt, "%s", timeShowNow ()) ;
  vtxtBreak (txt) ;

  crExportData (ba, txt) ;

  if (ba->html)
    crHtmlClose (txt) ;
  
  aceOut (ao, vtxtPtr (txt)) ;

  ac_free (ao) ;
  ac_free (h) ;
} /* crExport */

/*************************************************************************************/
/*************************************************************************************/

static void crReport (BA *ba)
{
  crGetRuns (ba) ;
  crGetData (ba) ;
  crExport (ba) ;

  return ;
} /* crReport */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// bestali.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   to analyse the output of the AceView aligner:\n"
	    "// Database\n"
	    "//   -db ACEDB : acedb database holding the metaData and results of the experiments\n"
	    "// Output\n"
	    "//   -title \'text\' :\n"
	    "//      use this title in the output\n"
	    "//   -o out_file_name [-gzo] : default stdout\n"
	    "//      redirect the output, each option adds its own suffix to the file name\n"
	    "//      if -gzo is specified, invokes gzip and add a .gz suffix to the file name\n"
	    "//   -html\n"
	    "//      export the table in html format\n"
	    "//\n"      
	    "// Filters\n"
	    "//   -project project_name\n"
	    "//       Export a table for all groups and runs of the given projectt\n"
	    "//   -run run\n"
	    "//       Export the report for a single run\n"
            "//   -strategy [ Exome | Genome | RNA_seq ]\n"
	    "//        optional, should match the choice in clipalign\n"
	    "// Help\n"
	    "//   -help : export this on line help\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' in the input files are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  BA ba ;
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp, *dbName = 0 ;
  char *errors = 0 ;

  memset (&ba, 0, sizeof (BA)) ;
  ba.h = h ;

  if (argc == 1)
    usage (0) ;

  /* optional arguments */

  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "-o", &(ba.outFileName)) ;
  getCmdLineOption (&argc, argv, "-run", &(ba.run)) ;
  getCmdLineOption (&argc, argv, "-title", &(ba.title)) ;
  getCmdLineOption (&argc, argv, "-project", &(ba.project)) ;
  ba.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;
  ba.html = getCmdLineOption (&argc, argv, "-html", 0) ;
  if (getCmdLineOption (&argc, argv, "-strategy", &ccp))
    {
      if (! strcasecmp (ccp, "Exome") ||
	  ! strcasecmp (ccp, "Genome")
	  )
	{
	  ba.RNA_seq = FALSE ;
	}
      else if (! strcasecmp (ccp, "RNA_seq"))
	{
	  ba.RNA_seq = TRUE ;
	}
    }
  ba.target_classDict = dictHandleCreate (10000, ba.h) ;
    {
        /* synchronize with bin/target2target_class.txt and with baExportAliProfile */
      const char **cl ;
      const char *classes[] = { "0_SpikeIn", "A_mito", "B_rrna", "C_chloro"
			       , "ET_av", "FT_cloud", "KT_RefSeq", "DT_seqc", "LT_UCSC", "MT_EBI", "NT_HINV", "OT_rnaGene", "PT_tRNA"
			       , "S_est", "U_introns", "W_Line", "Z_genome", "z_gdecoy"
			       , 0 } ;
      for (cl = classes ; *cl ; cl++)
	dictAdd (ba.target_classDict, *cl, 0) ;
    }

  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  if (argc > 1)
    {
      fprintf (stderr, "Unknown argument %s, try -help", argv[argc-1]) ;
      exit (1) ;
    }

  /* actions */

  fprintf (stderr, "Start %s\n", timeShowNow ()) ;

  if (! ba.project)
    messcrash ("Please specify a project") ;
  if (dbName)
    {
      ba.db = ac_open_db (dbName, &errors);
      if (! ba.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }
  else
    messcrash ("Please specify a database, it should be the MetaDB of a CALI project") ;

  ba.cells = arrayHandleCreate (10000, CL, h) ;
  ba.s = stackHandleCreate (1000000, h) ; pushText (ba.s, "toto") ;
  stackTextOnly (ba.s) ;
  ba.rowNameStack = stackHandleCreate (3000, h) ;
  stackTextOnly (ba.rowNameStack) ;
  ba.rowNames = keySetHandleCreate (h) ;
  ba.runNonEmpty = keySetHandleCreate (h) ;

  crInit (&ba) ;
  crReport (&ba) ;
     
  /* clean up */
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  
  ac_free (ba.db) ;
  fprintf (stderr, "Done %s\n", timeShowNow ()) ;
  ac_free (h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

