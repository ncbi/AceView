/*  File: taxtree.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
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

#include "ac.h"
#include "keyset.h"
#include "dna.h"
#include "mytime.h"
#include "freeout.h"
#include "dict.h"
#include "mytime.h"
#include "waceview/cgi.h"
#include "aceio.h"

extern KEY keyGetKey (KEY key, KEY tag) ;  
typedef struct taxTreeStruct { AC_DB db ; AC_HANDLE h ; int depth, rr ; BOOL gb, tt, cluster ; } TAXTREE ;
typedef struct gbCmdStruct { int level, addTo ; const char *tag, *div, *qq ; } GBCMD ;
GBCMD allGbCmds [] =
  {
    { 1, 1, "Gb", "nuccore", "NOT srcdb_refseq[PROP]" } ,
      { 2, 2, "Gb_rna", "nuccore", "AND biomol_rna[PROP] NOT srcdb_refseq[PROP]" } ,
      { 2, 3, "Gb_dna", "nuccore", "NOT biomol_rna[PROP] NOT srcdb_refseq[PROP]" } ,
    { 1, 0, "Refseq", "nuccore", "srcdb_refseq[PROP]" } , 
      { 2, 0, "Refseq_rna", "nuccore", "AND biomol_rna[PROP] AND srcdb_refseq[PROP]" } ,
      { 2, 0, "Refseq_dna", "nuccore", "NOT biomol_rna[PROP] AND srcdb_refseq[PROP]" } ,
    { 1, 1, "dbEST", "nucest", 0} ,
      { 2, 2, "dbEST_rna", "nucest", "AND biomol_rna[PROP]" } ,
      { 2, 3, "dbEST_dna", "nucest", "NOT biomol_rna[PROP]" } ,
    { 1, 1, "GSS", "nucgss", 0 } ,
      { 2, 2, "GSS_rna", "nucgss", "AND biomol_rna[PROP]" } ,
      { 2, 3, "GSS_dna", "nucgss", "NOT biomol_rna[PROP]" } ,
    { 0, 0, 0 } ,
  } ;
/*
All_all = gb + NOTrefseq + dbEST + GSS + TraceDb_no_accession + TraceDb_454_no_accession
All_rna = idem RNA
All_dna = idem DNA
*/
typedef struct ttCmdStruct { int level, addTo ; const char *tag ;  char type ; int hasAccession ; int is454 ; } TTCMD ;
TTCMD allTtCmds [] =
  {
    { 1, 0, "TraceDb",                      'A',  0,  0} ,
      { 2, 0, "TraceDb_with_accession",     'A',  1, -1} ,
      { 2, 1, "TraceDb_no_accession",       'A', -1, -1} ,
      { 2, 0, "TraceDb_454_with_accession", 'A',  1,  1} ,
      { 2, 1, "TraceDb_454_no_accession",   'A', -1,  1} ,
    { 2, 2, "TraceDb_rna_no_accession",     'R', -1, -1} ,
    { 2, 3, "TraceDb_dna_no_accession",     'D', -1, -1} ,
    { 2, 1, "TraceDb_rna_454_no_accession", 'R', -1,  1} ,
    { 2, 1, "TraceDb_dna_454_no_accession", 'D', -1,  1} ,

    { 0, 0, 0 , 0} ,
  } ;

/*************************************************************************************/
/*************************************************************************************/

static AC_DB taxtreeInit (const char *dbName)
{
  AC_DB db = 0 ;
  const char *errMessage = 0 ;

  if (dbName)
    {
      if (!(db = ac_open_db (dbName, &errMessage)))
	{
	  fprintf (stderr, "// Cannot open database, sorry:: %s\n",
		   errMessage ? errMessage : "unkown reason"
		   ) ;
	  exit (1) ;
	}
    }
  else
    {
      fprintf (stderr, "missing parameter -db, run without argument to get the on line help") ;
      exit (1) ;
    }

  return db ;
} /* taxtreeInit */

/*************************************************************************************/
/*****************************  actual work ******************************************/
/*************************************************************************************/
/* genbank queries */
static int taxtreeGbGet (TAXTREE *pp, AC_OBJ taxId, const char *latin, GBCMD *cmd)
{
  AC_HANDLE h = ac_new_handle ()  ;
  const char *ccp, *ccq ;
  char *cp, *cr ;
  int ii, nn = 0 ;
  vTXT txt1 = vtxtHandleCreate (h), webPage = 0 ;
  vTXT txt2 = vtxtHandleCreate (h) ;

  vtxtPrint (txt1, "https://www.ncbi.nlm.nih.gov/sites/entrez?") ;
  cp = url_encode(latin) ;
  vtxtPrintf (txt1, "db=%s&cmd=search&term=%s[Organism]", cmd->div, cp) ; 
  free (cp) ;
  if (cmd->qq)
    {
      vtxtPrintf (txt2, " %s", cmd->qq) ;
      cp = url_encode(vtxtPtr (txt2)) ;
      vtxtPrint (txt1, cp) ;
      free (cp) ;
    }
  
  if (0) freeOutf ("// %s\n", vtxtPtr (txt1)) ; 
  webPage = vtxtHandleGetURL (vtxtPtr (txt1), 0, h) ;
  
  if (webPage && 
      (ccp = strstr (vtxtPtr(webPage), "Items")) &&
      (ccq = strstr (ccp, "of")) &&
      (cr = strstr (ccq+3, "<"))
      )
    {
      *cr = 0 ;
      if (sscanf (ccq+3, "%d", &ii) == 1 && ii >= 0)
	nn = ii ;
    }
  
  ac_free (h) ;
  return nn ;
} /* taxtreeGbGet  */

/*************************************************************************************/
/* trace queries */
static int taxtreeTtGet (TAXTREE *pp, AC_OBJ taxId, const char *latin, TTCMD *cmd)
{
  AC_HANDLE h = ac_new_handle ()  ;
  int ii, nn = 0 ;
  vTXT txt2 = vtxtHandleCreate (h) ;
  ACEIN ai ;

  char *isRNA = " and SOURCE_TYPE != 'SYNTHETIC' and  (SOURCE_TYPE = 'NON GENOMIC' or TRACE_TYPE_CODE='EST' or TRACE_TYPE_CODE='RT-PCR' or STRATEGY='EST' or STRATEGY='cDNA'  or STRATEGY='CCS' or STRATEGY='RT-PCR' or STRATEGY= 'MODEL VERIFY')" ;
  char *isDNA = " and SOURCE_TYPE != 'SYNTHETIC' and SOURCE_TYPE='GENOMIC' and TRACE_TYPE_CODE != 'RT-PCR' and TRACE_TYPE_CODE != 'EST'  and  (STRATEGY=NULL or (STRATEGY!='EST' and STRATEGY!='cDNA'  and STRATEGY!='CCS' and STRATEGY!='RT-PCR' and STRATEGY !='MODEL VERIFY'))" ;
  
  vtxtPrintf (txt2, "query_tracedb \"query count species_code='%s' ", latin) ;

  if (cmd->type == 'D') /* DNA/danielle would prefer gemomic, the problem is with RNA-viruses */
    {  vtxtPrint (txt2, isDNA) ; }
  else if (cmd->type == 'R') /* RNA */
    {  vtxtPrint (txt2, isRNA) ; }

  if (cmd->hasAccession == -1) 
    vtxtPrint (txt2, " and ACCESSION = NULL") ;
  else if (cmd->hasAccession == 1) 
    vtxtPrint (txt2, " and ACCESSION != NULL") ;

  if (cmd->is454 == 1)
    vtxtPrint (txt2, " and TRACE_TYPE_CODE = '454' ") ;
  else if (cmd->is454 == -1)
    vtxtPrint (txt2, " and (TRACE_TYPE_CODE = NULL or TRACE_TYPE_CODE != '454') ") ;

  vtxtPrint (txt2, " and LOAD_DATE > 'Jan 01 1990'\"") ;

  if (0)   fprintf (stderr, "%s\n", vtxtPtr (txt2)) ;
  ai = aceInCreateFromPipe (vtxtPtr (txt2), "r", 0, h) ;
  if (aceInCard(ai) && aceInInt (ai, &ii))
    nn = ii ;

  ac_free (h) ;
  return nn ;
} /* taxtreeTtGet  */

/*************************************************************************************/

static int taxtreeGbRun (TAXTREE *pp)
{
  AC_HANDLE h = 0 ;
  int dd, nSp = 0, nn, dateDiff ;
  AC_KEYSET ks = 0 ; /* current list of non empty parents */
  AC_ITER iter ;     /* iterator on their children */
  AC_OBJ taxId = 0 ; /* active taxId, if non empty, add it to future ks */
  GBCMD *cmd ;
  const char *latin ;
  mytime_t date1, date2 ;

  date1 = timeParse ("2008-04-01") ;
  for (dd = 0 ; dd <= pp->depth ; ac_free(h), dd++)
    {
      h = ac_new_handle () ;
      if (dd == 0)
	iter = ac_dbquery_iter (pp->db, "find taxid !parent ; latin == \"root\" ", pp->h) ;
      else
	iter = ac_ksquery_iter (ks, ">child", pp->h) ;
	
      ac_free (ks) ;
      ks = ac_new_keyset (pp->db, pp->h) ;
      while (ac_free (taxId), taxId = ac_iter_obj (iter))
	{
	  cmd = allGbCmds ;
	  latin = strnew (ac_tag_printable (taxId, "Latin", "xx"), h) ;
	  date2 = ac_tag_date (taxId, "Gb_date", 0) ;
	  if (date2 && 
	      timeDiffDays (date1, date2, &dateDiff) &&
	      dateDiff > 0
	      )
	    {
	      nn = ac_tag_int (taxId, "Gb_all", 0) ;
	      if (nn)
		ac_keyset_add (ks, taxId) ; /* usefull at depth */	
	      continue ;
	    }

	  nn = taxtreeGbGet (pp, taxId, latin, cmd) ;
	  freeOutf ("TaxId %s // %s\nGb_date \"%s\"\n-D Gb_all\n", ac_name (taxId), latin, timeShowNow()) ;  
	  for (cmd = allGbCmds + 1 ; cmd->level > 1 ; cmd++)
	    freeOutf ("-D %s\n", cmd->tag) ;
	  if (nn)
	    {
	      ac_keyset_add (ks, taxId) ; /* usefull at depth */
	      nSp++ ;
	      freeOutf ("Gb_all %d\n", nn) ;
	      for (cmd = allGbCmds + 1 ; cmd->level ; cmd++)
		if ((nn = taxtreeGbGet (pp, taxId, latin, cmd)))
		  freeOutf ("%s %d\n", cmd->tag, nn) ;  	
	    }
	  freeOutf ("\n") ;
	  fflush (stdout) ;
	}
    }

  ac_free(h) ;
  return nSp ;
} /* taxtreeGbRun */

/*************************************************************************************/

static int taxtreeTtRun (TAXTREE *pp)
{
  AC_HANDLE h = 0 ;
  int nSp = 0, nn, dateDiff ;
  AC_ITER iter ;     /* iterator on their children */
  AC_OBJ taxId = 0 ; /* active taxId, if non empty, add it to future ks */
  TTCMD *cmd ;       /* iterate on all the different queries */
  BOOL debug = FALSE ;
  const char *latin ;
  mytime_t date1, date2 ;

  date1 = timeParse ("2008-12-31") ;

  h = ac_new_handle () ;
  if (debug)
    iter = ac_dbquery_iter (pp->db, "find taxid  IS 91422 OR IS 9258", pp->h) ;
  else
    iter = ac_dbquery_iter (pp->db, "find taxid  TraceDb", pp->h) ;

  while (ac_free (taxId), taxId = ac_iter_obj (iter))
    {
      cmd = allTtCmds ;
      latin = strnew (ac_tag_printable (taxId, "Latin", "xx"), h) ;
      date2 = ac_tag_date (taxId, "Trace_date", 0) ;
      if (! debug && date2 && 
	  timeDiffDays (date1, date2, &dateDiff) &&
	  dateDiff > 0
	  )
	{
	  continue ;
	}

      nn = taxtreeTtGet (pp, taxId, latin, cmd) ;
      freeOutf ("TaxId %s // %s\nTrace_date \"%s\"\n-D TraceDb_all\n", ac_name (taxId), latin, timeShowNow()) ;  
      for (cmd = allTtCmds + 1 ; cmd->level > 1 ; cmd++)
	freeOutf ("-D %s\n", cmd->tag) ;
      if (nn)
	{
	  nSp++ ;
	  freeOutf ("TraceDb %d\n", nn) ;
	  for (cmd = allTtCmds + 1 ; cmd->level ; cmd++)
	    if ((nn = taxtreeTtGet (pp, taxId, latin, cmd)))
	      {
		freeOutf ("%s %d", cmd->tag, nn) ;  	
		if (cmd->is454)
		  freeOutf (" %g \"kb approximately\"", 150.0*nn/1000.0) ;  	
		else 
		  freeOutf (" %g \"kb approximately\"", 700.0*nn/1000.0) ;  
		freeOutf ("\n") ;
	      }
	}
      freeOutf ("\n") ;
      fflush (stdout) ;
    }

  ac_free(h) ;
  return nSp ;
} /* taxtreeTtRun */

/*************************************************************************************/
typedef struct segStruct { KEY key, parent ; 
  long int tt[20], ctt[20]
         , gb[20], ngb[20], cgb[20], cngb[20];
  float ntt[20], cntt[20] ; } SEG ;
static int taxtreeCluster (TAXTREE *pp)
{
  AC_HANDLE h = 0 ;
  int i, ii, ir, jj, nseg = 0, nseg0 ;
  AC_ITER iter ;     /* iterator on their children */
  AC_TABLE tbl ;
  AC_OBJ taxId = 0 ; /* active taxId, if non empty, add it to future ks */
  TTCMD *ttCmd ;       /* iterate on all the different queries */
  GBCMD *gbCmd ;       /* iterate on all the different queries */
  SEG *seg, *segp, *seg0 ;
  Associator ass = assHandleCreate (h) ;
  Array segs = arrayHandleCreate (500000, SEG, h) ;
  KEY parent ;
  void *vp ;
  KEY str2tag (const char *tagname) ;
  const char *myTag ;

  /* gather all the taxid that have some sequence data */
  h = ac_new_handle () ;
  iter = ac_dbquery_iter (pp->db, "find taxid Sequence", pp->h) ;
  while (ac_free (taxId), taxId = ac_iter_obj (iter))
    {
      seg = arrayp (segs, nseg++, SEG) ;
      seg->key = ac_obj_key (taxId) ;
      i = nseg - 1 ;
      assInsert (ass, assVoid(seg->key), assVoid (i)) ;
      seg->parent = ac_tag_key (taxId, "Parent", 0) ;
      if (seg->parent == seg->key)
	seg->parent = 0 ;
      tbl = ac_tag_table (taxId, "Sequence", 0) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  myTag = ac_table_printable (tbl, ir, 1, "x") ;
	  for (ii = 0, ttCmd = allTtCmds ; ttCmd->level ; ii++, ttCmd++)
	    if (! strcasecmp (ttCmd->tag, myTag))
	      {
		seg->ctt[ii] = seg->tt[ii] = ac_table_int (tbl, ir, 2, 0) ;
		seg->cntt[ii] = seg->ntt[ii] = ac_table_float (tbl, ir, 3, 0) ;
		break ;
	      }
	  for (ii = 0, gbCmd = allGbCmds ; gbCmd->level ; ii++, gbCmd++)
	    if (! strcasecmp (gbCmd->tag, myTag))
	      {
		seg->cgb[ii] = seg->gb[ii] = ac_table_int (tbl, ir, 2, 0) ;
		seg->cngb[ii] = seg->ngb[ii] = ac_table_int (tbl, ir, 3, 0) ;
		break ;
	      }
	}
      ac_free (tbl) ;
    }
  ac_free (iter) ;
  /* recurse up from those seg that do have basic data */
  nseg0 = nseg ;
  for (ii = 0 ; ii < nseg0 ; ii++)
    {
      seg = arrayp (segs, ii, SEG) ;
      if (0) fprintf(stderr, "key=%s\n", ac_key_name(seg->key)) ;
      jj = ii ; /* jj is the index of seg */
      if (seg->parent == seg->key)
	seg->parent = 0 ;
      while ((parent = seg->parent))
	{ 
	  if (0) fprintf(stderr, "  parent=%s\n", ac_key_name(seg->parent)) ;
	  if (!assFind (ass, assVoid(parent), &vp))
	    {
	      assInsert (ass, assVoid(parent), assVoid(nseg)) ;
	      segp = arrayp (segs, nseg++, SEG) ;
	      segp->key = parent ;
	      segp->parent = keyGetKey (parent, str2tag("Parent")) ;
	      if (segp->parent == segp->key)
		segp->parent = 0 ;
	      if (0) fprintf(stderr, "    grand-parent=%s\n", ac_key_name(segp->parent)) ;
	      seg = arrayp (segs, jj, SEG) ; /* may have been reallocated */
	      jj = nseg - 1 ; /* jj becomes the index of segp */
	    }
	  else
	    {
	      jj = assInt (vp) ; /* jj becomes the index of segp */
	      segp = arrp (segs, jj, SEG) ;
	      if (0) fprintf(stderr, "    grand-parent=%s\n", ac_key_name(segp->parent)) ;
	    }
	  seg0 = arrayp (segs, ii, SEG) ;
	  for (i = 0, ttCmd = allTtCmds ; ttCmd->level ; i++, ttCmd++)
	    {
	      segp->ctt[i] += seg0->tt[i] ;
	      segp->cntt[i] += seg0->ntt[i] ;
	    }
	  for (i = 0, gbCmd = allGbCmds ; gbCmd->level ; i++, gbCmd++)
	    {
	      segp->cgb[i] += seg0->gb[i] ;
	      segp->cngb[i] += seg0->ngb[i] ;
	    }
	  seg = segp ;  /* so jj is again the index of seg */
	}  
    }
  /* export */
  for (ii = 0 ; ii < nseg ; ii++)
    {
      seg = arrayp (segs, ii, SEG) ;
      freeOutf ("TaxId \"%s\" // ii=%d\n", ac_key_name (seg->key), ii) ;
      /* remove previous Cumulated data */
      for (i = 0, ttCmd = allTtCmds ; ttCmd->level ; i++, ttCmd++)
	freeOutf ("-D C_%s\n", ttCmd->tag) ;
      for (i = 0, gbCmd = allGbCmds ; gbCmd->level ; i++, gbCmd++)
	freeOutf ("-D C_%s\n", gbCmd->tag) ;
      /* export the correct values */
      for (i = 0, ttCmd = allTtCmds ; ttCmd->level ; i++, ttCmd++)
	if (seg->ctt[i])
	  freeOutf ("C_%s %ld\t%g \"kb approximately\"\n", ttCmd->tag, seg->ctt[i], seg->cntt[i]) ;
      for (i = 0, gbCmd = allGbCmds ; gbCmd->level ; i++, gbCmd++)
	if (seg->cgb[i])
	  freeOutf ("C_%s %ld\t%ld kb\n", gbCmd->tag, seg->cgb[i], seg->cngb[i]) ;

      {
	/* export the totals */
	float nt1, cnt1, nt2, nt3, cnt2, cnt3 ;
	long int t1, t2, t3 ;
	long int ct1, ct2, ct3 ;

	t1 = t2 = t3 = nt1 = nt2 = nt3 = 0 ;
	ct1 = ct2 = ct3 = cnt1 = cnt2 = cnt3 = 0 ;
	for (i = 0, ttCmd = allTtCmds ; ttCmd->level ; i++, ttCmd++)
	  {
	    if (ttCmd->addTo == 1) 
	      {
		t1 += seg->tt[i] ; nt1 += seg->ntt[i] ;
		ct1 += seg->ctt[i] ; cnt1 += seg->ntt[i] ;
	      }
	    if (ttCmd->addTo == 2) 
	      {
		t2 += seg->tt[i] ; nt2 +=  seg->ntt[i] ;
		ct2 += seg->ctt[i] ; cnt2 += seg->ntt[i] ;
	      }
	    if (ttCmd->addTo == 3) 
	      {
		t3 += seg->tt[i] ; nt3 +=  seg->ntt[i] ;
		ct3 += seg->ctt[i] ; cnt3 += seg->ntt[i] ;
	      }
	  }
	for (i = 0, gbCmd = allGbCmds ; gbCmd->level ; i++, gbCmd++)
	  {
	    if (gbCmd->addTo == 1) 
	      {
		t1 += seg->gb[i] ; nt1 += seg->ngb[i] ;
		ct1 += seg->cgb[i] ; cnt1 += seg->cngb[i] ;
	      }
	    if (gbCmd->addTo == 2) 
	      {
		t2 += seg->gb[i] ; nt2 += seg->ngb[i] ;
		ct2 += seg->cgb[i] ; cnt2 += seg->cngb[i] ;
	      }
	    if (gbCmd->addTo == 3) 
	      {
		t3 += seg->gb[i] ; nt3 += seg->ngb[i] ;
		ct3 += seg->cgb[i] ; cnt3 += seg->cngb[i] ;
	      }
	  }
	freeOutf ("All %ld %g  \"kb approximately\"\n", t1, nt1) ;
	freeOutf ("All_rna %ld %g  \"kb approximately\"\n", t2, nt2) ;
	freeOutf ("All_dna %ld %g  \"kb approximately\"\n", t3, nt3) ;

	freeOutf ("C_all %ld %g  \"kb approximately\"\n", ct1, cnt1) ;
	freeOutf ("C_all_rna %ld %g  \"kb approximately\"\n", ct2, cnt2) ;
	freeOutf ("C_all_dna %ld %g  \"kb approximately\"\n", ct3, cnt3) ;
      }
      freeOutf ("\n") ;
    }

  ac_free(h) ;
  return 0 ;
} /* taxtreeCluster */

/*************************************************************************************/

static void taxTreeDistantRelatives (AC_DB db, int NN)
{
  /*  // cathy please come over for ginger snaps */
  char *levels[] =
    {
      "superkingdom",
      "kingdom",
      "subkingdom",
      
      "superphylum",
      "phylum",
      "subphylum",
      
      "superclass",
      "class",
      "subclass",
      "infraclass",
      
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "parvorder",
      
      "superfamily",
      "family",
      "subfamily",
      
      "tribe",
      "subtribe",
      
      "genus",
      "subgenus",
      
      "species_group",
      "species_subgroup",
      "species",
      "subspecies",
      
      
      "varietas",
      
      "forma",
      
      0
    } ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl ;
  int ir, ii, n1, n2 ;
  vTXT qq = vtxtHandleCreate (h) ;
  const char *errors = 0 ;
  const char *ccp ;
  
  vtxtPrintf (qq, "Colonne 1\n"
	      " Class taxid\n"
	      " \n"
	      " Colonne 2\n"
	      " From 1\n"
	      " Tag level\n"
	      " Mandatory\n"
	      " Integer\n"
	      " \n"
	      " Colonne 3\n"
	      " From 1\n"
	      " Tag child\n"
	      " Mandatory\n"
	      " Class taxid\n"
	      " \n"
	      " Colonne 4\n"
	      " From 3\n"
	      " Tag level\n"
	      " Mandatory\n"
	      " Integer\n"
	      " \n"
	      " Colonne 5\n" 
	      " From 1\n" 
	      " Tag latin\n"
	      " Class Text\n"
	      " \n"
	      " Colonne 6\n" 
	      " From 3\n" 
	      " Tag latin\n"
	      " Class Text\n"
	      " \n"
	      ) ;

  tbl = ac_tablemaker_table (db, vtxtPtr (qq), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;

  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      ccp = ac_table_printable (tbl, ir, 1, 0) ;
      for (n1 = -1, ii = 0 ; ccp && levels[ii] ; ii++)
	if (! strcasecmp (levels[ii], ccp)) { n1 = ii ; break ; }
      ccp = ac_table_printable (tbl, ir, 3, 0) ;
      for (n2 = -1, ii = 0 ; ccp && levels[ii] ; ii++)
	if (! strcasecmp (levels[ii], ccp)) { n2 = ii ; break ; }
      if (n1 >= 0 && n2 >= 0 && n2 - n1 >= NN)
	{
	  freeOutf ("%d\t", n2 - n1) ;
	  freeOutf ("%s\t", ac_table_printable (tbl, ir, 1, "?")) ;
	  freeOutf ("%s\t", ac_table_printable (tbl, ir, 3, "?")) ;
	  freeOutf ("%s\t", ac_table_printable (tbl, ir, 0, "?")) ;
	  freeOutf ("%s\t", ac_table_printable (tbl, ir, 2, "?")) ;
	  freeOutf ("%s\t", ac_table_printable (tbl, ir, 4, "?")) ;
	  freeOutf ("%s\n", ac_table_printable (tbl, ir, 5, "?")) ;
	}
    }
  ac_free (h) ;
} /* taxTreeDistantRelatives */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr, "// Usage: taxtree -db . -... \n") ;
  fprintf (stderr, "// Example:  taxtree -db . -depth 3\n") ;  
  fprintf (stderr, "// -db ACEDB : name of the ACEDB taxtree database\n") ;
  fprintf (stderr, "// -o fileName : redirects stdout\n") ;
  fprintf (stderr, "// -depth number : (-gb case) maximal depth of eploration of the taxtree, 0->whole tree\n") ;
  fprintf (stderr, "//   -gb : explore GenBank unbelievably slow when depth > 3: days and weeks !\n") ;
  fprintf (stderr, "//   -tt : query tracedb on all taxid with tag tracedb_all (with date limit)\n") ;
  fprintf (stderr, "//   -cluster : cluster the sequence counts upwards the taxtree into class kindom etc\n") ;
  fprintf (stderr, "//   -distant number : export list of distant relatives\n") ;
  fprintf (stderr, "//   -?? : do ??\n") ;
  fprintf (stderr, "// You said: %s\n", commandBuf) ;
  
  fprintf (stderr, "// ########## ERROR: Unknown argument ") ;
  for (i = 1 ; i < argc ; i++)
    fprintf (stderr, "%s ", argv[i]) ;
  fprintf (stderr, "\n") ;
  
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  char *cp ;
  const char *ccp ;
  const char *outFileName = 0, *dbName = 0 ;
  int ix, nngb=0, nntt=0, outlevel = 0 ;
  TAXTREE taxtree ;
  AC_HANDLE h = 0 ;
  char commandBuf [1000] ;

  freeinit () ; 
  fprintf (stderr, "// %s start\n", timeShowNow()) ;

  h = ac_new_handle () ;
  memset (&taxtree, 0, sizeof (TAXTREE)) ;
  taxtree.h = h ;

  for (ix = 0, cp = commandBuf ; cp < commandBuf + 900 && ix < argc ; cp += strlen (cp), ix++)
    sprintf(cp, "%s ", argv[ix]) ;
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-o", &outFileName) ;
  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  taxtree.gb = getCmdLineOption (&argc, argv, "-gb", 0) ;
  taxtree.tt = getCmdLineOption (&argc, argv, "-tt", 0) ;
  taxtree.cluster = getCmdLineOption (&argc, argv, "-cluster", 0) ;
  if (getCmdLineOption (&argc, argv, "-distant", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >= 0)   
	taxtree.rr = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }
  if (!dbName ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, argc, argv) ;

  if (getCmdLineOption (&argc, argv, "-depth", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >= 0)   
	taxtree.depth = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }
  if (argc != 1)
    {
      usage (commandBuf, argc, argv) ;
    }

  if (outFileName)
    {
      f = filopen (outFileName, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  taxtree.db = taxtreeInit(dbName) ;
  if (taxtree.gb) 
    nngb = taxtreeGbRun (&taxtree) ;
  if (taxtree.tt) 
    nntt = taxtreeTtRun (&taxtree) ;
  if (taxtree.cluster)
    nntt = taxtreeCluster(&taxtree) ;
  if (taxtree.rr)
    taxTreeDistantRelatives (taxtree.db, taxtree.rr) ;
  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  /*
  if (taxtree.targetFile)
    filclose (taxtree.targetFile) ;
  */

  fprintf (stderr, "// %s done: studied %d Gb species, %d Trace species\n", timeShowNow(), nngb, nntt) ;

  ac_free (taxtree.h) ; /* will close taxtree.db */

  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

