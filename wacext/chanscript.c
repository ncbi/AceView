/*  File:  chanscript.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2017
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
 *
 * Purpose of this code
 *  If a script must be run on a large set of parameters, one way is to submit
 *  each instance, say 3000 cases, to the farm, waiting in the queue 
 *  to get access to a machine
 *  Another way used here is to pass the set to the present program.
 *  It will dispatch on the farm a number of copies of itself, say 100,
 *  open a remote channel and pass the parameters to its children
 *  as soon as a child finishes an instance it will start a new one
 *  To run the desired program, chanscript simply exe the scripts and waits for 
 *  completion.
 *  A further refinement could be to resubmit the missing cases after 95%
 *  have succuueded, and the delay exceeds the elapsed time of the first 19%
 *
 *  The public functions are
*/

#include "ac.h"
#include "remote_channel.h"
#include "wego.h"

typedef struct csStruct {
  AC_HANDLE h ;
  AC_DB db ;
  BOOL debug ;
  BOOL random, kallisto ;
  int nAction ;
  const char *dbName ;
  const char *dir, *target, *run, *project ;
  const char *inChanName, *outChanName ;
  const char *outFileName ;
  CHAN *inChan, *outChan ;
  CHAN *actionInChan, *actionOutChan, *collationDone ;
} CS ;

typedef struct buf64struct { char  t[64] ;} BUF64 ;
typedef struct buf128struct { char t[128] ;} BUF128 ;

typedef struct actionInStruct {AC_TABLE t ;} AIN ;
typedef struct actionOutStruct { BOOL done ; char run[64], sublib[64] ; double z ; int N0, N5, N10, N15, N20 ;} AOUT ;

typedef struct serverInStruct {char buf[128] ;} SIN ;
typedef struct serverOutStruct {char buf [128] ;} SOUT ;

typedef struct hitStruct { double z ; } HIT ;

/*************************************************************************************/
/********************************** action code **************************************/
/* the action code should ideally be linked in from another file 
 * with fixed function prototype
 * probably the stuct AIN AOUT CIN COUT
 * should be defined there
 */
/*************************************************************************************/
/* parse the kallisto read counts */
/*  tmp/Kallisto/DRR014223/FlyBase/abundance.tsv 
 * target_id	length	eff_length	est_counts	tpm
 * FBtr0081624|Gene|FBgn0000003	299	40.2138	11	1.61717
 * FBtr0071763|Gene|FBgn0000008	4847	4548	62.8265	0.0816693
 * FBtr0071764|Gene|FBgn0000008	5173	4874	6.63286	0.00804548
 */

static int klstParse (CS *cs, Array a0, Array a, double *totalp, DICT *dict, const char *nam)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, gene ;
  char *fNam = hprintf (h, "%s/Kallisto/%s/%s/abundance.tsv", cs->dir, nam, cs->target) ;
  ACEIN ai ;
  HIT *up, *vp ;
  
  if (! filCheckName (fNam, 0, "r"))
    {
      ac_free (h) ;
      return 0 ;
    }
  ai = aceInCreate (fNam, 0, h) ;
  aceInSpecial (ai, "\n") ;

  *totalp = 0 ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      float z = 0 ;

      if (!cp) continue ;
      dictAdd (dict, cp, &gene) ;

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;

      aceInStep (ai, '\t') ;
      if (! aceInFloat (ai, &z))
	continue ;

      up = arrayp (a0, gene, HIT) ;
      vp = arrayp (a, gene, HIT) ;
      up->z += z ;
      vp->z = z ;
      *totalp += z ;
      nn++ ;
    }

  ac_free (h) ;
  return nn ;
}

/*************************************************************************************/

static int geneGeneLanesParse (CS *cs, Array a0, Array a, double *totalp, DICT *dict, const char *nam)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, gene ;
  char *fNam1 = hprintf (h, "%s/GENELANES/%s/f.1.magic.GENE.geneSupport.u.gz", cs->dir, nam, nam, cs->target) ;
  char *fNam2 = hprintf (h, "%s/GENELANES/%s/f2.1.magic.GENE.geneSupport.u.gz", cs->dir, nam, nam, cs->target) ;
  ACEIN ai ;
  HIT *up, *vp ;
  char *qq = hprintf (h, "gunzip -c %s/GENELANES/%s/*.magic.GENE.geneSupport.u.gz", cs->dir, nam, nam, cs->target) ;

  if (! filCheckName (fNam1, 0, "r") && ! filCheckName (fNam2, 0, "r"))
    {
      ac_free (h) ;
      return 0 ;
    }
  ai = aceInCreateFromPipe (qq, "r", 0, h) ;
  aceInSpecial (ai, "\n") ;

  *totalp = 0 ;
  while (aceInCard (ai))
    {
      /*  awaiting: LT_magic	G_t_2L_0_1442_1	u	SRR1543246	40.0	51.00	3848 (bp) */
      char *cp = aceInWord (ai) ;
      float z = 0 ;

      if (!cp) continue ; /* target_class */

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;   /* gene name */
      if (!cp) continue ;
      dictAdd (dict, cp, &gene) ;

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* u, nu, ignore */

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* run name, ignore */

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* seq count, ignore */
      aceInStep (ai, '\t') ;
      if (! aceInFloat (ai, &z)) /* tag count */
	continue ;

      up = arrayp (a0, gene, HIT) ;
      vp = arrayp (a, gene, HIT) ;
      up->z += z ;
      vp->z = z ;
      *totalp += z ;
      nn++ ;
    }

  ac_free (h) ;
  return nn ;
} /* geneGeneLanesParse */

/*************************************************************************************/
/* parse the Magc pipeline GeneRun  read counts */
/*   tmp/GENERUNS/SRR1197399/SRR1197399.magic.GENE.u.geneSupport.ace.gz
 Gene G_t_2L_0_55
magic
Run_U "SRR1197399" 0.00 533.00 seqs 1068.00 tags 89.89 kb 6813.50 reads 1157.00 compRead  3203.00 err  0.00 a2g 1914.00 partial 3677.50 orphan 65.00 badTopo 23.50 Multi 0.00 AmbStrand
 */

static int geneGeneRunsParse (CS *cs, Array a0, Array a, double *totalp, DICT *dict, const char *nam)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, gene = 0 ;
  char *fNam = hprintf (h, "%s/GENERUNS/%s/%s.%s.GENE.u.geneSupport.ace.gz", cs->dir, nam, nam, cs->target) ;
  ACEIN ai ;
  HIT *up, *vp ;
  
  if (! filCheckName (fNam, 0, "r"))
    {
      if (0) fprintf (stderr, "Missing file %s\n", fNam) ;
      ac_free (h) ;
      return 0 ;
    }
  ai = aceInCreate (fNam, 0, h) ;

  *totalp = 0 ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      float z = 0 ;

            if (!cp) continue ;
      if (! strcmp (cp, "Gene"))
	{
	  gene = 0 ;
	  cp = aceInWord (ai) ;
	  if (!cp) continue ;
	  dictAdd (dict, cp, &gene) ;
	  continue ;
	}
      if (! gene)
	continue ;
      if (strcmp (cp, "Run_U"))
	continue ;
      /* we are on the Run_U line: Run_U "SRR1197399" 0.00 533.00 seqs 1068.00 tags */
      cp = aceInWord (ai) ; /* run name, ignore */
      if (! cp) continue ;
  
      cp = aceInWord (ai) ; /* score, ignore */
      if (! cp) continue ;

      cp = aceInWord (ai) ; /* sequence count, ignore */
      if (! cp) continue ;

      cp = aceInWord (ai) ; /* seqs, ignore */
      if (! cp) continue ;

      aceInStep (ai, '\t') ;
      if (! aceInFloat (ai, &z))  /* tag count */
	continue ;

      up = arrayp (a0, gene, HIT) ;
      vp = arrayp (a, gene, HIT) ;
      up->z += z ;
      vp->z = z ;
      *totalp += z ;
      nn++ ;
    }

  ac_free (h) ;
  return nn ;
} /* geneGeneRunsParse */

/*************************************************************************************/

static void action (void *vp)
{
  CS *cs = (CS *)vp ;
  char oldRun[64] ;
  AIN ain ;
  AC_HANDLE h = ac_new_handle () ;
  DICT *geneDict = dictHandleCreate (30000, h) ;
  double log2 = log(2.0) ;
  AOUT aout ;
  Array b0 = 0 ;

  oldRun[0] = 0 ;
  while (channelMultiGet (cs->actionInChan, &ain, 1, AIN))
    {
      Array a, a0, aa, cumuls ;
      const char *nam ;
      AC_HANDLE h1 = ac_new_handle () ;
      int ir, gene ;
      double z ;
      AC_TABLE tt = ain.t ;
      BOOL ok = TRUE ;

      cumuls = arrayHandleCreate (tt->rows, double, h1) ;
      a0 = arrayHandleCreate (30000, HIT, h1) ;
      aa = arrayHandleCreate (32, Array, h1) ;
      ok = TRUE ;
      for (ir = 0 ; ok && ir < tt->rows ; ir++)
	{
	  array (aa, ir, Array) = a = arrayHandleCreate (30000, HIT, h1) ;
	  nam = ac_table_printable (tt, ir, 1, 0) ;
	  if (cs->kallisto)
	    {
	      if (!klstParse (cs, a0, a, &z, geneDict, nam))
		ok = FALSE ;
	    }
	  else if (0) /* parse the ace file of a Runs */
	    {
	      if (!geneGeneRunsParse (cs, a0, a, &z, geneDict, nam))
		ok = FALSE ;
	    }
	  else    /* parse the tab file of the alnes of a sub lib */
	    {
	      if (!geneGeneLanesParse (cs, a0, a, &z, geneDict, nam))
		ok = FALSE ;
	    }
	  array (cumuls, ir, double) = z ;
	}
      if (!ok)
	continue ;
      for (ir = 0 ; ir < tt->rows ; ir++)
	{
	  double x, y, X = 0, X2 = 0, Y = 0, Y2 = 0, XY = 0 ;
	  int N = 0, N5 = 0, N10 = 0, N15 = 0, N20 = 0 ;
	  double index,  cumul = array (cumuls, ir, double) ;

	  a = array (aa, ir, Array) ;
	  if (! cs->random) b0 = a0 ;
	  for (gene = 0 ; b0 && gene < arrayMax (b0) ; gene++)
	    {
	      HIT *up = arrp (b0, gene, HIT) ;
	      HIT *vp = arrayp (a, gene, HIT) ;

	      x = up->z ; y = vp->z ;
	      if (x < 50)
		continue ;
	      index = 3 * log (y * 10000/cumul) / log2 ;
	      if (index > 5) 
		N5++ ;
	      if (index > 10) 
		N10++ ;
	      if (index > 15) 
		N15++ ;
	      if (index > 20) 
		N20++ ;
	      x = log(10+x) ; y = log(10+y) ;
	      X += x ; X2 += x * x ; Y += y ; Y2 += y * y ; XY += x*y ;
	      N++ ;
	    }   
	  x = y = z = 0 ;
	  if (N)
	    {
	      x = X2/N - X*X/(N * N) ;
	      y = Y2/N - Y*Y/(N * N) ;
	      z = XY/N - X*Y/(N * N) ;
	    }
	  if (x * y == 0) 
	    z = 0 ;
	  else
	    z = z / sqrt (x*y) ;

	  memset (&aout, 0, sizeof (AOUT)) ;
	  strncpy (aout.run, cs->random ? oldRun : ac_table_printable (tt, ir, 0, ""), 63) ;
	  strncpy (aout.sublib, ac_table_printable (tt, ir, 1, ""), 63) ;

	  aout.z = z ;
	  aout.N0 = N ;
	  aout.N5 = N5 ;
	  aout.N10 = N10 ;
	  aout.N15 = N15 ;
	  aout.N20 = N20 ;
	  channelPut (cs->actionOutChan, &aout, AOUT) ;
	}
      if (cs->random)
	{
	  ac_free (b0) ;
	  b0 = arrayCopy (a0) ;
	  strncpy (oldRun,  ac_table_printable (tt, 0, 0, ""), 63) ;
	}
      ac_free (h1) ;
      ac_free (tt) ; /* delegated from main task */
    }
  if (cs->random)
    ac_free (b0) ;
  if (1) /* action done signal */
    {
      memset (&aout, 0, sizeof (AOUT)) ;
      aout.done = TRUE ;
      channelPut (cs->actionOutChan, &aout, AOUT) ;
    }
  ac_free (h) ;
  return ;
} /* action */

/*************************************************************************************/

static void collation (void *vp)
{
  CS *cs = (CS *)vp ;
  int C = cs->nAction ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (cs->outFileName, ".sublib_correl.txt", 0, h) ;
  AOUT aout ;

  aceOutf (ao, "# Run\tSublib\tCorrel\tN0\tN5\tN10\tN15\tN20\n") ;
  /* collate the results */
  while (channelMultiGet (cs->actionOutChan, &aout, 1, AOUT))
    {
      if (aout.done) /* action done signal */
	{  
	  if (cs->debug)
	    wego_log (hprintf (h, "Action done signal received, remains %d : %s\n", C, timeShowNow ())) ;
	  if (--C <= 0) /* all actions are done */
	    {  
	      channelClose (cs->outChan) ;
	      break ;
	    }
	}
      else if (aout.N0)
	aceOutf (ao, "%s\t%s\t%g\t%d\t%d\t%d\t%d\t%d\n"
		 , aout.run, aout.sublib, aout.z
		 , aout.N0, aout.N5, aout.N10, aout.N15, aout.N20
		 ) ;
    }
  
  ac_free (h) ;
  channelPut (cs->collationDone, &C, int) ; 
  return ;
} /* collation */

/*************************************************************************************/
/*************************************************************************************/
static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// chanscript: run a script familly using remote channels\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, Mars 2017, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// When a given scripts must be run for a large number of parameters\n"
	    "// this code wil spawn a number of copy of itself on the farm and pass \n"
	    "// the parameters via a remote channel. Then each child will exe the script\n"
	    "// and report the exit value. In this way, we do not need ot submit as many\n"
	    "// independant programs to the farm\n"
	    "// Usage:\n"
	    "// chanscript [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "//   --db : host:port mandatory the acedb MetaDB server\n"
	    "//   --dir : directory to find the data\n"
	    "//   --target : needed to compose the data file name\n"
	    "//   --run : apply to a single run\n"
	    "//   --project : apply to T runs of a given project\n"
	    "//   --cin : name of the incoming remote channel\n"
	    "//   --cout : name of the outgoing remote channel\n"
	    "//   --debug : debug verbose \n"
	    "//   -h -help --help : this help\n"
	    "//   --CP common_parameters : passed to each script\n"
	    "//   --VP variable_parameters : coma delimited, one black for each instance\n"
	    "//   -C  <int> : number parallel action,  default 4"
	    "//   -T  <int> : number analyzed instances,  default 12"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  chanscript --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  CS cs ;
  TASK *task ;
  int T, C ;
  int i ;
  const char *errors = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  freeinit () ; 
  memset (&cs, 0, sizeof (CS)) ;

  /* optional arguments */
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  cs.debug = getCmdLineBool (&argc, argv, "--debug") ; 
  cs.random = getCmdLineBool (&argc, argv, "--random") ; 
  cs.kallisto = getCmdLineBool (&argc, argv, "--kallisto") ; 

  getCmdLineOption (&argc, argv, "--db", &cs.dbName) ; 
  getCmdLineOption (&argc, argv, "--dir", &cs.dir) ;
  getCmdLineOption (&argc, argv, "--target", &cs.target) ;
  getCmdLineOption (&argc, argv, "--project", &cs.project) ;
  getCmdLineOption (&argc, argv, "--run", &cs.run) ;
  getCmdLineOption (&argc, argv, "--cin", &cs.inChanName) ;
  getCmdLineOption (&argc, argv, "--cout", &cs.outChanName) ;
  getCmdLineOption (&argc, argv, "-o", &cs.outFileName) ;
  T = 100000 ;  C = 4 ; /* defaults */
  getCmdLineInt (&argc, argv, "-T", &T) ;
  getCmdLineInt (&argc, argv, "-C", &C) ;
  cs.nAction = C ;

  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[argc-1])) ;
 
  if (!cs.dir || !cs.target)
    messcrash ("missing arguments, --dir and --target are mandatory, try --help\n") ;


  fprintf (stderr, "//start: %s\n", timeShowNow()) ;
  if (cs.dbName[0])
    {
      errors = 0 ;
      if (! (cs.db = ac_open_db (cs.dbName, &errors)))
	messcrash ("cannor open database %s : %s\n"
		   , cs.dbName, errors ? errors : "") ;
    }

  /* create the communication channels with the action code */
  cs.collationDone = channelCreate (C, int, h) ;
  cs.actionInChan = channelCreate (T, AIN, h) ;
  cs.actionOutChan = channelCreate (T, AOUT, h) ;
  if (cs.debug)
    {
      channelDebug (cs.collationDone, 0, "collationDoneChan") ;
      channelDebug (cs.actionInChan, 0, "actionInChan") ;
      channelDebug (cs.actionOutChan, 0, "actionOutChan") ;
    }

  /* launch the action code */
  wego_max_threads (8) ;
  wego_go (collation, &cs, CS) ; 
  for (i = 0 ; i < C ; i++)
    { wego_go (action, &cs, CS) ; }

  /* create the communication to the server, or fake it */
  if (cs.inChanName && cs.outChanName)
    {
      task = taskClientInit (& argc, argv, cs.h) ;
      /* initialze the channels */
      
      cs.inChan = taskCreateReaderChannel (task, cs.inChanName, 100, SIN) ;
      cs.outChan = taskCreateReaderChannel (task, cs.outChanName, 1, SOUT) ;

      if (cs.debug)
	{
	  channelDebug (cs.inChan, 0, "inChan") ;
	  channelDebug (cs.outChan, 0, "outChan") ;
	}
    }
  else if (cs.run)
    {
      SIN sin ;
      cs.inChan = channelCreate (T, SIN, h) ;
      strncpy (sin.buf, cs.run, 127) ;
      channelPut (cs.inChan, &sin, SIN) ;
      channelClose (cs.inChan) ;
    }
  else if (cs.project)
    {
      char qq[1024] ;
      AC_TABLE t ; 
      int i ;
      
      sprintf (qq, "select r from p in class project where p == \"%s\", r in p->run where r#sublibraries", cs.project) ; 
      if (cs.debug)
	wego_log (hprintf (h,"query : %s\n", qq)) ;
      errors = 0 ;
      t = ac_bql_table (cs.db, qq, 0, 0, &errors, 0) ;
      if (errors)
	messcrash ("Bad query %s\n", qq) ;
      if (cs.debug)
	wego_log (hprintf (h,"query recovered a table with %d lines\n", t ? t->rows : 0)) ;

      cs.inChan = channelCreate (T, SIN, h) ;
      for (i = 0 ; i < t->rows && i < T ; i++)
	{ 
	  SIN sin ;
	  strncpy (sin.buf, ac_table_printable (t, i, 0, ""), 127) ;
	  channelPut (cs.inChan, &sin, SIN) ;
	}
      channelClose (cs.inChan) ;
    }

  if (cs.inChan)
    {
      SIN sin ;
      AIN ain ;

      while (channelMultiGet (cs.inChan, &sin, 1, SIN))
	{
	  char qq[1024] ;
	  
	  sprintf (qq, "select r,s from r in class run where r == \"%s\", s in r->sublibraries", sin.buf) ; 
	  if (cs.debug)
	    wego_log (hprintf (h,"query : %s\n", qq)) ;
	  errors = 0 ;
	  ain.t = ac_bql_table (cs.db, qq, 0, 0, &errors, 0) ;
	  if (errors)
	    messcrash ("Bad query %s\n", qq) ;
	  if (cs.debug)
	    wego_log (hprintf (h,"query recovered a table with %d lines\n", ain.t ? ain.t->rows : 0)) ;
	  if (ain.t && ain.t->rows)
	    channelPut (cs.actionInChan, &ain, AIN) ; /* will free tt */
	}
      channelClose (cs.actionInChan) ;
    }
  /* wait for colation */
  channelMultiGet (cs.collationDone, &i, 1, int) ;

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  wego_flush () ;
  ac_free (cs.h) ;

  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

