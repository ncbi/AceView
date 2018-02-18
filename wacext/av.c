/*  File: av.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2010
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
 *  This is the master program driving alignment and analysis
 *  of mRNA relative to a reference genome
 *  The input may be full length mRNAs (i.e. RefSeq from genbank)
 *  ESTs (i.e. from dbEST) or next-gen deep sequencing
 *  short (i.e. ABI/Solexa) or long (Roche) or LifeTech-SOLiD
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

#include "ac.h"
#include "aceio.h"
#include "drmaa.h"  /* sun grid engine Distributed Ressoucre management Application API */

static BOOL debug = FALSE ;
static BOOL is64 = FALSE ;

/*************************************************************************************/
/*********************************************************************/
typedef struct avStruct { 
  AC_HANDLE h ;
  ACEIN ai ;
  ACEOUT ao ;
  const char *target, *project ;
  const char *maxMem, *queue, *dbName ;
  const char *reportDir, *inFileName, *outFileName ;
  BOOL gzi, gzo, drmaa, blastp ;
  drmaa_job_template_t *jt ;
  AC_DB db ;
} AV ;

/*************************************************************************************/

static void avInit (AV *av)
{
  freeinit () ; 
  
  if (sizeof(void*) == 8)
    is64 = TRUE ;
  
  memset (av, 0, sizeof (AV)) ;
  av->h = ac_new_handle () ;
  av->reportDir = "/home/mieg/ACEVIEWHELP" ;
  av->maxMem = "8G" ;
  av->queue = "low" ;
  return ;
} /* avInit */

/*********************************************************************/

static void avOpenInputOutput (AV *av, char *type)
{
  if (av->outFileName)
    { 
      ac_free (av->ao) ;
      if (av->gzo)
	av->ao = aceOutCreateToPipe (hprintf (av->h, "gzip  > %s.%s.gz", av->outFileName, type), av->h) ;
      else
	av->ao = aceOutCreateToFile (hprintf (av->h, "%s.%s", av->outFileName, type), "wb", av->h) ;
      if (! av->ao)
	exit (1) ;
    }
  else if (! av->ao) /* Re-open the file only if it is named */
    {
      if (av->gzo)
	av->ao = aceOutCreateToPipe ("gzip ", av->h) ;
      else
	av->ao = aceOutCreateToStdout (av->h) ;
    }
  if (av->inFileName)
    { 
      int k = strlen (av->inFileName) ;
      const char *ccp = k > 3 ? av->inFileName + k - 3 : av->inFileName ;
      ac_free (av->ai) ;
      ccp = av->inFileName ;
      if (av->gzi || ! strcmp (ccp, ".gz"))
	av->ai = aceInCreateFromPipe (hprintf (av->h, "gzip  -dc %s", av->inFileName), "r", 0, av->h) ;
      else
	av->ai = aceInCreateFromFile (av->inFileName, "r", 0, av->h) ;
      if (! av->ai)
	exit (1) ;
    }
  else if (! av->ai) /* Re-open the file only if it is named */
    {
      if (av->gzi)
	av->ai = aceInCreateFromPipe ("gzip -dc -", "r", 0, av->h) ;
      else
	av->ai = aceInCreateFromStdin (0, 0, av->h) ; 
    }
  if (! av->ai)
    messcrash ("cannot read input file %s", av->inFileName ?  av->inFileName : "stdin") ;
  
} /* avOpenInputOutput */

/*************************************************************************************/
/*************************************************************************************/

static BOOL drmaaExit (AV *av)
{
  char error [DRMAA_ERROR_STRING_BUFFER] ;
  int errnum ;
  
  if (av->jt)
    {
      errnum = drmaa_delete_job_template (av->jt, error, DRMAA_ERROR_STRING_BUFFER);
      
      if (errnum != DRMAA_ERRNO_SUCCESS) 
	fprintf (stderr, "Could not delete job template: %s\n", error);
    }
  
  errnum = drmaa_exit (error, DRMAA_ERROR_STRING_BUFFER) ;
  
  if (errnum != DRMAA_ERRNO_SUCCESS)
    messcrash ("Could not shut down the DRMAA library: %s\n", error) ;
  
  fprintf (stderr, "DRMAA library was closed successfully\n");

  return TRUE ;
}

/*************************************************************************************/

static BOOL drmaaInit (AV *av)
{
  char error [DRMAA_ERROR_STRING_BUFFER] ;
  int errnum ;
  
   errnum = drmaa_init (NULL, error, DRMAA_ERROR_STRING_BUFFER) ;
   
   if (errnum != DRMAA_ERRNO_SUCCESS)
     messcrash ("error initialising drmaa library : ", error) ;
 
   fprintf (stderr, "DRMAA library was started successfully\n");
   
   errnum = drmaa_allocate_job_template (&(av->jt), error, DRMAA_ERROR_STRING_BUFFER);
   
  if (errnum != DRMAA_ERRNO_SUCCESS) 
    {
      drmaaExit (av) ;
      messcrash ("Could not create job template: %s\n", error);
    }
   return TRUE ;
}

/*************************************************************************************/

static int av_get_proteins (AV *av)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, ii, NN = 1000 ;
  vTXT txt = vtxtHandleCreate (h) ;
  AC_KEYSET ks, ks1 ;

  system ("mkdir Blastp") ;


  if (1) /* worm */
    vtxtPrint (txt, "find kantor COUNT { >product peptide } > 0 AND (peptide:2 > 50 || COUNT {>product (NH2_Complete && COOH_Complete) || Refseq } > 0) ") ;
  else /* human */
    vtxtPrint (txt, "Find mrna best_in_gene ; > product peptide ; > kantor ; peptide:2 > 100") ;
  
  vtxtPrint (txt, "  AND (NOT Blastp_date || Blastp_date < 2010-11-09 )") ;

  ks = ac_dbquery_keyset (av->db, vtxtPtr (txt), h) ;
  nn = ac_keyset_count (ks) ;

  for (ii = 0 ; ii * NN < nn ; ii++)
    {
        /* subset starting at x0 of length nx */
      ks1 = ac_subset_keyset (ks, ii*NN, NN, h) ;
      ac_command_keyset (av->db
		       , messprintf ("peptide Blastp/f.%d.fasta", ii)
		       , ks1
		       , h) ;
    }
  ac_free (h) ;
  return nn ;
} /* av_get_proteins */

/*************************************************************************************/

static int av_drmaa_run_blastp (AV *av)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, ii ;
  vTXT txt = vtxtHandleCreate (h) ;
  int errnum ;
  char error[DRMAA_ERROR_STRING_BUFFER];
  drmaa_job_ids_t *ids = NULL;

  if (! filCheckName ( "bin/run_blastp","tcsh", "x"))
    {
      fprintf (stderr, "missing executable bin/run_blastp.tcsh") ;
      return 0 ;
    }

  errnum = drmaa_set_attribute (av->jt
				, DRMAA_REMOTE_COMMAND, "bin/run_blastp.tcsh"
				, error, DRMAA_ERROR_STRING_BUFFER
				);
  if (errnum != DRMAA_ERRNO_SUCCESS) 
    {
      fprintf (stderr, "Could not set attribute \"%s\": %s\n",
	       DRMAA_REMOTE_COMMAND, error) ;
      return 0 ;
    }
  
  for (ii = 0 ; ii < 5 ; ii++)
    {
      /* check is the data file is ready, and not yet analysed */
      vtxtClear (txt) ;
      vtxtPrintf (txt, "Blastp/f.%d", ii) ;
      if (! filCheckName(vtxtPtr (txt), "fasta", "r"))
	continue ;
      if (filCheckName(vtxtPtr (txt), "blastp", "r"))
	continue ;
      
      /* submit the job to the farm */
       if (1)
	 {
	   char buf[64] ;
	   char jobid[DRMAA_JOBNAME_BUFFER];
	   const char *args[2] ;
	   args[0] = buf ;
	   args[1] = NULL ;

	   sprintf (buf, "%d", ii) ;
	   errnum = drmaa_set_vector_attribute (av->jt
						, DRMAA_V_ARGV, args
						, error, DRMAA_ERROR_STRING_BUFFER
						);
	   if (errnum != DRMAA_ERRNO_SUCCESS)
	     {
	       fprintf (stderr, "Could not set attribute \"%s\": %s\n",
			DRMAA_REMOTE_COMMAND, error);
	       continue ;
	     }
	   
	  
	   errnum = drmaa_run_job (jobid, DRMAA_JOBNAME_BUFFER
				   , av->jt
				   , error, DRMAA_ERROR_STRING_BUFFER
				   ) ;
	   
	   if (errnum != DRMAA_ERRNO_SUCCESS) 
	     {
	       fprintf (stderr, "Could not submit job: %s\n", error);
	       continue ;
	     }


	   printf ("Your job has been submitted with id %s\n", jobid);
	   
	   /*
	   if (0)
	     {
	       errnum = drmaa_synchronize (jobids, DRMAA_TIMEOUT_WAIT_FOREVER,
					   1, error, DRMAA_ERROR_STRING_BUFFER);
	       
	       if (errnum != DRMAA_ERRNO_SUCCESS) {
		 fprintf (stderr, "Could not wait for jobs: %s\n", error);
	       }
	       else {
		 printf ("All job tasks have finished.\n");
	       }
	     }
	   
	   drmaa_release_job_ids (ids);
	   */

	 }
    }

  ac_free (h) ;
  return nn ;
} /*  av_drmaa_run_blastp */

/*************************************************************************************/
/*****************************  actual work ******************************************/
/*************************************************************************************/

static BOOL av_submit (AV *av, vTXT txt, const char *command, const char *fnam)
{
  vtxtPrintf (txt, "scripts/submit %s \t%s\\n", fnam, command) ;
  
  return 1 ;
} /* av_submit */

/*************************************************************************************/

static int a2_align (AV *av)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT command = vtxtHandleCreate (h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  vTXT fnam = vtxtHandleCreate (h) ;
  vTXT mkdir = vtxtHandleCreate (h) ;
  int nn = 0 ;
  char *isSolid = "" ;
  char *species = "hs" ;
  char *target = "mito" ;
  char *clipPolyA = "-clipPolyA" ;
  char *clipPolyT = "-clipPolyT" ;
  int minEntropy4 = 16 ;
  int seedLength = 18 ;
  int minAli = 24 ;
  int intronMaxLength = 25000 ;
  char *isSplice = "" ;
  char *slBonus = "" ;
  char *targetBonus = "" ;
  char *Xintrons = "" ;
  char *isStranded = "" ;
  char **targetp ;
  char *targets[] = {"mito", 0} ;
  char **manipp ;
  char *manips[] = {"IBC", 0 } ;
  char **tissuep ;
  char *tissues[] = {"s930", 0 } ;
  char **lanep ;
  char *lanes[] = {"f.1", 0} ;

  /* create a tcsh command */
  vtxtPrintf (command, 
	      "#!/bin/tcsh -ef\n"
	      "#$ -S /bin/tcsh\n"
	      "#$ -P %s\n"
	      "#$ -j y\n"
	      "#$ -m e\n"
	      "#$ -v SGE_SUMMARY=\"stderr\"\n"
	      "#$ -v SGE_NOMAIL\n"
	      "\n"
	      , av->queue
	      ) ;


  /* mkdir the relevant directories */
  
  vtxtPrintf (mkdir, " if (! -d tmp) mkdir tmp\n") ;
  system (vtxtPtr (mkdir)) ;

  for (targetp = targets ; *targetp ;  targetp++)
    {
      vtxtClear (mkdir) ;
      vtxtPrintf (mkdir, " if (! -d tmp/PHITS_%s) mkdir tmp/PHITS_%s\n"
		  , *targetp
		  , *targetp
		  ) ;
      system (vtxtPtr (mkdir)) ;
      for (manipp = manips ; *manipp ; manipp++)
	{ 
	  vtxtClear (mkdir) ;
	  vtxtPrintf (mkdir, " if (! -d tmp/PHITS_%s/%s) mkdir tmp/PHITS_%s/%s\n"
		      , *targetp, *manipp
		      , *targetp, *manipp
		      ) ;
	  system (vtxtPtr (mkdir)) ;
	  for (tissuep = tissues ; *tissuep ; tissuep++)
	    { 
	      vtxtClear (mkdir) ;
	      vtxtPrintf (mkdir, " if (! -d tmp/PHITS_%s/%s/%s) mkdir tmp/PHITS_%s/%s/%s\n"
			  , *targetp, *manipp, *tissuep
			  , *targetp, *manipp, *tissuep
			  ) ;
	      system (vtxtPtr (mkdir)) ;

	      for (lanep = lanes ; *lanep ; lanep++)
		{ 
		  vtxtClear (fnam) ;		  
		  vtxtPrintf (fnam, "tmp/PHITS_%s/%s/%s"
			      , *targetp, *manipp, *tissuep, *lanep
			      ) ;

		  vtxtPrintf (txt, "clipalign -gzo -maxHit 10 %s", isSolid) ; 
		  vtxtPrintf (txt, " -p Fastc/%s/%s/%s.fastc.gz", *manipp, *tissuep, *lanep) ;
		  vtxtPrintf (txt, " -t Targets/%s.%s.fasta", species, target) ;
		  vtxtPrintf (txt, " %s %s" , clipPolyA, clipPolyT) ;
		  vtxtPrintf (txt, " -minEntropy ", minEntropy4) ;
		  vtxtPrintf (txt, " -seedOffset 1 -seedShift 5 -seedLength %d", seedLength) ;
		  vtxtPrintf (txt, " -minAli %d", minAli) ;
		  vtxtPrintf (txt, " %s %s %s ", isSplice,  slBonus, targetBonus) ;
		  vtxtPrintf (txt, " $Xintrons $isStranded") ;
		  vtxtPrintf (txt, " -intronMaxLength %d", intronMaxLength) ;
		  vtxtPrintf (txt, " -o %s", vtxtPtr(fnam)) ;
		  
		  av_submit (av, command, vtxtPtr (txt), vtxtPtr (fnam)) ;
		}
	    }
	}
    }

  if (1)
    fprintf (stderr, vtxtPtr(command)) ;
  system (vtxtPtr(command)) ;
  ac_free (h) ;
  return nn ;
}

/*************************************************************************************/
/* Report a synthesis of the SNP. i.e. the mito fingreprint of the sample
 *
 * Tissue external (not tag-count-count-count)
 * faux, tester sur la position 8565 dans IBC, faire une approche probabibliste, et les 3 dernieres colonnes sont fausses
 */

/*
      if (! -e tmp/SNP/ALL/$CALI.$target.snplist.5.20.*txt) continue
*/

static void av_snp_report (AV *av, int s1, int s2)
{
  AC_HANDLE h = ac_new_handle () ;

  /* Expected Input: target position strand manip tissue snp count BRS:bR,bS,bT */

  aceOutf (av->ao, "Exported %d SNPs in file  tmp/SNP/%s.%s.map.%d.%d.txt\n" , 0, av->project, av->target, s1, s2) ;
  if ( av->outFileName && av->reportDir)
    system (hprintf (h,"cp %s %s\n", av->outFileName, av->reportDir)) ;
  
  ac_free (h) ;
  return ;
}  /* av_snp_report  */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: av -p project ......\n"
	   "//    Master code for alignemnt and analysis of transcriptome data.\n"
	   "//    Author: Danielle et Jean Thierry-Mieg, mieg@ncbi.nlm.nih.gov\n"
	   "//                                    www.ncbi.nlm.nih.gov/AceView\n"
	   "// Example:  av -p body_map -align\n"
	   "// [-p project [-config config_file]] : project name\n"
	   "//    project : defaults to the environment variable AV_PROJECT\n"
	   "//    config_file : deafults to ./av.<project>.config\n"
	   "//    The configuration file contains the description of the project\n"
	   "//    See the documentation in README.project_config\n"
	   "//\n"
	   "// -q queue : [default unified], on of the SGE queues (low,unified... see pstat)\n"
	   "// -memory value : [default 8G] , memory requested on the SGE farm\n"
	   "// -i fileName : input file name, equivalent to redirecting stdin\n"
	   "// -gzi : gunzip the input file (default for .gz file)\n"
	   "// -o fileName : output file name prefix\n"
	   "//    .txt, .err or some other endding are added to the output file name\n"
	   "// -gzo : gzip the output file\n"
	   "//    further add .gz to all exportted files\n"
	   "// -drmaa : open the Sun Grid Engine interface\n"
	   "// -db ACEDB : open an acedb database\n"
	   "// -blastp : run blastp on all the relevant protein of the acedb database\n"
	   "// -silent : suppress title lines and status reports from the output\n"
	   
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

int main (int argc, const char **argv)
{
  char commandBuf [1000] ;
  char *error ;
  AV av ;
  
  avInit (&av) ;
  
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;
  
  getCmdLineOption (&argc, argv, "-o", &(av.outFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(av.inFileName)) ;
  getCmdLineOption (&argc, argv, "-q", &(av.queue)) ;
  getCmdLineOption (&argc, argv, "-memory", &(av.maxMem)) ;
  getCmdLineOption (&argc, argv, "-db", &(av.dbName)) ;
  av.gzi = getCmdLineOption (&argc, argv, "-gzi", 0) ;
  av.gzo = getCmdLineOption (&argc, argv, "-gzo", 0) ;
  av.drmaa = getCmdLineOption (&argc, argv, "-drmaa", 0) ;
  av.blastp = getCmdLineOption (&argc, argv, "-blastp", 0) ;
  
  if (!  getCmdLineOption (&argc, argv, "-p", &(av.project)))
    av.project = getenv ("AV") ;
  if (! av.project)
    messcrash ("Please specify the project as -p option, or as setenv AV") ;
  
  avOpenInputOutput (&av, "txt") ;
  
  
  if (av.dbName)
    {
      av.db = ac_open_db (av.dbName, &error) ;
      if (! av.db)
	messcrash ("cannot open acedb database %s", av.dbName) ;
      if (av.blastp)
	av_get_proteins (&av) ;
    }


  if (av.drmaa)
    {
      drmaaInit (&av) ;
      if (av.blastp)
	av_drmaa_run_blastp (&av) ;
    }

  /* quo vadis */
  if (getCmdLineOption (&argc, argv, "-s3", 0))
    av_snp_report (&av, 5, 20) ;

  if (av.drmaa)
    drmaaExit (&av) ;
 
  ac_db_close (av.db) ;
  ac_free (av.h) ;

  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
