/*  File: nnaligner.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2019
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
 *
 *
 * OBJECTIVE: 
 *   Train and run a DNA aligner using a neural-net as its core extension engine
 */

/*
 #define ARRAY_CHECK 
 #define MALLOC_CHECK 
*/
#include "ac.h"

#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

typedef struct nnaStruct {
  AC_HANDLE h ;
  BOOL checkSam ;
  BOOL bestSam ;
  BOOL score ;
  BOOL exportIntrons ;
  BOOL gzo, gzi ;
  const char *target ;
  const char *samFiles ;
  const char *genomeFile ;
  const char *outFileName ;
  Array targetDnas ;
} NNA ;

typedef struct itsStruct { int target, a1, a2, n ;} ITS ;

static DICT *targetDict = 0 ;

static int itsOrder (const void *va, const void *vb)
{
  const ITS *a = (const ITS *)va, *b = (const ITS *)vb ;
  int n ;
  
  n = a->target - b->target ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/* scan a collection of sam files and perform the action */

static void nnParseGenomeFile (NNA *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (pp->genomeFile, pp->gzi, h) ;
  const char *targetName = pp->target ;
  int target = 0 ;
  int state = 0 ;
  int x = 0, dx ;
  Array dna = 0 ;
  int nChroms = 0 ;
  long int nBases = 0 ;

  if (! pp->targetDnas)
    pp->targetDnas = arrayHandleCreate (64, Array, pp->h) ;

  if (! ai)
    messcrash ("Cannot open genome fasta file -g %s\n", pp->genomeFile) ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || ! *cp)
	continue ;
      switch ((int) (*cp))
	{
	case '>':
	  if (dna)
	    dnaEncodeArray (dna) ;
	  dna = 0 ;
	  if (!cp[1] || (targetName  && strcmp (targetName, cp + 1)))
	    state = 0 ; 
	  else
	    {
	      state = 1 ;
	      nChroms++ ;
	      x = 0 ;
	      dictAdd (targetDict, cp+1, &target) ;
	      dna = array (pp->targetDnas, target, Array) =
		arrayHandleCreate (1000000, char, pp->h) ;
	    }
	  break ;
	default:
	  if (state != 1)
	    break ;
	  dx = strlen (cp) ;
	  nBases += dx ;
	  array (dna, x + dx, char) = 0 ; /* add a double terminal zero */
	  memcpy (arrp (dna, x, char), cp, dx) ;
	  x += dx ;
	  break ;
	}
    }
  if (dna)
    dnaEncodeArray (dna) ;
  fprintf (stderr, "%s Parsed %d targets, %ld bases in %s\n"
	   , timeShowNow() 
	   , nChroms, nBases, pp->genomeFile) ;
 ac_free (h) ;
}  /*  nnParseGenomeFile */

/*************************************************************************************/
/*************************************************************************************/
/* scan a collection of sam files and perform the action */

static void nnSamScanFiles (NNA *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  ACEOUT ao = 0 ;
  int nFile = 0 ; 
  DICT *intronDict = 0 ;
  KEYSET intronSupport = 0 ;

  if (pp->score)
    ao = aceOutCreate (pp->outFileName, ".sam_score", pp->gzo, h) ; 
  if (pp->exportIntrons)
    {
      intronDict = dictHandleCreate (1000000, h) ;
      intronSupport = keySetHandleCreate (h) ;
    }

  while ((ai = aceInCreateFromList (ai, &nFile, pp->samFiles, FALSE, h)))
    {
      if (pp->score)
	samFileScore (ai, ao, pp->target, targetDict, pp->targetDnas, intronDict, intronSupport) ;
      if (pp->exportIntrons)
	samFileExportIntrons (ai, pp->target, targetDict, intronDict, intronSupport) ;
      else
	samFileCheckCigars (ai, pp->target)  ;
    }

  if (pp->exportIntrons && keySetMax (intronSupport))
    {
      int ni = 0, ns = 0, i, iMax = keySetMax (intronSupport) ;
      ITS *up ;
      Array hits = arrayHandleCreate (100000, ITS, h) ; 
      ACEOUT ao = aceOutCreate (pp->outFileName, ".introns.txt", pp->gzo, h)  ;
      aceOutDate (ao, "###", "Intron support found int the sam files") ;
      aceOutf (ao, "### %s\n", pp->samFiles) ;
      aceOutf (ao, "# Chromosome\tFrom\tTo\tNumber of supporting reads\n") ;

      for (i = 1 ; i < iMax ; i++)
	{
	  int n = keySet (intronSupport, i) ;
	  if (n)
	    {
	      ns += n ; ni++ ;
	      up = arrayp (hits, arrayMax (hits), ITS) ;
	      up->n = n ;
	      sscanf (dictName (intronDict, i), "%d\t%d\t%d", &(up->target), &(up->a1), &(up->a2)) ;
	    }
	}

      arraySort (hits, itsOrder) ;
      iMax = arrayMax (hits) ;
      for (i = 0, up = arrp (hits, 0, ITS) ; i < iMax ; i++)
	aceOutf (ao, "%s\t%d\t%d\t%d\n"
		 , dictName (targetDict, up->target)
		 , up->a1, up->a2, up->n
		 ) ;
      fprintf (stderr, "Exported %d introns with %d supports in file %s\n"
	       , ni, ns, aceOutFileName (ao)) ;
    }

  ac_free (h) ;
}  /*  nnSamScanFiles */

/*************************************************************************************/
/*************************************************************************************/

static void usage (const char *message)
{
  if (! message)
  fprintf (stderr,
	   "// nnaligner: data preparation for a neural-net aligner\n"
	   "// Authors: Danielle and Jean Thierry-Mieg, NCBI, 2019, mieg@ncbi.nlm.nih.gov\n"
	   "// PURPOSE\n"
	   "//   Read (colllections of) sam/bam files, check and score them\n"
	   "//   select the best alignemnt, export these in BAM format or prepare\n"
	   "//   supervised training data for a local DNA/RNA NN aligner\n"
	   "//\n"
	   "// USAGE :\n"
	   "//   nnaligner [options]\n"
	   "//   try: -h --help \n"
	   "// The options are all optional and may be specified in any order\n"
	   "//\n"
	   "// INPUT\n"
	   "//   -i fileList : [default stdin] \n"
	   "//      example:   A:f1.sam,B:f2.bam,f3.sam.gz\n"
	   "//      The labels (A: B: etc) are optional\n"
	   "//      Any file called *.gz is gunzipped\n"
	   "//      Other preporcessing can be specified by using < as prefix\n"
	   "//        example:  -i \'< obj_get f.fasta.gz | gunzip\'\n"
	   "//   --gzi\n"
	   "//      Forces gunzip, usefull if parsing from stdin\n"
	   "//   -I input_format:\n"
	   "//      SAM (used by default on files called *.sam or *.sam.gz\n"
	   "//      BAM (used by default on files called *.bam\n"
	   "//   -g fastaFile : optional\n"
	   "//      Genome of Genes target fasta file\n"
	   "//      Needed to score the substitutions or construct the seeds\n"
	   "//\n"
	   "// OUTPUT\n"
	   "//   -o outFilePrefix\n"
	   "//      example -o dd/f : all output files will be called dd/f.*\n"
	   "//   --gzo\n"
	   "//     all large output files will be gzipped\n"

	   "//\n"
	   "// ACTIONS\n"
	   "//   --checkSam\n"
	   "//     Check the format and the cigar strings\n"
	   "//     Export only error messages\n"
	   "//   --bestSam\n"
	   "//     For each read, export lines with best scores\n"
	   "//   -u --unique\n"
	   "//     For each read, export a single line with best score\n"
	   "//   --introns\n"
	   "//     Export the list of introns (N in theCigar string)\n"
	   "//   --score\n"
	   "//     Export the list of introns (N in theCigar string)\n"
	   ) ;


  if (message)
      fprintf (stderr, "// %s\nFor more information try:  nnaligner -h or --help\n", message) ;
   exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  NNA p ;
  AC_HANDLE h = 0 ;
  char commandBuf [32000] ;
  const char *error ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&p, 0, sizeof (NNA)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  p.h = h ;

  if (argc < 2)
    usage (0) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help") ||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (0) ;

  if (getCmdLineBool (&argc, argv, "-version") || 
      getCmdLineBool (&argc, argv, "--version"))
    {
      fprintf (stderr, "nnaligner version 0.0.1\n") ;
      exit (0) ;
    }

  getCmdLineOption (&argc, argv, "-i", &(p.samFiles)) ;
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  getCmdLineOption (&argc, argv, "-g", &(p.genomeFile)) ;
  getCmdLineOption (&argc, argv, "-t", &(p.target)) ;
  p.checkSam = getCmdLineBool (&argc, argv, "--checkSam") ;
  p.score = getCmdLineBool (&argc, argv, "--score") ;
  p.bestSam = getCmdLineBool (&argc, argv, "--bestSam") ;
  p.exportIntrons = getCmdLineBool (&argc, argv, "--introns") ;

  p.gzi = getCmdLineBool (&argc, argv, "--gzi") ;
  p.gzo = getCmdLineBool (&argc, argv, "--gzo") ;
  
  if (argc > 1)
    {
      int i ;
      
      i = dnaEncodeChar ['a'] ; /* to please the linker */
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\nTry -h or --help\n") ;
      exit (1) ;
    }
  
  targetDict = dictHandleCreate (256, h) ;

  if (p.samFiles)
    aceInCheckList (p.samFiles, &error, p.h) ;
  if (p.genomeFile)
    nnParseGenomeFile (&p) ;
  if (p.samFiles)  /* san the files and perform the requested work */
    nnSamScanFiles (&p) ;

  ac_free (p.h) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
    }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
 
