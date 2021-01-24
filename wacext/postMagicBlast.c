/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 *  Copyright (C) D and J Thierry-Mieg 2020
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
 * This file is part of the AceView/Magic project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov

 * 15feb 2020
 * postmagicBlast.c
 *   Input: Output from MagicBlast (sam or tabular)
 *   Output: 
 *      best ali, ali.ace statitistics, introns, SNPs 
 *      i.e. every interesting post-treatment we can think of
 *      using the additive .tsf format
 */
#define MALLOC_CHECK  
#define ARRAY_CHECK  

#include "../wac/ac.h"

typedef struct pmbStruct { 
  AC_HANDLE h ;
  
  const char *run ; 
  
  const char *inFileName ; 
  const char *outFileName ; 

  BOOL gzi, gzo ;
} PMB ;

/*************************************************************************************/

static void pmbIntrons (PMB *pmd)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (pmd->inFileName, pmd->gzi, h) ;

  while (aceInCard (ai))
    {
      const char *ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;
    }
  
  ac_free (h) ;
} /* pmbIntrons */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// postMagicBlast: post treatment of the MagicBLAST output\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, with help from Tom Madden and Greg Boratyn, NCBI, October 2020, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a fasta or fastq file, aligned using MagicBlast\n"
	    "//    import the results\n"
	    "//    export all useful post treatment we may want in .tsf format\n"
	    "//\n"
	    "// Syntax:\n"
	    "// postMagicBlast -run run -i f.tabular [options]\n"
	    "//   run MANDATORY name of the run being analyzed\n"
	    "// All other parameters are optional and may be specified in any order\n"
	    "//   -db <dir>: MetaDB directory (default MetaDB)\n" 
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action\n" 
            "//   -gzo : gzip all output files\n"
	    "//   -export [TPbafpmg...] : only export some groups of columns, in the requested order\n"
	    "//      default: if -export is not specified, export all columns in default order\n"
	    "//            E: main Results1\n"
	    "// Caveat:\n"
	    "//   Caption lines at the top start with a #\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  qcsummary -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  PMB pmb ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&pmb, 0, sizeof (PMB)) ;
  pmb.h = h ;
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  pmb.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  pmb.gzo = getCmdLineBool (&argc, argv, "-gzo") ;

  getCmdLineOption (&argc, argv, "-run", &(pmb.run)) ;
  getCmdLineOption (&argc, argv, "-i", &(pmb.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(pmb.outFileName)) ;
  if (! pmb.run)
    {
      fprintf (stderr, "Sorry, missing parameter -run, please try postMagicBlast -help\n") ;
      exit (1) ;
    }
  if (argc > 1) usage (messprintf ("Unknown parameters %s", argv[1])) ;

  /*
    pmbParse (&pmb) ;
  */
  pmbIntrons (&pmb) ;
  
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

