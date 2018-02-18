/*  File: aql.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: aqlmain.c,v 1.8 2017/01/19 21:52:45 mieg Exp $
 * Description: top main-function for stand-alone AQL application
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 24 10:32 1999 (fw)
 * * Dec  2 16:27 1998 (edgrif): Added code to record build time of this module.
 * * Aug 17 16:53 1998 (fw): changed tableOut style to 'a' (like tace now)
 * * Aug  3 14:53 1998 (fw): accept an optional debug level parameter for -d option
 * * Jul 30 09:26 1998 (fw): check for env-var $ACEDB
 * * Jul 29 10:09 1998 (fw): reintroduced the prompt for number of line to display
 * * Jul 28 17:42 1998 (fw): use of aqlQueryExecute(), testQuery removed
 * * Jul 23 14:53 1998 (fw): removed the single value return-value thing
                             as aqlEval now always returns a table
 * * Jul 16 16:55 1998 (fw): accounted for scalar return value instead of table
 * * Jul 16 16:35 1998 (fw): slightly changed the interface,
                             global variable isDebug, activated
			     by -d option, removed -u option instead
 * Created: Wed Jul 22 11:16:44 1998 (fw)
 *-------------------------------------------------------------------
 */
/***************************************************************/

#include "aql.h"
#include "table.h"
#include "freeout.h"
#include "session.h"
#include "version.h"

/***************************************************************/

static void printUsage (void)
{
  freeOut ("\n"
	  "  This is a simple command line interface for AQL.  \n");
  freeOut ("  Type in a query (multiple lines are fine), followed by a blank line.\n"
	   "  Syntax and examples on http://www.sanger.ac.uk/Software/Acedb/AQL/\n"
	   "  Please record crashes and errors for Fred Wobus (details please).\n\n") ;

  freeOut (" Usage:  \n");
  freeOut ("       taql -h       prints this info\n\n");

  freeOut ("       taql [ options ] [<database>]\n");
  freeOut ("    Where options may include\n");
  freeOut ("           -q        quiet run - no messages (prompts etc.) to terminal\n"
	   "                     useful if the program is run as 'taql -q < queryfile.txt'\n"
	   "                     e.g. from within a script.\n\n");

  freeOut ("           -a        output in Acedb format (strings are quoted etc.)\n"
	   "           -A        output in extended Acedb format (results are quoted with their class)\n"
	   "           -j        output in Java parsable format\n"
	   "           -J        output in 'Java Style' format\n"
	   "               the default output style is plain, i.e. just text\n\n");

  freeOut ("           -param <name> <type> <value>\n"
	   "                     pass a scalar context variable\n"
	   "               name  - variable name, the query may then refer to $name\n"
	   "               type  - one of Int, Float, Text, DateType\n"
	   "               value - the scalar value, text values may have to be quoted\n\n");

  freeOut ("           -d        gives all technical debug info\n"
	   "               level 0 - no debug info\n"
	   "                     1 - internal query structure (default if -d is set)\n"
	   "                     2 - all intermediate parsetrees\n"
	   "                     3 - some tracking infor during query exection\n"
	   "                     4 - all debug info during query execution\n\n") ;

  freeOut ("    The <database> parameter is the dircetory the ACeDB database to be accessed.\n"
	   "           if no database is specified on the command line the value\n"
	   "           of the environment variable $ACEDB is enquired.\n\n");
} /* printUsage */


/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()

/***************************************************************/
/* command-line interface for AQL in a stand-alone application */
/***************************************************************/

int main (int argc, char **argv)
{
  char *buf = malloc(10000) ;
  char *line = malloc (1000) ;
  char *dbDir;
  BOOL isQuiet = FALSE ;
  char outputStyle = 'h';
  int debugLevel=0;
  int paramValueInt;
  float paramValueFloat;
  mytime_t paramValueDate;
  char *paramValueText = NULL;
  char *paramName = NULL;
  char *paramType = NULL;

  messErrorInit(argv[0]) ;	/* Record program name 
				   for crash messages*/
  freeOutInit() ;

  ++argv ; --argc ;
  if (argc && (!strcmp("-h", *argv) || !strcmp("-help", *argv)))
    {
      printUsage();
      return EXIT_FAILURE;
    }


  while (argc)
    {
      if (strcmp("-d", *argv) == 0)
	{ 
	  debugLevel = 1;		/* default if -d is given */
	  ++argv ; --argc ;

	  if (argc && (*argv[0] >= '0' && *argv[0] <= '4'))
	    {
	      int l;
	      sscanf (*argv, "%d", &l);
	      debugLevel = l;
	      ++argv; --argc;
	    }
	}
      else if (strcmp("-q", *argv) == 0)
	{ 
	  isQuiet = TRUE ;
	  putenv ("ACEDB_NO_BANNER="); /* silence the kernel as well */
	  ++argv ; --argc ;
	}
      else if (strcmp("-a", *argv) == 0)
	{
	  outputStyle = 'a';
	  ++argv ; --argc ;
	}
      else if (strcmp("-A", *argv) == 0)
	{
	  outputStyle = 'A';
	  ++argv ; --argc ;
	}
      else if (strcmp("-j", *argv) == 0)
	{
	  outputStyle = 'j';
	  ++argv ; --argc ;	  
	}
      else if (strcmp("-J", *argv) == 0)
	{
	  outputStyle = 'J';
	  ++argv ; --argc ;	  
	}
      else if (strcmp("-param", *argv) == 0 && paramName == NULL)
	{
	  ++argv ; --argc ;	  

	  if (argc < 3)
	    {
	      printUsage();
	      return EXIT_FAILURE;
	    }
	  
	  paramName = strnew (*argv, 0);
	  ++argv ; --argc ;	  

	  paramType = strnew (*argv, 0);
	  ++argv ; --argc ;	  

	  if (strncasecmp(paramType, "Int", 3) == 0)
	    sscanf (*argv, "%d", &paramValueInt);
	  else if (strcasecmp(paramType, "Float") == 0)
	    sscanf (*argv, "%f", &paramValueFloat);
	  else if (strncasecmp(paramType, "Date", 4) == 0)
	    paramValueDate = timeParse(*argv);
	  else if (strcasecmp(paramType, "Text") == 0 || strcasecmp(paramType, "String") == 0)
	    paramValueText = strnew(*argv, 0);
	  else
	    {
	      printUsage();
	      return EXIT_FAILURE;
	    }

	  ++argv ; --argc ;	  
	}
      else
	break;
    }
  /* the next command line argument could be a database directory */

  if (argc == 0)
    /* no more command line args, check for $ACEDB */
    {
      if (!(dbDir = getenv("ACEDB")))
	/* no env-var set, so show usage */
	{
	  printUsage();
	  return EXIT_FAILURE;
	}
    }
  else
    /* we do have another command line argument, 
       so we take that to be the database directory */
    dbDir = *argv;

  aceInit (dbDir) ;

  if (!isQuiet) { freeOut ("\n<AQL> "); }

  *buf = 0 ;
  while (!feof (stdin))		/* break with CTRL-D */
    { 
      *line = 0;
      
      if (fgets(line,sizeof(line),stdin) == NULL)
	break;

      if (strlen(line) > 0)
	/* we typed a line */
	{ 
	  if (*buf) strcat (buf, " ") ; 
	  strcat (buf, line) ;
	}
      /* line was blank */
      else if (strlen(buf) > 0)
	/* we already have a line in the buffer */
	{ 
	  AQL    aql;
	  TABLE *result;

	  if (paramName != NULL)
	    {
	      aql = aqlCreate2(buf, 0, debugLevel, 0, 0, 1, messprintf("$%s", paramName), 0); 

	      if (strncasecmp(paramType, "Int", 3) == 0)
		result = aqlExecute(aql, 0, 0, 0, messprintf ("%s $%s", paramType, paramName), paramValueInt, 0);
	      else if (strcasecmp(paramType, "Float") == 0)
		result = aqlExecute(aql, 0, 0, 0, messprintf ("%s $%s", paramType, paramName), paramValueFloat, 0);
	      else if (strncasecmp(paramType, "Date", 4) == 0)
		result = aqlExecute(aql, 0, 0, 0, messprintf ("%s $%s", paramType, paramName), paramValueDate, 0);
	      else if (strcasecmp(paramType, "Text") == 0 || strcasecmp(paramType, "String") == 0)
		result = aqlExecute(aql, 0, 0, 0, messprintf ("%s $%s", paramType, paramName), paramValueText, 0);
	    }
	  else
	    {
	      aql = aqlCreate2 (buf, 0, debugLevel, 0, 0, 0, 0, 0); 
	      result = aqlExecute(aql, 0, 0, 0, 0);
	    }


	  if (aqlIsError(aql))
	    {
	      freeOutf ("%s", aqlGetErrorReport(aql));
	    }
	  else
	    {
	      /* output the results table */

	      if (!isQuiet)
		freeOutf ("Result table with %d line%c\n", tableMax(result), tableMax(result) == 1 ? ' ' : 's') ;

	      if (!isQuiet && tableMax(result) >= 50)
		{
		  /* too many lines in the result table, so ask how many to display if not in quiet mode */
		  int lineNum = tableMax(result);

		  freeOut ("Number of lines to display : ");
		  scanf ("%d", &lineNum);
		  fflush (stdin);
		  freeOut ("------------------------------\n");
		  tableSliceOut (0, lineNum, result, '\t', outputStyle);
		}
	      else
		{
		  /* output the table with selected style */
		  if (!isQuiet) freeOut ("------------------------------\n");
		  tableOut (result, '\t', outputStyle) ;
		}
	      if (!isQuiet) freeOut ("------------------------------\n");
	    }
	  aqlDestroy (aql);

	  tableDestroy (result);


	  if (!isQuiet) { freeOut ("\n<AQL> ");  }
	  *buf = 0 ;
	}
      else
	/* another blank line - just another prompt */
	{
	  if (!isQuiet) { freeOut ("\n<AQL> "); }
	  *buf = 0;
	}
    }

  free(buf);
  free(line);

  /*** finish ***/

  aceQuit (FALSE);		/* clean-up without save */

  return (EXIT_SUCCESS);
} /* main */

/**********************************************************************/
