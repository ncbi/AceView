/*  File: gnbkserver.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:     Genbank indexing server
 *     a simple rpc daemon example

 * Runs as a daemon if phase == 5
   The server itself will disconnect after serverTimeOut seconds.
 
   If a time_out is set to 0 (zero) expiration is deactivated
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 25 23:00 1996 (mieg)
 * Created: Wed Nov 25 12:26:06 1992 (mieg)
 *-------------------------------------------------------------------
 */

 /* $Id: gnbkserver.c,v 1.2 2006/12/16 05:05:28 mieg Exp $ */

static char* masterFile = 0 ;
#define DEFAULT_PORT 0x20000300

/* Time outs */
static int serverTimeOut = 600 ;

#define NON_GRAPHIC

#include "regular.h"
#include "array.h"
#include "call.h"
#include "mytime.h"
#include <signal.h>
#include "gnbk.h"

static FILE *out = 0 ;
static BOOL isDaemon = FALSE;

/******************* INITIALISATION ***************************/
/**************************************************************/

void tStatus (void) { return ; } /* missing in the client */

void serverClose (void)
{ extern void closePortMap(void);
  
 if (!out) out = stderr ; 
  if (!isDaemon) closePortMap();
  fprintf(out,"\n\n#### Server normal exit %s\n  A bientot \n\n",  timeShowNow()) ;
  if (out && out!= stderr)
    filclose(out) ;
  gnbkClose () ;
  exit(0) ;
}

/**************************************************************/
  /* checkLife is called by alarm, 
   * but also after every query including from unauthorised clients
   */  

static   mytime_t last_t = 0 ;
static void doTimeOut (int dummy)  /* dummy:  to please signal() prototypes */
{ 
  int dt ;
  mytime_t t = timeNow () ;

  if (!out) out = stderr ; 
  if (serverTimeOut &&
      timeDiffSecs (last_t, t, &dt) &&
      10*dt > 9*serverTimeOut) 
    { fprintf(out,"// Reached timeout, server will shut down\n\n") ;
      serverClose () ;
    }
}

/****************** MAIN LOOP **********************/

Stack processQuery (char *question)
{ Stack s = 0 ;
  static int nTransactions = 0 ;
  extern Stack  gnbkQuery (char *question) ;

  alarm((unsigned int)0);

  nTransactions++ ;
  last_t = timeNow() ;	 
  s = gnbkQuery (question) ;
  fprintf (out, "Query %6d: %s %s :: %s\n",  nTransactions, 
	   timeShow (last_t), timeShowNow(), question) ;

  alarm ((unsigned int)serverTimeOut) ;
  return s ;
}

/**************************************************/
/********* ACEDB non graphic Query Server  ********/
/**************************************************/

static void openOutFile (void)
{ char buf [12*1024] ;
  strcpy (buf, masterFile) ;
  strcat (buf, ".log") ;
        
  out = filopen (buf, 0, "a") ;
  if (!out)
    out = filopen (buf, 0, "w") ;
  if (!out)
    out = stderr ;
  setbuf (out, NULL) ;
  setbuf (stdout, NULL) ;
}

static void showBanner(u_long port)
{
  fprintf(out,"**** genbank network server: Version 0.0 ****") ;
  fprintf(out,"Authors: Jean Thierry-Mieg (CNRS, France) mieg@kaa.crbm.cnrs-mop.fr\n") ;
  fprintf(out,"         Ned Lamb (CNRS, France) ned@vega.crbm.cnrs-mop.fr\n") ;
  fprintf(out,"You may redistribute this program and database subject to the\n") ;
  fprintf(out,"conditions in the accompanying copyright file.  Anyone interested in\n") ;
  fprintf(out,"maintaining an up to date version should contact one of the authors\n") ;
  fprintf(out,"at the above email addresses.\n\n") ;

  
  fprintf(out," Master file: %s\n", masterFile) ;
  fprintf(out," Using port %lu to listen for clients\n",port);
  fprintf(out," The server will wait %d seconds for clients.\n",serverTimeOut);
}

static void usage(void)
{
  messout ("Usage: gnbk {1 | 2 | 3 | 4} file [...] [-port <port> [ -timeout <t>]]") ;
  messout( "Example: gnbk 5 f1 f2 f3 -port 20000333 -timeout 600") ;
  messout( "port and timeout (in seconds) are needed only for phase 4 and 5") ;
  messout( "   Phase 1: clean the masterFile, produces masterFile.1") ;
  messout( "   Phase 2: creates the index file masterFile.idx") ;
  messout( "   Phase 3: runs as a simple c code on stdin/out") ;
  messout( "   Phase 4: runs as an rpc service on port in foreground") ;
  messout( "   Phase 5: runs as an rpc daemon") ;
  exit (1) ;
}

/* main: process the command lines, 
   initialises the system
   launches simple run if phase < 4
   runs rpc foreground if phase = 4
   runs as a daemon if pahse = 5 
*/

int main (int argc, char **argv)
{
  char dummy[15] ;
  extern void wait_for_client (u_long port, BOOL isDaemon) ;
  u_long p, port = DEFAULT_PORT;
  int phase = 0 , n, t ;
  mytime_t t0 = timeNow () ;
  extern int gnbkMax(void) ;
  Array fileNames ;

  out = stderr ; 
  if (argc < 3)
    { messout ("Sorry, you gave only %d argument(s) on the command line.", argc - 1) ;
      usage() ;
    }

  /* phase directs the way the program is used */
  if (sscanf(argv[1], "%d", &phase) != 1 || phase < 1  || phase > 5)
    usage() ;

  /* masterFile is the name of the big file to be indexed */
  fileNames = arrayCreate (argc - 2, char*) ;
  for (n = 2 ; n < argc ; n++)
    { if (!strcmp(argv[n], "-port")) break ;
    array(fileNames, n-2, char*) = argv[n] ;
    }
  if (arrayMax(fileNames) < 1)
    usage() ;
  masterFile = array(fileNames,0,char*) ;
  if (!gnbkInit (phase, fileNames))
    usage () ;
  switch (phase)
    {
    case 4: case 5: break ;
    default:  /* run a  simple C code */    
      return 0 ;
    }


  /* else we need a port */
  if (argc > n && 
      strlen(argv[n]) < 12 &&
      sscanf(argv[n++],"%s",dummy) == 1 &&
      !strcmp(dummy,"-port") &&      
      sscanf(argv[n++],"%lu",&p) == 1)
    port = p;
  else 
    usage () ;

  /* we now parse the timeouts */
  if (argc > n && 
      strlen(argv[n]) < 12 &&
      sscanf(argv[n++],"%s",dummy) == 1 &&
      !strcmp(dummy,"-timeout") &&      
      (sscanf(argv[n++],"%d",&t) == 1) &&
      t >= 0)
    serverTimeOut = t;
  
  openOutFile() ;

  if (phase == 5)
    { isDaemon = TRUE;
      fprintf(out,"\n\n#### Server forked at %s,\n####         ready at %s\n#### indexing  %d entries\n", 
	      timeShow (t0), timeShowNow(), gnbkMax() ) ;
      fprintf(out,"#### TimeOut: %d seconds\n",serverTimeOut) ;
    }
  else
      showBanner(port) ;

  last_t = timeNow() ;
  if (serverTimeOut)
    {
      signal (SIGALRM, doTimeOut) ;  
      alarm ((unsigned int)serverTimeOut) ;
    }
  wait_for_client(port, isDaemon) ;
  messcrash("gnbk network server error: main loop returned") ;
  return 1 ;
}

/**************************************************/
/**************************************************/

