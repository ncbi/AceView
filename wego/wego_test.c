/*  File: wego_test.c
 *  Author: Daniellle Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg, 2015
 *
 * Set of test program implementing in C example provided in the documentation of the GO language
 * The purpose of this code is to test the acedb 'wego' implementation of multi-threading and the 
 * acedb 'channel' communication interface
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *-------------------------------------------------------------------
 */

#include "ac.h"
#include "channel.h"
#include "wego.h"

typedef struct pStrct { 
  AC_HANDLE h ;
  const char *outFileName ;
  BOOL gzi, gzo ;
  BOOL silent ;
  int test ;
  int max_threads ;
} PP ;

/*************************************************************************************/
/*************************************************************************************/
/* Test 1
 * selfDestruct: the code should selfDestruct after a fixed time, but is revived or killed from stdin
 */
typedef struct test1Struct { CHAN *count, *countAgain, *abort ; } TEST1 ;
static void wegoTest1_countDown (void *vp)
{
  TEST1 *tt = (TEST1 *)vp ;
  BOOL b ;
  int i = 10 ;

  while (i--)
    {						
      if (channelTryGet (tt->countAgain, &b, BOOL))
	{
	  if (b)
	    i = 12 ;
	  else
	    {
	      i = -1 ;
	      channelPut (tt->count, &i, int) ;
	      break ;
	    }
	}
      channelPut (tt->count, &i, int) ;
      sleep (1) ;
    }
} /* wegoTest1_countDown */

/**********************/

static void wegoTest1_abort_from_keyBoard (void *vp)
{
  TEST1 *tt = (TEST1 *)vp ;
  char cc = 0 ;
  BOOL b = TRUE ;

  printf ("This program will self destruct in 10 seconds, type a to abort x to survive\n") ;
  while (b)
    {
      scanf ("%c", &cc) ;
      printf ("..got %c\n", cc) ;
      b = (cc == 'a' ? FALSE : TRUE) ;
      channelPut (tt->abort, &b, BOOL) ;
      if (! b)
	{
	  printf ("abort required from keyboard\n") ;
	  exit (0) ;
	}
    } 
} /* wegoTest1_abort_from_keyBoard */

/**********************/

static void wegoTest1 (PP *pp)
{
  TEST1 t ;
  int i ;
  BOOL b ; 

  t.abort = channelCreate (100, BOOL, pp->h) ;
  t.count = channelCreate (100, int, pp->h) ;
  t.countAgain = channelCreate (100, BOOL, pp->h) ;

  /* start the 2 threads */
  wego_go (wegoTest1_countDown, &t, TEST1) ;
  wego_go (wegoTest1_abort_from_keyBoard, &t, TEST1) ;

  while (1)
    {
      if (channelTryGet (t.count, &i, int))   /* GO : case i := <- count ; */
	{
	  if (i)
	    printf ("%d seconds remaining\n", i) ;
	  else
	    {
	      printf ("Counted down to zero \n") ;
	      break ;
	    }
	}
      if (channelTryGet (t.abort, &b, BOOL))  /* GO : case a := <- abort ; */
	{
	  if (! b)
	    {
	       printf ("abort required from keyboard\n") ;
	       channelPut (t.countAgain, &b, BOOL) ; /* signal countAgain to die */
	       i = 0 ;
	       while (channelGive (t.count, &i, int) != -1) {}  ;
	       break ;
	    }
	  else
	    {
	      printf ("12s revived required from keyboard\n") ;
	      channelPut (t.countAgain, &b, BOOL) ;
	    }
	}
      sleep (1) ;
    }
 
  return ;
} /* wegoTest1 */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: wego_test -test <int>  -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  wego_test -test 1 -o myFile \n"
	   "// -test : the number of the test to be executed, 99->run all tests\n"
	   "// -silent : suppress title lines and status reports from the output, just report the hits\n"
	   "// -o name : all output files will be called name.*\n"
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
  char *cp ;
  PP p ;
  AC_HANDLE h = 0 ;
  char commandBuf [1000] ;
  int ix ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&p, 0, sizeof (PP)) ;
  p.h = h ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;

  for (ix = 0, cp = commandBuf ; cp < commandBuf + 900 && ix < argc ; cp += strlen (cp), ix++)
    sprintf(cp, "%s ", argv[ix]) ;



  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  p.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  p.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);
  p.silent = getCmdLineOption (&argc, argv, "-silent", 0) ;
  p.max_threads = 8 ;
  p.max_threads = getCmdLineOption (&argc, argv, "-silent", 0) ;
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  getCmdLineInt (&argc, argv, "-test", &p.test);
  getCmdLineInt (&argc, argv, "-max_threads", &p.max_threads);

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
 
  aceInWaitAndRetry (0) ; /* does nothink, but triggers the linker on LINUX */
  wego_max_threads (p.max_threads) ;
  if (p.test == 1)
    wegoTest1 (&p) ;

  ac_free (p.h) ;
  if (!p.silent)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// done: %s max mem %dMb\n", timeShowNow (), mx) ;
    }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
