#include "ac.h"
#include "remote_channel.h"
#include "wego.h"

typedef struct pStruct { int max_threads ; const char *host ; int port ; } PP ;

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{   
  int i ; 

  fprintf (stderr,
	   "// Usage: chanserver  [--max_threads <int>] \n"
	   "//      try: -h --help \n"
	   "// --max_threads n : run at most b threads in parallel, default = 1 \n"
	   "// --test integer  \n"
	   "//       0 : do not create a client: useful only when running under debugger\n"
	   "//       1 : (default value) start ../bin.$ACEDB_MACHINE/chanclient\n" 
	   "//  --local | --submit : how to launch chanclient\n"
	   "//     --local  : (default value) execute bin.$ACEDB_MACHINE/chanclient on the local machine\n"
	   "//     --submit : call  system (submit  bin.$ACEDB_MACHINE/chanclient)\n"
	   "//        the script submit should be executable, and configured for your hardware\n"
	   "//        for example, it may submit the job to an SGE farm\n"
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
#define NN 12
#define NMAX 100

typedef struct tStruct {CHAN *c ; BOOL exit ; TASK *taskS ;} TT ;
void myReader (const void *vp)
{
  int i, nn ;
  const TT *ttp = (TT *)vp ; 
  CHAN *r1 = (ttp->c) ;
  for (i = 1, nn = 0 ; i < NMAX ; i++)
    {  
      int n ;
      
      n = channelGive (r1, &n, int) ;
      nn++ ;
    
      if (nn % 10 == 1) fprintf (stdout, "***** green1 received %d, total %d/%d values\n", n, nn, NMAX) ;
    } 
  fprintf (stdout, "*** received %d values %s \n", nn, timeShowNow()) ;
  
  taskSwitchContext (ttp->taskS) ;
  if (1) sleep (15) ;
  fprintf (stdout, "*** myread exiting %s\n", timeShowNow()) ;
  if (ttp->exit) exit (0) ;
  return  ;
}

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  char commandBuf [1000] ;
  BOOL sge = FALSE ;
  PP p ;
  CHAN *r1 = 0, *w1 = 0 ;
  TASK *task0 = 0 ;
  int i, test = 1 ;
  double z ;  
  TT tt ;

  memset (&p, 0, sizeof(PP)) ;

  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  p.max_threads = 5 ;
  getCmdLineInt (&argc, argv, "--max_threads", &p.max_threads);
  getCmdLineInt (&argc, argv, "--test", &test) ;
  getCmdLineOption (&argc, argv, "--local", 0) ; /* default */
  sge = getCmdLineOption (&argc, argv, "--submit", 0) ;

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }

  wego_max_threads (p.max_threads) ;

  if (sge)
    task0 = taskServerInit (argv[0], TASK_SUBMIT, h) ;
  else
    task0 = taskServerInit (argv[0], TASK_LOCAL, h) ;

  if (0) goto done ;
 



  for (i = 0 ; i < NN ; i++)
    switch (test)
      {
      case 1:
      case 2:
	taskDispatch ("chanclient", messprintf ("--type %d  --rnam red1 --wnam green1", i), 0) ; 
	break ;
      case 0:
	taskDispatch ("echo", "hello", 0) ; 
	break ;
      }

  if (1)
    {
      r1 = taskCreateReaderChannel (task0, "green1", NMAX + 1, int) ;
      if (! r1)
	messcrash ("server cannot create r:green1 channel") ;
      
      w1 = taskCreateWriterChannel (task0, "red1", NMAX + 1, double) ;
      if (! w1)
	messcrash ("server cannot create w:red1 channel") ; 
    }
  if (1)
    {
      sleep (10) ;
      fprintf (stderr, "switching\n") ;

      channelClose (r1) ;
      channelClose (w1) ;
      taskSwitchContext (task0) ;
      for (i = 0 ; i < NN ; i++)
	taskDispatch ("chanclient", messprintf ("--type %d  --rnam red1 --wnam green1", i + 100), 0) ; 
      r1 = taskCreateReaderChannel (task0, "green1", NMAX + 1, int) ;
      w1 = taskCreateWriterChannel (task0, "red1", NMAX + 1, double) ;
    }
  if (0) channelDebug (r1, TRUE, "green1") ;



  if (0)
    {  /* test of closing a classic channel */
      z = 10 ;
      channelPut (w1, &z, double) ;
      channelPut (w1, &z, double) ;
      channelPut (w1, &z, double) ;
      taskSwitchContext (task0) ;
      channelClose (w1) ;
    }

  if (0)
    {
      sleep (300) ;
      exit (0) ;
    }
  for (i = 1 ; i < NMAX ; i++)
    { 
      z = i ;
      channelPut (w1, &z, double) ;
      if (1) fprintf (stdout, "*** red1 exported %g\n", z) ;
      z = -i ;
    }

  if (0)
    {
      fprintf (stdout, "*** closing red1\n") ;
      channelClose (w1) ;
    }

  tt.c = r1 ; tt.taskS = task0 ;
  if (1)
    {
       myReader (&tt) ;
    }
  else
    {
      tt.exit = TRUE ;
      wego_go ( myReader, &tt, TT) ;
    }
  sleep (10) ;
  fprintf (stderr, "########## switching \n") ;
  if (0) taskSwitchContext (task0) ;
   
  sleep (0) ;
  if (0) exit (0) ;

 done:
  ac_free (h) ;
  fprintf (stderr, "########## done\n") ;

  return 0 ;
}









