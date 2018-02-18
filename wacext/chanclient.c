#include "ac.h"
#include "remote_channel.h"
#include "wego.h"

typedef struct pStruct {int type ; } PP ;

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: chanclient --type <int>  [--max_threads <int>] +rChan_h <host_name> +rChan_p <int> +rChan_s <int> +rChan_c <int>\n"
	   "//      try: -h --help \n"
	   "// Example:   chanclient --type 1 --max_threads  +rChan_h my_host +rChan_p 123452 +rChan_s 65432 +rChan_c 0 \n"
	   "// Objective: This code is called automatically by the chanserver program to test the remote_channel library. \n"
	   "// The 4 +rChan paramaters are added automatically by the chanserver call to taskDispatch () \n"
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
  PP p ;
  CHAN *r1 = 0, *r2 = 0, *w1 = 0, *w2 = 0 ;
  TASK *task ;
  double z1 = 0, z2 = 0 ; 
  const char *wnam = 0, *rnam = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  
  printf ("Hello from chanclient\n") ;

  memset (&p, 0, sizeof(PP)) ; 
  strcpy  (commandBuf, "Hello from chanclient") ;
  if (0) usage (commandBuf, argc, argv) ;

  if (argc < 3)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;
  
  fprintf (stderr, "*** chanclient type %d start\n",  p.type) ;
  if (0) sleep (30) ;

  task = taskClientInit (&argc, argv, h) ;  /* opens the tcp connection */
  if (! task)
    usage ("missing or wrong +rChan arguments, or connetction to the server cannot be established", argc, argv)  ;
  getCmdLineInt (&argc, argv, "--type", &p.type);
  wnam = "green1" ;
  getCmdLineOption (&argc, argv, "--wnam", &wnam);
  rnam = "red1" ;
  getCmdLineOption (&argc, argv, "--rnam", &rnam);
  
  if (argc != 1)

    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }

  
  if (1) 
    {
      r1 = taskCreateReaderChannel (task, rnam, 10, double) ;
      if (! r1)
	messcrash ("client cannot create r:%s channel", rnam) ;
    }
  w1 = taskCreateWriterChannel (task, wnam, 1, int) ;
  if (! w1)
    messcrash ("client cannot create w:%s channel", wnam) ; 
  
  channelDebug (r1, TRUE, "r1") ;
  channelDebug (w1, TRUE, "w1") ;
  if (0)
    {
      int i = 7 ;  
      channelPut (w1, &i, int) ;
      
      sleep (10000) ;
      channelClose (w1) ;
      exit (0) ;
    }
  if (0)
    {
      r2 = taskCreateReaderChannel (task, "red2", 1, double) ;
      if (! r2)
	messcrash ("client cannot create r:red2 channel") ;
    }
  if (0)
    {
      w2 = taskCreateWriterChannel (task, "green2", 1, int) ;
      if (! w2)
	messcrash ("client cannot create w:green2 channel") ;
    }

  if (0) 
    {  /*  verify that we separate the 2 channels r1 should give >0 numbers and r2 negative, expect 1 2 / 3 -1 / 4 -2 / 5 -3 */
      
      channelMultiGet (r1, &z1, 1, double) ;
      channelMultiGet (r1, &z2, 1, double) ;
      fprintf (stderr, "*** chanclient type %d received %g %g\n", p.type, z1, z2) ;
      channelMultiGet (r1, &z1, 1, double) ;
      channelMultiGet (r2, &z2, 1, double) ;
      fprintf (stderr, "*** chanclient type %d received %g %g\n", p.type, z1, z2) ;
      channelMultiGet (r1, &z1, 1, double) ;
      channelMultiGet (r2, &z2, 1, double) ;
      fprintf (stderr, "*** chanclient type %d received %g %g\n", p.type, z1, z2) ;
      channelMultiGet (r1, &z1, 1, double) ;
      channelMultiGet (r2, &z2, 1, double) ;
      fprintf (stderr, "*** chanclient type %d received %g %g\n", p.type, z1, z2) ;
      
      exit (0) ;
    }

  if (0) sleep (300) ;
  while (channelMultiGet (r1, &z1, 1, double))
    {
      int nn1 = z1, nn2 ;
      z2 = -999 ;
      if (r2 && ! channelMultiGet (r2, &z2, 1, double))
	break ;
      nn2 = z2 ;

      nn1 = 10000 * nn1 +  p.type ;
      if (0)
	switch (p.type)
	  {
	  case 1: nn1 = 10000 * nn1 ; nn2 = 100 * nn2 ; break ;
	  case 2: nn1 = -10000 * nn1 ; nn2 = -100 * nn2 ; break ;
	  }
      
      if (1) fprintf (stderr, "*** chanclient type %d received %g %g\n", p.type, z1, z2) ;

      if (0) sleep (.1 * (1+p.type)) ;
      if (0) sleep (1) ;
      fprintf (stderr, "***** chanclient type %d received %g %g, returned %d : %s \n", p.type, z1, z2, nn1, timeShowNow ()) ;
      if (! channelPut (w1, &nn1, int))
	{
	  fprintf (stderr, "*** chanclient type %d A received r1=closed, exiting\n",  p.type) ;
	  return 103 ;
	  break ;
	}
      nn1 *= 100 ;
      if (w2 && ! channelPut (w2, &nn1, int))
	break ;
    }
  if (0) channelClose (w1) ;
  if (0) sleep (10) ;
  fprintf (stderr, "*** chanclient type %d B received r1=closed, exiting %s\n",  p.type, timeShowNow()) ;

  ac_free (h) ;
  if (0) sleep (1) ; 
  fprintf (stderr, "*** ok\n") ;
  return 0 ;
}


