#include "ac.h"
#include "channel.h"

typedef struct piStruct { int N, C, sleep ; CHAN *inChan, *outChan ; } PP ;

/*************************************************************************************/
/********************************** actual work **************************************/
/* Throw N darts, check if they fall in the unit circle, return the counts */

static void piTrial (void *vp)
{
  PP *pp = (PP *)vp ;
  long s = 0 ;
  int T ;
  AC_HANDLE h = ac_new_handle () ;

  while (channelMultiGet (pp->inChan, &T, 1, int))
    {
      double x, y, z ;
      int i = pp->N, n = 0 ;
      s = rand () ;
      srandom (s) ;
      while (i--)
	{
	  if (0)
	    {
	      x = 2.0 * random ()/ RAND_MAX - 1.0 ; 
	      y = 2.0 * random ()/ RAND_MAX - 1.0 ; 
	    }
	  else
	    {
	      x = 2 * randfloat () - 1.0 ;
	      y = 2 * randfloat () - 1.0 ;
	    }
	  z = x * x + y * y ;
	  if (z < 1.0) n++ ;
	  /* printf ("%.2f %.2f %.2f %d %.2f\n", x, y, z, n, 4*(float)n/(++j)) ;  */
	}
      z = 4*n ; z /= pp->N ;  /* circle divided by square */
      channelPut (pp->outChan, &z, double) ;
      /* printf ("%.2f %.2f %d %.2f\n", x, y, n, z) ; */
      wego_log (hprintf (h,"Thead %d Trial %d \n", pp->C, T)) ;
      if (pp->sleep) sleep (pp->sleep) ;
    }
  ac_free (h) ;
  return ;
} /* piTRial */

/*************************************************************************************/
/*************************************************************************************/
static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// pi_channel_test: Compute Pi = 3.14...\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Testing the channel library by computing Pi using a distributed Monte-Carlo algorithm\n"
	    "//   In T independant trials, throw N darts in a square and count how many land in the inscribed circle\n"
	    "//     i.e. count caes where rand-x^2 + rand-y^2 < 1\n"
	    "//   The T trials are handled by C client, each in his own thread\n"
	    "//  Increasing N ad K improves the accuracy\n"
	    "//  Increasing C up to a point should improve speed, but may increase overhaed\n"
	    "//  The limit may correspond to the number of cores of the computer, or be higher\n"
	    "//\n"
	    "// Usage:\n"
	    "// pi_channel_test [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "//   -T <int> : [default 24] number of trials\n"
	    "//   -N <int> : [default 100] number of darts thrown in each trial\n"
	    "//   -C <int> : [default 4] number of calculating threads\n"
	    "//   -sleep <int> : [default 0] number of seconds each trial sleeps\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  pi_channel_test -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  PP p ;
  AC_HANDLE h = 0 ;
  int T, N, C ;
  int i ;
  double pi = 0, err, x ;
  double truePi = 3.14159265358979323846 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&p, 0, sizeof (PP)) ;
  aceInWaitAndRetry (0) ; /* does nothink, but triggers the linker on LINUX */
  /* optional arguments */
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  T = 24 ; N = 100 ; C = 4 ; /* defaults */
  getCmdLineInt (&argc, argv, "-T", &T) ;
  getCmdLineInt (&argc, argv, "-N", &N) ;
  getCmdLineInt (&argc, argv, "-C", &C) ;
  getCmdLineInt (&argc, argv, "-sleep", &p.sleep) ;

  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[argc-1])) ;
 
  fprintf (stderr, "//start: %s\n", timeShowNow()) ;
  /* initialze the channels */
  p.N = N ;
  p.inChan = channelCreate (T, int, h) ;
  p.outChan = channelCreate (T, double, h) ;

  /* launch the C claculators */ 
  wego_max_threads (C + 2) ;
  for (i = 0 ; i < C ; i++)
    { p.C = i ; wego_go (piTrial, &p, PP) ; }

  channelDebug (p.inChan, 0, "inChan") ;
  channelDebug (p.outChan, 0, "outChan") ;
  /* fill the inChannel, requesting T trials, they will be consumed by the calculators */
  for (i = 0 ; i < T ; i++)
    {
      channelPut (p.inChan, &i, int) ;
    }
  
  /* collate the results and compute the running average */
   for (i = 0 ; i < T ; i++)
     {
       x = channelGive (p.outChan, &x,double) ;
       pi = (i * pi + x) /(i+1) ; 
       /* printf ("%d %g %g\n",i, x, pi) ; */
     }
   err = pi - truePi ;
   fprintf (stderr , "# %d trials\t%d threads \t%d throws per trial \tPi %.9g \trelative error %g\n", T, C, N, pi, err/truePi) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  ac_free (h) ;

  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

