/*
 * Test the channel system

 * Test 1:
 *  Given a ring of nodes with value x[i], i = 0,1,2,... N
 *  compute the finite average by importing K times the value u and v of the 2 nearest neighbours
 *  averaging and reexporting the updated value x[i] = 1/4 ( u + 2 * x[i] + v) to the neighbours

 * The algorithm is implemented as 2 rings of channels cu and cv linking the neighbouring nodes in up and down 
 * directions cx[i] is used to transmit x[i] to node i+1, cy[i] to node i - 1
 * The code ends when each node has been updated K times
 */

#include "ac.h"
#include "channel.h"

typedef struct nodeStruct { double x ; CHAN *cx, *cy, *cu, *cv, *done ; int k ; } NODE ;

static void update (void *vp)
{
  NODE *node = (NODE *)vp ;
  double u, v, z ;
  int i ;
  
  for (i = 0 ; i < node->k ; i++)
    {
      z = (channelGet (node->cu, &u, double) 
	   +  2 * node->x 
	   + channelGet (node->cv, &v, double)
	   )/4.0 ; 
      node->x = z ;
      channelMultiPut (node->cx, &z, 1, double) ; 
      channelMultiPut (node->cy, &z, 1, double) ; 
    }
  channelMultiPut (node->done, &(node->x), 1, double) ;
  return ;
}

static void ring (int N, int K, int debug)
{
  AC_HANDLE h = ac_new_handle () ;
  int i ;
  NODE *node, nodes[N] ;
  CHAN *done, *cu[N], *cv[N] ;

  /* create the 2 rings of communication channels */
  for (i = 0 ; i < N ; i++)
    {
      cu[i] = channelCreate (K, double, h) ;
      cv[i] = channelCreate (K, double, h) ;
      channelDebug (cu[i], debug, hprintf (h, "cu.%d", i)) ;
      channelDebug (cv[i], debug, hprintf (h, "cv.%d", i)) ;
    }
  done = channelCreate (N, double, h) ;
  /* hook the nodes to the rings */
 for (i = 0 ; i < N ; i++)
   {
     node = nodes + i ;
     node->cx = cu[i] ;
     node->cu = cv[(i+1) % N] ;
     node->cy = cv[i] ;
     node->cv = cu[(N + i - 1) % N] ;
     node->done = done ;
     node->k = K ;
   }
 /* initialize the values and export */
  for (i = 0 ; i < N ; i++)
    {
      double z =  100 * randfloat () ;
      nodes[i].x = z ;
      channelPut (cu[i], &z, double) ; 
      channelPut (cv[i], &z, double) ; 
    }

  /* report */
  printf ("Initialisation") ;
  for (i = 0 ; i < N ; i++)
    printf ("\t%f", nodes[i].x) ;
  printf ("\n") ;
 
  /* run the program in parallel */
   for (i = 0 ; i < N ; i++)
     wego_go (update, nodes + i, NODE) ;

  /* wait for the answers */
  for (i = 0 ; i < N ; i++)
    channelMultiGet (done, &(nodes[i].x) , 1, double) ;
 
   /* report */
  printf ("After K steps") ;
  for (i = 0 ; i < N ; i++)
    printf ("\t%g", nodes[i].x) ;
  printf ("\n") ;
  return ;
}

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: chantest  \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  cnv \n"
	   
	   "// --N  <int> : default 4, number of nodes\n"
	   "// --K  <int> : default 4, number of iterations\n"
 	   "// --max_threads <int>   : default N, maximal number of simultaneous threads\n"
	   "//   the algorithem will block if max_threads < N and K is larger than N/2\n"
  	   "// --verbose : all channel operations are reported\n"
	   "//\n"
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

int main (int argc, const char *argv[])
{
  int N = 4 ; 
  int K = 4 ;
  int max_threads ;
  int debug = 0 ;
  char commandBuf [4000] ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
 
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  { 
    int ix ; char *cp ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }
  
  getCmdLineInt (&argc, argv, "--N", &N) ;  
  getCmdLineInt (&argc, argv, "--K", &K) ;
  max_threads = N ;
  debug = getCmdLineOption (&argc, argv, "--verbose", 0) ;

  getCmdLineInt (&argc, argv, "--max_threads", &max_threads);
 
  if (argc > 1)
    usage (commandBuf, argc, argv) ;
  
  wego_max_threads (max_threads) ;

   ring (N, K, debug) ;
   messout ("done") ; /* needed, otherwise the linker does not nclude the mess modules */

  return 0 ;
}
