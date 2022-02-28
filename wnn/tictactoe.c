
#include <wnn/nn.h>
/*
 * Define the rules and the display
 * construct a NN for predicting the score
 * train the net
 */

/*
 * who : -1/1 
 * x : -1 0 1
 */

typedef struct tttStruct { 
  AC_HANDLE h ;
  MX wins, onePos, isWin ;
  int winScore, maxWins, maxPos ;
  int size, maxThreads ;
  float *scores, *Qscores ;
  int *pp ;
  NN *NN ;
} TTT ;

/**********************************************************************/
/* show a single position */
static void tttShowPos (TTT *ttt, int ii0)
{
  int i, j, k ;
  int who = ii0 % 2 ? -1 : 1 ;
  char *cc = "0.X";
  int size = ttt->size ;
  float *s = ttt->scores ;
  float *q = ttt->Qscores ;
  
  for (k = 0 ; k < size ; k++)
    {
      if (s[ttt->size * ii0 + k] < 2)
	printf ("%.0f:%.0f\t"
		, 100*s[ttt->size * ii0 + k] 
		, 100*q[ttt->size * ii0 + k]
		) ;
      else
	printf ("\t") ;
    } 
  printf ("\n") ;
  
  for (i = 0 ; i < 3 ; i++)
    {
      for (k = 0 ; k < size ; k++)
	{
	  int *zip = ttt->pp + ttt->maxPos * (ttt->size * ii0 + k) ;
	  for (j = 0 ; j < 3 ; j++)
	    {
	      int n = zip[3*i + j] ;
	      printf("%c", who ? cc[n+1] : '.') ;
	    }
	  printf ("\t") ;
	}
      printf ("\n") ;
    }
  printf ("\n") ;
} /* tttShowPos */

/**********************************************************************/
/* show a game, i.e. a set of position till the end of the game */
static void tttShowGames (TTT *ttt)
{
  int i ;
  for (i = 0 ; i < ttt->maxPos ; i++)
    tttShowPos (ttt, i) ;
} /* tttShowGame */

/**********************************************************************/
/* construct a netwrork to predict the score */
static NN *tttCreateNN (TTT *ttt)
{
  AC_HANDLE h = ac_new_handle () ;
  NN *nn ;
  int dimIn = ttt->maxPos ;
  int dimOut = 1 ;
  int size = ttt->maxPos ;
  float x[dimIn * size] ;
  float y[dimOut * size] ;
  /*   int layerDims[] = { dimIn, 200, 40, 20, dimOut, 0 } ; */
  int layerDims[] = { dimIn, 90,dimOut, 0 } ;
  int layerActivations[] = { 0, TANH, TANH, 0 } ;
  /* int layerActivations[] = { 0, SIGMOID, SIGMOID, SIGMOID, SOFTMAX, 0 } ; */ 
  
  nn = nnInitialize (ttt->maxThreads, dimIn, dimOut, size, layerDims, layerActivations, 2) ;
  nn->learningRate = .1 ;
  nn->bMax = 1 ;
  
  memset (x, 0, sizeof(x)) ;
  memset (y, 0, sizeof(y)) ;
  nnSetX (nn, x, FALSE) ;
  nnSetY (nn, y, FALSE) ;
  nnSetX (nn, x, TRUE) ;
  nnSetY (nn, y, TRUE) ;
  
  ac_free (h) ;
  return nn ;
} /* tttCreateNN  */

/**********************************************************************/

static void tttInit (TTT  *ttt)
{
  float ww[] =
    { 1,1,1, 0,0,0, 0,0,0,
      0,0,0, 1,1,1, 0,0,0,
      0,0,0, 0,0,0, 1,1,1,
      1,0,0, 1,0,0, 1,0,0,
      0,1,0, 0,1,0, 0,1,0,
      0,0,1, 0,0,1, 0,0,1,
      1,0,0, 0,1,0, 0,0,1,
      0,0,1, 0,1,0, 1,0,0
    } ;
  
  ttt->winScore = 3 ;
  ttt->maxWins = 8 ;
  ttt->maxPos = 9 ;
  ttt->wins = mxCreate (ttt->h,  "winPos", MX_FLOAT, ttt->maxPos, ttt->maxWins, 0) ;
  ttt->onePos = mxCreate (ttt->h,  "onePos", MX_INT, ttt->maxPos, 0) ;
  ttt->NN =tttCreateNN (ttt) ;
  ttt->scores = (float*) halloc (ttt->size * ttt->maxPos * sizeof(float), ttt->h) ;
  ttt->Qscores = (float*) halloc (ttt->size * ttt->maxPos * sizeof(float), ttt->h) ;
  ttt->pp = (int*) halloc (ttt->size * ttt->maxPos  * ttt->maxPos * sizeof(int), ttt->h) ;
  
  mxSet (ttt->wins, ww) ;
  mxShow (ttt->wins) ;
  

  return ;
} /* tttInit */

/**********************************************************************/

static void tttReinit (TTT *ttt)
{
  int i, iMax = ttt->size * ttt->maxPos ;
  float *s = ttt->scores ;
  float *q = ttt->Qscores ;
  
  for (i = 0 ; i < iMax ; i++)
    q[i] = s[i] = 2 ;
  
  memset (ttt->pp, 0, ttt->maxPos * ttt->maxPos * ttt->size * sizeof(int)) ;

  return ;
} /* tttReinit */

/**********************************************************************/

static BOOL tttIsWinningPosition (TTT *ttt, int ii0, int kk)
{
  int i ;
  int who = ii0 % 2 ? -1 : 1 ;
  int winScore = ttt->winScore * who ;
  BOOL isWin = FALSE ;
  const int *zip ;
  const float *zfp ;
  
  zip = ttt->pp + ttt->maxPos * (ttt->size * ii0 + kk) ;
  
  mxSet (ttt->onePos, zip) ;
  ttt->isWin = mxContractFirstIndex (ttt->isWin, ttt->wins, ttt->onePos, ttt->h) ;
  mxValues (ttt->isWin, 0, &zfp, 0) ;
  for (i = 0 ; i < ttt->maxWins ; i++)
    if (zfp[i] == winScore)
      return TRUE ;  
  return isWin ;
} /* tttIsWinningPosition  */

/**********************************************************************/
/* tttScore: predict the score associated to this position
 * using the NN (or random) or as 100 for a winning position
 */
static float tttScore (TTT *ttt, int ii0, int kk, BOOL *donep)
{
  if (tttIsWinningPosition (ttt, ii0, kk))
    {
      float who = ii0 % 2 ? -1 : 1 ;
      *donep = TRUE ;
      return who ; 
    }
  else 
    {
      const float *yp ;
      
      nnSetIntX (ttt->NN, ttt->pp, FALSE) ;
      ttt->NN->test = FALSE ;
      nnNetworkForward (ttt->NN) ;
      nnOutValues (ttt->NN, 0, &yp, 0) ;
      return yp[ttt->size * ii0 + kk] ;
    }
} /* tttScore */

/**********************************************************************/
/* tttMove : perform a move, return FALSEis move is not allowed
 */
static BOOL tttMove (TTT *ttt, int ii0, int kk, int ii)
{
  int *zip, *zjp ;
  int ii1 = ii0 - 1 ;
  int who = ii0 % 2 ? -1 : 1 ;
  
  zjp = ttt->pp + ttt->maxPos * (ttt->size * ii0 + kk) ;
  if (ii0)
    {
      zip = ttt->pp + ttt->maxPos * (ttt->size * ii1 + kk) ;
      if (zip[ii]) /* the position is already occupied */
	return FALSE ;
      /* perform the move */
      memcpy (zjp, zip,  ttt->maxPos * sizeof(float)) ;
    }
  else
    memset (zjp, 0, ttt->maxPos * sizeof(float)) ;
  zjp[ii] = who ;
  return TRUE ;
} /* tttMove */

/**********************************************************************/
/* tttPlay : plays from current position to end of game
 */
static BOOL tttPlay (TTT *ttt, int ii0, int kk)
{
  int ii1 = ii0 + 1, iMax = ttt->maxPos ;
  int who, play, bestPlay = -1 ;
  float s, s1, bestScore = -1000 ;
  BOOL done = FALSE ;
  
  who = ii0 % 2 ? -1 : 1 ;
  
  for (play = 0 ; ! done && play < iMax ; play++)
    if (tttMove (ttt, ii0, kk, play))
      {
	s = who * tttScore (ttt, ii0, kk, &done) ;
	s1 = s + .01 * (randfloat () - .5) ;
	if (ii0 < -3 && randfloat () > .9) s1 = who * 100 ;
	if (s1 > bestScore)
	  {
	    bestScore = s ;
	    bestPlay = play ;
	  }
      }
  if (bestPlay < 0)
    {
      if (ii0)
	ttt->scores [(ii0 - 1) * ttt->size + kk] = 0 ;
    }
  else
    {
      tttMove (ttt, ii0, kk, bestPlay) ;
      ttt->scores [ii0 * ttt->size + kk] = who * bestScore ;
      if (! done)
	tttPlay (ttt, ii1, kk) ;
    }
  
  return TRUE ;
} /* tttPlay */

/**********************************************************************/
/* transfer the scores back to the lower steps of the game
 */
static void tttQscores (TTT *ttt)
{
  int ii, j ;
  int iiMax = ttt->size * ttt->maxPos ;
  float *s = ttt->scores ;
  float *q = ttt->Qscores ;
  
  for (ii = 0 ; ii < iiMax ; ii++)
    {
      float z = 1, zz = 0 ;
      q[ii] = 0 ;
      for (j = 1 ; j < ttt->maxPos ; j++)
	{
	  int k = ii + j * ttt->size ;
	  if (k < iiMax && s[k] < 2)
	    {
	      z *= .8 ; zz += z ;
	      q[ii] += z * s[k] ;
	    }
	}
      if (zz)
	q[ii] /= zz ;
      else
	q[ii] = s[ii] ;
      if (s[ii] == 2)
	q[ii] = 0 ;
    }
  return ;
} /* tttQscores */

/**********************************************************************/
/* use the Qscores as truth and train the network
 */
static void tttTrain (TTT *ttt)
{
  nnSetIntX (ttt->NN, ttt->pp, FALSE) ;
  nnSetY (ttt->NN, ttt->Qscores, FALSE) ;
  ttt->NN->test = FALSE ;
  nnTrainNetwork (ttt->NN, 1) ;
} /* tttTrain */

/**********************************************************************/

int main (int argc, const char *argv[])
{
  AC_HANDLE h = ac_new_handle () ;
  int i, test = 2, nIter = 100 ;
  TTT ttt ;
  freeinit () ; 
  messErrorInit (argv[0]) ;
  
  memset (&ttt, 0, sizeof (ttt)) ;
  ttt.size = 7 ;  
  ttt.maxThreads = 1 ;
  
  getCmdLineInt (&argc, argv, "-t", &(test)) ;
  getCmdLineInt (&argc, argv, "-n", &(nIter)) ;
  getCmdLineInt (&argc, argv, "-s", &(ttt.size)) ;
  getCmdLineInt (&argc, argv, "-th", &(ttt.maxThreads)) ;
  
  aceInWaitAndRetry  (0) ; /* triggers the linker */
  
  ttt.h = h ;
  tttInit (&ttt) ;
  for (i = 0 ; i < nIter ; i++)
    {
      int kk ;
      
      tttReinit (&ttt) ;
      for (kk = 0 ; kk < ttt.size ; kk++)
	{
	  tttPlay (&ttt, 0, kk) ;
	}
      tttQscores (&ttt) ;
      tttTrain (&ttt) ; 
    }
  tttShowGames (&ttt) ;  
  
  ac_free (h) ;
  return 0 ; 
} /* main */

/**********************************************************************/
/**********************************************************************/


