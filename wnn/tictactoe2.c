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

typedef struct tttState { 
  char buf[10], next[10] ;
  int value, delta, nextState ;
  int state, future[9], fDone ;
  int past ;
  int mm ;
} TSTATE ;

typedef struct tttStruct { 
  AC_HANDLE h ;
  MX wins, onePos, isWin ;
  short winningSignature[9] ;
  TSTATE states[1<<19] ;
  int maxState, nA, nB ;
  int winScore, maxWins, maxPos ;
  int xWins, oWins ;
  int size, maxThreads, nIter ;
  float *scores, *Qscores ;
  int *pp ;
  NN *NN ;
  DICT *dict ;
} TTT ;


/**********************************************************************/

static void tttShowState (TTT *ttt, int i0)
{
  if (i0 > 0 && i0 < ttt->maxState)
    {
      int i, j, ii ;
      TSTATE *ts = ttt->states + i0 ;
      
      fprintf(stderr, "##### State %d v=%d buf=%s next=%s\n", i0, ts->value, ts->buf, ts->next) ;
      fprintf(stderr, "\n") ;
      for (ii = i = 0 ; i < 3 ; i++)
	{
	  for (j = 0 ; j < 3 ; j++, ii++)
	    fprintf(stderr, "%c", ts->buf[ii]) ;
	  fprintf(stderr, "\n") ;
	}
    }	      
} /* tttShowState */

/**********************************************************************/

static void tttFindState (TTT *ttt, char *buf)
{
  int i0 ;
  if (buf && *buf && dictFind (ttt->dict, buf, &i0))
    tttShowState (ttt, i0) ;
} /* tttFindState */

/**********************************************************************/
/* construct the 8 winning signatures as 'short' bitsets */
static void tttGetWinningSignatures (TTT *ttt)
{
  int ww[] =
    { 1,1,1, 0,0,0, 0,0,0,
      0,0,0, 1,1,1, 0,0,0,
      0,0,0, 0,0,0, 1,1,1,
      1,0,0, 1,0,0, 1,0,0,
      0,1,0, 0,1,0, 0,1,0,
      0,0,1, 0,0,1, 0,0,1,
      1,0,0, 0,1,0, 0,0,1,
      0,0,1, 0,1,0, 1,0,0
    } ;
  int i,j ;
  
  for (i = 0 ; i < 8 ; i++)
    {
      int s = 0 ;
      for (j = 0 ; j < 9 ; j++)
	if (ww[9*i + j])
	  s |= (1<<j) ;
      ttt->winningSignature[i] = s ;
    }
  return ;
} /* */

/**********************************************************************/
/* construct the 3^9 configs and their associated data 
 * value = +100 if X wins
 * value = -100 if O wins
 * delta = 0 if X should play the next move 
 * delta = 1 if O should play the next move 
 * only licit states are constructed, up to one player wins, or board is full
 */
static int tttGetOneState (TTT *ttt, int mm)
{
  int ii, jj = 0, m ;
  TSTATE *ts = 0 ;
  DICT *dict = ttt->dict ;
  char ss[9] ;
  char buf[10] ;
  
  memset (ss, 0, sizeof(ss)) ;
  memset (buf, 0, sizeof(buf)) ;
  if (1)
    {
      char cMax = 0 ;
      int a = 0, b = 0, winner ;
      short ua = 0, ub = 0 ;
      
      ii = mm ;
      for (m = 0; m < 9 ; m++)
	{
	  char cc =ii & 0x3 ;
	  ii >>= 2 ;
	  ss[m] = cc ;
	  if (cc > cMax)
	    cMax = cc ;
	  buf[m] = '.' ;
	  if (cc == 1) { ua |= (0x1 << m) ; a++ ; buf[m] = 'X' ; }
	  if (cc == 2) { ub |= (0x1 << m) ; b++ ; buf[m] = 'O' ; }
	}
      if (cMax == 3)
	return 0 ;   /* bad config */
      for (winner = 0, m = 0 ; m < 8 ; m++)
	{
	  if ((ua & ttt->winningSignature[m]) == ttt->winningSignature[m])
	    winner |= 0x1 ;
	  if ((ub & ttt->winningSignature[m]) == ttt->winningSignature[m])
	    winner |= 0x2 ;
	}
      if (winner == 3)
	return  0 ; /* both win, game already over */
      
      if (a == b || a == b+1)
	{
	  dictAdd (dict, buf, &jj) ;
	  ts = ttt->states + jj ;
	  memcpy (ts->buf, buf, sizeof(buf)) ;
	  
	  ts->state = jj ;
	  ts->mm = mm ;
	  ts->value = 0 ;
	  ts->delta = a -  b ;
	  if (winner == 1)
	    { ttt->nA++ ; ts->value = 100 ; strcpy (ts->next, "X wins") ; }
	  if (winner == 2)
	    { ttt->nB++ ; ts->value = -100 ; strcpy (ts->next, "O wins") ; }
	}
    }
  if (jj > ttt->maxState)
    ttt->maxState = dictMax (dict) + 1 ;
  
  return ts->state ;
} /* tttGetOneState  */

static void tttGetStates (TTT *ttt, BOOL doVal)
{
  int mm ;
  int mMax = doVal ?  (1<<19) : 1 ;
  for (mm = 0 ; mm < mMax ; mm++)
    tttGetOneState ( ttt, mm) ;
  fprintf (stderr, "# constructed %d states %d A_wins %d B_wins\n", ttt->maxState, ttt->nA, ttt->nB) ;
  
  return ;
} /* tttGetStates  */

/**********************************************************************/

static void tttInit (TTT *ttt, BOOL doVal)
{
  int i ;
  
  ttt->dict = dictHandleCreate (30000, ttt->h) ;
  tttGetWinningSignatures (ttt) ;
  tttGetStates (ttt, doVal) ;
  for (i = 0 ; i < 1 ; i++)
    if ((ttt->states + i)->value == 100)
      tttShowState (ttt, i) ;
  return ;
}

/**********************************************************************/
/**********************************************************************/
/* This state is not allready winning, check if 
 * moving leads to a winning situation 
 * or if all moves lead to a losing situation
 */

static BOOL tttImproveOne (TTT *ttt, TSTATE *ts, int i0)
{
  static int nOk = 0 ;
  BOOL ok = FALSE ;
  int bad = 0 ;
  DICT *dict = ttt->dict ;
  int m, jj ;
  char *buf = ts->buf ;
  int v = ts->delta ? -100 : 100 ;    /* delta == 0, try to see if X can win */
  char cc =  ts->delta ? 'O' : 'X' ;  /* delta == 1, try to see if O can win */

  if (! ts->fDone)
    for (m = 0 ; m < 9 ; m++)
      if (buf[m] == '.')
	{
	  buf[m] = cc ;
	  if (! dictFind (dict, buf, &jj)) 
	    {
	      int k, mm = 0, i ;
	      for (i = 8 ; i >=  0 ; i--)
		{
		  k = 0 ;
		  if (buf[i] == 'X') k = 1 ;
		  if (buf[i] == 'O') k = 2 ;
		  mm <<= 2 ;
		  mm |= k ;
		}
	      ts->future[m] = tttGetOneState (ttt, mm) ;
	    }
	  buf[m] = '.' ;
	}
  ts->fDone = 1 ;

  for (m = 0 ; ! ok && m < 9 ; m++)
    if ((jj = ts->future[m]))
      {
	TSTATE *tt = ttt->states + jj ;
	if (tt->value == v)
	  {
	    ts->value = v ;
	    ts->nextState = jj ;
	    bad = 2 ; 
	    ok = TRUE ;
	    memcpy (ts->next, buf, 10) ;
	  }
	if (bad == 0 && tt->value == -v)
	  bad = 1 ;
	if (tt->value != -v)
	  bad = 2 ; /* there is a way out */
      }
  
  if (bad == 1)
    {
      ts->value = -v ;
      strcpy (ts->next, "looser") ;
      ok = TRUE ;
    }
  if (0 && v == -100 && ok && nOk++ < 8)
    tttShowState (ttt, i0) ;
  return ok ;
} /* tttImproveOne  */

/**********************************************************************/
/* If a state is not allready winning, check if 
 * moving leads to a winning situation 
 * or if all moves lead to a losing situation
 */
static void tttImprove (TTT *ttt)
{
  BOOL ok = TRUE ;
  int nn = 0 ;
  TSTATE *ts ;
  
  while (ok)
    {
      int i, iMax = ttt->maxState   ;
      ok = FALSE ;
      
      for (i = iMax - 1, ts = ttt->states + i ; i >= 0 ; i--, ts--)
	{
	  if (ts->value == 100)
	    continue ;
	  if (ts->value == -100)
	    continue ;
	  if (tttImproveOne (ttt, ts, i))
	    {
	      ok = TRUE ;
	      nn++ ;
	    }	      
	}
    }
  fprintf (stderr, "Improved %d states\n", nn) ;
  return ;
}

/**********************************************************************/

static int tttPlayOneMove (TTT *ttt, int jj)
{
  int i, ii, jj1, jj2 ;
  TSTATE *ts1, *ts = ttt->states + jj ;
  char buf[10] ;
  int v = ts->delta ? -100 : 100 ;
  float x, bestX = 0 ;
  
  if (! strcmp (ts->next, "X wins"))
    {
      ttt->xWins++ ;
      return 0 ;
    }
  if (! strcmp (ts->next, "O wins"))
    {
      ttt->oWins++ ;
      return 0 ;
    }
  
  tttImproveOne (ttt, ts, jj) ;
  if (ts->nextState)
    return ts->nextState ;
  /* otherwise select the best follower */
  strcpy (buf, ts->buf) ;
  for (jj1 = jj2 = 0, i = 0 ; i < 9 ; i++)
    if ((ii = ts->future[i]))
      {
	ts1 = ttt->states + ii ;
	jj2 = ii ;
	if (ts1->value != -v)
	  {
	    x = ts1->value ;
	    if (! jj1)
	      bestX = x ;
	    if (! ts1->past || x >= bestX)
	      {
		bestX = x ;
		jj1 = ii ; /* play a neutral move */
	      }
	  }
      }
  if (! jj1) /* play a losing move */
    {
       ts->value = -v ;
       return jj2 ;
    }
  
  return jj1 ; /* play a neutral move */
} /* tttPlayOneMove */

/**********************************************************************/

static void tttPlayOneGame (TTT *ttt, BOOL show)
{
  int jj1 = 1, jj2 ;
  TSTATE *ts ;
  
  if (show)
    fprintf (stderr, "########### GAME START\n") ;
  
  while (jj1)
    {
      jj2 = tttPlayOneMove (ttt, jj1) ;
      if (jj2)
	{
	  ts = ttt->states + jj2 ;
	  ts->past = jj1 ;
	  if(show)
	    tttShowState (ttt, jj2) ;
	}
      jj1 = jj2 ;
    }
  
  return ;
} /* tttPlayOneGame  */

/**********************************************************************/

static void tttPlay (TTT *ttt)
{
  int i, n = ttt->nIter, n1 = n/10 ;
  
  n1 = n/20 ;
  for (i = 0 ; i < n ; i++)
    {
      tttPlayOneGame (ttt, i>= n - 1 ? TRUE : FALSE) ;
      if (i + 1 == n || (i+1) % n1 == 0)
	{
	  fprintf (stderr, "After %d games wins are %d:%d\n"
		   , i+1, ttt->xWins, ttt->oWins
		   ) ;
	  ttt->xWins = ttt->oWins = 0 ;
	}
      return ;
    }
} /* tttPlay */
  
/**********************************************************************/
/**********************************************************************/

int main (int argc, const char *argv[])
{
  AC_HANDLE h = ac_new_handle () ;
  int test = 2 ; /* , nIter = 100 ; */
  TTT ttt ;
  freeinit () ; 
  messErrorInit (argv[0]) ;
  
  memset (&ttt, 0, sizeof (ttt)) ;
  ttt.size = 7 ;  
  ttt.maxThreads = 1 ;
  
  getCmdLineInt (&argc, argv, "-t", &(test)) ;
  getCmdLineInt (&argc, argv, "-n", &(ttt.nIter)) ;
  getCmdLineInt (&argc, argv, "-s", &(ttt.size)) ;
  getCmdLineInt (&argc, argv, "-th", &(ttt.maxThreads)) ;
  
  aceInWaitAndRetry  (0) ; /* triggers the linker */
  
  ttt.h = h ;
  
  if (1) /* exhaust all position before playing */
    {
      tttInit (&ttt, TRUE) ;
      tttImprove (&ttt) ;
      tttFindState (&ttt, ".........") ;
      tttFindState (&ttt, "O...XO..X") ;
    }
  else
    {
      tttInit (&ttt, FALSE) ;
      tttPlay (&ttt) ;
    }

  ac_free (h) ;
  return 0 ; 
} /* main */

/**********************************************************************/
/**********************************************************************/


