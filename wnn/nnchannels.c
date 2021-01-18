/*
###########################################################
### author mieg@ncbi.nlm.nih.gov
### created 2018_03_09
### A transposition in C of th
### Ng coursera lectures, week 4    Novembre 2017
###      Building your Deep Neural Network - Step by Step v5
### multi layers network
###########################################################
*/

#include <wnn/nn.h>

/**********************************************************************/
/**********************************************************************/

typedef struct nn_action {
  int action ;
  int b ;
  int n ;
  int layer ;
} NN_ACTION ;

/**********************************************************************/

static void nnAction (void *vp)
{
  NN *nn, *nn0 = *(NN **) vp ;
  NN_ACTION a ;
  Array nnn = nn0->nnn ;

  while (channelGet (nn0->actionChannel, &a, NN_ACTION))
    {
      switch (a.action)
	{
	case 1:
	  nn = arrp (nnn, a.n, NN) ;
	  nnNetworkForward (nn) ;
	  nnCost (nn) ;
	  channelPut (nn0->doneChannel, &a, NN_ACTION) ;
	  break ;
	case 2:
	  nn = arrp (nnn, a.n, NN) ;
	  nnNetworkPullback (nn) ;
	  channelPut (nn0->doneChannel, &a, NN_ACTION) ;
	  break ;
	case -1:
	  break ;
	}
    }
} /* nnAction */	

/**********************************************************************/

void nnTrainNetwork (NN *nn0, int nIterations)
{
  int iter ;
  int b, bMax = nn0->bMax, bestB = 0 ;
  int n, threadMax = nn0->threadMax ;
  int nMax = arrayMax (nn0->nnn) ;
  int microMax = nn0->microMax ;
  int quintille = nIterations > 100 ? nIterations / 5 : nIterations + 1 ;
  NN_ACTION a ;

  nn0->test = FALSE ;
  for (b = 0 ; b < bMax ; b++)
    {
      nn0[b].reportingStep = nIterations > 100 ? nIterations/20 : 1 ;
      nn0[b].learningRate = nn0[0].learningRate ;
      nn0[b].iterMax = nIterations ;
      nn0[b].test = FALSE ;
    }
	
  array (nn0->costs, nIterations - 1, float)  = 0 ; /* make room */

  for (iter = 0 ; iter < nIterations ; iter++)
    {
      for (b = 0 ; b < bMax ; b++)
	{
	  nn0[b].iter = iter ;
	  if (0 && iter > 8000) nn0[b].maxLearningRate = .00001 ;
	}
      switch (threadMax)
	{
	case 0:
	case 1:
	  for (n = 0 ; n < nMax ; n++)
	    {
	      nnNetworkForward (nn0+n) ;
	      nnCost (nn0+n) ;
	    }
	  bestB = nnCompareCostNetwork (nn0) ; /* best branch */
	  for (n = 0 ; n < microMax ; n++)
	    nnNetworkPullback (nn0 + bestB + n * bMax) ;
	    
	  break ;
	default:
	  for (n = 0 ; n < nMax ; n++)
	    {
	      a.action = 1 ;
	      a.n = n ;
	      a.layer = 0 ;
	      channelPut (nn0->actionChannel, &a, NN_ACTION) ;
	    }
	  for (n = 0 ; n < nMax ; n++)
	    channelGet (nn0->doneChannel, &a, NN_ACTION) ;
	  
	  bestB = nnCompareCostNetwork (nn0) ; /* best branch */
	  
	  for (n = 0 ; n < microMax ; n++)
	    {
	      a.action = 2 ;
	      a.n = bestB + n * bMax ;
	      a.layer = 0 ;
	      channelPut (nn0->actionChannel, &a, NN_ACTION) ;
	    }
	  for (n = 0 ; n < microMax ; n++)
	    channelGet (nn0->doneChannel, &a, NN_ACTION) ;
	  
	  break ; 
	}

      nnUpdateNetwork (nn0, bestB) ; /* update relative to the best branch */
      if (nn0->learningRate < 1e-6 || nn0->cost < 1e-6 )
	break ;
      
      if (0 && (iter %  quintille  == 0))
	nn0->learningRate /= 1.5 ;
    } 
  for (b = 0 ; b < bMax ; b++)
    nn0[b].iter = iter ;
    
  return ;
} /* nnTrainNetwork */

/**********************************************************************/

void nnTestNetwork (NN *nn0)
{
  int bMax = nn0->bMax = 1 ;  

  nn0->bMax = 1 ;
  nn0->test = TRUE ;
  nnNetworkForward (nn0) ;
  nnCost (nn0) ;

  nnCompareCostNetwork (nn0) ; /* cost and accuracy */

  /* restore */
  nn0->test = FALSE ;
  nn0->bMax = bMax ;
  return ;
} /* nnTestNetwork */

/**********************************************************************/
/* INITIALISATION */
/**********************************************************************/

void nnMultiThread (Array nnn, int threadMax, int bMax, int microMax)
{
  int n, b ;
  NN *nn0 = arrp (nnn, 0, NN) ;
  int nBlocs = 2 + bMax * microMax ;
  
  if (threadMax < 2)
    return ;

  nn0->actionChannel = channelCreate (nBlocs, NN_ACTION, nn0->h) ;
  nn0->doneChannel = channelCreate (nBlocs, NN_ACTION, nn0->h) ;

  for (b = 1 ; b < bMax ; b++)
    {
      nn0[b].actionChannel = nn0[0].actionChannel ;
      nn0[b].doneChannel = nn0[0].doneChannel ;
    }
  for (n = 0 ; n < threadMax ; n++)
    wego_go (nnAction, &nn0, NN*) ; /* nn1 will be copied by value */

  return ;
} /* nnMultithread  */

/**********************************************************************/
/**********************************************************************/
