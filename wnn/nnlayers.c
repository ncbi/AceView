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
/* INITIALISATION */
/**********************************************************************/

void nnFinalize (void *va)
{
  NN *nn = (NN*) va ;

  int layer, nLayers = arrayMax (nn->layers) ;
   for (layer = 0 ; layer < nLayers ; layer++)
    {
      LAYER *ly = arrayp (nn->layers, layer, LAYER) ;
      ac_free (ly->h) ;
    }
   ac_free (nn->h) ;
}  
  
/**********************************************************************/
/* Initialze a multi-layer NN
 * nn is the newtork
 * layer_dims is a zero terminated list of positive integers,
 *   it gives the number of neurons in each layer
 *   since each neuron returns a scalar, this number is the dim 
 *   of the vector of activation values exported by the layer 
 */
static void nnInitializeOneBatch (Array nnn, int bMax, int miniMax, int microMax, int threadMax
				  , int b, int n
				  , float ratio
				  , int dimIn, int dimOut, MX_TYPE X_type				  
				  , int globalSize, int microSize
				  , int *layer_dims, int *layer_activations, BOOL isComplex)
{
  AC_HANDLE h = ac_new_handle () ;
  int layer , nLayers = 0, *dimp ;
  NN *nn = arrayp (nnn, n * bMax + b, NN) ;
  NN *nn0 = arrayp (nnn, b, NN) ;
  NN *nn1 = arrayp (nnn, n * bMax, NN) ;
  
  nn->h = h ;

  nn->nnn = nnn ;
  nn->nn0 = nn0 ;
  nn->b = b ;
  nn->n = n ;
  nn->bMax = bMax ;
  nn->ratio = ratio ;
  nn->miniMax = miniMax ;
  nn->microMax = microMax ;
  nn->threadMax = threadMax ;

  nn->globalSize = globalSize ;
  nn->microSize = microSize ; 
  nn->dimIn = dimIn ;
  nn->dimOut = dimOut ; 
  nn->learningRate = .01 ;
  nn->W2_regularization = 0.0001 ;
  nn->maxLearningRate = 5 ;
  if (b + n == 0)
    nn->costs = arrayHandleCreate (10000, double, nn->h) ; 
  if (! b)
    {
      nn->X_test = mxCreate (nn->h, "X", X_type, dimIn, microSize, 0) ;
      nn->Y_test = mxCreate (nn->h, "Y", MX_FLOAT, dimOut, microSize, 0) ;
      nn->X_train = mxCreate (nn->h, "X", X_type, dimIn, microSize, 0) ;
      nn->Y_train = mxCreate (nn->h, "Y", MX_FLOAT, dimOut, microSize, 0) ;
    }
  else
    {
      nn->X_train = nn1->X_train ;
      nn->Y_train = nn1->Y_train ;
      nn->X_test = nn1->X_test ;
      nn->Y_test = nn1->Y_test ;
    }
  /* count the layers, up to a null dimension */
  for (dimp = layer_dims  ; *dimp ; dimp++)
    nLayers ++ ;    ;

  /* initialze the layers */
  if (layer_dims[nLayers - 1] != nn->dimOut)
    messcrash ("nnInitialise, wrong geometry nn->dimOut = %d != lastLayer->dim = %d", nn->dimOut,layer_dims[nLayers - 1]) ;

    
  nn->layers = arrayHandleCreate (nLayers, LAYER, nn->h) ;
  for (layer = 0 ; layer < nLayers ; layer++)
    {
      LAYER *ly0 = arrayp (nn0->layers, layer, LAYER) ;
      LAYER *ly = arrayp (nn->layers, layer, LAYER) ;

      ly->h = ac_new_handle () ;
      ly->dim = layer_dims [layer] ; 
      ly->learningWillingness  = 1 ;
      ly->updateMethod = 0 ;
      /* by default, use Relu in hidden layers, softmax in top layer */
      ly->activationMethod = (layer < nLayers - 1 ? RELU : SOFTMAX) ;
      if (layer_activations)
	ly->activationMethod = layer_activations[layer] ;
      if (! layer)  /* the zeroth layer has no adjustable parameters */
	{ 
	  MX X = nn->nn0->test ? nn->X_test : nn->X_train ;

	  ly->A = X ;
	  continue ;
	}
      /* ATTENTION, the matrix maps ly[-1] to ly[0] and the dims are given in this order */
      if (!nn->n) /* master of the b branch */
	{
	  int w1 =  ly[-1].dim ;
	  int w2 =  ly[0].dim ;
	  BOOL useB = TRUE ;

	  switch (ly->type)
	    {
	    case CONV2D2:
	      w1 = w2 = 2 ;
	      useB = FALSE ;
	      break ;
	    default:
	      break ;
	    }
	  ly->W = mxCreate (ly->h, "W"
			    , isComplex ? MX_COMPLEX : MX_FLOAT
			    , w1, w2
			    , 0   /* zero terminate th list of dimensions */
			    ) ;
	  mxRandomInit (ly->W) ;
	  switch (ly->type)
	    {
	    case CONV2D2:
	      /* ly->Wt = mxFlipTranspose (ly->Wt, ly->W, 0, 1, ly->h) ; */
	      break ;
	    default:
	      ly->Wt = mxTranspose (ly->Wt, ly->W, 0, 1, ly->h) ;
	      break ;
	    }

	  if (useB)
	    {
	      ly->b = mxCreate (ly->h, "b"
				, isComplex ? MX_COMPLEX : MX_FLOAT
				, ly->dim
				, 0   /* zero terminate th list of dimensions */
				) ;
	      if (1)
		mxRandomInit (ly->b) ;
	    }
	}
      else
	{
	  ly->W = ly0->W ;
	  ly->Wt = ly0->Wt ;
	  ly->b = ly0->b ;
	}
      ly->dW = mxCreate (ly->h, "dW"
			, isComplex ? MX_COMPLEX : MX_FLOAT
			, ly[-1].dim, ly[0].dim
			, 0   /* zero terminate th list of dimensions */
			) ;
 
      if (ly->b)
	ly->db = mxCreate (ly->h, "db"
			   , isComplex ? MX_COMPLEX : MX_FLOAT
			   , ly->dim
			   , 0   /* zero terminate th list of dimensions */
			   ) ;

      ly->dX = mxCreate (ly->h, "dX"
			 , isComplex ? MX_COMPLEX : MX_FLOAT
			 , ly[-1].dim, nn->microSize
			 , 0   /* zero terminate th list of dimensions */
			 ) ;
       if  (ly->type == MAXPOOL2D2) /* no tunable W */
	 ly->cache = mxCreate (ly->h, "cache"
			       , isComplex ? MX_COMPLEX : MX_FLOAT
			       , ly->dim, nn->microSize
			       , 0   /* zero terminate th list of dimensions */
			       ) ;
       else
	 ly->cache = mxCreate (ly->h, "cache"
			       , isComplex ? MX_COMPLEX : MX_FLOAT
			       , ly->dim, nn->microSize
			       , 0   /* zero terminate th list of dimensions */
			       ) ;
       ly->A = mxCreate (ly->h, "A"
			 , isComplex ? MX_COMPLEX : MX_FLOAT
			 , ly->dim, nn->microSize
			 , 0   /* zero terminate th list of dimensions */
			 ) ;
    }
  return ;
} /* nnInitializeOneBatch */

/**********************************************************************/

NN *nnInitialize  (int threadMax, int dimIn, int dimOut, int size, int *layer_dims, int *layer_activations, int isComplex)
{
  AC_HANDLE h = ac_new_handle () ;
  Array nnn = 0 ;
  int n, b, bMax = 2 ;
  NN *nn0 ;
  int miniMax, miniSize, microMax, microSize ;
  float ratio = 1.4 ;
  MX_TYPE X_type = MX_FLOAT ;
  
  if (isComplex == 2)
    {
      isComplex = 0 ;
      X_type = MX_INT ;
    }

  if (threadMax < 1)
    threadMax = 1 ;

  if (size > 10000)
    {
      microSize = 256 ;
      miniSize = 256 * threadMax ;
      miniMax = size / miniSize ;
    }
  else
    {
      microSize = miniSize = size ;
      miniMax = 1 ;
    }

  microMax = size/microSize ;
  
  nnn = arrayHandleCreate (bMax * microMax, NN, h) ;

  for (n = 0 ; n < microMax ; n++)
    for (b = 0 ; b < bMax ; b++)
      {
	nnInitializeOneBatch  ( nnn
				, bMax, miniMax, microMax, threadMax
				, b, n
				, ratio
				, dimIn, dimOut, X_type
				, size, microSize
				, layer_dims, layer_activations
				, isComplex) ;
      }

  nn0 = arrp (nnn, 0, NN) ;
  nn0->clipGradient = 1.2 * nn0->ratio ;
  nnMultiThread (nnn, threadMax, bMax, microMax) ;
  
  return nn0 ;
} /* nnInitialize */

/**********************************************************************/
/**********************************************************************/

void nnSetX (NN *nn0, float *xx, BOOL test)
{
  Array nnn = nn0->nnn ;
  int n, nMax = nn0->microMax ;
  int bMax = nn0->bMax ;
  NN *nn ;

  for (n = 0 ; n < nMax ; n++)
    {
      nn = arrp (nnn, n * bMax, NN) ;
      if (test)
	mxSet (nn->X_test, xx + n * nn0->microSize * nn0->dimIn) ;
      else
	mxSet (nn->X_train, xx + n * nn0->microSize * nn0->dimIn) ;
    }
} /* nnSetX */

/**********************************************************************/

void nnSetIntX (NN *nn0, int *xx, BOOL test)
{
  Array nnn = nn0->nnn ;
  int n, nMax = nn0->microMax ;
  int bMax = nn0->bMax ;
  NN *nn ;

  for (n = 0 ; n < nMax ; n++)
    {
      nn = arrp (nnn, n * bMax, NN) ;
      if (test)
	mxSet (nn->X_test, xx + n * nn0->microSize * nn0->dimIn) ;
      else
	mxSet (nn->X_train, xx + n * nn0->microSize * nn0->dimIn) ;
    }
} /* nnSetIntX */

/**********************************************************************/

void nnSetY (NN *nn0, float *yy, BOOL test)
{
  Array nnn = nn0->nnn ;
  int n, nMax = nn0->microMax ;
  int bMax = nn0->bMax ;
  NN *nn ;

  for (n = 0 ; n < nMax ; n++)
    {
      nn = arrp (nnn, n * bMax, NN) ;
      if (test)
	mxSet (nn->Y_test, yy + n * nn0->microSize * nn0->dimOut) ;
      else
	mxSet (nn->Y_train, yy + n * nn0->microSize * nn0->dimOut) ;
    }
} /* nnSetY */

/**********************************************************************/
/**********************************************************************/
/* FORWARD ACTION */
/**********************************************************************/
/* Implement the forward action of a layer:
 *   Z = WX + b
 *   A = activation (Z)
 *      X is either nn->X the input data of the network nn->X
 *        or the activation vector of the previous layer
 *      W and b are the auto-adjustable parameters of the layer
 *      Z is the linear part of the action of a layer
 *      The activation function is a Relu, a sigmoid, a cardioid etc.
 */

static void nnLayerForward (NN *nn, int layer)
{
  LAYER *ly = arrp (nn->layers, layer, LAYER) ;
  MX X0 = nn->nn0->test ? nn->X_test : nn->X_train ;
  MX X = layer > 1 ? ly[-1].A : X0 ;
  MX (*activation) (MX zOut, MX zIn, MX cache, AC_HANDLE h) ; 

  switch (ly->type)
    {
    case MAXPOOL2D2: 
      ly->A = mxMaxPooling (ly->A, X, ly->W, ly->h) ;
      break ;
    case CONV2D2:
      ly->A = mxConvolution (ly->A, X, ly->W, ly->h) ;
      ly->A = mxAdd (ly->A, ly->A, ly->b, ly->h) ;
      break ;
    default:
      ly->A = mxDot (ly->A, ly->W, X, ly->h) ;
      ly->A = mxAdd (ly->A, ly->A, ly->b, ly->h) ;
      break ;
    }
  switch (ly->activationMethod)
    {
    case NONE:
      activation = NULL ;
      break ;
    case RELU: 
      activation = mxRelu ;
      break ;
    case ELU: 
      activation = mxElu ;
      break ;
    case SOFT_RELU: 
      activation = mxSoftRelu ;
      break ;
    case TANH: 
      activation = mxTanh ;
      break ;
    case SIGMOID:
      activation = mxSigmoid ;
      break ;
    case SOFTMAX:
      activation = mxSoftMax ;
      break ;
    }
  if (activation)
    ly->A = activation(ly->A, ly->A, ly->cache, ly->h) ;

} /* nnLayerForward */

/**********************************************************************/
/**********************************************************************/
/* Implement the forward propagation of the complete network */
void nnNetworkForward (NN *nn)
{
  int layer, max = arrayMax (nn->layers) ;

  for (layer = 1 ; layer < max ; layer++)
    nnLayerForward (nn, layer) ;      
} /* nnNetworkForward */

/**********************************************************************/
/**********************************************************************/
/* LOSS FUNCTION */
/**********************************************************************/
/* Loss function, the Shannon entropy comparing observed to predicted values
 *
 * We assume that the last layer actvation vector
 * was computed using the softmax
 * which is directly interpretable as a theoretical (predicted) probability
 *
 * The loss function is the cross Shannon cross entropy
 *        Loss = - sum_over_all_classes { p_observed   log(p_theoretical) }
 * The loss function is similar to the Lagrangian, summimg over all gauge index
 *
 * The cost function is similar to the action
 * Rather than integratiing the lagrangain over all positions in space-time
 * We sum over all training set configurations
 *        Cost = sum_over_all_examples (Loss)
 *
 * p_theoretical = ly->A, the activation. i.e. the softmax of the last layer
 * p_observed = nn->Y, the truth vector of the training set
 */

double nnCost (NN *nn)
{
  AC_HANDLE h = ac_new_handle () ;
  LAYER *ly = arrp (nn->layers, arrayMax (nn->layers) - 1 , LAYER) ;  /* last layer */
  MX Log_a1 = 0 ;
  MX Y1 = 0 ;
  MX Y = nn->nn0->test ? nn->Y_test : nn->Y_train ;
  double cost= 0 ;
 
  /* assume that a1 is the softmax or a number */
  if (ly->activationMethod == SOFTMAX)
    {  /*  a1 is the softmax, compute sum {p lpg(p) } */
      Log_a1 = ly->cache ;
      Log_a1 = ly->cache = mxLog (Log_a1, ly->A, FALSE, ly->h) ;
      Y1 = Y ;
      cost = - creal (mxFullContraction (Y1, Log_a1, &(nn->nGood), &(nn->nBad))) ; /* trace over all classes and all examples */
      cost /= nn->microSize ;
    }
  else
    { /*  a1 is continuous value, not a class label, compute sum {(X-Y)^2 } */
      Y1 = mxCreate (h, "Y1", ly->A->type, nn->dimOut, nn->microSize, 0) ;
      Y1 = mxLinearCombine (Y1, 1.0, Y, - 1.0, ly->A, h) ;   /*Y1 = Y - ly->A :  substract */
      cost = creal (mxFullContraction (Y1, Y1, 0, 0)) ; /* trace over all classes and all examples */
      cost /= nn->microSize ;
     }
  ac_free (h) ;


  nn->cost = cost ;
  return cost ;
} /* nnCost */  

/**********************************************************************/
/* return the top corner float value of the output layer, set the **z */
float nnOutValues (NN *nn, const int **zip, const float **zfp, const complex float **zcp) 
{
  LAYER *ly = arrp (nn->layers, arrayMax (nn->layers) - 1 , LAYER) ;  /* last layer */
  mxValues (ly->A, zip, zfp, zcp) ;
  
  return  (zfp && *zfp) ? **zfp : 0 ;
} /* nnOutValues */

/**********************************************************************/
/**********************************************************************/
/* UPDATE PARAMETERS */
/**********************************************************************/

static float nnLayerUpdateBasic (const NN *nn0, int b, int bestB, LAYER *ly, LAYER *ly1)
{
  float newRate, rate = nn0->learningRate ;
  float W2 = nn0->W2_regularization ;

  /* uud learning rate modifier */
  if (nn0->bMax == 2 
      && nn0->iter >= 0  /* would work at iter = 0, jump it to be compatible with tensorflow uud.py */
      )
    {
      float ratio = nn0->ratio ;

      if (1 && bestB == 1)       /* disk brakes */
	{
	  rate /= ratio ;
	  if (0 &&    /* 4 wheels disk breaks */
	      nn0[bestB].cost > 1.5 * nn0->oldCost)
	    {
	      rate /= (ratio * ratio) ;
	    }
	}
      if (b == 1)
	{
	  rate /= ratio ;
	}
      else
	{
	  rate *= ratio ;
	}
      if (rate && rate < 1e-6) 
	rate = 1e-6 ;
      if (0 && rate && rate < .01)
	rate = .01 ;
      if (rate > nn0->maxLearningRate)
	rate = nn0->maxLearningRate ;
    }

  newRate = rate ;
  
  rate = rate * ly1->learningWillingness / (nn0->microSize * nn0->microMax) ;
  if (0 && W2 > rate/10)
    W2 = rate/10 ;
  
  /* the ly1 matrices are transferred to ly
   * but the learning rate depends on the branch 
   */
  ly->W = mxLinearCombine (ly->W, 1.0 - W2, ly1->W, -rate, ly1->dW, ly->h) ;
  if (ly->b)
    ly->b = mxLinearCombine (ly->b, 1.0 - W2, ly1->b, -rate, ly1->db, ly->h) ;

  return newRate ;
}  /* nnLayerUpdateBasic */

/**********************************************************************/
/* Select the update method
 * This is where we distinguish between the simple method
 * or we add intertia or Adam terms
 */
static float nnLayerUpdateParameters (const NN *nn0, int b, int bestB, LAYER *ly, LAYER *ly1)
{
  float newRate = 0 ;
  
  switch (ly1->updateMethod)
    {
    case 0:
    default:
      newRate = nnLayerUpdateBasic (nn0, b, bestB, ly, ly1) ;
    }
  ly->Wt = mxTranspose (ly->Wt, ly->W, 0, 1, ly->h) ;
  return newRate ;
}  /* nnLayerUpdateParameters */

/**********************************************************************/

void nnUpdateNetwork (NN *nn0, int bestB)
{
  NN *nn1, *nn ;
  Array nnn = nn0->nnn ;     
  int b, bMax = nn0->bMax ;
  int n, nMax = arrayMax (nnn) ;
  int layer ;
  float newRate = nn0->learningRate ;
  
  /* in best branch cumulate for each layer the parallelized batches given by nnn */
  nn1 = nn0 + bestB ;
  /* gradient clipping */

  if (nn0->clipGradient > 0)
    {
      double dw2 = 0 ;
      double ratio = 1 ;
      
      for (layer = 1 ; layer < arrayMax (nn1->layers) ; layer++)
	{
	  LAYER *ly = arrp (nn1->layers, layer, LAYER) ;
	  dw2 += mxFNorm (ly->dW) ;
	  if (ly->db)
	    dw2 += mxFNorm (ly->db) ; 
	}
      dw2 = sqrt (dw2) * nn0->learningRate ;
      if (nn0->iter && nn0->dw2 > 0)
	ratio = dw2/nn0->dw2 ;
      nn0->dw2 = dw2 ;
      if (ratio > nn0->clipGradient)
	nn0->learningRate *= nn0->clipGradient/ratio ;
    }
      
  for (layer = 1 ; layer < arrayMax (nn0->layers) ; layer++)
    {
      LAYER *ly, *ly1 = arrp (nn1->layers, layer, LAYER) ;
      for (n = 0 ; n < nMax ; n++)
	{
	  nn = nn0 + n ;
	  if (nn->b == bestB && nn != nn1)
	    {
	      ly = arrp (nn->layers, layer, LAYER) ;

	      ly1->dW = mxAdd (ly1->dW, ly1->dW, ly->dW, ly1->h) ;
	      ly1->db = mxAdd (ly1->db, ly1->db, ly->db, ly1->h) ;
	    }
	}
      
      /* update all branches except best, since we need its original value */
      for (b = 0 ; b < bMax ; b++)
	if (b != bestB)
	  {
	    nn = nn0 + b ;
	    ly = arrp (nn->layers, layer, LAYER) ;
	    nnLayerUpdateParameters (nn0, b, bestB, ly, ly1) ;
	  }
      /* update the best branch */
      {
	float newRate2 ;
	newRate2 = nnLayerUpdateParameters (nn0, bestB, bestB, ly1, ly1) ;
	if (newRate2 > 0)
	  nn0->oldCost = nn0[bestB].cost ;
	else
	  newRate2 = - newRate2 ;
	newRate = newRate2 ;
      }
   }

  nn0->learningRate = newRate ;  
  if (0 && nn0->iter % nn0->reportingStep == 0)
    {
      LAYER *ly = arrayp (nn0->layers, 1, LAYER) ;
      fprintf (stderr, "***iter %d\tb0 %f\tb1 %f\tbest %d %f\trate = %f\tW=%f B=%f\n"
	       , nn0->iter
	       , nn0[0].cost,  nn0[1].cost, bestB, nn0[bestB].cost
	       , nn0->learningRate
	       , creal(mxTopCorner (ly->W)) 
	       , creal(mxTopCorner (ly->b)) 
	       ) ;
    }

  nn0->learningRate = newRate ;
  return ;
} /* nnUpdateNetwork */

/**********************************************************************/

int nnCompareCostNetwork (NN *nn0)
{
  int b, bestB, bMax = nn0->bMax ;
  int n, nMax = arrayMax (nn0->nnn) ;
  double bestC, cc[bMax+1] ;
  int nGood[bMax+1] ;
  int nBad[bMax+1] ;

  memset (cc, 0, sizeof (cc)) ;
  memset (nGood, 0, sizeof (nGood)) ;
  memset (nBad, 0, sizeof (nBad)) ;

  for (n = 0 ; n < nMax ; n++) 
    {
      cc[nn0[n].b] += nn0[n].cost ;
      nGood[nn0[n].b] += nn0[n].nGood ;
      nBad[nn0[n].b] += nn0[n].nBad ;
    }
  bestB = 0 ; /* nn0->iter % bMax ; */
  bestC = cc[bestB] ; 
  for (b = 0 ; b < bMax ; b++)
    if (1.01 * cc[b] < bestC /* 2018_06_09: favor acceleration by 1% of cost, otherwise XOR stalls */
	|| (b == 0 && nn0->iter == 0) /* for compatibility with TensorFlow uud.py */
	)
      { bestC = cc[b] ; bestB = b ; nn0[b].cost = cc[b] ; }

  if (nn0->test || (nn0->iter >= nn0->iterMax - 1) || (nn0->iter % nn0->reportingStep == 0))
	{
	  LAYER *ly = arrayp (nn0->layers, 1, LAYER) ;
	  if (0) 
	    fprintf (stderr, "iter %d\tb0 %f\tb1 %f\tbest %d %f\trate = %f\tnGood = %d\tnBad = %d W=%f B=%f\n"
		     , nn0->iter
		     , cc[0]/nMax, cc[1]/nMax, bestB, cc[bestB]/nMax
		     , nn0->learningRate
		     , nGood[bestB], nBad[bestB]
		     , creal(mxTopCorner (ly->W)) 
		     , creal(mxTopCorner (ly->b)) 
		     ) ;
	}

  return bestB ;
} /* nnCompareCostNetwork */

/**********************************************************************/
/**********************************************************************/
/* PULL BACK */
/**********************************************************************/
/* Implement the backward action of a layer:
 *   Z = WX + b
 *   A = activation (Z)
 *
 * During the pull back
 *   We import dA = ly[1].dX, i.e. the dX of the layer above
 *   The we compute dZ, using the derivative of the activation function
 *   which was cached during the forward action
 *
 * Then from dZ, 
 *    we compute dW and dB and immediatly update the parameters
 *    finally we compute dX, which will be used by the lower layer
 */
static void nnLayerPullback (NN *nn, int layer)
{
  AC_HANDLE h = ac_new_handle () ;
  LAYER *ly = arrp (nn->layers, layer, LAYER) ;
  MX dZ = 0, dZ2 = 0 ;
  MX dA_up = 0 ;

  /* compute dZ, the gradient of Z
   * where Z is the linear action of the layer
   *     Z = W X + b
   */
  if (layer == arrayMax (nn->layers) - 1)
    { /* Top layer, using the softmax, the gradient dZ is deduced from the truth vector
       * thanks to the magic property of the softmax formula 
       */ 
      MX Y = nn->nn0->test ? nn->Y_test : nn->Y_train ;
      dZ = mxCreate (h, "dZ", ly->A->type, nn->dimOut, nn->microSize, 0) ;
      dZ = mxLinearCombine (dZ, 1.0, ly->A, -1.0, Y, h) ;  /* dZ = ly->A - Y */
    }
  else
    {  /* Hidden layer, the dA_up gradient is inherited from the layer above */
      dA_up = ly[1].dX ;
      /* Since 
       *    A = activation (Z)
       * we have (with L the Loss fucntion, i.e. the 'Lagrangian' of the network 
       *   curly L/curly Z = cutly L/curly A * curly A/curly Z
       * or in computer notation, where d denotes the 'gradient' in nn terminology
       *   dZ = activation_prime (Z) dA
       * where activation_prime is total derivative of the activation function
       * this derivative was stored in cache when the activation was computed
       */ 
      switch (ly->activationMethod)
	{
	case 0:
	  dZ = dA_up ;
	  break ;
	case SOFTMAX:
	  if (0)
	    {dZ = mxElementWiseMultiplication (dZ, dA_up, ly->cache, h) ;
	      dZ2 = mxElementWiseMultiplication (dZ2, ly->A, ly->cache, h) ;
	      dZ = mxLinearCombine (dZ, 1, dZ, 1, dZ2, h) ;
	      break ;
	    }
	default:
	  dZ = mxElementWiseMultiplication (dZ, dA_up, ly->cache, h) ;
	  break ;
	}

    }
  
  /* compute the gradient of the parameters W b  of the linear action 
   *      Z = W X + b  => dW = X_transposed . dZ,   db = (sum) dZ
   */
  
  switch (ly->type)
    {
    case MAXPOOL2D2: /* no tunable W */
      break ;
    case CONV2D2: /* we should do something */
      /* ly->dW = mxContractLastIndex (ly->dW, ly[-1].A, dZ, 0) ; */
      break ;
    case FULL:
    default:
      ly->dW = mxContractLastIndex (ly->dW, ly[-1].A, dZ, 0) ;
      if (ly->db)
	ly->db = mxSumLastIndex (ly->db, dZ, 0) ;
      break ;
    }
  
  /* compute the gradient of the input of the current layer
   * which will be passed down and used by the lower layer
   *      Z = W X + b  => dX = W_transposed dZ
   * where X, the input of the present layer is the A output of the previous layer
   */
  if (layer > 0)
    {
      switch (ly->type)
	{
	case MAXPOOL2D2: 
	  ly->dX = mxMaxPoolingBack (ly->dX, dZ, ly[-1].A, ly->W, ly->h) ;
	  break ;
	case CONV2D2:
	  ly->dX = mxConvolution (ly->dX, ly->Wt, dZ, ly->h) ;
	  break ;
	default:
	  ly->dX = mxDot (ly->dX, ly->Wt, dZ, ly->h) ;
	  break ;
	}
    }

  /* recover the memory */
  ac_free (h) ;
  
  return ;
} /* nnLayerPullback */

/**********************************************************************/
/* Implement the backward propagation of the complete network */
void nnNetworkPullback (NN *nn)
{
  int layer, max = arrayMax (nn->layers) ;

  for (layer = max - 1 ; layer > 0 ; layer--)
    nnLayerPullback (nn, layer) ;    
  
} /* nnNetworkPullback */

/**********************************************************************/
/**********************************************************************/
/* FULL NETWORK */
/**********************************************************************/

/**********************************************************************/
/**********************************************************************/

void nnCheckNetwork (NN *nn0)
{
  LAYER *ly ;
  const int *zip ;
  const float *zfpX, *zfpA, *zfpY ;
  const complex float *zcp ;
  int i, j ;
  int n, nMax = arrayMax (nn0->nnn) ;
  
  for (n = 0 ; n < nMax ; n++)
    {  
      NN *nn = arrp (nn0->nnn, n, NN) ;
      MX X = nn->nn0->test ? nn->X_test : nn->X_train ;
      MX Y = nn->nn0->test ? nn->Y_test : nn->Y_train ;

      
      ly = arrp (nn->layers, arrayMax (nn->layers) - 1, LAYER) ;
      if (! ly->A)
	continue ;
      mxValues (ly->A, &zip, &zfpA, &zcp) ;
      mxValues (X, &zip, &zfpX, &zcp) ;
      mxValues (Y, &zip, &zfpY, &zcp) ;
      for (i = 0 ; i  + n * nn->microSize < 20 && i < nn->microSize ; i++)
	{
	  fprintf (stderr, "%d :", i + n * nn->microSize) ;
	  fprintf (stderr, "\tInput") ;
	  for (j = 0 ; j < nn->dimOut ; j++)
	    fprintf (stderr, "  %g", zfpX[nn->dimOut * i + j]) ;

	  fprintf (stderr, "\tTruth") ; 
	  for (j = 0 ; j < nn->dimOut ; j++)
	    fprintf (stderr, "  %f", zfpY[nn->dimOut * i + j]) ;
	  
	  fprintf (stderr, "\tPrediction") ;
	  for (j = 0 ; j < nn->dimOut ; j++)
	    fprintf (stderr, "  %g", zfpA[nn->dimOut * i + j]) ;
	  fprintf (stderr, "\n") ;
	}
      
    }
} /* nnCheckNetwork */

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

