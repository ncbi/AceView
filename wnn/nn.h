
#ifndef NEURAL_NET_DEF
#define NEURAL_NET_DEF 1

#include <ac.h>
#include <matrix.h>
#include <channel.h>

/**********************************************************************/
/* DELARE STRUTURES */
/**********************************************************************/

typedef struct nnStruct {
  AC_HANDLE h ;
  const struct nnStruct *nn0 ; /* self reference, used when multithreading */
  Array costs ;
  Array layers ;
  Array nnn ;    /* array of NN, one per batch */
  MX X_train, X_test ;         /* the input data with shapes (dimIn, size, 0)  */
  MX Y_train, Y_test ;         /* the truth vector with shapes (dimOut, size, 0) */
  CHAN *doneChannel ;
  CHAN *actionChannel ;
  BOOL test ;
  int batch ;
  int dimIn ;    /* the number of values in each input vector */
  int dimOut ;   /* the number of values in each truth vector */
  int globalSize ; /* number of X_all entry vectors */
  int miniMax, miniSize ;   /* number of entry vectors in a mini batch */
  int microMax, microSize ;  /* number of entry vectors in a micro batch (micro * nthreads = mini) */
  int epoch, epochMax ;
  int iter, iterMax ;
  int reportingStep ;
  int n, b, bMax, threadMax ;
  int nGood, nBad ;
  int threadId ;
  float cost, oldCost ;
  float learningRate, ratio ;
  float maxLearningRate ;
  float W2_regularization ;
  float clipGradient ;
  double dw2 ;
} NN ;

/**********************************************************************/

typedef enum {NO_METHOD=0, FULL, CONV2D2, MAXPOOL2D2} NEURON_TYPE ;
typedef enum {NONE=0, SIGMOID, TANH, ELU, RELU, SOFT_RELU, SOFTMAX} ACTIVATION_METHOD ;

typedef struct layerStruct {
  AC_HANDLE h ;
  int dim ;   /* the number of neurons of the layer, hence A->shapes[0] */
  float learningWillingness  ;
  int updateMethod ; /* 0: basic, 1: inertia, 2: Adam */
  NEURON_TYPE type ;
  ACTIVATION_METHOD activationMethod ; 
  int nCumul ;
  MX A ;      /* the activation vector of this layer, exported upwards*/
  MX W, b ;   /* the auto adaptative parameters of the linear action of the layer */
  MX oldW, oldB ;   /* the auto adaptative parameters of the linear action of the layer */
  MX Wt ;
  MX cache ;  /* cache the derivative of the activation function */
  MX dX ;     /* the gradient of input to this layer, i.e. the activation of the lower layer */
  MX dW ;     /* the gradient of the parameters */
  MX db ;
} LAYER ;

/**********************************************************************/
/**********************************************************************/

NN *nnInitialize (int maxThreads, int dimIn, int dimOut, int size, int *layer_dims, int *layer_activations, int isIntOrComplex) ;
void nnSetX (NN *nn0, float *xx, BOOL test) ;
void nnSetIntX (NN *nn0, int *xx, BOOL test) ;
void nnSetY (NN *nn0, float *yy, BOOL test) ;

void nnMultiThread (Array nnn, int threadMax, int bMax, int microMax) ;  /* start the synchronization channels */
void nnCheckNetwork (NN *nn) ;
void nnTrainNetwork (NN *nn, int nIterations) ;
void nnTestNetwork (NN *nn0) ;
double nnCost (NN *nn) ;
float nnOutValues (NN *nn, const int **zip, const float **zfp, const complex float **zcp) ;
int nnCompareCostNetwork (NN *nn) ;
void nnUpdateNetwork (NN *nn, int bestB) ;
void nnNetworkPullback (NN *nn) ;
void nnNetworkForward (NN *nn) ;

/**********************************************************************/
/**********************************************************************/

#endif
