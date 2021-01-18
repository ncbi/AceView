/*
 * Create a model unsing python3 and tensorflow
 * Run it in C
 * Original code from Shankar, may 2019
 */

#include "ac.h"
#include "../wtf/tf2c.h"

TF_Buffer* ReadFile(const char* filename);
TFXY makeTrainingBatch ;
TFXY makeTestBatch ;

/***********************************************************/

int main (int argc, char** argv) 
{
  int p = 5 ;
  int batchSize = 8 ;
  int nTrainingSteps = 20000 ;
  float accuracy = 0 ;
  int xRank = 2, yRank = 1 ;
  int xDims[] = {5, 2} ;
  int yDims[] = {5} ;
  TF2C *tf = 0 ;

  if (1)
    { /* optional, if you are not using the acedb libs in your program */
      freeinit () ;
      messErrorInit ("tfxor5") ;
    }

  tf = tf2cCreate ("tfxor5") ;
  tf2cSetShapeX (tf, xRank, xDims) ;
  tf2cSetShapeY (tf, yRank, yDims) ;
  printf("Initial predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 1) ;
  printf("Initial predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 1) ;

  
  printf("acc= %f\n\nTraining for a few steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  printf("Updated predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, p*p) ;

  printf ("Saving checkpoint\n");
  tf2cCheckPoint (tf) ;

  if (1) tf2cDestroy(tf);
  return 0 ;
} /* main */

/***********************************************************/
/***********************************************************/

void makeTrainingBatch (float *xp,
		    float *yp,
		    int64_t xShape[],
		    int64_t yShape[]
		    )
 
{
  int i, j, batchSize = xShape[0] ;
  int a, b , c ;
  int nn = batchSize * sizeof (float) ;

  memset (xp, 0, 10 * nn) ;
  memset (yp, 0, 5 * nn) ;

  for (i = 0; i < batchSize; ++i) 
    {
      a = 5 * (float)rand() / (float)RAND_MAX;
      b = 5 * (float)rand() / (float)RAND_MAX;
      c = (a + b) % 5 ;
      xp[10*i + a] = 1 ;
      xp[10*i + 5 + b] = 1 ;
      yp[5*i + c] = 1 ;
    }
} /* makeTrainingBatch */

/***********************************************************/

void makeTestBatch (float *xp,
		    float *yp,
		    int64_t xShape[],
		    int64_t yShape[]
		    )
{
  int i, batchSize = xShape[0], a, b, c ;
  int nn = batchSize * sizeof (float) ;

  memset (xp, 0, 10 * nn) ;
  memset (yp, 0, 5 * nn) ;

 
  for (i = 0 ; i < batchSize ; ++i) 
    {
      a = i/5 ; b = i % 5 ; c = (a + b) % 5 ;
      xp[10*i + a] = 1 ;
      xp[10*i + 5 + b] = 1 ;
      yp[5*i + c] = 1 ;
     }
} /* makeTestBatch */

/***********************************************************/
/***********************************************************/




