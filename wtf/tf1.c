/*
 * Create a model unsing python3 and tensorflow
 *   using: tf1.py
 * Run it in C
 *   using tf1.c and the lib tf2c.[ch]
 * Original python code from Shankar, may 2019
 */

#include "ac.h"
#include "../wtf/tf2c.h"

TF_Buffer* ReadFile(const char* filename);
TFXY makeTrainingBatch ;
TFXY makeTestBatch ;

/***********************************************************/

int main (int argc, char** argv) 
{
  int batchSize = 32 ;
  int nTrainingSteps = 50 ;
  float accuracy = 0 ;
  int xRank = 2, yRank = 2 ;
  int xDims[] = {1, 1} ;
  int yDims[] = {1, 1} ;
  TF2C *tf = 0 ;

  if (1)
    { /* optional, if you are not using the acedb libs in your program */
      freeinit () ;
      messErrorInit ("tf1") ;
    }

  tf = tf2cCreate ("tf1") ;
  tf2cSetShapeX (tf, xRank, xDims) ;
  tf2cSetShapeY (tf, yRank, yDims) ;
  printf("Initial predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 3) ;
  printf("Initial predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 3) ;

  
  printf("acc= %f\n\nTraining for a few steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  if (1) tf2cCheckPoint (tf) ;

  printf("acc= %f\n\nTraining for a few more steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  if (1) tf2cCheckPoint (tf) ;

  printf("acc= %f\n\nTraining for a few more steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  printf("Updated predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 3) ;


  printf("Updated predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 3) ;

  printf ("Saving checkpoint\n");
  if (1) tf2cCheckPoint (tf) ;

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
  int i, batchSize = xShape[0] ;
  
  for (i = 0; i < batchSize; ++i) 
    {
      xp[i] = (float)rand() / (float)RAND_MAX;
      yp[i] = 3.0 * xp[i] + 2.0;
    }
} /* makeTrainingBatch */

/***********************************************************/

void makeTestBatch (float *xp,
		    float *yp,
		    int64_t xShape[],
		    int64_t yShape[]
		    )
{
  int i, batchSize = xShape[0] ;
  
  for (i = 0 ; i < batchSize ; ++i) 
    {
      xp[i] = (i + 1)/10.0 ;
      yp[i] = 3.0 * xp[i] + 2.0;
    }
} /* makeTestBatch */

/***********************************************************/
/***********************************************************/




