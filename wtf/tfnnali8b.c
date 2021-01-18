/*
 * Create a model unsing python3 and tensorflow
 * Run it in C
 * Original code from Shankar, may 2019a simple dna align predictor
 * starts with B=8 bases read + B=8 bases genome encoded as 1-hot (4x2xB =64 values) export one of 3 values (sub ins del)
 */

#include "ac.h"
#include "../wtf/tf2c.h"

TF_Buffer* ReadFile(const char* filename);
TFXY makeTrainingBatch ;
TFXY makeTestBatch ;

#define BB 8

/***********************************************************/

int main (int argc, char** argv) 
{
  int batchSize = 1024 ;
  int nTrainingSteps = 10000 ;
  float accuracy = 0 ;
  int xRank = 1, yRank = 1 ;
  int xDims[] = {8 * BB} ;
  int yDims[] = {3} ;
  TF2C *tf = 0 ;

  if (1)
    { /* optional, if you are not using the acedb libs in your program */
      freeinit () ;
      messErrorInit ("tf1") ;
    }

  tf = tf2cCreate ("nnali") ;
  tf2cSetShapeX (tf, xRank, xDims) ;
  tf2cSetShapeY (tf, yRank, yDims) ;
  printf("Initial predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 12) ;

  
  printf("acc= %f\n\nTraining for a few steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  printf("acc= %f\n\nTraining for a few more steps\n", accuracy);
  tf2cTrain (tf, makeTrainingBatch, batchSize, nTrainingSteps) ;

  printf("Updated predictions\n");
  accuracy = tf2cTest (tf, makeTestBatch, 12) ;

  printf ("Saving checkpoint\n");
  if (0) tf2cCheckPoint (tf) ;

  if (0) tf2cDestroy(tf);
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
  
  memset (xp, 0, 8 * BB * sizeof(float)) ;
  memset (yp, 0, 3 * sizeof(float)) ;
  for (i = 0 ; i < batchSize ; ++i) 
    {
      int c =  randint () % 3 ;
      yp[3*i + c] = 1 ;
      for (j = 0 ; j < BB ; j++)
	{
	  int b = randint () % 4 ;
	  xp[8*i + b] = 1 ;
	  if (j < 2)
	    {
	      xp[64*i + 8*j + 4 + b] = 1 ; /* copy the first 2 bases */
	      continue ;
	    }
	  if (j == 2 && c == 0)
	    {
	      int b1 = (b + (randint () % 3)) % 4 ;
	      xp[64*i + 8*j + 4 + b1] = 1 ; /* a substitution */
	      continue ;
	    }
	  if (c == 0)
	    {
	      xp[64*i + 8*j + 4 + b] = 1 ; /* copy the first other bases */
	      continue ;
	    }
	  if (c == 1 && j == 2)   /* insert a base */
	    {
	      int b1 = randint () % 4 ;
	      xp[64*i + 8*j + 4 + b1] = 1 ; /* an insertion */
	      continue ;
	    }
	  if (c == 1 && j > 2 && j < BB - 1)   /* insert a base */
	    {
	      xp[64*i + 8*(j+1) + 4 + b] = 1 ; /* copy behind the insert */
	      continue ;
	    }
	  if (c == 2 && j > 2)   /* delete a base */
	    {
	      xp[64*i + 8*(j-1) + 4 + b] = 1 ; /* copy one base backwards */
	      continue ;
	    }
	}
      if (c == 2) 
	{
	  int b1 = randint () % 4 ;
	  xp[64*i + 8*(BB-1) + 4 + b1] = 1 ; /* a new base at the end of the strecth */
	  continue ;
	}
    }
} /* makeTrainingBatch */

/***********************************************************/

void makeTestBatch (float *xp,
		    float *yp,
		    int64_t xShape[],
		    int64_t yShape[]
		    )
 
{
  int i, j, batchSize = xShape[0] ;
  
  memset (xp, 0, 8 * BB * sizeof(float)) ;
  memset (yp, 0, 3 * sizeof(float)) ;
  for (i = 0 ; i < batchSize ; ++i) 
    {
      int c =  randint () % 3 ;

      yp[3*i + c] = 1 ;
      for (j = 0 ; j < BB ; j++)
	{
	  int b = randint () % 4 ;
	  xp[8*i + b] = 1 ;
	  if (j < 2)
	    {
	      xp[64*i + 8*j + 4 + b] = 1 ; /* copy the first 2 bases */
	      continue ;
	    }
	  if (j == 2 && c == 0)
	    {
	      int b1 = (b + (randint () % 3)) % 4 ;
	      xp[64*i + 8*j + 4 + b1] = 1 ; /* a substitution */
	      continue ;
	    }
	  if (c == 0)
	    {
	      xp[64*i + 8*j + 4 + b] = 1 ; /* copy the first other bases */
	      continue ;
	    }
	  if (c == 1 && j == 2)   /* insert a base */
	    {
	      int b1 = randint () % 4 ;
	      xp[64*i + 8*j + 4 + b1] = 1 ; /* an insertion */
	      continue ;
	    }
	  if (c == 1 && j > 2 && j < BB - 1)   /* insert a base */
	    {
	      xp[64*i + 8*(j+1) + 4 + b] = 1 ; /* copy behind the insert */
	      continue ;
	    }
	  if (c == 2 && j > 2)   /* delete a base */
	    {
	      xp[64*i + 8*(j-1) + 4 + b] = 1 ; /* copy one base backwards */
	      continue ;
	    }
	}
      if (c == 2) 
	{
	  int b1 = randint () % 4 ;
	  xp[64*i + 8*(BB-1) + 4 + b1] = 1 ; /* a new base at the end of the strecth */
	  continue ;
	}
    }
} /* makeTrainingBatch */

/***********************************************************/
/***********************************************************/




