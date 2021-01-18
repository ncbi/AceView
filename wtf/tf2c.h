/*
 * This very simple interface allows to run in C
 * a tensorflow graph prepared in python
 * The graph can be trained and tested in C.
 * The advantages are that it becomes possible:
 *   a) To prepare in C the training data, possibly in a separate thread 
 *   b) To run the graph on C-code production data without latency.  
 *
 * This interface is derived from a short example written by Shankar
 * available at
 *    .......
 *
 * The code expects to find in the directory project
 * a tensor flow graph exported in python by the command
 *     ......
 * 
 * Adjustable parameters defined as tensorflow place-holders
 * in python, can be accessed and modified from C.
 */

#ifndef TF2C_H
#define TF2C_H

#include <tensorflow/c/c_api.h>

typedef struct tf2cStruct TF2C ;

  /* A shape is a zero terminated list of dimensions
   * read by tfCreate(project) from the file  <project>/graph
   *
   * For example a 7x12 pixels image with 3 colors
   *         to be assigned to 55 classes via softMax
   *         is represented by tensors of shapes
   *    xRank=3, xShape = {*, 7,12, 3, 0, 0, 0, 0}
   *    yRank=1, yShape = {*,55, 0, 0, 0, 0, 0, 0}
   * the first number, here denoted *, is the batch-size
   *
   * The role of the TFXY functions to be provided by the user
   * of this interface is to fill the x and y tensors with
   * float numbers, while respecting the given shapes.
   */

typedef void (TFXY) (float *xp, float *yp
		     , int64_t xShape[]
		     , int64_t yShape[]) ;


TF2C  *tf2cCreate  (const char *project) ;
void  tf2cDestroy (TF2C *tf) ;

void tf2cSetShapeX (TF2C *tf, int rank, int shape[]) ;
void tf2cSetShapeY (TF2C *tf, int rank, int shape[]) ;

float tf2cPlaceHolder (TF2C *tf, const char *name, float value) ; 

BOOL  tf2cTrain (TF2C *tf, TFXY *f, int batchSize, int nSteps) ;
float tf2cTest (TF2C *tf, TFXY *f, int batchSize) ;
float *tf2cRun (TF2C *tf, float *productionData, int batchSize) ;

int  tf2cCheckPoint (TF2C* tf) ;

/* Example:
 * Provide a TFXY fuction to fill the tensoirs as explained above. 
 * Then point to a tensor flow project directory using
     tf = tf2cCreate ("my_project") ;
 
 * tf is then passed to the training and testing operations.
 * one can train the program with 1000 steps and batch size 32
     tf2cTrain (tf, makeBatchData, 32, 1000) ;
 * save the weigths in a checkPoint
     tf2cCheckPoint (tf) ;
 * and exit.  
     tf2cDestroy (tf) ;
 *
 * In a later session
     tf = tf2cCreate ("my_project") ;
 * will recover the trained weights, and the network can be
 * run on a vector of 12 production data using
     float *results = tf2cRun (tf, productionData, 12) ;
 * where productionData respects xShape and results yShape.
 * At the end, please call 
      tf2cDestroy (tf)
 * to recover the memory
 * Any later access to tf is invalid 
 *
 * If say the learning rate is a place-holder called "lr" in python 
 * one can increase its default value by 10% by calling
      float lr = tfPlaceHolder (tf, "lr", 0) ;
      tfPlaceHolder (tf, "lr", 1.1 * lr) ;     /* not yet implemented */
 */
#endif

