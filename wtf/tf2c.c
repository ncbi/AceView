/*
 * Create a model unsing python3 and tensorflow
 * Run it in C
 * Original code from Shankar, may 2019
 */

#include "ac.h"
#include "../wtf/tf2c.h"

typedef struct tf2cStruct {
  void *magic ;
  AC_HANDLE h ;
	
  const char *nam ;
  int trainSize, testSize ;

  /* Shapes
   * example a 7x12 pixels image with 3 colors
   *         to be assigned to 55 classes via softMax
   *    xRank=3, xShape = {*, 7,12, 3, 0, 0, 0, 0}
   *    yRank=1, yShape = {*,55, 0, 0, 0, 0, 0, 0}
   * the first number is the batch-size
   */
  int xRank, yRank ;
  int64_t xShape[8], yShape[8] ;
  
  TF_Tensor *xTrain, *yTrain ;
  TF_Tensor *xTest,  *yTest, *zTest ;
  
  TF_Graph   *graph;
  TF_Session *session;
  TF_Status  *status;

  TF_Output input, target, output, loss, checkpoint ;

  TF_Operation *init_op, *train_op, *loss_op, *save_op, *restore_op ;
} TF2C ;

static int TFMAGIC = 10203040 ;

static BOOL tf2cWeightInit (TF2C *tf) ;
static TF_Tensor *ScalarStringTensor(const char *str, TF_Status *status) ;
static TF_Buffer *ReadFile (const char* filename) ;
static int ModelPredict(TF2C* tf, float* batch, int batchSize) ;

/***********************************************************/

static void tf2cCrash (TF2C *tf, const char *message)
{
  fprintf (stderr, "FATAL ERROR in tf2c (project=%s) : %s\n"
	   , tf->nam ? tf->nam : "tf->nam == NULL"
	   , message ? message : "Sorry"
	   ) ;
    
  exit (1) ; 
} /* tf2cCrash */

  /***********************************************************/
  
BOOL tf2cOK (TF2C *tf, const char *message)
{
  if (! tf->magic ||
      tf->magic != &TFMAGIC
      )
    messcrash ("tf2cOK received a bad tf, please fix the code %s", message) ;
  if (TF_GetCode(tf->status) != TF_OK) 
    {
      fprintf (stderr, "ERROR: %s\n%s\n", message, TF_Message(tf->status));
      exit (1) ;
    }
  return TRUE ;
} /* tf2cOK */

/***********************************************************/
/* Create the TF graph created in python */
static void tf2cLoadGraph (TF2C *tf, const char *graph_def_filename)
{
   TF_Graph* g = tf->graph = TF_NewGraph() ;

   TF_Buffer* graph_def = ReadFile(graph_def_filename);
   if (graph_def == NULL) tf2cCrash (tf, "tf2cLoadGraph could not acces the graph file") ;
   fprintf (stderr, "Read GraphDef of %zu bytes\n", graph_def->length);
   TF_ImportGraphDefOptions* opts = TF_NewImportGraphDefOptions();
   TF_GraphImportGraphDef(g, graph_def, opts, tf->status);
   TF_DeleteImportGraphDefOptions(opts);
   TF_DeleteBuffer(graph_def);

   tf2cOK (tf, "tf2cLoadGraph") ;
   fprintf (stderr, " tf2cLoadGraph (%s) success\n", graph_def_filename) ;
} /* tf2cLoadGraph */

/***********************************************************/
/* Create the session */
static void tf2cCreateSession (TF2C *tf)
{
  TF_SessionOptions* opts = TF_NewSessionOptions ();
  tf->session = TF_NewSession(tf->graph, opts, tf->status);
  TF_DeleteSessionOptions(opts);

  tf2cOK (tf, "tf2cCreateSession") ;
} /* tf2cCreateSession */

/***********************************************************/
/* Hook the main operators */
static void tf2cHookOperators (TF2C *tf)
{
  TF_Graph *g = tf->graph ;

  tf->input.oper = TF_GraphOperationByName(g, "input");
  tf->input.index = 0;
  tf->target.oper = TF_GraphOperationByName(g, "target");
  tf->target.index = 0;
  tf->output.oper = TF_GraphOperationByName(g, "output");
  tf->output.index = 0;
  tf->checkpoint.oper = TF_GraphOperationByName(g, "save/Const");
  tf->checkpoint.index = 0;
  
  tf->init_op = TF_GraphOperationByName(g, "init");
  tf->train_op = TF_GraphOperationByName(g, "train");
  tf->loss_op = TF_GraphOperationByName(g, "loss");
  tf->save_op = TF_GraphOperationByName(g, "save/control_dependency");
  tf->restore_op = TF_GraphOperationByName(g, "save/restore_all");

  tf2cOK (tf, "tf2cHookOperators") ;
} /* tf2cHookOperators */

/***********************************************************/

static BOOL  existsDir (const char* dirname) 
{
  struct stat buf ;
  return stat(dirname, &buf) == 0 ? TRUE : FALSE ;
}  /* existstDir */

/***********************************************************/

TF2C *tf2cCreate (const char *nam)
{
  char graph_def_filename[1024] ;
  BOOL ok ;

  TF2C *tf = (TF2C *)malloc (sizeof (TF2C)) ;
  memset (tf, 0, sizeof (TF2C)) ;
  
  tf->magic = &TFMAGIC ;
  tf->nam = nam ;
  if (! tf->nam)
    {
      tf2cCrash (tf, "tf2cCreate : missing the name of the directory holding the graph and checkpoints") ;
    }
  ok = existsDir(tf->nam) ;
  if (! ok)
    {
      tf2cCrash (tf, "the project directory is not accessible") ;
    }
  
  sprintf (graph_def_filename, "%s/graph.pb", tf->nam) ;
  tf->status = TF_NewStatus();
  tf2cLoadGraph (tf, graph_def_filename) ;
  tf2cCreateSession (tf) ; 
  tf2cHookOperators (tf) ;
  tf2cWeightInit (tf) ; /* initialize or restore the weights */
  
  
  return tf ;
} /* tf2cCreate */

/***********************************************************/

void tf2cSetShapeX (TF2C *tf, int rank, int shape[])
{
  int i, n = 1 ;
  tf->xRank = rank ;
  if (rank > 7)
    tf2cCrash (tf, "The rank of the X tensor must not exceed 7, sorry") ;
  for (i = 0 ; i < 8 ; i++)
    tf->xShape[i] = 0 ;
  for (i = 0 ; i < rank ; i++) 
    {
      tf->xShape[i+1] = shape[i] ;  /* xShape[0] reserved for the batchSize */ 
      n *= shape[i] ;
    }
  if (n <= 0)
    {
      fprintf (stderr, "tf2cSetShapeY received  a bad shape, containing dim <= 0") ;
      for (i = 0 ; i < rank ; i++)
	fprintf (stderr, "i=%d xShape=%d\n", i, shape[i]) ;
      exit (1) ;
    }
}

/***********************************************************/

void tf2cSetShapeY (TF2C *tf, int rank, int shape[])
{
  int i, n = 1 ;
  
  if (rank > 7)
    tf2cCrash (tf, "The rank of the Y tensor must not exceed 7, sorry") ;
  tf->yRank = rank ;
  for (i = 0 ; i < 8 ; i++)
    tf->yShape[i] = 0 ;
  for (i = 0 ; i < rank ; i++)
    {
      tf->yShape[i+1] = shape[i] ;  /* yShape[0] reserved for the batchSize */
      n *= shape[i] ;
    }
  if (n <= 0)
    {
      fprintf (stderr, "tf2cSetShapeY received  a bad shape, containing dim <= 0") ;
      for (i = 0 ; i < rank ; i++)
	fprintf (stderr, "i=%d yShape=%d\n", i, shape[i]) ;
      exit (1) ;
    }
} /* tf2cSetShapeY */

/***********************************************************/

void tf2cDestroy (TF2C* tf) 
{
  if (tf && tf->magic == &TFMAGIC)
    {
      if (1)
	{
	  TF_DeleteTensor (tf->xTrain) ; /* input */
	  TF_DeleteTensor (tf->yTrain) ; /* desired output */
	  
	  TF_DeleteTensor (tf->xTest) ;  /* input */
	  TF_DeleteTensor (tf->yTest) ;  /* desired output */

	}
      if (1)
	{
	  TF_DeleteGraph  (tf->graph) ;
	  TF_DeleteSession(tf->session, tf->status);
	  TF_DeleteStatus (tf->status) ;
	  
	  TF_DeleteTensor (tf->zTest) ;  /* actual output */
	}
      free (tf) ;
    }
} /* tf2cDestroy */

/***********************************************************/

static int tf2cDoCheckPoint (TF2C* tf, BOOL doSave)
{
  char fNam[1024] ;
  TF_Tensor    *t ;
  TF_Output    inputs[1] ;
  TF_Tensor    *inTensors[1] ;
  TF_Status  *status;
  const TF_Operation *ops[1] ;

  status = TF_NewStatus();
  inputs[0] = tf->checkpoint ;
  sprintf (fNam, "%s/checkPoint", tf->nam) ;
  inTensors[0] = t = ScalarStringTensor (fNam, status) ;
  ops[0] = (doSave ? tf->save_op : tf->restore_op) ;
  
  TF_SessionRun (tf->session
		 , NULL, inputs  /* inputs */
		 , inTensors, 1  /* input tensors */
		 , NULL, NULL, 0 /* outputs */
		 , ops, 1        /* operations */
		 , NULL, tf->status
		 );

  if (1) TF_DeleteTensor(t); 
  if (1) TF_DeleteStatus (status) ;

  tf2cOK (tf, "tf2cDoCheckPoint") ;

  return tf->checkpoint.index ;
} /* tf2DoCheckPoint */

/***********************************************************/

static BOOL tf2cDoRunTrainingStep (TF2C *tf)
{
  TF_Output inputs[2] = {tf->input, tf->target};
  TF_Tensor* inTensors[2] = {tf->xTrain, tf->yTrain} ;
  const TF_Operation *ops[2] ;

  ops[0] = tf->train_op ;
  ops[1] = tf->loss_op ;
  TF_SessionRun (tf->session
		 , NULL, inputs  /* inputs */
		 , inTensors, 2  /* input tensors */
		 , NULL, NULL, 0 /* outputs */
		 , ops, 2        /* operations */
		 , NULL, tf->status
		 );

  return  tf2cOK (tf, "tf2cDoRunTrainingStep") ;
}

/***********************************************************/
/***********************************************************/

static BOOL tf2cWeightDoRestore (TF2C *tf)
{
  BOOL ok = FALSE ;
  char buf[1024] ;

  fprintf (stderr, "Restoring weights from %s/checkPoint, (remove the checkpoints directory to reset)\n"
	       , tf->nam
	       ) ;
 
  sprintf (buf, "%s/checkPoint.index", tf->nam) ;
  if (filCheckName (buf, 0, "r"))
    {
     tf2cDoCheckPoint (tf,  FALSE) ;
    }
  return ok ;
} /* tf2cWeightDoRestore  */

/***********************************************************/

int tf2cCheckPoint (TF2C* tf)
{
  return tf2cDoCheckPoint (tf, TRUE) ;
} /* tf2Checkpoint */

/***********************************************************/

static BOOL tf2cWeightDoInit (TF2C *tf)
{
  const TF_Operation *ops[1] ;

  ops[0] = tf->init_op ;
  TF_SessionRun(tf->session, NULL,
                
                NULL, NULL, 0, /* No inputs */
                NULL, NULL, 0, /* No outputs */
                ops, 1, /* single operation */
                NULL, tf->status /* No metadata */
		) ;
  
  tf2cOK (tf, "tf2cWeightDoInit") ;
  fprintf (stderr, "tf2c: Initialising the weigth using %s/checkPoint.%d\n"
	   , tf->nam
	   , tf->checkpoint.index
	   ) ;
  return TRUE ;
} /* tf2cWeightDoInit */

/***********************************************************/

static BOOL tf2cWeightInit (TF2C *tf)
{
  if (! tf->nam)
    tf->nam = "Untitled graph" ;
  if (tf2cWeightDoRestore (tf))
    return TRUE ;
  if (tf2cWeightDoInit (tf))
    return TRUE ;
      
  tf2cCrash (tf, "tf2cWeightInit Cannot initialize the weights") ;
  
  return FALSE ;
} /* tf2cWeightInit */

/***********************************************************/

BOOL tf2cTrain (TF2C *tf, TFXY *xyLoad, int batchSize, int nTrainingSteps)
{
  int i ;

  tf->xShape[0] = batchSize ;
  tf->yShape[0] = batchSize ;

  if (tf->xTrain && batchSize != tf->trainSize)
    {
      TF_DeleteTensor (tf->xTrain) ; tf->xTrain = 0 ;
      TF_DeleteTensor (tf->yTrain) ; tf->yTrain = 0 ;
    }
  if (! tf->xTrain ||  batchSize != tf->trainSize)
    {
      int i ; 
      unsigned long nbytes ;
      
      tf->xShape[0] = batchSize ;
      tf->yShape[0] = batchSize ;
      
      for (i = 0, nbytes = sizeof(float) ; i < tf->xRank + 1 ; i++)
	nbytes *= tf->xShape [i] ;
      printf ("aaa xTrain %lu\n", nbytes) ;
      tf->xTrain = TF_AllocateTensor (TF_FLOAT, tf->xShape, tf->xRank+1, nbytes);
      for (i = 0, nbytes = sizeof(float) ; i < tf->yRank + 1 ; i++)
	nbytes *= tf->yShape [i] ;
      printf ("aaa xTrain %lu\n", nbytes) ;
      tf->yTrain = TF_AllocateTensor (TF_FLOAT, tf->yShape, tf->yRank+1, nbytes);

      tf->trainSize = batchSize ;
    }

  if (! xyLoad)
    tf2cCrash (tf, "Null dataLoader in call to tf2tcTrain") ;
  
  for (i = 0 ; i< nTrainingSteps ; i++)
    {
      float *xTrainp = TF_TensorData (tf->xTrain) ;
      float *yTrainp = TF_TensorData (tf->yTrain) ;
      
      xyLoad (xTrainp, yTrainp, tf->xShape, tf->yShape) ;
      tf2cDoRunTrainingStep (tf) ;
      if (10*i %  nTrainingSteps  == 0)
	fprintf (stderr, "trainig step %d\n", i) ;
    }
  return TRUE ;
} /* tf2cTrain */

/***********************************************************/

float tf2cTest (TF2C *tf, TFXY *xyLoad, int batchSize)
{
  float accuracy = 0 ;
  if (0)printf ("aaa\n") ;
  if (tf->xTest && batchSize != tf->testSize)
    {
      TF_DeleteTensor (tf->xTest) ; tf->xTest = 0 ;
      TF_DeleteTensor (tf->yTest) ; tf->yTest = 0 ;
    }
  if (0) printf ("aaa2\n") ;

  tf->xShape[0] = batchSize ;
  tf->yShape[0] = batchSize ;

  if (! tf->xTest ||  batchSize != tf->testSize)
    {
      int i ;
      unsigned long nbytes = sizeof(float) ;
      
      for (i = 0, nbytes = sizeof(float) ; i < tf->xRank + 1 ; i++)
	  nbytes *= tf->xShape [i] ;
      printf ("aaa xTest %lu\n", nbytes) ;
      tf->xTest = TF_AllocateTensor (TF_FLOAT, tf->xShape, tf->xRank+1, nbytes);
      for (i = 0, nbytes = sizeof(float) ; i < tf->yRank + 1 ; i++)
	  nbytes *= tf->yShape [i] ;
      printf ("aaa yTest %lu\n", nbytes) ;
      tf->yTest = TF_AllocateTensor (TF_FLOAT, tf->yShape, tf->yRank+1, nbytes);

      tf->testSize = batchSize ;
    }
  if (0)
    printf ("aaa3\n") ;
  if (! xyLoad)
    tf2cCrash (tf, "Null dataLoader in call to tf2tcTest") ;
  if (0)
    printf ("aaa4\n") ;
  if (1) /* run the test just once */
    {
      float *xTestp = TF_TensorData (tf->xTest) ;
      float *yTestp = TF_TensorData (tf->yTest) ;
      
      xyLoad (xTestp, yTestp, tf->xShape, tf->yShape) ;
      printf ("aaa5\n") ;   
      accuracy = ModelPredict (tf, xTestp, batchSize) ;
      printf ("aaa6\n") ;
    }
  return accuracy ;
}
/***********************************************************/

static TF_Tensor *ScalarStringTensor(const char* str, TF_Status* status) 
{
  size_t nbytes = 8 + TF_StringEncodedSize (strlen(str)) ;
  TF_Tensor *t = TF_AllocateTensor(TF_STRING, NULL, 0, nbytes) ;
  void *data = TF_TensorData(t);
  memset (data, 0, 8) ;  /* 8-byte offset of first string. */
  TF_StringEncode (str, strlen(str), data + 8, nbytes - 8, status) ;
  return t ;
} /* ScalarStringTensor */

/***********************************************************/

static TF_Buffer *ReadFile (const char* filename) 
{
  int fd = open(filename, 0);
  if (fd < 0) {
    perror("failed to open file: ");
    return NULL;
  }
  struct stat stat;
  if (fstat(fd, &stat) != 0) {
    perror("failed to read file: ");
    return NULL;
  }
  char* data = (char*)malloc(stat.st_size);
  ssize_t nread = read(fd, data, stat.st_size);
  if (nread < 0)
    {
      perror("failed to read file: ");
      free(data);
      return NULL;
    }
  if (nread != stat.st_size)
    {
      fprintf(stderr, "read %zu bytes, expected to read %ld\n"
	      , nread, stat.st_size
	      ) ;
      free(data);
      return NULL;
    }
  TF_Buffer* ret = TF_NewBufferFromString(data, stat.st_size);
  free(data);
  return ret;
} /* ReadFile */

/***********************************************************/
/***********************************************************/
/***********************************************************/

static int ModelPredict(TF2C* tf, float* batch, int batchSize)
{
  // batch consists of 1x1 matrices.
  int i ;
  unsigned long nbytes = sizeof(float);
  
  TF_Output inputs[1] = {tf->input};
  TF_Tensor* inTensors[1] ;
  TF_Output outputs[1] = {tf->output};
  TF_Tensor* outTensors[1] ;

  inTensors[0]  = tf->xTest ;
  outTensors[0] = tf->yTest ; 
  printf ("bbb0 \n") ;
  printf ("bbb0 s=%ldz\n", TF_TensorByteSize(outTensors[0])) ;   
  TF_SessionRun(tf->session, NULL
		, inputs, inTensors, 1
		, outputs, outTensors, 1
                , NULL, 0, NULL   /* no target */
		, tf->status
		);
   printf ("bbb000\n") ;   
   if (! outTensors[0])
     printf ("bbb000 ERROR no yTest\n") ;
    if (0)  tf2cOK (tf, "Modelpredict") ;
   printf ("bbb001\n") ;   

  for (i = 0 ; i < tf->yRank + 1 ; i++)
    {
      if (1) fprintf (stderr, "i=%d yShape=%d\n", i,(int) tf->yShape[i]) ;
      nbytes *= tf->yShape [i] ;
    }   
  printf ("bbb1 nbytes=%lu\n", nbytes) ;   
  if (! outTensors[0])
    printf ("bbb1 ERROR no yTest\n") ;
  if (TF_TensorByteSize(outTensors[0]) != nbytes)
    { 
      printf ("bbb2ERROR\n") ;   
      fprintf(stderr,
	      "ERROR: Expected predictions tensor to have %zu bytes, has %zu\n",
	      nbytes, TF_TensorByteSize(outTensors[0]));
      if (0) TF_DeleteTensor(outTensors[0]);
      exit (1) ;
    }
   printf ("bbb2\n") ;   
   float predictions[nbytes/sizeof(float)] ;
   memset (predictions, 0, sizeof (predictions)) ;
   memcpy (predictions, TF_TensorData(outTensors[0]), nbytes);
   printf ("bbb3\n") ;   
   printf ("Predictions:\n");
   for (i = 0; i < batchSize; ++i)
     fprintf(stdout, "\t x = %f, predicted y = %f, expected %f\n"
	     , batch[i], predictions[i]
	     , 3 * batch[i] + 2 
	     );
   if (0) TF_DeleteTensor(outTensors[0]);
   
   return 1;
}  /* ModelPredict */

