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
/**********************************************************************/

static void uudTest (void)
{
  float a, b, c, x, y, z, t, x1, y1, z1, x2, y2, z2, t1, t2, dx, dy , dz, eta, eta0, eta1, eta2, rho, rho2 ;
  float ddx, ddx1, ddy, ddy1, ddz, ddz1 ;
  int i, j, nnn[6], nn[6], nMax, n, h, hMax ;
  int d = 3 ;
  float u = 0 ;

  hMax =  1000 ;
  nMax = 1000000 ;
  eta0 = .0001 ;
  rho = 2 ;
  
  eta = eta0 ;
  rho = sqrt (rho) ;

  for (rho2 = 1 ; rho2 < 6 ; rho2 += .2)
    {  
      rho = sqrt (rho2) ;
      memset (nnn, 0, sizeof (nnn)) ;
      ddx1 = ddy1 = ddz1 = 0 ;
      for (h = 0 ; h < hMax ; h++) /* loop om n examples */
	{
	  a = 5 * randfloat ()  ;
	  b = 5 * randfloat ()  ;
	  c = 5 * randfloat ()  ;
	  x = y = z = 5 ;
	  memset (nn, 0, sizeof (nn)) ;
	  eta = eta0 ;
	  for (n = 0 ; n < nMax && ! nn[4] ; n++)
	    {
	      t = 0 ;
	      t += a * x * x /2 ;
	      t += b * y * y /2 ;
	      t += c * z * z /2 ;
	      ddx = a * x ;
	      ddy = b * y ;
	      ddz = c * z ;
	      dx = (1-u) * ddx + u * ddx1 ;
	      dy = (1-u) * ddy + u * ddy1 ;
	      dz = (1-u) * ddz + u * ddz1 ;
	      ddx1 = ddx ;
	      ddy1 = ddy ;
	      ddz1 = ddz ;

	      eta1 = eta / rho ;
	      eta2 = eta * rho ;
	      x1 = x - eta1 * dx ; 
	      y1 = y - eta1 * dy ; 
	      z1 = z - eta1 * dz ; 
	      x2 = x - eta2 * dx ; 
	      y2 = y - eta2 * dy ; 
	      z2 = z - eta2 * dz ; 
	      t1 = 0 ;
	      t1 += a * x1 * x1 /2 ;
	      if (d > 1) t1 += b * y1 * y1 /2 ;
	      if (d > 2) t1 += c * z1 * z1 /2 ;
	      t2 = 0 ;
	      t2 += a * x2 * x2 /2 ;
	      if (d > 1) t2 += b * y2 * y2 /2 ;
	      if (d > 2) t2 += c * z2 * z2 /2 ;
	      if (t1 < t2)
		{ t = t1 ; x = x1 ; y = y1 ; z = z1 ; if (eta > 1e-7) eta = eta1 ; }
	      else
		{ t = t2 ; x = x2 ; y = y2 ; z = z2 ; if (eta < 20) eta = eta2 ; }
	      for (i = 0, j = 10 ; i < 6 ; i++, j *= 10)
		if (! nn[i] && j * t < 1)
		  nn[i] = n ; 
	      if (0) printf ("n=%d\trho=%.2f\teta0=%f\teta=%f\tt=%f\n", n, rho, eta0, eta, t) ;
	    }
	  if (! nn[4]) nn[4] = n ;
	  for (i = 0 ; i < 5 ; i++)
	    nnn[i] += nn[i] ;
	}
      printf ("rho=%.2f\teta0=%f\teta=%f", rho2, eta0, eta) ;
      for (i = 0 ; i < 5 ; i++)
	printf ("\t10-%d::%d", i+1, nnn[i]/hMax) ;
      printf ("\n") ;
    }
}

/**********************************************************************/
/* APPLICATION */
/**********************************************************************/
			      
NN *nnCreateDna1()
{
  /* learning_rate = 0.10
     # Train the network a large learning rate does not work for the cplx case
   
     # n=1000 ln=2 m=400 (ln,8,4,1)  cplx at 10^-5 rapidly much better than real
     # n=1000 ln=20 m=400 (ln,8,4,1) cplx marche pas
     # n=1000 ln=2 m=40 (2,4,1) cplx relu, rate=.2 converge puis remonte
     # n=1000 ln=2 m=40 (2,4,1) cplx sigmoid, rate=.2 converge toujours
     # n=1000 ln=50 m=4000 (50,8,4,1) cplx sigmoid, rate=.1 NO CHEAT converge toujours plateau long a .36 entre 2000 et 600 iterations puis descend  vers .16 a 1000, mais pas fini
     # n=1000 ln=50 m=4000 (50,8,4,1) cplx sigmoid, rate=.1 NO CHEAT converge toujours plateau a .15 aprtir de 3000 
  */
  
  int size = 1000 ;
  NN *nn ;
  int ln = 50 ;
  int layerDims[] = { ln, 8, 4, 1, 0} ;
  nn = nnInitialize (1, 1, 1, size, layerDims, 0, FALSE) ;
  
  return nn ;
}

/* find the equation of  astraight line, understand the instability seen in tensorflow */
NN *nnCreateLineTest (int size, int maxThreads)
{
  NN *nn ;
  int b ;
  int dimIn = 1 ;
  int dimOut = 1 ;
  int layerDims[] = { dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, 0, 0 } ;

  nn = nnInitialize (maxThreads, 1, 1, size, layerDims, layerActivations, FALSE) ;
  size = 101 ;
  if (1)
    {
      int i ;
      float xx[size], yy[size] ;
     for (i = 0 ; i < size ; i++)
	{
	  xx[i] = i/100.0 ;
	  yy[i] = 1 * (30 * xx[i] + 40) ;
	}  


     nnSetX (nn, xx, FALSE) ;
     nnSetY (nn, yy, FALSE) ;
    }

  for (b = 0 ; b < nn->bMax ; b++)
    {
      float w0 = -2 ; float b0 = 0 ;
      LAYER *ly = arrp (nn[b].layers, 1, LAYER) ;
      mxSet (ly->W, &w0) ;
      mxSet (ly->b, &b0) ; 
      nn[b].W2_regularization = 0 ;
    } 
  nn->learningRate = 0.1 ; /* .34 ; */
  return nn ;
}

/* try to discover the equation y = 3x + 4
 * stops when the quadratic error is low enough 10^-6 so the coeff are 10^-3 off
 * works well with uud method: bMax =2, very slow with bMax = 1
 */
NN *nnCreateNumericTest1 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 1 ;
  int dimOut = 1 ;
  int layerDims[] = { dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, 0, 0 } ;

  nn = nnInitialize (maxThreads, 1, 1, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .01 ;
  size = (9+size)/10 ; size *= 10 ;
  if (1)
    {
      int i ;
      float xx[size], yy[size] ; 
      
      for (i = 0 ; i < 10 ; i++)
	{
	  xx[i] = 1  * randfloat () ;
	  yy[i] = 3 * xx[i] + 4 ;
	}  

      if (size > 10 && (size % 10) == 0)
	for (i = 1 ; i < size/10 ; i++)
	  {
	    memcpy (xx + 10 * i * dimIn, xx, 10 * dimIn * sizeof (float)) ;
	    memcpy (yy + 10 * i * dimOut, yy, 10 * dimOut * sizeof (float)) ;
	  }

      nnSetX (nn, xx, FALSE) ;
      nnSetY (nn, yy, FALSE) ;
    }

  return nn ;
}

/* ATTENTION , if for each dimIn sum(xi) != 0, then bi diverges and the NN does not converge
 * there is no issue with sum (yi) which is taken care of by bj
*/

/* try to discover the equation y = u + 4 * v - 6.5 * w + 2
 * stops when the quadratic error is low enough 10^-6 so the coeff are 10^-3 off
 * works well with uud method: bMax =2, very slow with bMax = 1
 */
NN *nnCreateNumericTest2 (int size, int maxThreads)
{
  NN *nn ;
  int i ;
  int dimIn = 3 ;
  int dimOut =1 ;
  int layerDims[] = { dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, 0, 0 } ;
  float x[size*dimIn] ;
  float u, v, w, y[size*dimOut] ;
  
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .01 ;

  memset (x, 0, sizeof(x)) ;
  memset (y, 0, sizeof(y)) ;
  for (i = 0 ; i < size ; i++)
    {
      u = x[0 + dimIn * i] = 2 * (randfloat () - .5) ;
      v = x[1 + dimIn * i] = 2 * (randfloat () - .5) ;
      w = x[2 + dimIn * i] = 2 * (randfloat () - .5) ;

      y[i] = u + 4*v - 6.5 * w + 2 ;
    }
  
  nnSetX (nn, x, FALSE) ;
  nnSetY (nn, y, FALSE) ;

  return nn ;
}

NN *nnCreateNumericTest3 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 2 ;
  int layerDims[] = { dimIn, 1, 0 } ;
  int layerActivations[] = { 0, 0, 0 } ;

  size = 4 ;
  nn = nnInitialize (maxThreads, dimIn, 1, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .01 ;

  if (1)
    {   /* hard code  4 training vectors : { 1,1} and {1,-1} with truth x1 + x2 
	 * if we impose to strat with w1,w2,b = 1/2 we are bad and stable ?
	 */
      int i ;
      float x[] = {1,1,1,-1,1,1,1,-1} ;
      float y[] = {2,0,2,0};   /* x1 + x2 */

      if (1)
	for (i = 0 ; i < 8 ; i++)
	  x[i] += 0.0001 * randfloat() ;
      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }
  if (0)
    {
      LAYER *ly = arrp (nn->layers, 1, LAYER) ;
      float w[] = {.5, .5,.5, .5};       
      mxSet (ly->W, w) ;
      if (1) mxSet (ly->b, w) ;
    }

  return nn ;
}
NN *nnCreateNumericTest4 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 2 ;
  int layerDims[] = { dimIn, 1, 0 } ;
  int layerActivations[] = { 0, 0, 0 } ;
  size = 4 ;
  
  nn = nnInitialize (maxThreads, dimIn, 1, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .1 ;

  if (1)
    {   /* hard code  4 training vectors : { 1,1} and {2,-2} with truth x1 + x2 */
      float x[] = {1,1,2,-2,-1,-2,-2,3} ;
      float y[] = {2,0,-3,1};   /* x1 + x2 */

      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }
  if (1)
    {
      LAYER *ly = arrp (nn->layers, 1, LAYER) ;
      float w[] = {.5, .5, .5, .5};       
      mxSet (ly->W, w) ;
      if (0) mxSet (ly->b, w) ;
    }

  return nn ;
}


NN *nnCreateNumericTest5 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 2 ;
  int layerDims[] = { dimIn, 1, 0, 0 } ;
  int layerActivations[] = { 0, 0, 0, 0 } ;

  nn = nnInitialize (maxThreads, dimIn, 1, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .1 ;

 if (1)
    {
      int i, j ;
      float s, x[size * dimIn], y[size] ; 

      memset (x, 0, sizeof(x)) ;
      memset (y, 0, sizeof(y)) ;
      for (i = 0 ; i < size ; i++)
	{
	  s = 30 ;
	  for (j = 0 ; j < dimIn ; j++)
	    {
	      x[j + dimIn * i] = 1 * (randfloat () - .5) ;
	      s += (j+1) * x[j + dimIn * i] ;
	    }
	  y[i] = s ;
	}
      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }

 
  return nn ;
}

/* 2 hidden layers */

static void nnNormalizeY (float *y, int iMax, int jMax)
{
  double Y[iMax], YY[iMax] ;
  int i, j ;
  
  memset (Y, 0, sizeof (Y)) ;
  memset (YY, 0, sizeof (YY)) ;
  /* compute Y = <y> */
  for (i = 0 ; i < iMax ; i++)
    for (j = 0 ; j < jMax ; j++)
      Y[i] += y[i + iMax * j] ;
  for (i = 0 ; i < iMax ; i++)
    Y[i] /= jMax ;
  
  /* centralize y: substract the average */
  for (i = 0 ; i < iMax ; i++)
    for (j = 0 ; j < jMax ; j++)
      y[i + iMax * j] -= Y[i] ;

  /* YY = sum (y ^ 2) */
  for (i = 0 ; i < iMax ; i++)
    for (j = 0 ; j < jMax ; j++)
      YY[i] += y[i + iMax * j] * y[i + iMax * j] ;

  /* sigma = sqrt (<YY>) */
  for (i = 0 ; i < iMax ; i++)
    {
      YY[i] /= jMax ;
      YY[i] = sqrt (YY[i]) ;
    }
  /* y -> y/sigma */
   for (i = 0 ; i < iMax ; i++)
    for (j = 0 ; YY[i]>0 && j < jMax ; j++)
      y[i + iMax * j] /= YY[i] ;
   return ;
}

NN *nnCreateNumericTest6 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 1 ;
  int layerDims[] = { dimIn, 2, dimOut, 0 } ;
  int layerActivations[] = { 0, SIGMOID, 0, 0 } ;

 if (1)
    {
      int i, j ;
      float s, x[size * dimIn], y[size * dimOut] ; 

      memset (x, 0, sizeof(x)) ;
      memset (y, 0, sizeof(y)) ;
      for (i = 0 ; i < size ; i++)
	{
	  s = 0 ;
	  for (j = 0 ; j < dimIn ; j++)
	    {
	      x[j + dimIn * i] = 1 * (randfloat () - .5) ;
	      s += (j+1) * x[j + dimIn * i] ;
	    }
	  y[i] = s + 0.001 * randfloat () ;
	}
      if (1) nnNormalizeY (y, dimOut, size) ;	

      nn = nnInitialize (maxThreads, dimIn, 1, size, layerDims, layerActivations, FALSE) ;
      nn->learningRate = .1 ;

    
      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }

 
  return nn ;
}
/* a simple test with 2 qualitative classes separated by a linear equation in 2 dimemsions */
NN *nnCreateClassifierTest11 (int size, int maxThreads)
{
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 2 ;
  int layerDims[] = { dimIn,  dimIn,  dimOut, 0 } ;
  int layerActivations[] = { 0, SIGMOID, SOFTMAX, 0 } ; /* 99 : softmax */

  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .01 ;

  if (1)
    {
      int i,j ;
      float x[size * dimIn], y[size * dimOut] ; 

      for (i = 0 ; i < size ; i++)
	{
	  for (j = 0 ; j < dimIn ; j++)
	    {
	      float z = randfloat () - .5 ;
	      z = z > 0 ? 1 : -1 ;
	      x[dimIn * i + j] = z ;
	    }
	  if (x[dimIn * i + 0] * x[dimIn * i + 1] < 0)  /* x + 2y < 1, an oblique straight line */
	    { y[dimOut * i + 0] = 1 ; y[dimOut * i + 1] = 0 ; }
	  else
	    { y[dimOut * i + 0] = 0 ; y[dimOut * i + 1] = 1 ; }
	  
	  if (i < 20)
	    fprintf (stderr, "%d:\t%.2f %.2f\t%f %f\n"
		     , i
		     , x[dimIn * i + 0], x[dimIn * i + 1]
		     , y[dimOut * i + 0], y[dimOut * i + 1]
		     ) ;
	}
      nnNormalizeY (x, dimIn, size) ;	
      if (0) nnNormalizeY (y, dimOut, size) ;	

      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }

  return nn ;
}
/* a simple XOR test in 2 dimension */
NN *nnCreateClassifierTest12 (int size, int maxThreads)
{ 
  int H = size ;
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 2 ;
  int layerDims[] = { dimIn, dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, SIGMOID, 0, 0 } ; /* 99 : softmax */
  
  size = 4*H ;
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .8 ;

  if (1)
    { int i, h ;
      LAYER *ly = arrp (nn->layers, 1, LAYER) ;
      float xx[8*H] ;
      float yy[8*H] ;
      float x[] = {1,1,1,-1,-1,1,-1,-1} ;
      float y[] = { 1,-1, -1,1,  -1,1, 1,-1} ;
      float w1[] = {1,1,-1,-1} ;
      float w2[] = {4,4,-4,-4} ;
      float b1[] ={-1.5, -1.5} ;
      float b2[] ={-1, 1} ;

      /* without the noise, the algorithm does NOT converge
       * with randfloat()/10, it does converge sudently after around 1000 iterations
       *
       * RELU does NOT converge          BUG ?
       * Sigmoid does converge
       * cost L2 does converge
       * cost Softmax does not converge  BUG ?
       */
      for (h = 0 ; h < H ; h++)
	for (i = 0 ; i < 8 ; i++)
	  {
	    xx[8*h+i] = x[i] ;
	    yy[8*h+i] = y[i] ;
	  }
      for (i = 0 ; i < size * dimIn ; i++)
	xx[i] += randfloat ()/10 ;

      nnSetX (nn, xx, FALSE) ;
      nnSetY (nn, yy, FALSE) ;
      if (0)
	{  /* this is the solution */
	  mxSet (ly->W, w1) ;
	  mxSet (ly->b, b1) ;
	  ly++ ;
	  mxSet (ly->W, w2) ;
	  mxSet (ly->b, b2) ;
	}
    }

  return nn ;
}
  
NN *nnCreateClassifierTest13 (int size, int maxThreads)
{ 
  int H = 1 ;
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 2 ;
  int layerDims[] = { dimIn, dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, RELU, SOFTMAX, 0 } ; /* 99 : softmax */

  size = 4*H ;
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .8 ;

  if (1)
    { int i, h ;
      LAYER *ly = arrp (nn->layers, 1, LAYER) ;
      float xx[8*H] ;
      float yy[8*H] ;
      float x[] = {1,1,1,-1,-1,1,-1,-1} ;
      float y[] = { 1,0, 0,1,  0,1, 1,0} ;
      float w1[] = {1,1,-1,-1} ;
      float w2[] = {20,20,-20,-20} ;
      float b1[] ={-1.5, -1.5} ;
      float b2[] ={0, 10} ;

      /* without the noise, the algorithm does NOT converge
       * with randfloat()/10, it does converge sudently after around 1000 iterations
       *
       * RELU does NOT converge          BUG ?
       * Sigmoid does converge
       * cost L2 does converge
       * cost Softmax does not converge  BUG ?
       */
      for (h = 0 ; h < H ; h++)
	for (i = 0 ; i < 8 ; i++)
	  {
	    xx[8*h+i] = x[i] ;
	    yy[8*h+i] = y[i] ;
	  }
      if (0)
	for (i = 0 ; i < size * dimIn ; i++)
	  xx[i] += randfloat ()/10 ;
      
      nnSetX (nn, xx, FALSE) ;
      nnSetY (nn, yy, FALSE) ;
      if (1)
	{  /* this is the solution */
	  mxSet (ly->W, w1) ;
	  mxSet (ly->b, b1) ;
	  ly++ ;
	  mxSet (ly->W, w2) ;
	  mxSet (ly->b, b2) ;
	}
    }

  return nn ;
}
NN *nnCreateClassifierTest22 (int size, int maxThreads)
{ 
  int H = size ;
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 2 ;
  int layerDims[] = { dimIn, dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, SIGMOID, 0, 0 } ; /* 99 : softmax */

  size = 4*H ;
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .001 ;
  size = (size+99)/100 ; size *= 100 ; H = size ;
  if (1)
    { int i, h ;
      LAYER *ly = arrp (nn->layers, 1, LAYER) ;
      float xx[8*H] ;
      float yy[8*H] ;
      float x[] = {1,1,1,-1,-1,1,-1,-1} ;
      float y[] = { 1,-1, -1,1,  -1,1, 1,-1} ;
      float w1[] = {1,1,-1,-1} ;
      float w2[] = {4,4,-4,-4} ;
      float b1[] ={-1.5, -1.5} ;
      float b2[] ={-1, 1} ;

      /* without the noise, the algorithm does NOT converge
       * with randfloat()/10, it does converge sudently after around 1000 iterations
       *
       * RELU does NOT converge          BUG ?
       * Sigmoid does converge
       * cost L2 does converge
       * cost Softmax does not converge  BUG ?
       */
      if (0) H = 200 ;
      if (H>size) H=size ;

      for (h = 0 ; h < H ; h++)
	for (i = 0 ; i < 8 ; i++)
	  {
	    xx[8*h+i] = x[i] ;
	    yy[8*h+i] = y[i] ;
	  }
      if (1)
	for (i = 0 ; i < H * dimIn ; i++)
	  xx[i] += randfloat ()/10 ;

      if (size > H && (size % H) == 0)
	for (i = 1 ; i < size/H ; i++)
	  {
	    memcpy (xx + H * i * dimIn, xx,  H * dimIn * sizeof (float)) ;
	    memcpy (yy + H * i * dimOut, yy, H * dimOut * sizeof (float)) ;
	  }

      if (0) nnNormalizeY (xx, dimIn, size) ;
      nnSetX (nn, xx, FALSE) ;
      nnSetY (nn, yy, FALSE) ; 
      if (0)
	{  /* this is the solution */
	  mxSet (ly->W, w1) ;
	  mxSet (ly->b, b1) ;
	  ly++ ;
	  mxSet (ly->W, w2) ;
	  mxSet (ly->b, b2) ;
	}
    }

  return nn ;
}
  
/* a general  XOR test in 2 dimension */
NN *nnCreateClassifierTest23 (int size, int maxThreads)
{ 
  int H = 2000 ;
  NN *nn ;
  int dimIn = 2 ;
  int dimOut = 2 ;
  int layerDims[] = { dimIn, dimIn, dimOut, 0 } ;
  int layerActivations[] = { 0, SIGMOID, 0, 0 } ; /* 99 : softmax */

  size = H ;
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .8 ;

  if (1)
    { 
      int  i, h ;
      float a, b ;
      float x[2*H], y[2*H] ;

      float xx[] = {1,1,1,-1,-1,1,-1,-1} ;
      float yy[] = { 1,-1, -1,1,  -1,1, 1,-1} ;

      for (i = 0 ; i < 8 ; i++)
	{
	  y[i] = yy[i] ;
	  x[i] = xx[i] + randfloat ()/10 ;
	}
      
      for (h = 4 ; h < H ; h++)
	{
	  a = x[2*h + 0] = randfloat () ;
	  b = x[2*h + 1] = randfloat () ;
	  if (
	      (a > b && a > -b) ||
	      (-a > b && -a > -b)
	      )
	    { y[2*h+0] = 1 ; y[2*h+1]=-1 ;}
	  else
	    { y[2*h+0] = -1 ; y[2*h+1]=1 ;}
	}

      nnSetX (nn, x, FALSE) ;
      nnSetY (nn, y, FALSE) ;
    }
  
  return nn ;
}
  
/* A general XOR[prime] */

NN *nnCreateXOR (int size0, int pp, int method, float lr, int maxThreads)
{ 
  NN *nn ;
  int i, pass, size ;
  int layerDims[pp+2] ;
  int layerActivations[pp+2] ;


  memset (layerDims, 0, sizeof (layerDims)) ;
  memset (layerActivations, 0, sizeof (layerActivations)) ;

  size = size0 ;
  size = 10 * pp * pp ;
  switch (method)
    {
    case 0:
      fprintf (stderr, "XOR[%d} using RELU size=%d\n", pp, size) ;
      layerDims[0] = 2 * pp ;
      layerDims[1] = pp ;
      layerDims[2] = pp ;
      
      layerActivations[1] = RELU ;
      layerActivations[2] = SOFTMAX ;
      break ;

    case 1:
      fprintf (stderr, "XOR[%d} using SOFTRELU size=%d\n", pp, size) ;
      layerDims[0] = 2 * pp ;
      layerDims[1] = pp ;
      layerDims[2] = pp ;
      
      layerActivations[1] = SOFT_RELU ;
      layerActivations[2] = SOFTMAX ;
      break ;

    case 2:
      fprintf (stderr, "XOR[%d} using SIGMOID size=%d\n", pp, size) ;
      layerDims[0] = 2 * pp ;
      layerDims[1] = pp ;
      layerDims[2] = pp ;
      
      layerActivations[1] = SIGMOID ;
      layerActivations[2] = SOFTMAX ;
      break ;

    case 3: 
      fprintf (stderr, "XOR_2p[%d} using RELU size=%d\n", pp, size) ;
      layerDims[0] = 2 * pp ;
      layerDims[1] = 2 * pp ;
      layerDims[2] = pp ;
      
      layerActivations[1] = RELU ;
      layerActivations[2] = SOFTMAX ;
      break ;

    case 4:
      fprintf (stderr, "XOR[%d} using ELU size=%d\n", pp, size) ;
      layerDims[0] = 2 * pp ;
      layerDims[1] = pp ;
      layerDims[2] = pp ;
      
      layerActivations[1] = ELU ;
      layerActivations[2] = SOFTMAX ;
      break ;

    }
  nn = nnInitialize (maxThreads, 2*pp, pp, size, layerDims, layerActivations, FALSE) ;  
  nn->learningRate = lr ;
  nn->maxLearningRate = 10 ;
  nn->W2_regularization = .001 ;

  for (pass = 0 ; pass < 2 ; pass++)
    { 
      float x[2*pp*size] ;
      float y[pp*size] ;
      int k1, k2, k3 ;

      if (pass) size = pp * pp ;
      for (i = 0 ; i < 2*pp*size ; i++)
	x [i] = 0 ;
      for (i = 0 ; i < pp*size ; i++)
	y [i] = 0 ;

      for (i = 0 ; i < size ; i++)
	{
	  if (0 &&  i > pp * pp)
	    {
	      k1 = randint () ;
	      k2 = randint () ;
	      k3 = k1 + k2 ;
	    }
	  else
	    {
	      k1 = i/pp ;
	      k2 = i ;
	      k3 = k1 + k2 ;
	    }

	  if (k1 < 0) k1 = -k1 ;
	  if (k2 < 0) k2 = -k2 ;
	  k1 = k1 % pp ;
	  k2 = k2 % pp ;
	  k3 = (k1 + k2) % pp ;


	  x [pp * (2*i) + k1] = 1 ;     /* 1-hot k1 pp vector */
	  x [pp * (2*i + 1) + k2] = 1 ; /* 1-hot k2 pp vector */
	  y [pp * i + k3] = 1 ;         /* 1-hot k3 pp vector */
	}
      if (0)
	{
	  for (i = 0 ; i < 2*pp*size ; i++)
	    x [i] += 0.1 * randfloat()  ;
	  for (i = 0 ; i < pp*size ; i++)
	    y [i] += 0.1 * randfloat()  ;
	}
      if (0) nnNormalizeY (y, pp, size) ;	
      if (0) nnNormalizeY (x, 2*pp, size) ;
      nnSetX (nn, x, pass) ;
      nnSetY (nn, y, pass) ;
    }
   
  return nn ;
}
  
static int nnParseMNIST (float *xx, float *yy, int dimIn, int dimOut, int size)
{
  int jj, n = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate ("/home/mieg/CODE/DARKNET/MNIST/mnist_train.csv", 0, h) ;
  aceInSpecial (ai, "\n") ;

  jj = -1 ;
  while (jj < size - 1 && aceInCard (ai))
    {
      char cutter ;
      int i, z ;
      const char *ccp = aceInWordCut (ai, ",", &cutter) ;

      if (! ccp)
	break ;
      z = atoi (ccp) ;
      if (z < 0 || z >= dimOut)
	continue ;
      jj++ ;
      yy [dimOut * jj + z] = 1 ;
      for (i = 0 ; i < dimIn ; i++)  /* 784 = 28 * 28 */
	{
	  ccp = aceInWordCut (ai, ",", &cutter) ;
	  if (! ccp)
	    break ;
	  z = atoi (ccp) ;
	  if (z < 0 || z > 255)
	    break ;
	  xx[784 * jj + i] = z ;
	}
    }      
  
  ac_free (h) ;
  return n ;
}

static NN *nnCreateMNIST (int size, int maxThreads)
{
  AC_HANDLE h = ac_new_handle () ;
  NN *nn ;
  int dimIn = 784 ;
  int dimOut = 8 ;
  float x[dimIn * size] ;
  float y[dimOut * size] ;
  /*   int layerDims[] = { dimIn, 200, 40, 20, dimOut, 0 } ; */
  int layerDims[] = { dimIn, 256,256, dimOut, 0 } ;
  int layerActivations[] = { 0, 0, 0, SOFTMAX, 0 } ; 
  /* int layerActivations[] = { 0, SIGMOID, SIGMOID, SIGMOID, SOFTMAX, 0 } ; */
  
  nn = nnInitialize (maxThreads, dimIn, dimOut, size, layerDims, layerActivations, FALSE) ;
  nn->learningRate = .01 ;
  
  nnParseMNIST (x, y, dimIn, dimOut, size) ; 
  nnNormalizeY (x, dimIn, size) ;

  nnSetX (nn, x, FALSE) ;
  nnSetY (nn, y, FALSE) ;
  
  ac_free (h) ;
  return nn ;
}

/**********************************************************************/
/**********************************************************************/

int main (int argc, const char *argv[])
{
  int test = 2 ;
  int nIter = 100 ;
  int size = 100 ;
  int maxThreads = 0 ;
  int bMax = 1 ;
  int method = 0 ;
  float lr = .1 ;
  
  NN *nn = 0, *nn2 = 0 ;
  /*  randseed (1) ; */
 freeinit () ; 
 messErrorInit (argv[0]) ;

  if (0)
    { /* continuous form of XOR, verification that the distrib is uniform */
      double x, y ;
      int i, k, iMax = 1000000, p = 7 ;
      int nn[p] ;

      for (k = 0 ; k < p ; k++)
	nn[k] = 0 ;

      for (i = 0 ; i < iMax ; i++)
	{
	  x = randfloat() ;
	  y = randfloat() ;
	  k = (p - 0) * (x + y) + .5 ;
	  k = k % p ;
	  nn[k]++ ;
	  if (0) printf ("%g\t%g\t%d\n", x, y, k) ;
	}
      for (k = 0 ; k < p ; k++)
	printf ("%d\t%d\n", k, nn[k]) ;
      exit (0) ;
   }
  if (0)
    { /* optimal linearized approximation of ELU using 3 straight lines */
      double a, b, x, y1, y2, s, bests = 100000000, besta = 0, bestb = 0 ;

      for (a = 0 ; a < 1 ; a+= .001)
	for (b = -1 ; b < 0 ; b+= .001)
	  {
	    s = 0 ;
	    for (x = -5 ; x < 3 ; x += .001)
	      {
		if (x < 0)
		  y1 = exp(x) - 1 ;
		else
		  y1 = x ;
		y2 = a * (x+0) + b ;
		if (y2 < x) y2 = x ;
		if (y2 < -1) y2 = -1 ;
		s += (y2 - y1)*(y2 - y1) ;
	      }
	    if (s < bests)
	      {
		bests = s ;
		besta = a ; bestb = b ;
	      }
	  }
      printf ("a=%f b =%f s = %f\n", besta, bestb, bests) ;
      exit (0) ;

    }


 getCmdLineInt (&argc, argv, "-t", &(test)) ;
 getCmdLineInt (&argc, argv, "-n", &(nIter)) ;
 getCmdLineInt (&argc, argv, "-s", &(size)) ;
 getCmdLineInt (&argc, argv, "-m", &(method)) ;
 getCmdLineInt (&argc, argv, "-th", &(maxThreads)) ;
 getCmdLineFloat (&argc, argv, "-lr", &(lr)) ;
 getCmdLineInt (&argc, argv, "-bMax", &(bMax)) ;
 switch (test)
   {
   case 200:
     {
       uudTest () ;
       exit (0) ;
     }
   case 97:   /* 1D convolution */
     {
       AC_HANDLE h = ac_new_handle () ;
       int dx = 4, w = 2, k = 2 ;
       MX a = mxCreate (h, "A", MX_FLOAT, dx, k, 0) ;
       MX b = mxCreate (h, "B", MX_FLOAT, dx, k, 0) ;
       MX dX = mxCreate (h, "B", MX_FLOAT, w*dx, k, 0) ;
       MX W = mxCreate (h, "W", MX_FLOAT, w, 0, 0) ;
       int i,  size = k * dx ;
       float x[size] ;

       for (i = 0 ; i < size ; i++)
	 x[i] = 100 * randfloat() ;
       mxSet (b, x) ;

       for (i = 0 ; i < size ; i++)
	 x[i] = 1 ; /* pure sum */
       mxSet (W, x) ;

       mxConvolution (a, b, W,h) ;
       mxShow (a) ;
       mxShow (b) ;

       mxMaxPoolingBack (dX, a, b, W, h) ;
       mxShow (dX) ;

       return 0 ;
     }
   case 98:   /* 1D max pooling */
     {
       AC_HANDLE h = ac_new_handle () ;
       int dx = 2, w = 3, k = 5 ;
       MX a = mxCreate (h, "A", MX_FLOAT, dx, k, 0) ;
       MX b = mxCreate (h, "B", MX_FLOAT, w*dx, k, 0) ;
       MX dX = mxCreate (h, "B", MX_FLOAT, w*dx, k, 0) ;
       MX W = mxCreate (h, "W", MX_FLOAT, w, 0, 0) ;
       int i,  size = k * w * dx ;
       float x[size] ;

       for (i = 0 ; i < size ; i++)
	 x[i] = 100 * randfloat() ;

       mxSet (b, x) ;
       mxMaxPooling (a, b, W, h) ;
       mxShow (a) ;
       mxShow (b) ;

       mxMaxPoolingBack (dX, a, b, W, h) ;
       mxShow (dX) ;

       return 0 ;
     }
   case 99: /* 2D max pooling */
     {
       AC_HANDLE h = ac_new_handle () ;
       int dx = 2, w = 3 ;
       MX a = mxCreate (h, "A", MX_FLOAT, dx, dx, 0) ;
       MX b = mxCreate (h, "B", MX_FLOAT, w*dx, w*dx, 0) ;
       MX dX = mxCreate (h, "B", MX_FLOAT, w*dx, w*dx, 0) ;
       MX W = mxCreate (h, "W", MX_FLOAT, w, w, 0) ;
       int i,  size = w * w * dx * dx ;
       float x[size] ;

       for (i = 0 ; i < size ; i++)
	 x[i] = 100 * randfloat() ;

       mxSet (b, x) ;
       mxMaxPooling (a, b, W, h) ;
       mxShow (a) ;
       mxShow (b) ;

       mxMaxPoolingBack (dX, a, b, W, h) ;
       mxShow (dX) ;

       return 0 ;
     }

     /* NUMERIC FIT */
     /* Test the recovery of a linear equation, no activation, no classes */


   case 1:
     nn = nnCreateNumericTest1 (size, maxThreads) ; /* ok 2018_04_09  y = 3 * x + 4.  Fast if bMax=2, SLOW if bMax = 1 */
     break ;
   case 2:
     nn = nnCreateNumericTest2 (size, maxThreads) ; /* ok 2018_04_09 y = u + 4 *v - 6.5 * w + 2.  Fast if bMax=2, SLOW if bMax = 1 */
     break ; 
   case 3:
     nn = nnCreateNumericTest3 (size, maxThreads) ; /* ok 2018_04_09 y = u + v.  Fast if  bMax=2, FAILS with bMax = 1 */
     break ;
   case 4:
     nn = nnCreateNumericTest4 (size, maxThreads) ; /* ok 2018_04_09 y = u + v.  Fast if  bMax=2, FAILS with bMax = 1 */
     break ;
   case 5:
     nn = nnCreateNumericTest5 (size, maxThreads) ; /* ok 2018_04_09, y = u + 2 * v + 3 * w ..., Fast if bMax = 2, FAILS if bMax = 1 */
     break ;
   case 6:
     nn = nnCreateNumericTest6 (size, maxThreads) ; /* ok 2018_04_09 after we normalize the Y , pas evedient 2-18_08 idem 2 layersY */
     break ;
 
   case 8:
     nn = nnCreateLineTest (size, maxThreads) ; /* ok 2018_04_09 y=30*x+40, works 5548 bMax=2, FAILS bMax=1*/
     break ;

     /* SIMPLE CLASSIFIER */
     /* Test the recovery of a linear classification */
   case 11:
     nn = nnCreateClassifierTest11 (size, maxThreads) ; /* not ok 2018_04_09 */
     break ;
   case 12:
     nn = nnCreateClassifierTest12 (size, maxThreads) ; /* XOR using square distance, works 2018_04_12 */
     break ;
  case 13:
    nn = nnCreateClassifierTest13 (size, maxThreads) ; /* XOR using softmax, does not work */
     break ;
  case 22:
     nn = nnCreateClassifierTest22 (size, maxThreads) ; /* XOR works 2018_04_12 */
     break ;
   case 23:
     nn = nnCreateClassifierTest22 (size, maxThreads) ; /* XOR works 2018_04_12 */
     nnTrainNetwork (nn, nIter) ;
     nn2 = nnCreateClassifierTest23 (size, maxThreads) ; /* XOR works 2018_04_12 */
     /* transfer knowledge */
     {
       int i ;
       LAYER *ly1 = arrp (nn->layers, 0, LAYER) ;
       LAYER *ly2 = arrp (nn2->layers, 0, LAYER) ;
       for (i = 1 ; i < arrayMax (nn->layers) ; i++)
	 {
	   ly2[i].W = ly1[i].W ; 
	   ly2[i].b = ly1[i].b ; 
	 }
       nIter = 100 ;  
       nn2->learningRate = nn->learningRate ;
       nn = nn2 ;
     }
     break ;
   case 400:
     nn = nnCreateMNIST (size, maxThreads) ;
     break ;
   default:
     if (test > 1000 && test < 2000)
       nn = nnCreateXOR (size, test - 1000, method, lr, maxThreads) ;
     else
       messcrash ("Unknown -t %d\n", test) ;
     break ;
   }
 nn->bMax = bMax ;
  printf ("Train\n") ;
 nnTrainNetwork (nn, nIter);
 printf ("Check\n") ;
 nnCheckNetwork (nn) ;
  printf ("Test\n") ;
  nnTestNetwork (nn) ;

  printf ("Done method=%d test=%d pp=%d lr=%.3f\n", method, test, test-1000, lr ) ;
 /*
    train
    test
  */
 ac_free (nn) ;
  return 0 ; 
} /* main */

/**********************************************************************/
/**********************************************************************/
/* APPLICATION */
/**********************************************************************/
#ifdef JUNK


def getFastc (fileName, ln=50, m=10000):
    np.random.seed (4)
    b2n = {'a':1, 't':-1, 'g':2, 'c':-2, 'n':0}
    b2c = {'t':1, 'a':-1, 'c':2, 'g':-2, 'n':0}
    X = np.zeros((ln,m))
    Y = np.zeros((1,m))
    n = -1
    kk=0
    with open(fileName) as f:
        for seq in f:
            if seq[0] != '>':
                kk += 1
                if (0*kk % 100) > 0:
                    continue
                n += 1
                if n >= m:
                    break
                i = -1
                k = [0,0,0,0,0,0]
                kk = 0
                y = (np.random.randn () < .3)
                Y[0,n] = y
                # print ("n=", n, "seq=", seq)
                for c in seq:
                    i += 1
                    kk += 1
                    if 0*i==21:
                        print (c,seq)
                    if i >= 20 and i < 20 + ln:
                        if y == 1:
                            X[i-20,n] = b2n[c]
                            k[3+b2n[c]] += 1
                        else:
                            X[i-20,n] = b2c[c]
                            k[3+b2c[c]] += 1 
                #X[ln-1, n] = 16*k[3-3]/kk
                #X[ln-2, n] = 16*k[3-2]/kk
                #X[ln-3, n] = 16*k[3-1]/kk
                #X[ln-4, n] = 16*k[3+1]/kk
                #X[ln-5, n] = 16*k[3+2]/kk
                # X[ln-5, n] = y
                #X[ln-1,n] = y # this is cheating, but now the NN converges
                #X[0,n] = y # this is cheating, but now the NN converges
    #X /= 4 
    return X, Y


#endif  
/**********************************************************************/
/**********************************************************************/
