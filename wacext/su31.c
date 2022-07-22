#include "ac.h"
#include <complex.h>
#include "matrix.h"

/* Create june 2022
 * consruct the superalgebra sl(3/1)
 */

#include "ac.h"
#include <complex.h>
#include "matrix.h"

/*************************************************************************************/


/*****  SU(2/1) representation theory. This is used by, but does not depend on the analysis above of the Feynman diagrams **********************/
/*****  Casimir studies with Pater Jarvis, mars 2021 **********************************************/
/***** Scalar anomaly paper is below **********************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

typedef struct kasStruct { MX *mu, *QQ, gg, GG, kas2, CHI, ccc, CCC, cccGhost, CCCGhost, c4, c5, C4, C5, kas3 ; int a1, a2, b, d, d1, d2, d3, d4, d5, d6,d7, d8, chi, scale, NN ; BOOL show ; float zc4, zC4 ; AC_HANDLE h ; } KAS ;
typedef struct comtpStruct { int a, b, c, n, s ; } LC ;
static MX KasCommut (MX a, MX b, int sign, KAS *kas) ;
  
static MX *constructSU3Matrices (KAS *kas)
{
  MX muE1, muE2, muE3 ;
  MX muF1, muF2, muF3 ;
  MX muH1, muH2, muH3  ; /* the 8 generators of SU(3) in the Chevalley basis */

  MX muY, muY1, muY2, muY3 ;
  MX muU1, muU2, muU3 ;
  MX muV1, muV2, muV3 ;

  MX CHI ;
  
  MX *mu ;
  int d ;
  int i ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (60 * sizeof (MX), kas->h) ;

  kas->chi = 1 ;
  
  kas->a1 = 0 ;
  kas->a2 = 0 ;
  int b = kas->b ;
  
  kas->d = d = 8 ;

#ifdef JUNK
  int d0, d1, d2, d3, d4, d5, d6, d7 ;
  d0 = 1 ; d1 = 3 ; d2 = d3 = 0 ; d4 = 3 ; d5 = d6 = 0 ; d7 = 1 ; 
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = d3 ;
  kas->d4 = d4 ;
  kas->d5 = d5 ;
  kas->d6 = d6 ;
  kas->d7 = d7 ;
#endif
  
  mu[0] = CHI = mxCreate (h,  "chi", MX_INT, d, d, 0) ;
  mu[1] = muH1 = mxCreate (h,  "muH1", MX_INT, d, d, 0) ;
  mu[2] = muH2 = mxCreate (h,  "muH2", MX_INT, d, d, 0) ;
  mu[3] = muH3 = mxCreate (h,  "muH3", MX_INT, d, d, 0) ;

  mu[4] = muY  = mxCreate (h,  "muY", MX_INT, d, d, 0) ;
  mu[5] = muY1 = mxCreate (h, "muY1", MX_INT, d, d, 0) ;
  mu[6] = muY2 = mxCreate (h, "muY2", MX_INT, d, d, 0) ;
  mu[7] = muY3 = mxCreate (h, "muY3", MX_INT, d, d, 0) ;

  
  mu[11] = muE1 = mxCreate (h,  "muE1", MX_INT, d, d, 0) ;
  mu[12] = muE2 = mxCreate (h,  "muE2", MX_INT, d, d, 0) ;
  mu[13] = muE3 = mxCreate (h,  "muE3", MX_INT, d, d, 0) ;

  mu[21] = muF1 = mxCreate (h,  "muF1", MX_INT, d, d, 0) ;
  mu[22] = muF2 = mxCreate (h,  "muF2", MX_INT, d, d, 0) ;
  mu[23] = muF3 = mxCreate (h,  "muF3", MX_INT, d, d, 0) ;

  mu[31] = muU1 = mxCreate (h,  "muU1", MX_INT, d, d, 0) ;
  mu[32] = muU2 = mxCreate (h,  "muU2", MX_INT, d, d, 0) ;
  mu[33] = muU3 = mxCreate (h,  "muU3", MX_INT, d, d, 0) ;

  mu[41] = muV1 = mxCreate (h,  "muV1", MX_INT, d, d, 0) ;
  mu[42] = muV2 = mxCreate (h,  "muV2", MX_INT, d, d, 0) ;
  mu[43] = muV3 = mxCreate (h,  "muV3", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  int chi [] = {
		1,  0, 0, 0,    0, 0, 0,   0,
		
		0, -1, 0, 0,    0, 0, 0,   0,
		0,  0,-1, 0,    0, 0, 0,   0,
		0,  0, 0,-1,    0, 0, 0,   0,
		
		0,  0, 0, 0,    1, 0, 0,   0,
		0,  0, 0, 0,    0, 1, 0,   0,
		0,  0, 0, 0,    0, 0, 1,   0,

		0,  0, 0, 0,    0, 0, 0,  -1
  };
  mxSet (CHI, chi) ;
  mxShow (CHI) ;
 
  int xx1 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  1, 0, 0,    0, 0, 0,   0,
		0,  0,-1, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,   -1, 0, 0,   0,
		0,  0, 0, 0,    0, 1, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muH1, xx1) ;
  mxShow (muH1) ;
 
  int xx2 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 1, 0,    0, 0, 0,   0,
		0,  0, 0,-1,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0,-1, 0,   0,
		0,  0, 0, 0,    0, 0, 1,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muH2, xx2) ;
  mxShow (muH2) ;
 
  int xx3 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  1, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0,-1,    0, 0, 0,   0,
		
		0,  0, 0, 0,   -1, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 1,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muH3, xx3) ;
  mxShow (muH3) ;


  int xxy [] = {
		b,  0, 0, 0,    0, 0, 0,   0,
		
		0,b-1, 0, 0,    0, 0, 0,   0,
		0,  0,b-1,0,    0, 0, 0,   0,
		0,  0,0,b-1,    0, 0, 0,   0,
		
		0,  0, 0, 0,  b-2, 0, 0,   0,
		0,  0, 0, 0,    0,b-2,0,   0,
		0,  0, 0, 0,    0,0,b-2,   0,

		0,  0, 0, 0,    0, 0, 0, b-3
  };
  mxSet (muY, xxy) ;
  mxShow (muY) ;

  int B = 2*b ;
  int xxy1 [] = {
		B,  0, 0, 0,    0, 0, 0,   0,
		
		0,  B, 0, 0,    0, 0, 0,   0,
		0,  0,B-3,0,    0, 0, 0,   0,
		0,  0,0,B-3,    0, 0, 0,   0,
		
		0,  0, 0, 0,  B-6, 0, 0,   0,
		0,  0, 0, 0,    0,B-3,0,   0,
		0,  0, 0, 0,    0,0,B-3,   0,

		0,  0, 0, 0,    0, 0, 0, B-6
  };
  mxSet (muY1, xxy1) ;
  mxShow (muY1) ;
 
  int xxy2 [] = {
		B,  0, 0, 0,    0, 0, 0,   0,
		
		0,B-3, 0, 0,    0, 0, 0,   0,
		0,  0, B, 0,    0, 0, 0,   0,
		0,  0,0,B-3,    0, 0, 0,   0,
		
		0,  0, 0, 0,  B-3, 0, 0,   0,
		0,  0, 0, 0,    0,B-6,0,   0,
		0,  0, 0, 0,    0,0,B-3,   0,

		0,  0, 0, 0,    0, 0, 0, B-6
  };
  mxSet (muY2, xxy2) ;
  mxShow (muY2) ;
 
  int xxy3 [] = {
		B,  0, 0, 0,    0, 0, 0,   0,
		
		0,B-3, 0, 0,    0, 0, 0,   0,
		0,  0,B-3,0,    0, 0, 0,   0,
		0,  0,0,  B,    0, 0, 0,   0,
		
		0,  0, 0, 0,  B-3, 0, 0,   0,
		0,  0, 0, 0,    0,B-3,0,   0,
		0,  0, 0, 0,    0,0,B-6,   0,

		0,  0, 0, 0,    0, 0, 0, B-6
  };
  mxSet (muY3, xxy3) ;
  mxShow (muY3) ;
 
   /* even lowering operator */
  int xf1 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 1, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    1, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muF1, xf1) ;
  mxShow (muF1) ;
 
  int xf2 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 1,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 1, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muF2, xf2) ;
  mxShow (muF2) ;

  int xf3 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 1,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,   -1, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  mxSet (muF3, xf3) ;
  mxShow (muF3) ;

  /* even raising operator */
  int xe1 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  1, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 1, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };
  mxSet (muE1, xe1) ;
  mxShow (muE1) ;

  int xe2 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 1, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 1,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  memset (xx, 0, sizeof (xx)) ;
  i = 2 ; xx[d * i + 1] = 1 ;
  mxSet (muE2, xe2) ;
  mxShow (muE2) ;

  int xe3 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  1, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0,-1,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  mxSet (muE3, xe3) ;
  mxShow (muE3) ;

  
  int xv1 [] = {
		0,  1, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 1,   0,
		0,  0, 0, 0,    0, 1, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   1,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  mxSet (muV1, xv1) ;
  mxShow (muV1) ;

  int xv2 [] = {
		0,  0, 1, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0,-1,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    1, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,  -1,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  mxSet (muV2, xv2) ;
  mxShow (muV2) ;

  int xv3 [] = {
		0,  0, 0, 1,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0,-1, 0,   0,
		0,  0, 0, 0,   -1, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   1,

		0,  0, 0, 0,    0, 0, 0,   0
  };

  mxSet (muV3, xv3) ;
  mxShow (muV3) ;

    int xu1 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		B,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0,B-3,    0, 0, 0,   0,
		0,  0,B-3, 0,    0, 0, 0,   0,

		0,  0, 0, 0,  B-6, 0, 0,   0
  };

  mxSet (muU1, xu1) ;
  mxShow (muU1) ;

  int xu2 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		B,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0,B-3,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		0,3-B, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0,-B+6,0,  0
  };

  mxSet (muU2, xu2) ;
  mxShow (muU2) ;

  int xu3 [] = {
		0,  0, 0, 0,    0, 0, 0,   0,
		
		0,  0, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,
		B,  0, 0, 0,    0, 0, 0,   0,
		
		0, 0,3-B, 0,    0, 0, 0,   0,
		0,  3-B, 0, 0,    0, 0, 0,   0,
		0,  0, 0, 0,    0, 0, 0,   0,

		0,  0, 0, 0,    0,0,B-6,   0
  };

  mxSet (muU3, xu3) ;
  mxShow (muU3) ;


  return mu ;
} /* constructSU3Matrices */

/*************************************************************************************/
/*************************************************************************************/

static MX KasCommut (MX a, MX b, int sign, KAS *kas)
{
  MX p = mxMatMult (a, b, kas->h) ;
  MX q = mxMatMult (b, a, kas->h) ;
  MX r = mxCreate (kas->h, "r", MX_INT, kas->d, kas->d, 0) ;

  r = sign == 1 ? mxAdd (r, p, q, kas->h) : mxSubstract (r, p, q, kas->h) ;
  
  return r ;
}

/*************************************************************************************/

static MX checkOneCommutator (LC *up, KAS *kas)
{
  int d = kas->d ;
  int dd = kas->d * kas->d ;
  
  MX a = kas->mu[up->a] ;
  MX b = kas->mu[up->b] ;
  MX c = kas->mu[up->c] ;
  MX r = mxCreate (kas->h, "r", MX_INT, d,d,0) ;
  MX s = mxCreate (kas->h, "s", MX_INT, d,d,0) ;
  MX t = mxCreate (kas->h, "t", MX_INT, d,d,0) ;
  const int *xx ;
  int yy [dd] ;
  int i, k ;

  MX ab = KasCommut (a, b, up->s, kas) ;
  int scale = (up->s == 1 ? (kas->scale ? kas->scale : 1) : 1) ;

  mxValues (c, &xx, 0, 0) ;
  for (i = 0 ; i < dd ; i++)
    yy[i] = up->n * scale * xx[i] ;
  mxSet (s, yy) ;
  t = mxSubstract (t, ab, s, kas->h) ;
  mxValues (t, &xx, 0, 0) ;
  for (i = k = 0 ; i < dd ; i++)
    k += xx[i] * xx[i] ;
  if (k > 0)
    {
      mxShow (ab) ;
      mxShow (c) ;
      if (up->s == -1)
	messcrash ("\nKasChect Failed [%s,%s] = %d * %s\n", a->name, b->name, up->n, c->name) ;
      else
	messcrash ("\nKasChect Failed {%s,%s} = %d * %s\n", a->name, b->name, up->n, c->name) ;
    }
  else if (0)
    {
      if (up->s == -1 && up->n == 0)
	printf ("[%s,%s] = 0\n",  a->name, b->name) ;
      else if (up->s == 1 && up->n == 0)
	printf ("{%s,%s} = 0\n",  a->name, b->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1)
	printf ("[%s,%s] = %d %s\n",  a->name, b->name, up->n, c->name) ;

      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1)
	printf ("{%s,%s} = %d %s\n",  a->name, b->name, up->n, c->name) ;
    }
  return r ;
}

/*************************************************************************************/

static void checkCommutators (KAS *kas)
{
  LC *up, *XXX ;
  LC Su3XXX[] = {
		 /* Check the SU(3) algebra */
		 
		 {1,2,2,0,-1},     /* [h1,h2] = 0 */
				 
		 {11,12,13,1,-1},  /* [e1,e2] = e3 */
		 {11,13,12,0,-1},  /* [e1,e3] = 0 */
		 {12,13,11,0,-1},  /* [e2,e3] = 0 */

		 {21,22,23,-1,-1},  /* [f1,f2] = -f3 */
		 {21,23,22,0,-1},  /* [f1,f3] = 0 */
		 {22,23,21,0,-1},  /* [f2,f3] = 0 */



		 
		 {1,11,11,2,-1},   /* [h1,e1] = 2 e1 */
		 {1,12,12,-1,-1},  /* [h1,e2] = -e2 */
		 {1,13,13,1,-1},   /* [h1,e3] = e3 */


		 {1,21,21,-2,-1},  /* [h1,f1] = -2 f1 */
		 {1,22,22,1,-1},   /* [h1,f2] = f2 */
		 {1,23,23,-1,-1},  /* [h1,f3] = -f3 */

		 {2,11,11,-1,-1},  /* [h2,e1] = -e1 */
		 {2,12,12,2,-1},   /* [h2,e2] = 2e2 */
		 {2,13,13,1,-1},   /* [h2,e3] = e3 */


		 {2,21,21,1,-1},   /* [h2,f1] = f1 */
		 {2,22,22,-2,-1},  /* [h2,f2] = -2 f2 */
		 {2,23,23,-1,-1},  /* [h2,f3] = -f3 */

		 {11,21,1,1,-1},   /* [e1,f1] = h1 */
		 {11,22,23,0,-1},  /* [e1,f2] = 0 */
		 {11,23,22,-1,-1},  /* [e1,f3] = -f2 */

		 {12,21,1,0,-1},   /* [e2,f1] = 0 */
		 {12,22,2,1,-1},   /* [e2,f2] = h2 */
		 {12,23,21,1,-1},  /* [e2,f3] = 0 */

		 {13,21,12,-1,-1}, /* [e3,f1] = -e2 */
		 {13,22,11,1,-1},  /* [e3,f2] = e1 */
		 {13,23,3,1,-1},   /* [e3,f3] = h1+h2 */

		 /* Check that the v[123] form an  SU(3) module */

		 {0,41,41, 0, 1},   /* {chi,v1} = 0 */
		 {1,41,41, 1,-1},   /* [h1,v1] = v1 */
		 {2,41,41, 0,-1},   /* [h2,v1] = 0 */
		 {4,41,41, -1,-1},   /* [Y,v1] = -v1 */

		 {0,42,42, 0, 1},   /* {chi,v1} = 0 */
		 {1,42,42,-1,-1},   /* [h1,v2] = -v2 */
		 {2,42,42, 1,-1},   /* [h2,v2] = v2 */
		 {4,42,42, -1,-1},   /* [Y,v2] = -v2 */

		 
		 {0,43,43, 0, 1},   /* {chi,v1} = 0 */
		 {1,43,43, 0,-1},   /* [h1,v3] = 0 */
		 {2,43,43,-1,-1},   /* [h2,v3] = -v3 */
		 {4,43,43, -1,-1},   /* [Y,v3] = -v3 */

		 {21,41,42, 1,-1},   /* [f1,v1] = v2 */
		 {21,42,42, 0,-1},   /* [f1,v2] = 0 */
		 {21,43,42, 0,-1},   /* [f1,v3] = 0 */
		 
		 {22,41,43, 0,-1},   /* [f2,v1] = 0 */
		 {22,42,43, 1,-1},   /* [f2,v2] = v3 */
		 {22,43,43, 0,-1},   /* [f2,v3] = 0 */
		 
		 {23,41,43, 1,-1},   /* [f3,v1] = v3 */
		 {23,42,43, 0,-1},   /* [f3,v2] = 0 */
		 {23,43,43, 0,-1},   /* [f3,v3] = 0 */
		 
		 {11,41,42, 0,-1},   /* [e1,v1] = 0 */
		 {11,42,41, 1,-1},   /* [e1,v2] = v1 */
		 {11,43,42, 0,-1},   /* [e1,v3] = 0 */
		 
		 {12,41,43, 0,-1},   /* [e2,v1] = 0 */
		 {12,42,43, 0,-1},   /* [e2,v2] = 0 */
		 {12,43,42, 1,-1},   /* [e2,v3] = v2 */
		 
		 {13,41,43, 0,-1},   /* [e3,v1] = 0 */
		 {13,42,43, 0,-1},   /* [e3,v2] = 0 */
		 {13,43,41, 1,-1},   /* [e3,v3] = v1 */
		 
		 /* Check that the u[123] form an  SU(3) module */

		 {0,31,31, 0, 1},   /* {chi,u1} = 0 */
		 {1,31,31,-1,-1},   /* [h1,u1] = -u1 */
		 {2,31,31, 0,-1},   /* [h2,u1] = 0 */
		 {4,31,31, 1,-1},   /* [Y,u1] = u1 */

		 {0,31,31, 0, 1},   /* {chi,u1} = 0 */
		 {1,32,32, 1,-1},   /* [h1,u2] = u2 */
		 {2,32,32, -1,-1},   /* [h2,u2] = -u2 */
		 {4,31,31, 1,-1},   /* [Y,u1] = u1 */
		 
		 {0,31,31, 0, 1},   /* {chi,u1} = 0 */
		 {1,33,33, 0,-1},   /* [h1,u3] = 0 */
		 {2,33,33, 1,-1},   /* [h2,u3] = u3 */
		 {4,31,31, 1,-1},   /* [Y,u1] = u1 */


		 {22,31,33, 0,-1},   /* [f2,u1] = 0 */
		 {23,32,33, 0,-1},   /* [f3,u2] = 0 */



		 {21,31,32, 0,-1},   /* [f1,u1] = 0 */
		 {21,33,32, 0,-1},   /* [f1,u3] = 0 */
		 
		 {22,31,33, 0,-1},   /* [f2,u1] = 0 */
		 {22,32,33, 0,-1},   /* [f2,u2] = 0 */
		 
		 {23,31,33, 0,-1},   /* [f3,u1] = 0 */
		 {23,32,33, 0,-1},   /* [f3,u2] = 0 */

		 
		 {11,32,31, 0,-1},   /* [e1,u2] = 0 */
		 {11,33,32, 0,-1},   /* [e1,u3] = 0 */
		 
		 {12,31,33, 0,-1},   /* [e2,u1] = 0 */
		 {12,33,32, 0,-1},   /* [e2,u3] = 0 */
		 
		 {13,32,33, 0,-1},   /* [e3,u2] = 0 */
		 {13,33,31, 0,-1},   /* [e3,u3] = 0 */

		 
		 {11,31,32,-1,-1},   /* [e1,u1] = -u2 */
		 {12,32,33,-1,-1},   /* [e2,u2] = -u3 */

		 {21,32,31, -1,-1},  /* [f1,u2] = -u1 */
		 {22,33,32,-1,-1},   /* [f2,u3] = -u2 */
		 {23,33,31,-1,-1},   /* [f3,u3] = -u1 */
		 {11,31,32,-1,-1},   /* [e1,u1] = -u2 */
		 {12,32,33,-1,-1},   /* [e2,u2] = -u3 */
		 {13,31,33,-1,-1},   /* [e3,u1] = -u3 */

		 
		 {31,41,5, 1,1},     /* {u1,v1} = Y1 */
		 {32,42,6, 1,1},     /* {u2,v2} = Y2 */
		 {33,43,7, 1,1},     /* {u2,v2} = Y2 */

		 {31,42,21,3,1},     /* {u1,v2} = 3 f1 */
		 {32,41,11,3,1},     /* {u2,v1} = 3 e1 */
		 
		 {32,43,22,3,1},     /* {u2,v3} = 3 f2 */
		 {33,42,12,3,1},     /* {u3,v2} = 3 e2 */

		 {31,43,23,3,1},     /* {u1,v3} = 3 f3 */
		 {33,41,13,3,1},     /* {u3,v1} = 3 e3 */

		 {31,31,31,0,1},     /* {u1,u1} = 0 */
		 {31,32,31,0,1},     /* {u1,u2} = 0 */
		 {31,33,31,0,1},     /* {u1,u3} = 0 */
		 
		 {32,31,31,0,1},     /* {u2,u1} = 0 */
		 {32,32,31,0,1},     /* {u2,u2} = 0 */
		 {32,33,31,0,1},     /* {u2,u3} = 0 */
		 
		 {33,31,31,0,1},     /* {u3,u1} = 0 */
		 {33,32,31,0,1},     /* {u3,u2} = 0 */
		 {33,33,31,0,1},     /* {u3,u3} = 0 */
		 
		 {41,41,41,0,1},     /* {v1,v1} = 0 */
		 {41,42,41,0,1},     /* {v1,v2} = 0 */
		 {41,43,41,0,1},     /* {v1,v4} = 0 */
		 
		 {42,41,41,0,1},     /* {v2,v1} = 0 */
		 {42,42,41,0,1},     /* {v2,v2} = 0 */
		 {42,43,41,0,1},     /* {v2,v4} = 0 */
		 
		 {43,41,41,0,1},     /* {v4,v1} = 0 */
		 {43,42,41,0,1},     /* {v4,v2} = 0 */
		 {43,43,41,0,1},     /* {v3,v3} = 0 */
		 
		 {0,0,0,0,0}
  } ;



  /* check that K1 = Y+H */
  if (0)
    {
      int i, d = kas->d ;
      const int *xxY = messalloc (d*d*sizeof(int)) ;
      const int *xxH = messalloc (d*d*sizeof(int)) ;
      const int *xxK1 = messalloc (d*d*sizeof(int)) ;
      const int *xxK2 = messalloc (d*d*sizeof(int)) ;
      
      mxValues (kas->mu[0], &xxY, 0, 0) ;
      mxValues (kas->mu[3], &xxH, 0, 0) ;
      mxValues (kas->mu[8], &xxK1, 0, 0) ;
      mxValues (kas->mu[9], &xxK2, 0, 0) ;
	    
      for (i = 0 ; i < d * d ; i++)
	if (2*xxK1[i] !=  xxY[i] + xxH[i])
	  messcrash ("K1=(Y+H)/2 failed for i=%d\n",i) ;
      for (i = 0 ; i < d * d ; i++)
	if (2*xxK2[i] !=  xxY[i] - xxH[i])
	  messcrash ("K2=(Y-H)/2 failed for i=%d\n",i) ;
    }

  XXX = Su3XXX ;

  
  for (up = XXX ; up->s ; up++)
    checkOneCommutator (up, kas) ;
  printf ("SUCCESS (a1=%d a2=%d, 0) all comutators have been verified\n", kas->a1, kas->a2) ;
  return ;
} /* checkCommutators */

/*************************************************************************************/
/*************************************************************************************/

static int mxSuperTrace (KAS *kas, MX m)
{
  AC_HANDLE h = ac_new_handle () ;
  
  MX M = mxMatMult (kas->mu[0], m, kas->h) ;
  int i, t = 0, d = kas->d  ;
  const int *xx ;
  mxValues (M, &xx, 0, 0) ;
  for (i = 0 ; i < d ; i++)
    t += xx[d*i + i] ;

  ac_free (h) ;
  return t ;
} /* mxSuperTrace */

/*************************************************************************************/

static void lowerMetric (KAS *kas)
{
  int i,j ;
  int xx [2500] ;

  memset (xx, 0, sizeof (xx)) ;
  /* even sector */
  for (i = 1 ; i < 30 ; i++)
    for (j = 1 ; j < 30 ; j++)
      {
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;

	if (a && b)
	  {
	    MX m = mxMatMult (a,  b, kas->h) ;
	    int  t = mxSuperTrace (kas, m) ;
	    xx[50*i+j] = t ;
	    xx[50*j+i] = -t ;
	    
	    if (t)
	      printf ("gg[%d,%d] = %d\n", i,j,t) ;
	  }
      }

  /* odd sector */
  for (i = 31 ; i < 50 ; i++)
    for (j = 31 ; j < 50 ; j++)
      {
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;

	if (a && b)
	  {
	    MX m = mxMatMult (a, b, kas->h) ;
	    int  t = mxSuperTrace (kas, m) ;
	    xx[50*i+j] = t ;
	    xx[50*j+i] = -t ;
	    
	    if (t)
	      printf ("gg[%d,%d] = %d\n", i,j,t) ;
	  }
      }
  
  return ;
} /* lowerMetric */

/*************************************************************************************/

static void upperMetric (KAS *kas)
{
    
  return ;
} /* upperMetric */

/*************************************************************************************/
/*************************************************************************************/

static void KillingCasimir (KAS *kas)
{
    
  return ;
} /* KillingCasimir */

/*************************************************************************************/

static void ghostCasimirs (KAS *kas)
{
  MX uv1 = mxMatMult (kas->mu[31],kas->mu[41], kas->h) ;
  MX uv2 = mxMatMult (kas->mu[32],kas->mu[42], kas->h) ;
  MX uv3 = mxMatMult (kas->mu[33],kas->mu[43], kas->h) ;
  MX muY = kas->mu[4] ;

  MX vu1 = mxMatMult (kas->mu[41],kas->mu[31], kas->h) ;
  MX vu2 = mxMatMult (kas->mu[42],kas->mu[32], kas->h) ;
  MX vu3 = mxMatMult (kas->mu[43],kas->mu[33], kas->h) ;

  MX uv12 = mxCreate (kas->h, "uv12", MX_INT, kas->d, kas->d, 0) ;
  MX uv123 = mxCreate (kas->h, "uv123", MX_INT, kas->d, kas->d, 0) ;
  
  MX vu12 = mxCreate (kas->h, "vu12", MX_INT, kas->d, kas->d, 0) ;
  MX vu123 = mxCreate (kas->h, "vu123", MX_INT, kas->d, kas->d, 0) ;

  uv12 = mxAdd (uv12, uv1, uv2, kas->h) ;
  uv123 = mxAdd (uv123, uv12, uv3, kas->h) ;
  vu12 = mxAdd (vu12, vu1, vu2, kas->h) ;
  vu123 = mxAdd (vu123, vu12, vu3, kas->h) ;


  MX w1 = mxMatMult (uv123, vu123, kas->h) ;
  MX w2 = mxMatMult (vu123,uv123, kas->h) ;
  int a1 = kas->a1 ;
  int a2 = kas->a2 ;
  int b = kas->b ;
  
  printf ("##### GHOST CASIMIR\n") ;
  printf ("##### expect %d\n",
	  (2*b - (2*a1 + a2))*(2*b -( 3+a1-a2))*(2*b -(6 + a1 + 2*a2))
	  ) ;
  mxShow (uv123) ;
  mxShow (vu123) ;
  mxShow (muY) ;
  MX w12 = mxCreate (kas->h, "w12", MX_INT, kas->d, kas->d, 0) ;
  w12 = mxAdd (w12, w1, w2, kas->h) ;

  
  mxShow (w1) ;
  mxShow (w2) ;
  printf ("##### GHOST CASIMIR w12  of course uv123 + vu123 = 6 Y by construction, not interesting\n") ;
  mxShow (w12) ;


  /* Gorelik ghost Casimir */
  MX u12 = mxMatMult (kas->mu[31],kas->mu[32], kas->h) ;
  MX u123 = mxMatMult (u12,kas->mu[33], kas->h) ;
  MX u1234 = mxMatMult (u123,kas->mu[41], kas->h) ;
  MX u12345 = mxMatMult (u1234,kas->mu[42], kas->h) ;
  MX u123456 = mxMatMult (u12345,kas->mu[43], kas->h) ;

  mxShow (u123456) ;

  
  /* on the h.w. we get
   *  1/8 (2b - (2a1 + a2))(2b -( 3+a1-a2)) (2b -(6 + a1 + 2a2))
   * the other combinations give the same values on the other weights
   * up to lower terms
   */
    
  return ;
} /* checkCasimirs */


/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// su21: Construction of su(2/1) representations and Feynman diagrams\n"
	    "// Authors: Jean Thierry-Mieg, NCBI, 2020-, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Construct the matrices of irreducible and indecomposable representations\n"
	    "// Construct the Casimirs, super Casimirs, Gorelik ghost Casimir\n"
	    "// \n"
	    "// Also compute the anomalies and the feynman diagrams supporting my JHEP su(2/1) papers\n"
	    "//\n"
	    "// Syntax:\n"
	    "// su31 [options]\n"
	    "//   [] [-h] [-help] [--help] : this message\n"
	    "// A: Representations\n"
	    "//   su21 -a1 <int> -a2 <int> -b <int> [-N <int>]\n"
	    "//     export the matrices, Casimirs and verifications for the module with \n"
	    "//     Dynkin lables (a1,a2,b), a1,a2 positive integers, b signed integer\n"
	    "//     Number of generations N (N >= 2)\n"
	    "//       In theory, b can be any complex number,\n"
	    "//     for numerical convenience, we restrict here to signed integers\n"
	    "//     but the formulas like the Casimir eigen values are anlytic in b\n"
	    "//       When a or N are large, many outputs are suppressed, try first a<=3, N<=3\n"
	    "//\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  KAS kas ;
  
  freeinit () ;

  memset (&kas, 0, sizeof(kas)) ;
  
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  /*   BOOL SU3 = getCmdLineBool (&argc, argv, "-su3") ; */
  int NN = 0 ;
    
  getCmdLineInt (&argc, argv, "-N", &NN) ; /* Number of generations >= 2 */
  getCmdLineInt (&argc, argv, "-NN", &NN) ; /* synonim */

  int a1 = 0, a2 = 0, b = 0 ;

  getCmdLineInt (&argc, argv, "-a1", &a1) ;
  getCmdLineInt (&argc, argv, "-a2", &a2) ;
  getCmdLineInt (&argc, argv, "-b", &b) ;

  if (a1 < 0)
    usage ("SU(3) Dynkin weigth a should be a positiver integer") ; 
  if (a2 < 0)
    usage ("SU(3) Dynkin weigth a should be a positiver integer") ; 
  if (NN != 0 && NN < 2)
    usage ("The number of generations N should be an integer >= 2") ;

  if (!NN && (a1 + a2 || b))
    {
      /* 2021_03_18 
       * construct the 15 matrices for the generic irreps of su(3/1) with h.w. (a1,a2,b)
       * verify all commutations relations
       * compute the casimir tensors and operators 

       */

      kas.b = b ;
      kas.a1 = a1 ;
      kas.a2 = a2 ;
      constructSU3Matrices (&kas) ;
      checkCommutators (&kas) ;
      ghostCasimirs (&kas) ;
      lowerMetric (&kas) ;
      upperMetric (&kas) ;
      KillingCasimir (&kas) ;
#ifdef JUNK
      Kasimirs (1,0,1, FALSE) ;  /* adjopint */
      Kasimirs (1,0,0, FALSE) ;  /* fundamental */
      Kasimirs (b,a1,a2, TRUE) ;  
  if (0) muInit (h) ;   /* init the 4x4 matrices */
  if (0) muInit2 (h) ;  /* init the 2-families 8x8 rotated matrices */
  if (NN >= 2) muInitNMarcu (a,b, NN) ;  /* init the 2-families 8x8 marcu indecomposable matrices */

#endif
    }
  /* always init, otherwise the gcc linker is unhappy */

  ac_free (h) ;   
  return 0 ;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
#ifdef JUNK

Cher Jerome,

  Tout d'abord, j'espere que bebe va bien ainsi que les parents !
  
  Ceci dit, je suis bien ennuye de ton silence sur le draft slmn

  A:
    De mon point de vue, le resultat est tyres interessant, je me pose cette question activemnent
  depuis l'apparition des articles de Coquereaux ... en 89
  a) peut on construire des indecomposables a N>3 generations dans sl(2/1)
  b) peut on construire3 des indecomposables a3 3 geenrations dans sl(5/1)

    ce qui se trauit en physique par:
    a') sl(2/1) impose t'elle qu'il n'exsite pas plus de 3 generation (ce qui semble experimentallement demontre)
    b') le modle super-grand-unifie sl(5/1) qui unifie su(3) et su(2/1) ademet il les 3 generations

      D'autre part notre article recent avec Peter sur la cohomologie ne traite que le cas N=2, donc ne repond pas aux 2 questions a' et b'

  B:
    La solution proposee, utiliser la derivee u' des matrices impaires par rapport a l'hypercharge est tres simp-le et elegante

  C:
    Je pense que ces resultats sont nouveaux


  J'ai donc redige un brouillon dans mon style de physicien en deux jours et depuis j'edite des phrases par ci par la, sans trouver de defaut majeur.

 ====

  Je trouve que tu merites tout a fait de signer car je n'aurais jamais trouve l'idee de u' si tu ne m'avais pas envoye les matrices P
  qui dans le cas de sl(2/1) conduisait a v' (partie non bloc diagonale de v) valait zero
  Donc bien il me faut savoir si tu veux signer
  Dans les deux cas j'aimerais tes commentaires
      tu peux trouver que c'est deja connu, c'est faux, c'est mal regige c'est ok
      bien sur si ce'st faux ou deja connu, moi non plus je ne veux pas signer !

 ===

  avec mes remerciement pour le beau changement de variable P et toutes mes amities

jean
	 
#endif
