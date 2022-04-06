/*  File: matrix.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg yand R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Extension of Array, to define typed matrix and tensors
 *              Bool, Int, Float, Complex
 *              Bool is the same as Int modulo 2
 *
 *    matrices are used in particular in the Deep learning library
 *
 * Exported functions:
 *              See Header file: matrix.h
 * Public read-only subfields (can only be altered at creation time)
 *    rank : 1 = vector, 2 = matrix, 3... = tensor
 *    shape : int[], a list of rank positive integers 
 *    mxAdd, mxTranspose, mxDot, mxPointMultiply
 *    
 * HISTORY:
 * Last edited:
 * Created: Feb 3, 2018
 *-------------------------------------------------------------------
 */

#define MAXRANK 64

#include "ac.h"
#include <complex.h>

/* SIYNCHRONIZE VERY CAREFULLY WITH matrix.h matrix.c matrix_slices.c */
/* Type of a matrix */
typedef enum {
  MX_NULL = 0, MX_BOOL, MX_INT, MX_FLOAT, MX_COMPLEX
} MX_TYPE ;

typedef struct mxStruct *MX ;
typedef struct mxStruct { 
  const char *name ;
  const MX_TYPE type ;
  const int rank ;
  const int shapes[MAXRANK] ;
  const int slice0, slice1 ;

  const void *magic ;
  complex float *zc ;
  float *zf ;
  int *zi ;
  const int size ;
  const int typeSize ;
  const int deltas[MAXRANK] ; 
} MM ;

static const char *typeName[] = { "NULL", "BOOL", "INT", "FLOAT", "COMPLEX", "overflow"} ;

#define MATRIX_IS_DEFINED
#include "matrix.h"

/**********************************************************************/
/**********************************************************************/

MX mxSliceLastIndex (MX a, MX b, int s0, int s1)
{ 
  int i, delta, rank, gaMax, gbMax ;
  mxCheck (a, "mxSliceLastIndex (a,...)") ;
  mxCheck (b, "mxSliceLastIndex (...,b,...)") ;
  
  rank = b->rank ;
  if (a->type != b->type)
    messcrash ("a->tye != b->type in  mxSliceLastIndex") ;
  if (a->rank != b->rank)
    messcrash ("a->tye != b->type in  mxSliceLastIndex") ;
  gaMax = a->shapes [rank - 1] ;
  gbMax = b->shapes [rank - 1] ;
  if (s1 < s0)
    messcrash ("(s1 = %d) < (s0 = %d) in mxSliceLastIndex", s1, s0) ;
  if (s0 < 0)
    messcrash ("(s0 = %d) < 0 ) in mxSliceLastIndex", s0) ;
  if (s0 >= gbMax)
    messcrash ("(s0 = %d) >= b->shapes[lastSize] = %d)) in mxSliceLastIndex", s0, gbMax) ;
  if (s0 >= gbMax)
    messcrash ("(s0 = %d) >= b->shapes[lastSize] = %d)) in mxSliceLastIndex", s0, gbMax) ;

  for (i = 0 ; i < rank - 1 ; i++)
    if (a->shapes[i] != b->shapes[i])
      messcrash ("a->shapes[%d] != b->shapes[%d] in  mxSliceLastIndex", i, i) ;
  if (gaMax != s1 - s0)
    messcrash ("a->shapes[rak-1] != s1 - s0 in  mxSliceLastIndex") ;

  delta = b->deltas[b->rank - 1] ;
  if (s1 > b->shapes[rank - 1])
    s1 = b->shapes[rank - 1] ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
	{
	  int * Restrict zia = a->zi ;
	  int * Restrict zib = b->zi + s0 * delta ;
	  memcpy (zia, zib, (s1 - s0) * delta * b->typeSize) ;
	}
      break ;
    case MX_FLOAT:
	{
	  float * Restrict zfa = a->zf ;
	  float * Restrict zfb = b->zf + s0 * delta ;
	  memcpy (zfa, zfb, (s1 - s0) * delta * b->typeSize) ;
	}
      break ;
    case MX_COMPLEX:
	{
	  complex float * Restrict zca = a->zc ;
	  complex float * Restrict zcb = b->zc + s0 * delta ;
	  memcpy (zca, zcb, (s1 - s0) * delta * b->typeSize) ;
	}
      break ;
    }  
  return a ;
} /* mxSliceLastIndex */

/**********************************************************************/
/* Matrix addition */
MX mxLinearCombine (MX a, complex float betaC, MX b, complex float gammaC, MX c, AC_HANDLE h)
{
  int i, ka, kb, kc, kaMax, kbMax, kcMax ;
  float betaR = creal (betaC) ;
  float gammaR = creal (gammaC) ;

  int * Restrict zia, * Restrict zib, * Restrict zic ;
  float * Restrict zfa, * Restrict zfb, * Restrict zfc ;
  complex float * Restrict zca, * Restrict zcb, * Restrict zcc ;

  if (a) mxCheck (a, "mxAdd(a,...)") ;
  mxCheck (b, "mxAdd(...,b,...)") ;
  mxCheck (c, "mxAdd(...,c)") ;

  if (b->rank < c->rank) /* switch, addition is Abelian */
    { MX m = b ; b = c ; c = m ; }

  
  if (!a)
    {
      a = mxCreate (h, hprintf (h, "%s +- %s", b->name, c->name), b->type, -999, b->shapes) ;
    }


  /* check the types */
  if (a->type < b->type)
     messcrash ("mxAdd A = B + C :: %s = %s + %s\nMatrix type mismatch, matrix A was received predeclared and cannot hold the results since  A->type = %s < B->type = %s"
	       , a->name
	       , b->name
	       , c->name
	       , typeName[a->type]
	       , typeName[b->type]
		) ;
  if (a->type < c->type)
     messcrash ("mxAdd A = B + C :: %s = %s + %s\nMatrix type mismatch, matrix A was received predeclared and cannot hold the results since  A->type = %s < C->type = %s"
	       , a->name
	       , b->name
	       , c->name
	       , typeName[a->type]
	       , typeName[c->type]
		) ;
  /* check b can broadcast into a */
  if (b->rank > a->rank)
    messcrash ("mxAdd A = B + C :: %s = %s + %s, A->rank = %d < B->rank = %d"
	       , a->name
	       , b->name
	       , c->name
	       , a->rank
	       , b->rank
	       ) ;
  if (c->rank > a->rank)
    messcrash ("mxAdd A = B + C :: %s = %s + %s, A->rank = %d < C->rank = %d"
	       , a->name
	       , b->name
	       , c->name
	       , a->rank
	       , c->rank
	       ) ;
  for (i = 0 ; i < b->rank ; i++)
    if (a->shapes[i] != b->shapes[i])
      messcrash ("mxAdd A = B + C :: %s = %s + %s, A->shapes[%d] = %d < B->shapes[%d] = %d"
		 , a->name
		 , b->name
		 , c->name
		 , i
		 , a->shapes[i]
		 , i
		 , b->shapes[i]
		 ) ;
      
  for (i = 0 ; i < c->rank - 1 ; i++)
    if (a->shapes[i] != c->shapes[i])
      messcrash ("mxAdd A = B + C :: %s = %s + %s, A->shapes[%d] = %d != C->shapes[%d] = %d"
		 , a->name
		 , b->name
		 , c->name
		 , i
		 , a->shapes[i]
		 , i
		 , c->shapes[i]
		 ) ;
      
  i = c->rank - 1 ;
    if (a->shapes[i] != c->shapes[i] && c->shapes[i] != 1)
      messcrash ("mxAdd A = B + C :: %s = %s + %s, A->shapes[%d] = %d != C->shapes[%d] = %d"
		 , a->name
		 , b->name
		 , c->name
		 , i
		 , a->shapes[i]
		 , i
		 , c->shapes[i]
		 ) ;
      
  kaMax = a->size / b->size ;
  kbMax = b->size / c->size ;
  kcMax = c->size ;

  zia = a->zi ;
  zib = b->zi ;
  zic = c->zi ;
 
  zfa = a->zf ;
  zfb = b->zf ;
  zfc = c->zf ;

  zca = a->zc ;
  zcb = b->zc ;
  zcc = c->zc ;


  for (ka = 0 ; ka < kaMax ; ka++)
    {
      int iia = ka * b->size ;
      for (kb = 0 ; kb < kbMax ; kb++)
	{
	  int iib = kb * c->size ;

	  switch (a->type)
	    {
	    case MX_NULL:
	    case MX_BOOL:
	    case MX_INT:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zia = a->zi + iia + iib ;
		      zib = b->zi + iib ;
		      zic = c->zi ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zia[kc] = zib[kc] * betaR + zic[kc] * gammaR  ;
		      break ;
		    case MX_FLOAT:
		      break ;
		    case MX_COMPLEX:
		      break ;
		    }
		  break ;
		default:
		  break ;
		}
	      break ;
	    case MX_FLOAT:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zfa = a->zf + iia + iib ;
		      zib = b->zi + iib ;
		      zic = c->zi ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zfa[kc] = zib[kc] * betaR + zic[kc] * gammaR  ;
		      break ;
		    case MX_FLOAT:
		      zfa = a->zf + iia + iib ;
		      zib = b->zi + iib ;
		      zfc = c->zf ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zfa[kc] = zib[kc] * betaR + zfc[kc] * gammaR  ;
		      break ;
		    case MX_COMPLEX:
		      break ;
		    default:
		      break ;
		    }
		  
		  break ;
		case MX_FLOAT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zfa = a->zf + iia + iib ;
		      zfb = b->zf + iib ;
		      zic = c->zi ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zfa[kc] = zfb[kc] * betaR + zic[kc] * gammaR  ;
		      break ;
		    case MX_FLOAT:
		      zfa = a->zf + iia + iib ;
		      zfb = b->zf + iib ;
		      zfc = c->zf ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zfa[kc] = zfb[kc] * betaR + zfc[kc] * gammaR  ;
		      break ;
		    case MX_COMPLEX:
		      break ;
		    default:
		      break ;
		    }
		  break ;
		default:
		  break ;
		}
	      break ;
	    case MX_COMPLEX:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zca = a->zc + iia + iib ; 
		      zib = b->zi + iib ;
		      zic = c->zi ;
 		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zib[kc] * betaC + zic[kc] * gammaC  ;
		      break ;
		    case MX_FLOAT:
		      zca = a->zc + iia + iib ;
		      zib = b->zi + iib ;
		      zfc = c->zf ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zib[kc] * betaC + zfc[kc] * gammaC  ;
		      break ;
		    case MX_COMPLEX:
		      zca = a->zc + iia + iib ;
		      zib = b->zi + iib ;
		      zcc = c->zc ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zib[kc] * betaC + zcc[kc] * gammaC  ;
		    }
		  
		  break ;
		case MX_FLOAT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zca = a->zc + iia + iib ;
		      zfb = b->zf + iib ;
		      zic = c->zi ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zfb[kc] * betaC + zic[kc] * gammaC  ;
		      break ;
		    case MX_FLOAT:
		      zca = a->zc + iia + iib ;
		      zfb = b->zf + iib ;
		      zfc = c->zf ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zfb[kc] * betaC + zfc[kc] * gammaC  ;
		      break ;
		    case MX_COMPLEX:
		      zca = a->zc + iia + iib ;
		      zfb = b->zf + iib ;
		      zcc = c->zc ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zfb[kc] * betaC + zcc[kc] * gammaC  ;
		    }
		  break ;
		case MX_COMPLEX:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zca = a->zc + iia + iib ;
		      zcb = b->zc + iib ;
		      zic = c->zi ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zcb[kc] * betaC + zic[kc] * gammaC  ;
		      break ;
		    case MX_FLOAT:
		      zca = a->zc + iia + iib ; 
		      zcb = b->zc + iib ;
		      zfc = c->zf ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zcb[kc] * betaC + zfc[kc] * gammaC  ;
		      break ;
		    case MX_COMPLEX:
		      zca = a->zc  + iia + iib ;
		      zcb = b->zc + iib ;
		      zcc = c->zc ;
		      for (kc = 0 ; kc < kcMax ; kc++)
			zca[kc] = zcb[kc] * betaC + zcc[kc] * gammaC  ;
		    }
		  break ;
		}
	      break ;
	    }
	}
    }
  return a ;
} /* mxDoAdd */

/**********************************************************************/
