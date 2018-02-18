#include "regular.h"
#include "ac.h"

typedef struct cxStruct { double x , y ;} CXs ;
typedef CXs *CX ;

/*******************************************************/
/* return a new complex c=0 */
CX cNew (AC_HANDLE h)
{
  CX c = (CX) halloc (sizeof(CXs), h) ;
  return c ;
}

/*******************************************************/
/* return a static complex c=0 */
CX cXY (double x, double y, AC_HANDLE h)
{
  CX z = cNew (h) ;
  z->x = x ; z->y = y ;
  return z ;
}

/*******************************************************/
/* kill the sixth decimal, so that sqrt(3)^2 = 3 exactly */
double cRound (double x)
{
  double y = 1 + x < 0 ? -x : x ;
  int n     = 1000 * y + .1 ;
  double dx = 1000 * y - n ;
  if (0 && dx < .01)
    {
      y = n ; y /= 1000 ;
      y-- ;
      y = x < 0 ? -y : y ;
    }
  else
    y = x ;
  return y ;    
}

/*******************************************************/

char *cShow (CX z)
{
  static char buf[256] ;
  double x = cRound (z->x) ;
  double y = cRound (z->y) ;
  if (x*x < .00000001) x = 0 ;
  if (y*y < .00000001) y = 0 ;
  if (x == 0 && y == 0) sprintf (buf, ".") ;
  else if (y == 0)  sprintf (buf, "%g", x) ;
  else if (x == 0 && y == 1)  sprintf (buf, "i") ;
  else if (x == 0 && y == -1)  sprintf (buf, "-i") ;
  else if (x == 0)  sprintf (buf, "%gi", y) ;
  else if (y == 1)  sprintf (buf, "%g+i", x) ;
  else if (y == -1)  sprintf (buf, "%g-i", x) ;
  else sprintf (buf, "%g+%gi", x,y) ;
  return buf ;
}

/*******************************************************/
CX cCopy (CX a, AC_HANDLE h)
{ 
  CX b = cNew (h) ; *b = *a ; 
  return b ;
}
     
/*******************************************************/
/* Real part */
double cR (CX a)
{ return a->x ; }

/*******************************************************/
/* Imaginary part */
double cIm (CX a)
{ return a->y ; }

/*******************************************************/
/* Norm */
double cNorm (CX a)
{ return a->x * a->x + a->y * a->y ; }

/*******************************************************/
/* Conjugates en place:  a conj= b */
CX cJJ (CX a)
{ a->y = - a->y ; return a ; }

/*******************************************************/
/* Add en place:  a += b */
CX cAA (CX a, CX b)
{ a->x += b->x ; a->y += b->y ; return a ; }

/*******************************************************/
/* Multiply en place:  a *= b */
CX cMM (CX a, CX b)
{ 
  double x, y ;
  x = a->x * b->x - a->y * b->y ;
  y = a->x * b->y + a->y * b->x ;
  a->x = x ; a->y = y ;
  return a ;
}

/*******************************************************/
/*******************************************************/
/* a matrix is made of CX numbers
   we can create, add, mutiply, conjugate, scale, trace
*/

typedef struct matrixStruct { int dim ; CX z ; } MMs ;
typedef MMs *MM ;

/*******************************************************/
/*  matrix element: returns an address IN the matrix */
#define mE(_a,_i,_j) ((CX)((char*)((_a)->z) + sizeof(CXs)*((_a)->dim * (_i) + (_j))))

/*******************************************************/
/* return a new complex m=0 */
MM mNew (int dim, AC_HANDLE h)
{
  MM m = (MM) halloc (sizeof(MMs), h) ;
  m->dim = dim ;
  m->z = (CX) halloc (dim * dim * sizeof(CXs), h) ;
  return m ;
}

/*******************************************************/
/* sets the REAL and IMAGINARY part of a matrix */
void mSet (MM a, double *x, double *y)
{
  int i, j, n = a->dim ;
  CX z ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      {
	z = mE(a,i,j) ;
	if (x) z->x = x[i * a->dim + j] ;
	if (y) z->y = y[i * a->dim + j] ;
      }
  return ;
}

/*******************************************************/
MM mCopy (MM a, AC_HANDLE h)
{ 
  MM b = mNew (a->dim, h) ;

  memcpy (b->z, a->z, a->dim * a->dim * sizeof(CXs)) ;
  return b ;
}
     
/*********************/
/* show the matrix */
static void mShow (MM a, char *title)
{
  int i, j ;
  
  printf ("\n matrix %s=\n", title) ;
  for (i = 0 ; i < 8 && i < a->dim ; i++)
    {
      for (j = 0 ; j < 8 && j < a->dim ; j++)
	{
	  printf ("\t%s", cShow(mE(a,i,j))) ;
	}
      printf ("\n") ;
    }
} /* matrixShow */

/*******************************************************/
/* Conjugates en place:  a conj= b */
void mJ (MM a)
{
  CX z ;
  int i, j, n = a->dim ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      {
	z = mE(a,i,j) ;
	cJJ (z) ;
      }
  return ;
}

/*******************************************************/
/* Scale en place:  a *= z */
MM mScale (double x, double y, MM a)
{
  CXs cc ;
  CX za, z = &cc ;
  z->x = x ; z->y = y ;
  int i, j, n = a->dim ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      {
	za = mE(a,i,j) ;
	cMM (za, z) ;
      }
  return a ;
}

/*******************************************************/
/* Trace */
CX mTrace (MM a, AC_HANDLE h)
{
  CX za, t = cNew (h) ;
  int i, n = a->dim ;
  for (i = 0 ; i < n ; i++)
    {
      za = mE(a,i,i) ;
      cAA (t, za) ;
    }
  return t ;
}

/*******************************************************/
/* Add en place:  a += b */
MM mAA (MM a, MM b)
{
  CX za, zb ;
  int i, j, n = a->dim ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      {
	za = mE(a,i,j) ;
	zb = mE(b,i,j) ;
	cAA (za, zb) ;
      }
  return a ;
}

/*******************************************************/
/* Multiply en place:  a *= b */
MM mMM (MM a, MM b)
{
  CXs cc ;
  CX za, zb, zc, z = &cc ;
  int i, j, k, n = a->dim ;
  MM c = mNew (a->dim, 0) ; 
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      {
	zc = mE(c,i,j) ;
	for (k = 0 ; k < n ; k++)
	{
	  
	  za = mE(a,i,k) ;
	  zb = mE(b,k,j) ;
	  *z = *za ;
	  cMM (z, zb) ;
	  cAA (zc, z) ;
	}
      }
  memcpy (a->z, c->z, a->dim * a->dim * sizeof(CXs)) ;
  ac_free (c->z) ;
  ac_free (c) ;
  return a ;
}

/*******************************************************/
/* Real part */
MM mR (MM a, AC_HANDLE h)
{ 
  MM b = mCopy (a, h) ;
  mJ (b) ;
  mAA (b, a) ;
  mScale (.5, 0, b) ;
  return b ;
}

/*******************************************************/
/* Imaginary part */
MM mI (MM a, AC_HANDLE h)
{ 
  MM b = mCopy (a, h) ;
  mJ (b) ;
  mScale (-1, 0, b) ;
  mAA (b, a) ;
  mScale (.5, -1, b) ;
  return b ;
}

/*******************************************************/
/*******************************************************/
/* a 4 tensor is made of CX numbers
   we can create, set parts, matrix multiply and change the rank, contract
*/
/* dim will always be 8, names wil be 1238, or 4567 */
typedef struct tensorStruct3 { int rank, dim ; char *names ; CX z ; } TT3s ;
typedef TT3s *TT3 ;
typedef struct tensorStruct4 { int rank, dim ; char *names ; CX z ; } TT4s ;
typedef TT4s *TT4 ;

/*******************************************************/
/* tensor element: returns an address IN the ensor */

#define t3E(_t,_i,_j,_k) ((_t)->z+((_t)->dim * (_t)->dim * (_i) + (_t)->dim * (_j) + (_k)))
#define t4E(_t,_i,_j,_k,_l) ((_t)->z+((_t)->dim * (_t)->dim * (_t)->dim * (_i) + (_t)->dim * (_t)->dim * (_j) + (_t)->dim * (_k) + (_l)))

static void t3Show (TT3 t3)
{
  unsigned int show, shown ;
  int ii, jj, kk, nn = t3->dim ;
  CX z ;
  
  for (ii = 0 ; ii < nn ; ii++)
    {
      shown = show = 0 ;
      if (ii == 4) show |= 0x1 ;
      for (jj = 0 ; jj < nn ; jj++)
	{
	  if (jj == 4) show |= 0x2 ;
	  else show &= ~(0x2) ;
	  for (kk = 0 ; kk < nn ; kk++)
	    {
	      if (kk == 4) show |= 0x4 ;
	      else show &= ~(0x4) ;
	      z = t3E (t3,ii,jj,kk) ;
	      if (cNorm (z) > .00001)
		{
		  printf ("\t%d:%d:%d -> %s", ii, jj, kk, cShow (z)) ;
		  shown |= 0x7 ;
		}		      
	      if (shown & 0x4) printf ("\n") ;
	      shown &= ~(0x4) ;
	    }
	  if (shown & 0x2) printf ("\n") ; 
	  shown &= ~(0x2) ;	  
	}
      if (shown) printf ("\n\n") ;
    }
} /* t3Show */

/*******************************************************/

static void t4Show (TT4 t4)
{
  unsigned int show, shown ;
  int ii, jj, kk, ll, nn = t4->dim ;
  CX z ;
  
  for (ii = 0 ; ii < nn ; ii++)
    {
      shown = show = 0 ;
      if (ii == 4) show |= 0x1 ;
      for (jj = 0 ; jj < nn ; jj++)
	{
	  if (jj == 4) show |= 0x2 ;
	  else show &= ~(0x2) ;
	  for (kk = 0 ; kk < nn ; kk++)
	    {
	      if (kk == 4) show |= 0x4 ;
	      else show &= ~(0x4) ;
	      for (ll = 0 ; ll < nn ; ll++)
		{
		  if (ll == 4) show |= 0x8 ;
		  else show &= ~(0x8) ;
		  z = t4E (t4,ii,jj,kk,ll) ;
		  if (show && cNorm (z) > .00001)
		    {
		      printf ("\t%d:%d:%d:%d -> %s", ii, jj, kk, ll, cShow (z)) ;
		      shown |= 0xf ;
		    }		      
		}
	      if (shown & 0x4) printf ("\n") ;
	      shown &= ~(0x4) ;
	    }
	  if (shown & 0x2) printf ("\n") ; 
	  shown &= ~(0x2) ;	  
	}
      if (shown) printf ("\n\n") ;
    }
} /* t4Show */

/*******************************************************/
/* Scale en place:  a *= z */
TT3 t3Scale (double x, double y, TT3 a)
{
  CXs cc ;
  CX za, z = &cc ;
  z->x = x ; z->y = y ;
  int i, j, k, n = a->dim ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	{
	  za = t3E(a,i,j,k) ;
	  cMM (za, z) ;
	}
  return a ;
}

/*******************************************************/
/* Scale en place:  a *= z */
TT4 t4Scale (double x, double y, TT4 a)
{
  CXs cc ;
  CX za, z = &cc ;
  z->x = x ; z->y = y ;
  int i, j, k, l, n = a->dim ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    za = t4E(a,i,j,k,l) ;
	    cMM (za, z) ;
	  }
  return a ;
}

/*******************************************************/
/* return a new complex null 3 tensor */
TT3 t3New (int dim, char *names, AC_HANDLE h)
{
  TT3 t = (TT3) halloc (sizeof(TT3s), h) ;
  t->rank = 3 ;
  t->dim = dim++ ;
  t->z = halloc (dim * dim * dim * sizeof(CXs), h) ;
  if (names) t->names = strnew (names, h) ;
  return t ;
}

/*******************************************************/
/* return a new complex null 4 tensor */
TT4 t4New (int dim, char *names, AC_HANDLE h)
{
  TT4 t = (TT4) halloc (sizeof(TT4s), h) ;
  t->rank = 4 ;
  t->dim = dim++ ;
  t->z = halloc (dim * dim * dim * dim * sizeof(CXs), h) ;
  if (names) t->names = strnew (names, h) ;
  return t ;
}

/*******************************************************/

TT3 t3Copy (TT3 a, char *names, AC_HANDLE h)
{ 
  TT3 b = t3New (a->dim, names, h) ;

  memcpy (b->z, a->z, a->dim * a->dim * a->dim * sizeof(CXs)) ;
  return b ;
}
     
/*******************************************************/

TT4 t4Copy (TT4 a, char *names, AC_HANDLE h)
{ 
  TT4 b = t4New (a->dim, names, h) ;

  memcpy (b->z, a->z, a->dim * a->dim * a->dim * a->dim * sizeof(CXs)) ;
  return b ;
}
     
/*******************************************************/
/* Add en place:  a += b */
TT3 t3AA (TT3 a, TT3 b)
{
  CX za, zb ;
  int i, j, k, n = a->dim ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	{
	  za = t3E(a,i,j,k) ;
	  zb = t3E(b,i,j,k) ;
	  cAA (za, zb) ;
	}
  return a ;
}

/*******************************************************/
/* Add en place:  a += b */
TT4 t4AA (TT4 a, TT4 b)
{
  CX za, zb ;
  int i, j, k, l, n = a->dim ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    za = t4E(a,i,j,k,l) ;
	    zb = t4E(b,i,j,k,l) ;
	    cAA (za, zb) ;
	  }
  return a ;
}

/*******************************************************/
/* Roll the index order left  b(jkli) = a(ijkl) */
TT4 t4RollLeft (TT4 a, char *names, AC_HANDLE h)
{
  CX za, zb ;
  TT4 b = t4New (a->dim, names, h) ;
  int i, j, k, l, n = a->dim ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    za = t4E(a,i,j,k,l) ;
	    zb = t4E(b,j,k,l,i) ;
	    *zb = *za ;
	  }
  return b ;
}

/*******************************************************/
/* Symmetrize the 2 pairs of indices */
void t4DoubleSym (TT4 a)
{
  AC_HANDLE h = ac_new_handle () ;
  CX z, z0 ;
  int ii,jj,kk,ll, n = a->dim ;
  TT4 tt = t4New (n, "ijkl", h) ;

  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z0 = t4E (a,ii,jj,kk,ll) ;

	    z = t4E(tt,ii,jj,kk,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,jj,ll,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,ii,kk,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,ii,ll,kk) ;
	    cAA (z, z0) ;
	  }
  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z = t4E(tt,ii,jj,kk,ll) ;
	    z0 = t4E(a,ii,jj,kk,ll) ;
	    z0->x = z->x/4 ; z0->y = z->y/4 ;
	  }
  ac_free (h) ;
  return ;
} /*t4DoubleSym */

/*******************************************************/
/* Symmetrize the LRLR indices in LL and RR */
void t4HalfSym (TT4 a)
{
  AC_HANDLE h = ac_new_handle () ;
  CX z, z0 ;
  int ii,jj,kk,ll, n = a->dim ;
  TT4 tt = t4New (n, "ijkl", h) ;

  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z0 = t4E (a,ii,jj,kk,ll) ;

	    z = t4E(tt,ii,jj,kk,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,ll,kk,jj) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,jj,ii,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,ll,ii,jj) ;
	    cAA (z, z0) ;
	  }
  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z = t4E(tt,ii,jj,kk,ll) ;
	    z0 = t4E(a,ii,jj,kk,ll) ;
	    z0->x = z->x/4 ; z0->y = z->y/4 ;
	  }
  ac_free (h) ;
  return ;
} /*t4HalfSym */

/*******************************************************/
/* Fully symmetrize the 4 indices */
void t4FullSym (TT4 a)
{
  AC_HANDLE h = ac_new_handle () ;
  CX z, z0 ;
  int ii,jj,kk,ll, n = a->dim ;
  TT4 tt = t4New (n, "ijkl", h) ;

  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z0 = t4E (a,ii,jj,kk,ll) ;
	    z = t4E(tt,ii,jj,kk,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,jj,ll,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,kk,jj,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,kk,ll,jj) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,ll,jj,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ii,ll,kk,jj) ;
	    cAA (z, z0) ;
	    
	    z = t4E(tt,jj,ii,kk,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,ii,ll,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,kk,ii,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,kk,ll,ii) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,ll,ii,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,jj,ll,kk,ii) ;
	    cAA (z, z0) ;
	    
	    z = t4E(tt,kk,jj,ii,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,jj,ll,ii) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,ii,jj,ll) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,ii,ll,jj) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,ll,jj,ii) ;
	    cAA (z, z0) ;
	    z = t4E(tt,kk,ll,ii,jj) ;
	    cAA (z, z0) ;
	    
	    z = t4E(tt,ll,jj,kk,ii) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ll,jj,ii,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ll,kk,jj,ii) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ll,kk,ii,jj) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ll,ii,jj,kk) ;
	    cAA (z, z0) ;
	    z = t4E(tt,ll,ii,kk,jj) ;
	    cAA (z, z0) ;
	  }
  for (ii = 0 ; ii < n ; ii++)
    for (jj = 0 ; jj < n ; jj++)
      for (kk = 0 ; kk < n ; kk++) 
	for (ll = 0 ; ll < n ; ll++)
	  {
	    z = t4E(tt,ll,jj,kk,ii) ;
	    z0 = t4E(a,ii,jj,kk,ll) ;
	    z0->x = z->x/24 ; z0->y = z->y/24 ;
	  }
  ac_free (h) ;
  return ;
} /*t4FullSym */

/*******************************************************/
/* Returns a 4 tensor by contracting one leg with a metric
   by contracting the columns a1<->b1
   names are the name sof the 4 remaining indices 
*/
TT4 t4t2contract (TT4 a, MM b, char *names, int a1, int b1, AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, z = &cc ;
  int i, j, k, l, i1, n = a->dim ;
  TT4 c = t4New (a->dim, names, h) ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    zc = t4E (c,i,j,k,l) ;
	    for (i1 = 0 ; i1 < n ; i1++)
	      {
		if (a1== 4) za = t4E(a,i,j,k,i1) ;
		else if (a1== 3) za = t4E(a,i,j,i1,l) ;
		else if (a1== 2) za = t4E(a,i,i1,j,l) ;
		else if (a1== 1) za = t4E(a,i1,j,k,l) ;

		if (b1==1) zb = mE(b,i1,j) ;
		else if (b1==2) zb = mE(b,i,i1) ;

		*z = *za ;
		cMM (z, zb) ;
		cAA (zc, z) ;
	      }
	  }
  return c ;
}

/*******************************************************/
/* returns a double contracted 4 tensor 
   by contracting the columns a1<->b1 and a2<->b2
   names are the name sof the 4 remaining indices 
*/
TT4 t4t4doubleContract (TT4 a, TT4 b, char *names, int a1, int a2, int b1, int b2, AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, z = &cc ;
  int i, j, k, l, i1, j1, n = a->dim ;
  TT4 c = t4New (a->dim, names, h) ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    zc = t4E (c,i,j,k,l) ;
	    for (i1 = 0 ; i1 < n ; i1++)
	      for ( j1 = 0 ; j1 < n ; j1++)
		{
		  if (a1==3 && a2==4) za = t4E(a,i,j,i1,j1) ;
		  else if (a1==4 && a2==3) za = t4E(a,i,j,j1,i1) ;
		  else if (a1==2 && a2==4) za = t4E(a,i,i1,j,j1) ;
		  else if (a1==4 && a2==2) za = t4E(a,i,j1,j,i1) ;
		  else if (a1==1 && a2==4) za = t4E(a,i1,i,j,j1) ;
		  else if (a1==4 && a2==1) za = t4E(a,j1,i,j,i1) ;
		  else if (a1==2 && a2==3) za = t4E(a,i,i1,j1,j) ;
		  else if (a1==3 && a2==2) za = t4E(a,i,j1,i1,j) ;
		  else if (a1==1 && a2==3) za = t4E(a,i1,i,j1,j) ;
		  else if (a1==3 && a2==1) za = t4E(a,j1,i,i1,j) ;
		  else if (a1==1 && a2==2) za = t4E(a,i1,j1,i,j) ;
		  else if (a1==2 && a2==1) za = t4E(a,j1,i1,i,j) ;
	
		  if (b1==3 && b2==4) zb = t4E(b,k,l,i1,j1) ;
		  else if (b1==4 && b2==3) zb = t4E(b,k,l,j1,i1) ;
		  else if (b1==2 && b2==4) zb = t4E(b,k,i1,l,j1) ;
		  else if (b1==4 && b2==2) zb = t4E(b,k,j1,l,i1) ;
		  else if (b1==1 && b2==4) zb = t4E(b,i1,k,l,j1) ;
		  else if (b1==4 && b2==1) zb = t4E(b,j1,k,l,i1) ;
		  else if (b1==2 && b2==3) zb = t4E(b,k,i1,j1,l) ;
		  else if (b1==3 && b2==2) zb = t4E(b,k,j1,i1,l) ;
		  else if (b1==1 && b2==3) zb = t4E(b,i1,k,j1,l) ;
		  else if (b1==3 && b2==1) zb = t4E(b,j1,k,i1,l) ;
		  else if (b1==1 && b2==2) zb = t4E(b,i1,j1,k,l) ;
		  else if (b1==2 && b2==1) zb = t4E(b,j1,i1,k,l) ;

		  *z = *za ;
		  cMM (z, zb) ;
		  cAA (zc, z) ;
		}
	  }
  return c ;
}

/*******************************************************/
/* returns a 4 tensor by contracting 2 3-tensors 
   by contracting the columns a1<->b1
   names are the name sof the 4 remaining indices 
*/
TT4 t3t3contract (TT3 a, TT3 b, char *names, int a1, int b1, AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, z = &cc ;
  int i, j, k, l, i1, n = a->dim ;
  TT4 c = t4New (a->dim, names, h) ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    zc = t3E (c,i,j,k) ;
	    for (i1 = 0 ; i1 < n ; i1++)
	      {
		if (a1==3) za = t3E(a,i,j,i1) ;
		else if (a1==2) za = t3E(a,i,i1,j) ;
		else if (a1==1) za = t3E(a,i1,i,j) ;
		
		if (b1==3) zb = t3E(b,k,l,i1) ;
		else if (b1==2) zb = t3E(b,k,i1,l) ;
		else if (b1==1) zb = t3E(b,i1,k,l) ;
		
		*z = *za ;
		cMM (z, zb) ;
		cAA (zc, z) ;
	      }
	  }
  return c ;
}

/*******************************************************/
/* returns a 3 tensor by contracting 3 3-tensors 
   by contracting the columns a2<->b1, b2<->c1, c2<->a1 
   names are the name sof the 4 remaining indices 
*/
TT3 t3t3t3tripleContract (TT3 a, TT3 b, TT3 c, char *names
			  , int a2, int b1, int b2, int c1, int c2, int a1
			  , AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, zd = 0, z = &cc ;
  int i, j, k, i1, j1, k1, n = a->dim ;
  TT3 d = t3New (a->dim, names, h) ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	{
	  zd = t3E (d,i,j,k) ;
	  for (i1 = 0 ; i1 < n ; i1++) /* a2<->b1 */
	    for (j1 = 0 ; j1 < n ; j1++) /* b2<->c1 */
	      for (k1 = 0 ; k1 < n ; k1++) /* c2<->a1 */
	    {
	      if (a1==2 && a2==3) za = t3E(a,i,k1,i1) ;
	      else if (a1==3 && a2==2) za = t3E(a,i,i1,k1) ;
	      else if (a1==1 && a2==3) za = t3E(a,k1,i,i1);
	      else if (a1==3 && a2==1) za = t3E(a,i1,i,k1) ;
	      else if (a1==1 && a2==2) za = t3E(a,k1,i1,i) ;
	      else if (a1==2 && a2==1) za = t3E(a,i1,k1,i) ;
	
	      if (b1==2 && b2==3) zb = t3E(b,j,i1,j1) ;
	      else if (b1==3 && b2==2) zb = t3E(b,j,j1,i1) ;
	      else if (b1==1 && b2==3) zb = t3E(b,i1,j,j1);
	      else if (b1==3 && b2==1) zb = t3E(b,j1,j,i1) ;
	      else if (b1==1 && b2==2) zb = t3E(b,i1,j1,j) ;
	      else if (b1==2 && b2==1) zb = t3E(b,j1,i1,j) ;
	
	      if (c1==2 && c2==3) zc = t3E(c,k,j1,k1) ;
	      else if (c1==3 && c2==2) zc = t3E(c,k,k1,j1) ;
	      else if (c1==1 && c2==3) zc = t3E(c,j1,k,k1);
	      else if (c1==3 && c2==1) zc = t3E(c,k1,k,j1) ;
	      else if (c1==1 && c2==2) zc = t3E(c,j1,k1,k) ;
	      else if (c1==2 && c2==1) zc = t3E(c,k1,j1,k) ;

	      *z = *za ;
	      cMM (z, zb) ;
	      cMM (z, zc) ;
	      cAA (zd, z) ;
	    }
	}
  return d ;
}

/*******************************************************/
/* returns a 4 tensor by contracting 2 3-tensors and a 4 tensor 
   by contracting the columns a2<->b1, b2<->c1, c2<->a1 
   names are the name sof the 4 remaining indices 
*/
TT4 t3t3t4tripleContract (TT3 a, TT3 b, TT4 c, char *names
			  , int a2, int b1, int b2, int c1, int c2, int a1
			  , AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, zd = 0, z = &cc ;
  int i, j, k, l, i1, j1, k1, n = a->dim ;
  TT4 d = t4New (a->dim, names, h) ;
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    zd = t4E (d,i,j,k,l) ;
	    for (i1 = 0 ; i1 < n ; i1++) /* a2<->b1 */
	      for (j1 = 0 ; j1 < n ; j1++) /* b2<->c1 */
		for (k1 = 0 ; k1 < n ; k1++) /* c2<->a1 */
		  {
		    if (a1==2 && a2==3) za = t3E(a,i,k1,i1) ;
		    else if (a1==3 && a2==2) za = t3E(a,i,i1,k1) ;
		    else if (a1==1 && a2==3) za = t3E(a,k1,i,i1);
		    else if (a1==3 && a2==1) za = t3E(a,i1,i,k1) ;
		    else if (a1==1 && a2==2) za = t3E(a,k1,i1,i) ;
		    else if (a1==2 && a2==1) za = t3E(a,i1,k1,i) ;
		    
		    if (b1==2 && b2==3) zb = t3E(b,j,i1,j1) ;
		    else if (b1==3 && b2==2) zb = t3E(b,j,j1,i1) ;
		    else if (b1==1 && b2==3) zb = t3E(b,i1,j,j1);
		    else if (b1==3 && b2==1) zb = t3E(b,j1,j,i1) ;
		    else if (b1==1 && b2==2) zb = t3E(b,i1,j1,j) ;
		    else if (b1==2 && b2==1) zb = t3E(b,j1,i1,j) ;
		    
		    if (c1==3 && c2==4) zc = t4E(c,k,l,j1,k1) ;
		    else if (c1==4 && c2==3) zc = t4E(c,k,l,k1,j1) ;
		    else if (c1==2 && c2==4) zc = t4E(c,k,j1,l,k1) ;
		    else if (c1==4 && c2==2) zc = t4E(c,k,k1,l,j1) ;
		    else if (c1==1 && c2==4) zc = t4E(c,j1,k,l,k1) ;
		    else if (c1==4 && c2==1) zc = t4E(c,k1,k,l,j1) ;
		    else if (c1==2 && c2==3) zc = t4E(c,k,j1,k1,l) ;
		    else if (c1==3 && c2==2) zc = t4E(c,k,k1,j1,l) ;
		    else if (c1==1 && c2==3) zc = t4E(c,j1,k,k1,l) ;
		    else if (c1==3 && c2==1) zc = t4E(c,k1,k,j1,l) ;
		    else if (c1==1 && c2==2) zc = t4E(c,j1,k1,k,l) ;
		    else if (c1==2 && c2==1) zc = t4E(c,k1,j1,k,l) ;
		    
		    *z = *za ;
		    cMM (z, zb) ;
		    cMM (z, zc) ;
		    cAA (zd, z) ;
		  }
	  }
  return d ;
}

/*******************************************************/
/* returns a 4 tensor by contracting 4 3-tensors
   by contracting the columns a2<->b1, b2<->c1, c2<->d1, d2<->a1
   names are the name sof the 4 remaining indices 
*/
TT4 t3quadrupleContract (TT3 a, TT3 b, TT3 c, TT3 d, char *names
			 , int a2, int b1, int b2, int c1, int c2, int d1, int d2, int a1
			 , AC_HANDLE h)
{
  CXs cc ;
  CX  za = 0, zb = 0, zc = 0, zd = 0, ze = 0, z = &cc ;
  int i, j, k, l, i1, j1, k1, l1, n = a->dim ;
  TT4 e = t4New (a->dim, names, h) ;

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      for (k = 0 ; k < n ; k++)
	for (l = 0 ; l < n ; l++)
	  {
	    zd = t4E (d,i,j,k,l) ;
	    for (i1 = 0 ; i1 < n ; i1++) /* a2<->b1 */
	      for (j1 = 0 ; j1 < n ; j1++) /* b2<->c1 */
		for (k1 = 0 ; k1 < n ; k1++) /* c2<->d1 */
		  for (l1 = 0 ; l1 < n ; l1++) /* d2<->a1 */
		    {
		      if (a1==2 && a2==3) za = t3E(a,i,l1,i1) ;
		      else if (a1==3 && a2==2) za = t3E(a,i,i1,l1) ;
		      else if (a1==1 && a2==3) za = t3E(a,l1,i,i1);
		      else if (a1==3 && a2==1) za = t3E(a,i1,i,l1) ;
		      else if (a1==1 && a2==2) za = t3E(a,l1,i1,i) ;
		      else if (a1==2 && a2==1) za = t3E(a,i1,l1,i) ;
		      
		      if (b1==2 && b2==3) zb = t3E(b,j,i1,j1) ;
		      else if (b1==3 && b2==2) zb = t3E(b,j,j1,i1) ;
		      else if (b1==1 && b2==3) zb = t3E(b,i1,j,j1);
		      else if (b1==3 && b2==1) zb = t3E(b,j1,j,i1) ;
		      else if (b1==1 && b2==2) zb = t3E(b,i1,j1,j) ;
		      else if (b1==2 && b2==1) zb = t3E(b,j1,i1,j) ;
		    
		      if (c1==2 && c2==3) zc = t3E(c,k,j1,k1) ;
		      else if (c1==3 && c2==2) zc = t3E(c,k,k1,j1) ;
		      else if (c1==1 && c2==3) zc = t3E(c,j1,k,k1);
		      else if (c1==3 && c2==1) zc = t3E(c,k1,k,j1) ;
		      else if (c1==1 && c2==2) zc = t3E(c,j1,k1,j) ;
		      else if (c1==2 && c2==1) zc = t3E(c,k1,j1,j) ;

		      if (d1==2 && d2==3) zd = t3E(d,k,k1,l1) ;
		      else if (d1==3 && d2==2) zd = t3E(d,k,l1,k1) ;
		      else if (d1==1 && d2==3) zd = t3E(d,k1,k,l1);
		      else if (d1==3 && d2==1) zd = t3E(d,l1,k,k1) ;
		      else if (d1==1 && d2==2) zd = t3E(d,k1,l1,j) ;
		      else if (d1==2 && d2==1) zd = t3E(d,l1,k1,j) ;

		      *z = *za ;
		      cMM (z, zb) ;
		      cMM (z, zc) ;
		      cMM (z, zd) ;
		      cAA (ze, z) ;
		    }
	  }
  return e ;
} /* t3quadrupleContract */

/*******************************************************/
/*******************************************************/


/* calcul matriciel du vertex des scalaires de su(2/1) */
/****************************************************************************/
/*************************  utilities ***************************************/

/* in C, accessing matrices is confusing, i create a macro valid
 * for matrices of doubles of size NMMU by NMMU
 */

/*********************/
/* create the su(2/1) matrices
 * dim 9, because labelled [1,8] as usual not as in C
 * but the lines and cols are labelled as usual [0,3[, [0,4], [0,15[
 */
MM M3[9], M3B[9], M4[9], M15[9], MQQ[9], MEE[9], PLQQ, PRQQ, PL15, PR15, M15L[9], M15R[9], MA[9], MB[9], MH[9], MK[9] ;
TT3 daij, faij ;
TT4 phi4, phi4a, phi4b, phi4f, phi4d, phi4c, phi4fd, H4, K4, HHKK, HKKH, KKHH, HKHK, ZHHKK ;
TT4 H4a, K4a, HHKKa, HKKHa, KKHHa, HKHKa ;
TT4 fHHKK, fKKKK, fHHHH ;

static BOOL su21ConstructFamillyMatrix (AC_HANDLE h)
{
  int ii, i, j, jj ;
  MM m, m15 ;
  double s2 = sqrt(2.0), s3 = sqrt(3.0) ;
  CXs un = {1,0} ;
  CXs moinsUn = {-1.0,0} ;
  CXs deux = {2.0,0} ;
  CXs moinsDeux = {-2.0,0} ;
  CXs plusI = {0,1.0};
  CXs moinsI = {0,-1.0};

  PL15 = mNew (15, h) ;
  m = PL15 ;
  *mE(m , 0, 0) = un ; *mE (m , 1, 1) = un ;
  for (jj = 3 ; jj < 15 ; jj += 4) 
    {
      *mE (m , jj + 1, jj + 1) = un ;   
      *mE (m , jj + 2, jj + 2) = un ;   
    }

  PR15 = mNew (15, h) ;
  m = PR15 ;
  *mE (m , 2, 2) = un ;   
  for (jj = 3 ; jj < 15 ; jj += 4) 
    {
      *mE (m , jj + 0, jj + 0) = un ;   
      *mE (m , jj + 3, jj + 3) = un ;   
    }

  for (ii = 1 ; ii <= 8 ; ii++)
    {
      M3[ii] = mNew (4,h) ;
      M3B[ii] = mNew (4,h) ;
      M4[ii] = mNew (4,h) ;
      M15[ii] = mNew (15,h) ;
     }

  m = M3[1] ; *mE (m , 0, 1) = *mE (m ,1,0) = un ; 
  m = M3[2] ; *mE (m , 0, 1) = moinsI ; *mE (m ,1,0) = plusI ;
  m = M3[3] ; *mE (m , 0, 0) = un ; *mE (m , 1, 1) = moinsUn ; 
  m = M3[8] ; *mE (m , 0, 0) = un ; *mE (m , 1, 1) = un ; *mE (m , 2, 2) = deux ;

  m = M3[4] ; *mE (m , 0, 2) = *mE (m ,2,0) = un ; 
  m = M3[5] ; *mE (m , 0, 2) = moinsI ; *mE (m ,2,0) = plusI ;
  m = M3[6] ; *mE (m , 1, 2) = un ; *mE (m , 2, 1) = un ; 
  m = M3[7] ; *mE (m , 1, 2) = moinsI ; *mE (m , 2, 1) = plusI ;

  m = M3B[1] ; *mE (m , 1, 2) = *mE (m ,2,1) = un ; 
  m = M3B[2] ; *mE (m , 1, 2) = moinsI ; *mE (m ,2,1) = plusI ;
  m = M3B[3] ; *mE (m , 2, 2) = un ; *mE (m , 3, 3) = moinsUn ; 
  m = M3B[8] ; *mE (m , 0, 0) = moinsDeux ; *mE (m , 1, 1) = moinsUn ; *mE (m , 2, 2) = moinsUn ;

  m = M3B[4] ; *mE (m , 0, 1) = *mE (m ,1,0) = un ; 
  m = M3B[5] ; *mE (m , 0, 1) = moinsI ; *mE (m ,1,0) = plusI ;
  m = M3B[6] ; *mE (m , 0, 2) = un ; *mE (m , 2, 0) = un ; 
  m = M3B[7] ; *mE (m , 0, 2) = moinsI ; *mE (m , 2, 0) = plusI ;


  m = M4[1] ; *mE (m , 1, 2) = *mE (m ,2,1) = un ;
  m = M4[2] ; *mE (m , 1, 2) = moinsI ; *mE (m ,2,1) = plusI ;
  m = M4[3] ; *mE (m , 1, 1) = un ; *mE (m , 2, 2) = moinsUn ;
  m = M4[8] ; 
  *mE (m , 0, 0) = *cXY(-4.0/3.0,0,h) ; 
  *mE (m , 1, 1) = *cXY(-1.0/3.0,0,h) ; 
  *mE (m , 2, 2) = *cXY(-1.0/3.0, 0, h) ; *mE (m , 3, 3) = *cXY(2.0/3.0, 0, h) ;

  m = M4[4] ; 
  *mE (m , 0, 2) = *cXY(-s2/s3,0,h) ;  *mE (m , 2, 0) = *cXY(s2/s3,0,h) ; 
  *mE (m , 1, 3) = *cXY(1.0/s3,0,h) ;  *mE (m , 3, 1) = *cXY(1.0/s3,0,h) ;
  m = M4[5] ; 
  *mE (m , 0, 2) = *cXY(0,s2/s3,h) ;  *mE (m , 2, 0) =  *cXY(0,s2/s3,h) ; 
  *mE (m , 1, 3) =  *cXY(0,-1.0/s3,h) ;  *mE (m , 3, 1) =  *cXY(0,1.0/s3,h) ;

  m = M4[6] ; 
  *mE (m , 0, 1) = *cXY(s2/s3,0,h) ;  *mE (m , 1, 0) = *cXY(-s2/s3,0,h) ; 
  *mE (m , 2, 3) = *cXY(1.0/s3,0,h) ;  *mE (m , 3, 2) = *cXY(1.0/s3,0,h) ;
  m = M4[7] ; 
  *mE (m , 0, 1) = *cXY(0,-s2/s3,h) ;  *mE (m , 1, 0) = *cXY(0,-s2/s3,h) ; 
  *mE (m , 2, 3) = *cXY(0,-1.0/s3,h) ;  *mE (m , 3, 2) = *cXY(0,1.0/s3,h) ;

  for (ii = 1 ; ii <= 8 ; ii++)
    {
      m15 = M15[ii] ;
      m = M3[ii] ; 
      for (i = 0 ; i < 3 ; i++)
	for (j = 0 ; j < 3 ; j++)
	  *mE (m15, i, j) = *mE (m, i,j) ;

      m = M4[ii] ;
      for (jj = 3 ; jj < 15 ; jj += 4) 
	for (i = 0 ; i < 4 ; i++)
	  for (j = 0 ; j < 4 ; j++)
	    {
	      *mE (m15, jj + i, jj + j) = *mE (m, i,j) ;
	    }
    }
  return TRUE ;
} /* su21ConstructFamillyMatrix */

/* construct an indecomposable representation with 2 quark families */
static BOOL su21Construct2quarkMatrix (AC_HANDLE h)
{
  int ii, i, j ;
  MM m, mqq ;
  CXs un = {1,0} ;
  CXs moinsUn = {-1.0,0} ;
  /*
  CXs moinsI = {0,-1.0};
  CXs deux = {2.0,0} ;
  double s2 = sqrt(2.0), s3 = sqrt(3.0) ;
  */
  PLQQ = mNew (8, h) ;
  m = PLQQ ;
  *mE(m , 1, 1) = un ; *mE (m , 2, 2) = un ;
  *mE(m , 5, 5) = un ; *mE (m , 6, 6) = un ;

  PRQQ = mNew (8, h) ;
  m = PRQQ ;
  *mE(m , 0, 0) = un ; *mE (m , 3, 3) = un ;
  *mE(m , 4, 4) = un ; *mE (m , 7, 7) = un ;

  for (ii = 1 ; ii <= 8 ; ii++)
    {
      mqq = MEE[ii] = mNew (8,h) ;
      m = M3[ii] ;
      for (i = 0 ; i < 4 ; i++)
	for (j = 0 ; j < 4 ; j++)
	  {
	    *mE (mqq , i, j) = *mE (m , i, j)  ;
	    *mE (mqq , i + 4, j + 4) = *mE (m , i, j)  ;
	  }
      mqq = MQQ[ii] = mNew (8,h) ;
      m = M4[ii] ;
      for (i = 0 ; i < 4 ; i++)
	for (j = 0 ; j < 4 ; j++)
	  {
	    *mE (mqq , i, j) = *mE (m , i, j)  ;
	    *mE (mqq , i + 4, j + 4) = *mE (m , i, j)  ;
	  }
    }
  mqq = MEE[4] ;
  m = M3[5] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
      }
  mqq = MEE[5] ;
  m = M3[4] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
	cMM (mE (mqq , i + 4, j), &moinsUn) ;
      }
  mqq = MEE[6] ;
  m = M3[7] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
      }
  mqq = MEE[7] ;
  m = M3[6] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
	cMM (mE (mqq , i + 4, j), &moinsUn) ;
      }


  mqq = MQQ[4] ;
  m = M4[5] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
      }
  mqq = MQQ[5] ;
  m = M4[4] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
	cMM (mE (mqq , i + 4, j), &moinsUn) ;
      }
  mqq = MQQ[6] ;
  m = M4[7] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
      }
  mqq = MQQ[7] ;
  m = M4[6] ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	*mE (mqq , i + 4, j) = *mE (m , i, j)  ;
	cMM (mE (mqq , i + 4, j), &moinsUn) ;
      }


  return TRUE ;
} /* su21Construct2quarkMatrix */

/****************************************************************************/

static CX su21SuperTrace (MM mm, AC_HANDLE h)
{
  int ii, jj ; 
  double p ;
  CX z = cNew (h), za ;
  
  if (mm->dim == 15) jj = 0 ;
  else if (mm->dim == 3) jj = 0 ;
  else if (mm->dim == 4) jj = 3 ;
  
  for (ii = 0 ; ii < mm->dim ; ii++)
    {
      za = mE (mm, ii, ii) ;
      p = mE(PL15, jj+ii, jj+ii)->x ;
      if (p == 1)
	{ z->x += za->x ; z->y += za->y ;}
      else 
	{ z->x -= za->x ; z->y -= za->y ;}
    }
  z->x = cRound (z->x) ;
  z->y = cRound (z->y) ;
  return z ;
}

/****************************************************************************/

static CX su21LeftTrace (MM mm, AC_HANDLE h)
{
  int ii, jj ; 
  double p ;
  CX z = cNew (h), za ;
  
  if (mm->dim == 15) jj = 0 ;
  else if (mm->dim == 3) jj = 0 ;
  else if (mm->dim == 4) jj = 3 ;
  
  for (ii = 0 ; ii < mm->dim ; ii++)
    {
      za = mE (mm, ii, ii) ;
      p = mE(PL15, jj+ii, jj+ii)->x ;
      if (p == 1)
	{ z->x += za->x ; z->y += za->y ;}
    }
  z->x = cRound (z->x) ;
  z->y = cRound (z->y) ;
  return z ;
}

/****************************************************************************/

static CX su21RightTrace (MM mm, AC_HANDLE h)
{
  int ii, jj ; 
  double p ;
  CX z = cNew (h), za ;
  
  if (mm->dim == 15) jj = 0 ;
  else if (mm->dim == 3) jj = 0 ;
  else if (mm->dim == 4) jj = 3 ;
  
  for (ii = 0 ; ii < mm->dim ; ii++)
    {
      za = mE (mm, ii, ii) ;
      p = mE(PR15, jj+ii, jj+ii)->x ;
      if (p == 1)
	{ z->x += za->x ; z->y += za->y ;}
    }
  z->x = cRound (z->x) ;
  z->y = cRound (z->y) ;
  return z ;
}

/****************************************************************************/
/****************************************************************************/
/* this code computes the numbe of configs of 4 numbers a,b,c,d
 * such that a+b+c+d = n
 */ 
static void titrationConfig (void)
{
  int a, b, c, d, n1, n2, n3, n4, nn1, nn2, nn3, nn4,  nn ;

  for (nn = 1 ; nn <= 100 ; nn++)
    {
      nn1 = nn2 = nn3 = nn4 = 0 ;
      for (a = 1 ; a <= nn ; a++)
	for (b = 0 ; b <= a && a + b <= nn ; b++)
	  for (c = 0 ; c <= b  && a + b + c <= nn ; c++)
	    for (d = 0 ; d <= c  && a + b + c + d <= nn ; d++)
	      {
		if (a + b + c + d != nn) continue ;
		if (a == b && a == c && a == d)
		  { n1 = 1 ; n2 = 0 ; n3 = 0 ; n4 = 1 ; }
		else if (
			 (a == b && a == c) ||
			 (b == c && b == d) 
			 )
		  { n1 = 2 ; n2 = 2 ; n3 = 0 ; n4 = 4 ; }
		else if (
			 (a == b && c == d)
			 )
		  { n1 = 2 ; n2 = 4 ; n3 = 0 ; n4 = 6 ;}
		else if ( a== b || b == c || c == d)
		  { n1 = 2 ; n2 = 10 ; n3 = 2 ; n4 = 10 ; }
		else
		  { n1 = 2 ; n2 = 22 ; n3 = 2 ; n4 = 22 ; }
		if (nn <= 10) printf ("N=%d\t(%d,%d,%d,%d)\t %d;%d\t %d;%d\n", nn, a,b,c,d,n1, n2, n3, n4) ;
		nn1 += n1 ; nn2 += n2 ;	nn3 += n3 ; nn4 += n4 ;
	      }
      printf ("N=%d any\t%d;%d\tFDR\t%.5f\t\t%d;%d\tFDR\t%.5f\n", nn, nn1, nn2,nn1/((double)nn1+nn2), nn3, nn4,nn3/((double)nn3+nn4)) ;
    }
  exit (0) ;
} /* titrationConfig */

/****************************************************************************/
/*************************  real work ***************************************/

int main (int argc, char **argv)
{
  int ii, jj, kk, ll ;
  CX z ;
  MM p ;
  AC_HANDLE h = 0 ;

  freeinit () ;
  h = ac_new_handle () ;


  if (0) titrationConfig() ;


  su21ConstructFamillyMatrix (h) ; 
  if(1)su21Construct2quarkMatrix (h) ; 
  for (ii = 1 ; ii <= 8 ; ii++)
    {
      if (1)
	{
	  if (0) mShow (M3[ii], hprintf (h, "Lepton %d", ii)) ;
	  if (0) mShow (M4[ii], hprintf (h, "Quark %d", ii)) ;
	  if (1) mShow (MEE[ii], hprintf (h, "TwoLeptons %d", ii)) ;
	  if (0) mShow (MQQ[ii], hprintf (h, "TwoQuark %d", ii)) ;
	}
      if (0) mShow (M15[ii], hprintf (h, "Familly %d", ii)) ;
    }

   /* construct the faij coefficient */
   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       faij = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (0 && aa > 3 && aa < 8)
	     continue ;
	   for (ii = 1 ; ii < 9 ; ii++)
	     for (jj = 1 ; jj < 9 ; jj++)
	       {
		 p1 =  mCopy (M3[aa], h) ;
		 mMM (p1, M3[ii]) ; 
		 p2 =  mCopy (M3[ii], h) ;
		 mMM (p2, M3[aa]) ;
		 mScale (-1, 0, p2) ;
		 mAA (p2, p1) ;
		 mMM (p2, M3[jj]) ;
		 if(aa==8) mScale (0,1/2.0,p2) ;
		 else mScale (0,-1/2.0,p2) ;
		 z = mTrace (p2, h) ;
		 *(t3E(faij,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("f[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (0) 
	 {
	   printf ("f[aij] leptons\n") ;
	   t3Show (faij) ;
	 }
     }

   /* construct the faij coefficient */
   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 faijEE = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (0 && aa > 3 && aa < 8)
	     continue ;
	   for (ii = 1 ; ii < 9 ; ii++)
	     for (jj = 1 ; jj < 9 ; jj++)
	       {
		 p1 =  mCopy (MEE[aa], h) ;
		 mMM (p1, MEE[ii]) ; 
		 p2 =  mCopy (MEE[ii], h) ;
		 mMM (p2, MEE[aa]) ;
		 mScale (-1, 0, p2) ;
		 mAA (p2, p1) ;
		 mMM (p2, MEE[jj]) ;
		 if(aa==8) mScale (0,1/2.0,p2) ;
		 else mScale (0,-1/2.0,p2) ;
		 z = mTrace (p2, h) ;
		 *(t3E(faijEE,aa,ii,jj)) = *z ;
	       }
	 }
       if (0) 
	 {
	   printf ("f[aij]EE double-leptons\n") ;
	   t3Show (faijEE) ;
	 }
       t3Scale (-1.0/2.0,0,faijEE) ;
       t3AA (faijEE, faij) ;
       if (1) 
	 {
	   printf ("f-leptons - f-double-lepton, should be empty\n") ;
	   t3Show (faijEE) ;
	 }

     }

   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 faijq = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   for (ii = 1 ; ii < 9 ; ii++)
	     for (jj = 1 ; jj < 9 ; jj++)
	       {
		 p1 =  mCopy (M4[aa], h) ;
		 mMM (p1, M4[ii]) ; 
		 p2 =  mCopy (M4[ii], h) ;
		 mMM (p2, M4[aa]) ;
		 mScale (-1, 0, p2) ;
		 mAA (p2, p1) ;
		 mMM (p2, M4[jj]) ;
		 if(aa==8) mScale (0,1/2.0,p2) ;
		 else mScale (0,-1/2.0,p2) ;
		 z = mTrace (p2, h) ;
		 *(t3E(faijq,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("f[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (0) 
	 {
	   printf ("f[aij] quarks\n") ;
	   t3Show (faijq) ;
	 }
       t3Scale (3.0,0,faijq) ;
       t3AA (faijq, faij) ;
       if (1) 
	 {
	   printf ("f-leptons + 3-f-quarks, should be empty\n") ;
	   t3Show (faijq) ;
	 }
     }

   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 faijf = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (aa > 3 && aa < 8)
	     continue ;
	   for (ii = 4 ; ii < 8 ; ii++)
	     for (jj = 4 ; jj < 8 ; jj++)
	       {
		 p1 =  mCopy (M15[aa], h) ;
		 mMM (p1, M15[ii]) ; 
		 p2 =  mCopy (M15[ii], h) ;
		 mMM (p2, M15[aa]) ;
		 mScale (-1, 0, p2) ;
		 mAA (p2, p1) ;
		 mMM (p2, M15[jj]) ;
		 z = mTrace (p2, h) ;
		 { z->x /= 2 ; z->y /= 2 ; }
		 *(t3E(faijf,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("f[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (1) 
	 {
	   printf ("f[aij] famille\n") ;
	   t3Show (faijf) ;
	 }
       t3Scale (1.0,0,faijf) ;
       if (1) 
	 {
	   printf ("f-famille, should be empty\n") ;
	   t3Show (faijf) ;
	 }
     }

   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 faijf = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (aa > 3 && aa < 8)
	     continue ;
	   for (ii = 4 ; ii < 8 ; ii++)
	     for (jj = 4 ; jj < 8 ; jj++)
	       {
		 p1 =  mCopy (MQQ[aa], h) ;
		 mMM (p1, MQQ[ii]) ; 
		 p2 =  mCopy (MQQ[ii], h) ;
		 mMM (p2, MQQ[aa]) ;
		 mScale (-1, 0, p2) ;
		 mAA (p2, p1) ;
		 mMM (p2, MQQ[jj]) ;
		 z = mTrace (p2, h) ;
		 { z->x /= 2 ; z->y /= 2 ; }
		 *(t3E(faijf,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("f[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (1) 
	 {
	   printf ("f[aij] 2quarks\n") ;
	   t3Show (faijf) ;
	 }
       t3Scale (1.0,0,faijf) ;
       if (1) 
	 {
	   printf ("f-2quarks, should be empty\n") ;
	   t3Show (faijf) ;
	 }
     }

   /* construct the daij coefficient */
   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       daij = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (aa > 3 && aa < 8)
	     continue ;
	   for (ii = 4 ; ii < 8 ; ii++)
	     for (jj = 4 ; jj < 8 ; jj++)
	       {
		 p1 =  mCopy (M3[ii], h) ;
		 mMM (p1, M3[jj]) ; 
		 p2 =  mCopy (M3[jj], h) ;
		 mMM (p2, M3[ii]) ; 
		 mAA (p2, p1) ;
		 mMM (p2, M3[aa]) ;
		 z = su21SuperTrace (p2, h) ;
		 if (aa==8) { z->x /= -2 ; z->y /= -2 ; }
		 else  { z->x /= 2 ; z->y /= 2 ; }
		 *(t3E(daij,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("d[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (1) 
	 {
	   printf ("d[aij] leptons\n") ;
	   t3Show (daij) ;
	 }
     }
   
   return 0 ;

   /* construct the daij coefficient */
   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 daijq = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (aa > 3 && aa < 8)
	     continue ;
	   for (ii = 4 ; ii < 8 ; ii++)
	     for (jj = 4 ; jj < 8 ; jj++)
	       {
		 p1 =  mCopy (M4[ii], h) ;
		 mMM (p1, M4[jj]) ; 
		 p2 =  mCopy (M4[jj], h) ;
		 mMM (p2, M4[ii]) ; 
		 mAA (p2, p1) ;
		 mMM (p2, M4[aa]) ;
		 z = su21SuperTrace (p2, h) ;
		 if (aa==8) { z->x /= -2 ; z->y /= -2 ; }
		 else  { z->x /= 2 ; z->y /= 2 ; }
		 *(t3E(daijq,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("d[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (0) 
	 {
	   printf ("d[aij] quark\n") ;
	   t3Show (daij) ;
	 }
       t3Scale (-1.0,0,daijq) ;
       t3AA (daijq, daij) ;
       if (1) 
	 {
	   printf ("d-leptons - d-quarks, should be empty\n") ;
	   t3Show (daijq) ;
	 }
     }
   /* construct the daij coefficient */
   if (1)
     {
       CX z ;
       int aa, ii, jj ;
       MM p1, p2 ;
       
       TT3 daijf = t3New (9,"aij", h) ;
       for (aa = 1 ; aa < 9 ; aa++)
	 {
	   if (aa > 3 && aa < 8)
	     continue ;
	   for (ii = 4 ; ii < 8 ; ii++)
	     for (jj = 4 ; jj < 8 ; jj++)
	       {
		 p1 =  mCopy (M15[ii], h) ;
		 mMM (p1, M15[jj]) ; 
		 p2 =  mCopy (M15[jj], h) ;
		 mMM (p2, M15[ii]) ; 
		 mAA (p2, p1) ;
		 mMM (p2, M15[aa]) ;
		 z = su21SuperTrace (p2, h) ;
		 if (aa==8) { z->x /= -2 ; z->y /= -2 ; }
		 else  { z->x /= 2 ; z->y /= 2 ; }
		 *(t3E(daijf,aa,ii,jj)) = *z ;
		 if (0 && aa==8)
		   printf ("d[8;%d,%d]=%s\n", ii, jj, cShow(z)) ;
	       }
	 }
       if (0) 
	 {
	   printf ("d[aij] famille\n") ;
	   t3Show (daijf) ;
	 }
       t3Scale (-.25,0,daijf) ;
       t3AA (daijf, daij) ;
       if (1) 
	 {
	   printf ("d-leptons - 1/4 d-famille, should be empty\n") ;
	   t3Show (daijf) ;
	 }
     }

   /***************************************************/
   /* JACOBI */
   if (1)
     {
       CXs cc0, cc ;
       CX z0 = &cc0, z=&cc, z1, z2 ;
       int a,b,c, a1,b1,c1, cycle;
       int d; /*free index */
       int e ; /* dummy index */
       int bose[10] = {0,1,1,1,0,0,0,0,1,0} ;

       for (d = 1 ; d < 4 ; d++)
	 {
	   if (!bose[d]) continue ;

	   for (a = 1 ; a < 4 ; a++)
	     for (b = 1 ; b < 4 ; b++)
	       for (c = 1 ; c < 4 ; c++)
		 {
		   if (!bose[a]) continue ;
		   if (!bose[b]) continue ;
		   if (!bose[c]) continue ;
		   z->x = z->y = 0 ;
		   for (cycle = 0 ; cycle < 3 ; cycle++)
		     {
		       switch (cycle)
			 {
			 case 0:
			   a1 = a ; b1 = b ; c1 = c ;
			   break ;
			 case 1:
			   a1 = b ; b1 = c ; c1 = a ;
			   break ;
			 case 2:
			   a1 = c ; b1 = a ; c1 = b ;
			   break ;
			 }		
		       for (e= 1 ; e < 9 ; e++)	 
			 {
			   z1 = t3E(faij,d,e,a1) ;
			   *z0 = *z1 ;
			   z2 = t3E(faij,e,b1,c1) ;
			   cMM (z0, z2) ;
			   cAA (z, z0) ;
			 }
		     }
		   if (0) printf ("%d; %d,%d,%d  z=%s\n",d,a1,b1,c1,cShow(z)) ;
		   if (cNorm(z) > .01)
		     messcrash ("Jacobi(%d ; (%d,%d,%d)) = %s not null\n", 
				d,a,b,c,cShow(z)) ;
		 }
	 }
     }
   printf ("Jacobi fdea*febc ok\n") ;
   /***************************************************/
   /* JACOBI */
   if (1)
     {
       CXs cc0, cc ;
       CX z0 = &cc0, z=&cc, z1, z2 ;
       int a,b,c, a1,b1,c1, cycle;
       int d; /*free index */
       int e ; /* dummy index */
       TT3 jabc = 0 ;
       int bose[10] = {0,1,1,1,0,0,0,0,1,0} ;
       for (d = 1 ; d < 9 ; d++)
	 {
	   if (bose[d]) continue ;
	   ac_free (jabc) ;
	   jabc = t3New (9,"abc", h) ;
	   for (a = 1 ; a < 9 ; a++)
	     for (b = 1 ; b < 9 ; b++)
	       for (c = 1 ; c < 9 ; c++)
		 {
		   if (bose[a]) continue ;
		   if (bose[b]) continue ;
		   if (bose[c]) continue ;
		   z->x = z->y = 0 ;
		   for (cycle = 0 ; cycle < 3 ; cycle++)
		     {
		       switch (cycle)
			 {
			 case 0:
			   a1 = a ; b1 = b ; c1 = c ;
			   break ;
			 case 1:
			   a1 = b ; b1 = c ; c1 = a ;
			   break ;
			 case 2:
			   a1 = c ; b1 = a ; c1 = b ;
			   break ;
			 }		
		       for (e= 1 ; e < 9 ; e++)	 
			 {
			   z1 = t3E(faij,d,e,a1) ;
			   *z0 = *z1 ;
			   z2 = t3E(daij,e,b1,c1) ;
			   cMM (z0, z2) ;
			   cAA (z, z0) ;
			 }
		       if (0) printf ("%d; %d,%d,%d  z=%s\n",d,a1,b1,c1,cShow(z)) ;
		     }
		   if (cNorm(z) > .01)
		     messcrash ("Jacobi(%d ; (%d,%d,%d)) = %s not null\n", 
				d,a,b,c,cShow(z)) ;
		 }
	 }
     }
   printf ("Jacobi fdei*dejk ok\n") ;

   /***************************************************/
   /* JACOBI */
   if (0)
     {
       CXs cc0, cc ;
       CX z0 = &cc0, z=&cc, z1, z2, zi ;
       int a,b,c, a1,b1,c1, cycle;
       int d; /*free index */
       int e ; /* dummy index */
       TT3 jabc = 0 ;
       int bose[10] = {0,1,1,1,0,0,0,0,1,0} ;
       zi = cXY (0,1,h) ;

       for (d = 1 ; d < 9 ; d++)
	 {
	   if (bose[d]) continue ;
	   ac_free (jabc) ;
	   jabc = t3New (9,"abc", h) ;
	   for (a = 1 ; a < 9 ; a++)
	     for (b = 1 ; b < 9 ; b++)
	       for (c = 1 ; c < 9 ; c++)
		 {
		   if (bose[a]) continue ;
		   if (bose[b]) continue ;
		   if (bose[c]) continue ;
		   z->x = z->y = 0 ;
		   for (cycle = 0 ; cycle < 2 ; cycle++)
		     {
		       switch (cycle)
			 {
			 case 0:
			   a1 = a ; b1 = b ; c1 = c ;
			   break ;
			 case 1:
			   a1 = b ; b1 = a ; c1 = c ;
			   break ;
			 }		
		       for (e= 1 ; e < 9 ; e++)	 
			 {
			   z1 = t3E(faij,e,b1,c1) ;
			   *z0 = *z1 ;
			   cMM (z0,zi) ;
			   z1 = t3E(daij,e,b1,c1) ;
			   cAA (z0,z1) ;
			   z2 = t3E(faij,d,e,a1) ;
			   cMM (z0, z2) ;
			   cAA (z, z0) ;
			 }
		     }
		   if (cNorm(z) > .01)
		     {
		       invokeDebugger() ;
		       messcrash ("Jacobi(%d ; (%d,%d,%d)) = %s not null\n", 
				  d,a,b,c,cShow(z)) ;
		     }
		 }
	 }
     }
   printf ("Jacobi fdei*(dejk + i fejk) ok\n") ;

   /***************************************************/
   /* propagateur */
  printf ("\nPropagateur\n") ;
  for (ii = 4 ; ii <= 7 ; ii++)
    {
      for (jj = 4 ; jj <= 7 ; jj++)
	{
	  p = mCopy (M15[ii], h) ;
	  mMM (p, M15[jj]) ;
	  z = su21LeftTrace (p, h) ;
	  z = su21RightTrace (p, h) ;
	  printf ("\t%d:%d:%s", ii, jj, cShow (z)) ;
	}
      printf ("\n") ;
    }
  if (0) printf ("m15(3,3) = %s\n", cShow (mE(M15[8],3,3))) ;

  if (0) /* daij */
    {
      p = mCopy (M15[4], h) ;
      mMM (p, M15[4]) ;
      mScale (-2,0,p) ;
      mAA (p, M15[3]) ;
      mAA (p, M15[8]) ;
      mShow (p, "M4 M4 - ( M3 + M8) ") ;
    }

   /* vertex 4 scalaires */
   printf ("\nVertex in LRLR, normalization H mu psiL-psiR\n") ;
   printf ("\nscale=3/16=1/4(stats)2(loops)2(trace epsison)1(use H), then divide by 16/3\n") ;
  phi4 = t4New (8,"iJkL",h) ;
   {
     for (ii = 4 ; ii <= 7 ; ii++)
       for (jj = 4 ; jj <= 7 ; jj++)
	 for (kk = 4 ; kk <= 7 ; kk++) 
	   for (ll = 4 ; ll <= 7 ; ll++)
	     {
	       p = mCopy (M15[ii], h) ;
	       mMM (p, M15[jj]) ;
	       mMM (p, M15[kk]) ;
	       mMM (p, M15[ll]) ;
	       z = su21RightTrace (p, h) ;
	       *t4E(phi4,ii,jj,kk,ll) = *z ;
	     }
     t4HalfSym (phi4) ;
     t4Scale (3.0/16.0, 0, phi4) ;
     printf ("phi4 computed as Tr() normalized as 1, divided by 16/3\n") ;
     if (0) 
       {
	 TT4 tt ;

	 tt = t4Copy (phi4, "ijkl", h) ;
	 t4Scale (1, 0, tt) ;
	 t4Show (tt) ;
       }
   }


   /* construction directe du terme LRLR = -1/8 (faij LR) (1,1,1/3) (fbkl LR) */
   if (1)
     {
       int ii, jj, kk, ll ;
       CXs cc ;
       CX z, z0 = &cc, z1, z2, z3, z8 ;

       phi4a = t4New (9, "LRLR", h) ;
       phi4f = t4New (9, "LRLR", h) ;
       for (ii = 4 ; ii <= 7 ; ii++)
	 for (jj = 4 ; jj <= 7 ; jj++)
	   for (kk = 4 ; kk <= 7 ; kk++) 
	     for (ll = 4 ; ll <= 7 ; ll++)
	       {
		 z1 = cCopy(t3E (faij, 1, ii, jj), h) ;
		 cMM (z1, t3E (faij, 1, kk, ll)) ;
		 z2 = cCopy(t3E (faij, 2, ii, jj), h) ;
		 cMM (z2, t3E (faij, 2, kk, ll)) ;
		 z3 = cCopy(t3E (faij, 3, ii, jj), h) ;
		 cMM (z3, t3E (faij, 3, kk, ll)) ;
		 z8 = cCopy(t3E (faij, 8, ii, jj), h) ;
		 cMM (z8, t3E (faij, 8, kk, ll)) ;

		 /* f (super-g) f */
		 z8->x *= -1 ; z8->y *= -1 ;
		 z0->x = z0->y = 0 ;
		 if(1)
		   {
		     cAA (z0, z1) ;
		     cAA (z0, z2) ;
		     cAA (z0, z3) ;
		   }
		 cAA (z0, z8) ;

		 z = t4E(phi4f,ii,jj,kk,ll) ;
		 cAA (z, z0) ;

		 /* f (su(3) g) f */ 
		 z8->x *= -3 ; z8->y *= -3 ; /* restore the sign and scale */
		 z0->x = z0->y = 0 ;
		 if(1)
		   {
		     cAA (z0, z1) ;
		     cAA (z0, z2) ;
		     cAA (z0, z3) ;
		   }
		 cAA (z0, z8) ;
		 z0->x *= -1 ; z0->y *= -1 ;
		 z = t4E(phi4a,ii,jj,kk,ll) ;
		 cAA (z, z0) ;

	       }
       t4Scale (1.0/8.0, 0, phi4a) ;
       t4HalfSym (phi4a) ;
       t4HalfSym (phi4f) ;
     }
   if (1)
     {
       TT4 tt ;

       printf ("phi4a computed as  -1/8(faij LR) (1,1,1/3) (fbkl LR)\n") ;
       if (0)
	 {
	   tt = t4Copy (phi4a, "ijkl", h) ;
	   t4Scale (1.0, 0, tt) ;
	   t4Show (tt) ;
	 }
     }

   /* construction directe du terme LRLR = 1/8 (daij LR) (3,3,3/1) (dbkl LR) */
   if (1)
     {
       int ii, jj, kk, ll ;
       CXs cc ;
       CX z, z0 = &cc, z1, z2, z3, z8 ;

       phi4b = t4New (9, "LRLR", h) ;
       phi4d = t4New (9, "LRLR", h) ;
       for (ii = 4 ; ii <= 7 ; ii++)
	 for (jj = 4 ; jj <= 7 ; jj++)
	   for (kk = 4 ; kk <= 7 ; kk++) 
	     for (ll = 4 ; ll <= 7 ; ll++)
	       {
		 z1 = cCopy(t3E (daij, 1, ii, jj), h) ;
		 cMM (z1, t3E (daij, 1, kk, ll)) ;
		 z2 = cCopy(t3E (daij, 2, ii, jj), h) ;
		 cMM (z2, t3E (daij, 2, kk, ll)) ;
		 z3 = cCopy(t3E (daij, 3, ii, jj), h) ;
		 cMM (z3, t3E (daij, 3, kk, ll)) ;
		 z8 = cCopy(t3E (daij, 8, ii, jj), h) ;
		 cMM (z8, t3E (daij, 8, kk, ll)) ;

		 /* d super-g d */
		 z8->x *= -1 ; z8->y *= -1 ;
		 z0->x = z0->y = 0 ;
		 if (1)
		   {
		     cAA (z0, z1) ;
		     cAA (z0, z2) ;
		     cAA (z0, z3) ;
		   }
		 z0->x *= 1 ; z0->y *= 1 ;
		 cAA (z0, z8) ;

		 z = t4E(phi4d,ii,jj,kk,ll) ;
		 cAA (z, z0) ;

		 /* d (3 3 3 / 1) d */
		 z8->x *= -1 ; z8->y *= -1 ; /* restore the sign */
		 z0->x = z0->y = 0 ;
		 if (1)
		   {
		     cAA (z0, z1) ;
		     cAA (z0, z2) ;
		     cAA (z0, z3) ;
		   }
		 z0->x *= 3 ; z0->y *= 3 ; /* scale */
		 cAA (z0, z8) ;

		 z = t4E(phi4b,ii,jj,kk,ll) ;
		 cAA (z, z0) ;
	       }
       t4Scale (1.0/8.0, 0, phi4b) ;
       t4HalfSym (phi4b) ;
       t4HalfSym (phi4d) ;
     }
   if (1)
     {
       TT4 tt ;

       printf ("phi4b computed as  1/8(daij LR) (3,3,3/1) (dbkl LR)\n") ;
       if (0)
	 {
	   tt = t4Copy (phi4b, "ijkl", h) ;
	   t4Scale (1.0, 0, tt) ;
	   t4Show (tt) ;
	 }
     }

   /* construction directe du terme LRLR = (faij LR) (1,1,1/-1) (dbkl LR) */
   if (1)
     {
       int ii, jj, kk, ll ;
       CXs cc ;
       CX z, z0 = &cc, z1, z2, z3, z8 ;

       phi4fd = t4New (9, "LRLR", h) ;
       for (ii = 4 ; ii <= 7 ; ii++)
	 for (jj = 4 ; jj <= 7 ; jj++)
	   for (kk = 4 ; kk <= 7 ; kk++) 
	     for (ll = 4 ; ll <= 7 ; ll++)
	       {
		 z1 = cCopy(t3E (faij, 1, ii, jj), h) ;
		 cMM (z1, t3E (daij, 1, kk, ll)) ;
		 z2 = cCopy(t3E (faij, 2, ii, jj), h) ;
		 cMM (z2, t3E (daij, 2, kk, ll)) ;
		 z3 = cCopy(t3E (faij, 3, ii, jj), h) ;
		 cMM (z3, t3E (daij, 3, kk, ll)) ;
		 z8 = cCopy(t3E (faij, 8, ii, jj), h) ;
		 cMM (z8, t3E (daij, 8, kk, ll)) ;

		 /* d (1,1,1 / 1) d */
		 z8->x *= -1 ; z8->y *= -1 ; /* scale */
		 z0->x = z0->y = 0 ;
		 if (1)
		   {
		     cAA (z0, z1) ;
		     cAA (z0, z2) ;
		     cAA (z0, z3) ;
		   }
		 z0->x *= 1 ; z0->y *= 1 ; /* scale */
		 cAA (z0, z8) ;

		 z = t4E(phi4fd,ii,jj,kk,ll) ;
		 cAA (z, z0) ;
	       }
       t4Scale (1.0/1.0, 0, phi4fd) ;
       t4HalfSym (phi4fd) ;
     }
   if (1)
     {
       TT4 tt ;

       printf ("(faij LR) (1,1,1/-1) (dbkl LR) should be empty\n") ;
       if (1)
	 {
	   tt = t4Copy (phi4fd, "ijkl", h) ;
	   t4Scale (1.0, 0, tt) ;
	   t4Show (tt) ;
	 }
     }

   if (1)
     {
       TT4 tt ;

       if (1)
	 {
	   printf ("verif phi4 - (phi4a + phi4b))  (should be empty) \n") ;
	   tt = t4Copy (phi4b, "ijkl", h) ;
	   t4AA (tt, phi4a) ;
	   t4Scale (-1, 0, tt) ;
	   t4AA (tt, phi4) ;
	   t4Show (tt) ;
	 }
       if (1)
	 {
	   printf ("verif -f g f + d g d == 0  (should be empty) \n") ;
	   tt = t4Copy (phi4f, "ijkl", h) ;
	   t4Scale (-1, 0, tt) ;
	   t4AA (tt, phi4d) ;
	   t4Show (tt) ;
	 }
     }
   printf ("Using:\n") ;
   printf ("f3_45=%s\t", cShow (t3E(faij,3,4,5))) ;
   printf ("d3_44=%s\n", cShow (t3E(daij,3,4,4))) ;
   printf ("f8_45=%s\t", cShow (t3E(faij,8,4,5))) ;
   printf ("d8_44=%s\n", cShow (t3E(daij,8,4,4))) ;
   printf ("f3_67=%s\t", cShow (t3E(faij,3,6,7))) ;
   printf ("d3_66=%s\n", cShow (t3E(daij,3,6,6))) ;
   printf ("f8_67=%s\t", cShow (t3E(faij,8,6,7))) ;
   printf ("d8_66=%s\n", cShow (t3E(daij,8,6,6))) ;

   if (0)
   {
     /* faux parce que le propagateur du phi doit etre ajoute */
     TT4 zphi4 = t4t4doubleContract(phi4,phi4,"iJkL",3,4,1,2,h);
     if (1) t4Show (zphi4) ;
   }

   /* matrices left et right */
   for (jj = 4 ; jj <= 7 ; jj++)
     {
       M15L[jj] = mCopy (M15[jj], h) ;
       M15R[jj] = mCopy (M15[jj], h) ;
       mMM (M15L[jj], PR15) ; /* reversed, becuase i multiply on the right */
       mMM (M15R[jj], PL15) ;
     }

   /* matrices A et B */
   for (jj = 4 ; jj <= 7 ; jj++)
     { 
       MA[jj] = mCopy (M15L[jj], h) ;
       mAA (MA[jj], M15R[jj]) ;
       mScale (.5,0,MA[jj]) ;
       if (0)
	 mShow (MA[jj], messprintf ("MA[%d]", jj)) ;
     }

   for (jj = 4 ; jj <= 7 ; jj += 2)
     { 
       MB[jj] = mCopy (M15R[jj+1], h) ;
       mScale (-1,0,MB[jj]) ;
       mAA (MB[jj], M15L[jj+1]) ;
       mScale (0,.5,MB[jj]) ;
       
       MB[jj+1] = mCopy (M15R[jj], h) ;
       mScale (-1,0,MB[jj+1]) ;
       mAA (MB[jj+1], M15L[jj]) ;
       mScale (0,-.5,MB[jj+1]) ;
       
       if (0)
	 {
	   mShow (MB[jj], messprintf ("MB[%d]", jj)) ;
	   mShow (MB[jj], messprintf ("MB[%d]", jj+1)) ;
	 }
     }
   
   
   /* matrices H et K */ 
   for (jj = 4 ; jj <= 7 ; jj++)
     {
       MH[jj] = mCopy (MA[jj], h) ;
       mAA (MH[jj], MB[jj]) ;
       MK[jj] = mCopy (MA[jj], h) ;
       mScale (-1,0,MK[jj]) ;
       mAA (MK[jj], MB[jj]) ;
       mScale (0,1,MK[jj]) ;

       if (1)
	 {
	   mShow (MH[jj], messprintf ("MH[%d]", jj)) ;
	   mShow (MK[jj], messprintf ("MK[%d]", jj)) ;
	 }
     }
   
   if (0)
     {
       jj = 6 ;
       mShow (M15[jj], messprintf ("M15[%d]", jj)) ;
       mShow (M15L[jj], messprintf ("M15L[%d]", jj)) ;
       mShow (M15R[jj], messprintf ("M15R[%d]", jj)) ;
       mShow (MA[jj], messprintf ("MA[%d]", jj)) ;
       mShow (MB[jj], messprintf ("MB[%d]", jj)) ;
       mShow (MH[jj], messprintf ("MH[%d]", jj)) ;
       mShow (MK[jj], messprintf ("MK[%d]", jj)) ;
       jj = 7 ;
       mShow (MH[jj], messprintf ("MH[%d]", jj)) ;
       mShow (MK[jj], messprintf ("MK[%d]", jj)) ;
     }

  /* vertex */
   printf ("\nVertex in 4H, normalization H mu psiL-psiR\n") ;
   printf ("\nscale=3/32=1/4(stats)2(trace epsison)1(use H), then divide by 16/3\n") ;
   {
     TT4 tt = t4New (8,"iJkL",h) ;
     for (ii = 4 ; ii <= 7 ; ii++)
       for (jj = 4 ; jj <= 7 ; jj++)
	 for (kk = 4 ; kk <= 7 ; kk++) 
	   for (ll = 4 ; ll <= 7 ; ll++)
	     {
	       p = mCopy (MH[ii], h) ;
	       mMM (p, MH[jj]) ;
	       mMM (p, MH[kk]) ;
	       mMM (p, MH[ll]) ;
	       z = mTrace (p, h) ;
	       cAA (t4E(tt,ii,jj,kk,ll), z) ;
	     }
     t4Scale (3.0/32.0, 0, tt) ;
     t4FullSym (tt) ;
     H4 = t4Copy (tt, "HHHH", h) ;
     if (1)
       {
	 printf ("\n16*Vertex HHHH \n") ;
	 tt = t4Copy (H4, "HHHH", h) ;
	 t4Scale (16.0, 0, tt) ;
	 t4Show (tt) ;
       }
   }


   printf ("\nertex in 4K, normalization H mu psiL-psiR\n") ;
   printf ("\nscale=3/32=1/4(stats)2(trace epsison)1(use H), then divide by 16/3\n") ;
   {
     TT4 tt = t4New (8,"iJkL",h) ;
     for (ii = 4 ; ii <= 7 ; ii++)
       for (jj = 4 ; jj <= 7 ; jj++)
	 for (kk = 4 ; kk <= 7 ; kk++) 
	   for (ll = 4 ; ll <= 7 ; ll++)
	     {
	       p = mCopy (MK[ii], h) ;
	       mMM (p, MK[jj]) ;
	       mMM (p, MK[kk]) ;
	       mMM (p, MK[ll]) ;
	       z = mTrace (p, h) ;
	       cAA (t4E(tt,ii,jj,kk,ll), z) ;
	     }
     t4FullSym (tt) ;
     t4Scale (3.0/32.0, 0, tt) ;
     K4 = t4Copy (tt, "KKKK", h) ;
     if (0)
       {
	 printf ("\n16*Vertex KKKK \n") ;
	 tt = t4Copy (K4, "KKKK", h) ;
	 t4Scale (16.0, 0, tt) ;
	 t4Show (tt) ;
       }
   }


   printf ("\n Vertex mixte HHKK, normalization H mu psiL-psiR\n") ;
   printf ("\nscale=2*3/16=1(stats)2(trace epsison)16(use 2H), then divide by 16/3\n") ;
   {
     TT4 tt = t4New (8,"iJkL",h) ;
     for (ii = 4 ; ii <= 7 ; ii++)
       for (jj = 4 ; jj <= 7 ; jj++)
	 for (kk = 4 ; kk <= 7 ; kk++) 
	   for (ll = 4 ; ll <= 7 ; ll++)
	     {
	       p = mCopy (MH[ii], h) ;
	       mMM (p, MH[jj]) ;
	       mMM (p, MK[kk]) ;
	       mMM (p, MK[ll]) ;
	       z = mTrace (p, h) ;
	       cAA (t4E(tt,ii,jj,kk,ll), z) ;
	     }
     t4DoubleSym (tt) ;
     t4Scale (2.0*3.0/16.0, 0, tt) ;
     HHKK = t4Copy (tt,"HHKK", h) ;
     HKKH = t4RollLeft (HHKK, "HKKH", h) ;
     KKHH = t4RollLeft (HKKH, "KKHH", h) ;
     if (1)
       {
	 printf ("\n16*Vertex mixte HHKK \n") ;
	 tt = t4Copy (HHKK, "HHKK", h) ;
	 t4Scale (16.0, 0, tt) ;
	 t4Show (tt) ;	
	 if (0) t4Show (HKKH) ;
       }
   }
   printf ("\nVertex mixte HKHK units, (should be empty)\n") ;
   {
     TT4 tt = t4New (8,"iJkL",h) ;
     for (ii = 4 ; ii <= 7 ; ii++)
       for (jj = 4 ; jj <= 7 ; jj++)
	 for (kk = 4 ; kk <= 7 ; kk++) 
	   for (ll = 4 ; ll <= 7 ; ll++)
	     {
	       p = mCopy (MH[ii], h) ;
	       mMM (p, MK[jj]) ;
	       mMM (p, MH[kk]) ;
	       mMM (p, MK[ll]) ;
	       z = mTrace (p, h) ;
	       cAA (t4E(tt,ii,jj,kk,ll), z) ;
	       cAA(t4E(tt,ii,ll,kk,jj), z) ;

	       cAA(t4E(tt,kk,jj,ii,ll), z) ;
	       cAA(t4E(tt,kk,ll,ii,jj), z) ;
	     }
     HKHK = tt ;
     t4Scale (4.0, 0, HKHK) ;
     if (1) t4Show (HKHK) ;
   }

   /* construction directe du terme H4 = 1/8 (daij HH) gab (dbkl HH) */
   if (1)
     {
       int ii, jj, kk, ll ;
       CXs cc ;
       CX z, z0 = &cc, z1, z2, z3, z8 ;

       H4a = t4New (9, "ijkl", h) ;
       for (ii = 4 ; ii <= 7 ; ii++)
	 for (jj = 4 ; jj <= 7 ; jj++)
	   for (kk = 4 ; kk <= 7 ; kk++) 
	     for (ll = 4 ; ll <= 7 ; ll++)
	       {
		 z1 = cCopy(t3E (daij, 1, ii, jj), h) ;
		 cMM (z1, t3E (daij, 1, kk, ll)) ;
		 z2 = cCopy(t3E (daij, 2, ii, jj), h) ;
		 cMM (z2, t3E (daij, 2, kk, ll)) ;
		 z3 = cCopy(t3E (daij, 3, ii, jj), h) ;
		 cMM (z3, t3E (daij, 3, kk, ll)) ;
		 z8 = cCopy(t3E (daij, 8, ii, jj), h) ;
		 cMM (z8, t3E (daij, 8, kk, ll)) ;
		 z8->x *= 1 ; z8->y *= 1 ;
		 z0->x = z0->y = 0 ;
		 cAA (z0, z1) ;
		 cAA (z0, z2) ;
		 cAA (z0, z3) ;
		 cAA (z0, z8) ;

		 z = t4E(H4a,ii,jj,kk,ll) ;
		 cAA (z, z0) ;
	       }
       t4Scale (1.0/8.0, 0, H4a) ;
       t4FullSym (H4a) ;
     }
   if (1)
     {
       TT4 tt ;

       if (1)
	 {
	   printf ("16*H4 computed as  1/8(daij HH) gab (dbkl HH)\n") ;
	   tt = t4Copy (H4a, "ijkl", h) ;
	   t4Scale (16.0, 0, tt) ;
	   t4Show (tt) ;
	 }
       if (1)
	 {
	   printf ("verif H4 == H4a  (should be empty) \n") ;
	   tt = t4Copy (H4a, "ijkl", h) ;
	   t4Scale (-1, 0, tt) ;
	   t4AA (tt, H4) ;
	   t4Show (tt) ;
	 }
     }
   /* construction directe du terme H2K2 = -1/8(daij HH) super-gab (dbkl KK) */
   if (1)
     {
       int ii, jj, kk, ll ;
       CXs cc ;
       CX z, z0 = &cc, z1, z2, z3, z8 ;

       HHKKa = t4New (9, "ijkl", h) ;
       for (ii = 4 ; ii <= 7 ; ii++)
	 for (jj = 4 ; jj <= 7 ; jj++)
	   for (kk = 4 ; kk <= 7 ; kk++) 
	     for (ll = 4 ; ll <= 7 ; ll++)
	       {
		 z1 = cCopy(t3E (daij, 1, ii, jj), h) ;
		 cMM (z1, t3E (daij, 1, kk, ll)) ;
		 z2 = cCopy(t3E (daij, 2, ii, jj), h) ;
		 cMM (z2, t3E (daij, 2, kk, ll)) ;
		 z3 = cCopy(t3E (daij, 3, ii, jj), h) ;
		 cMM (z3, t3E (daij, 3, kk, ll)) ;
		 z8 = cCopy(t3E (daij, 8, ii, jj), h) ;
		 cMM (z8, t3E (daij, 8, kk, ll)) ;
		 z8->x *= -1 ; z8->y *= -1 ;
		 z0->x = z0->y = 0 ;
		 cAA (z0, z1) ;
		 cAA (z0, z2) ;
		 cAA (z0, z3) ;
		 cAA (z0, z8) ;

		 z = t4E(HHKKa,ii,jj,kk,ll) ;
		 cAA (z, z0) ;
	       }
       t4DoubleSym (HHKKa) ;
       t4Scale (-1.0/8.0,0,HHKKa) ;
     }
   if (1)
     {
       TT4 tt, tt1 ;

       if (1)
	 {
	   printf ("16*HHKK computed as  -1/8 (daij HH) super-gab (dbkl KK)\n") ;
	   tt = t4Copy (HHKKa, "ijkl", h) ;
	   t4Scale (16.0,0,tt) ;
	   t4Show (tt) ;
	 }
       
       if (1)
	 {
	   printf ("verif: HHKK == HHKKa   (should be empty)\n") ;
	   tt = t4Copy (HHKKa, "ijkl", h) ; 
	   t4Scale (-1,0,tt) ;
	   tt1 = t4Copy (HHKK, "ijkl", h) ; 
	   t4Scale (1.0,0,tt1) ;
	   t4AA (tt, tt1) ;
	   t4Show (tt) ;
	 }
     }

   /****************************************************************************/
   /******************* contre termes scalaires ********************************/
   /****************************************************************************/

   /* pseudo vertex fHHHH fKKKK fHHKK induit par le propagateur du vecteur */
   {
     CXs cc ;
     CX za, zf, z = &cc ;
     int a, i,j,k,l, n = H4->dim ;
     TT4 fHHKKa = t4New (n, "ijKL", h) ;
     TT4 fHHKK8 = t4New (n, "ijKL", h) ;
     fHHKK = t4New (n, "ijKL", h) ;

     for (i = 4 ; i < 8 ; i++)
       for (j = 4 ; j < 8 ; j++)
	 for (k = 4 ; k < 8 ; k++)
	   for (l = 4 ; l < 8 ; l++)
	  {
	    zf = t4E (fHHKKa, i, j, k, l) ;
	    for (a = 1 ; a <= 3 ; a++)
	      {
		z->x = z->y = 0 ;
		za = t3E(faij, a, i, j) ; *z = *za ;
		za = t3E(faij, a, k, l) ; cMM (z, za) ;
		cAA (zf, z) ;
	      }

	    zf = t4E (fHHKK8, i, j, k, l) ;
	    for (a = 8 ; a <= 8 ; a++)
	      {
		z->x = z->y = 0 ;
		za = t3E(faij, a, i, j) ; *z = *za ;
		za = t3E(faij, a, k, l) ; cMM (z, za) ;
		cAA (zf, z) ;
	      }
	  }
     if (1)
       {
	 printf ("Pseudo vertex fHHKKa") ;
	 if (1) t4Show (fHHKKa) ;
	 printf ("Pseudo vertex fHHKK8") ;
	 if (1) t4Show (fHHKK8) ;
       }
     t4Scale (1, 0, fHHKK8) ;
     fHHKK = t4Copy (fHHKKa, "hhkk", h) ;
     t4AA (fHHKK, fHHKK8) ;
     fHHHH = t4Copy (fHHKK, "HHHH", h) ;
     fKKKK = t4Copy (fHHKK, "HHHH", h) ;
   }


   /* contre terme H4  scalaires */
   if (0)
     {
       CXs cc ;
       CX za, z = &cc ;
       int i,j,k,l, n = H4->dim ;
       
       TT4 ZH4 = t4New (n, "ijkl", h) ;
       TT4 ZH4A = t4New (n, "ijkl", h) ;
       TT4 ZH4B = t4New (n, "ijkl", h) ;
       TT4 tt1 = t4t4doubleContract(H4,H4,"iJkL",3,4,1,2,h);
       TT4 tt2 = t4t4doubleContract(HHKK,KKHH,"iJkL",3,4,1,2,h);
       
       for (i = 0 ; i < n ; i++)
	 for (j = 0 ; j < n ; j++)
	   for (k = 0 ; k < n ; k++)
	   for (l = 0 ; l < n ; l++)
	     {
	       z = t4E(ZH4A, i,j,k,l) ;
	       za = t4E (tt1,i,j,k,l) ; cAA (z, za) ;
	       za = t4E (tt1,i,k,l,j) ; cAA (z, za) ;
	       za = t4E (tt1,i,l,j,k) ; cAA (z, za) ;
	       z->x /= 2.0 ; z->y /= 2.0 ;
	       
	       z = t4E(ZH4B, i,j,k,l) ;
	       za = t4E (tt2,i,j,k,l) ; cAA (z, za) ;
	       za = t4E (tt2,i,k,l,j) ; cAA (z, za) ;
	       za = t4E (tt2,i,l,j,k) ; cAA (z, za) ;
	       z->x /= 2.0 ; z->y /= 2.0 ;
	       
	       z = t4E(ZH4, i,j,k,l) ;
	       za = t4E(ZH4A, i,j,k,l) ; cAA (z, za) ;
	       za = t4E(ZH4B, i,j,k,l) ; cAA (z, za) ;
	     }
       
       if (1)
	 {
	   printf ("Contre terme a 4 scalaires ZH4A") ;
	   if (1) t4Show (ZH4A) ;
	   printf ("Contre terme a 4 scalaires ZH4B") ;
	   if (1) t4Show (ZH4B) ;
	   printf ("Contre terme a 4 scalaires ZH4") ;
	   if (1) t4Show (ZH4) ;
	 }
     }
   /* contre terme H2K2 scalaires */
   if (1)
     {
       CX za, z1, z2, z3 ;
       int i,j,k,l, n = H4->dim ;
       TT4 ZHHKKA = t4New (n, "ijKL", h) ;
       TT4 ZHHKKB = t4New (n, "ijKL", h) ;
       TT4 ZHHKKf = t4New (n, "ijKL", h) ;
       TT4 ZHHKK = t4New (n, "ijKL", h) ;
       TT4 tt1 = t4t4doubleContract(HHKK, K4,"iJkL",3,4,1,2,h) ;
       TT4 tt2 = t4t4doubleContract(H4, HHKK,"iJkL",3,4,1,2,h) ;
       TT4 tt3 = t4t4doubleContract(HKKH,HKKH,"iJkL",3,4,2,1,h) ;
       TT4 ttf1 = t4t4doubleContract(HHKK, fKKKK,"iJkL",3,4,3,2,h) ;
       TT4 ttf2 = t4t4doubleContract(fHHHH, HHKK,"iJkL",2,3,1,2,h) ;
       TT4 ttf3 = t4t4doubleContract(fHHHH, HHKK,"iJkL",2,3,2,3,h) ;
       TT4 ttf4 = t4t4doubleContract(fHHHH, HHKK,"iJkL",2,3,2,4,h) ;
       TT4 ttf5 = t4t4doubleContract(fHHHH, HHKK,"iJkL",2,3,1,3,h) ;
       TT4 ttf6 = t4t4doubleContract(fHHHH, HHKK,"iJkL",2,3,2,4,h) ;
       

       for (i = 4 ; i < 8 ; i++)
	 for (j = 4 ; j < 8 ; j++)
	   for (k = 4 ; k < 8 ; k++)
	     for (l = 4 ; l < 8 ; l++)
	       {
		 z1 = t4E(ZHHKKA, i,j,k,l) ;
		 za = t4E(tt1,i,j,k,l) ; cAA (z1, za) ;
		 za = t4E(tt2,i,j,k,l) ; cAA (z1, za) ;
		 z1->x /= 2.0 ; z1->y /= 2.0 ;
		 
		 z2 = t4E(ZHHKKB, i,j,k,l) ;
		 za = t4E (tt3,j,k,l,i) ; cAA (z2, za) ;
		 za = t4E (tt3,j,l,k,i) ; cAA (z2, za) ;
		 
		 z3 = t4E(ZHHKKf, i,j,k,l) ;
		 za = t4E (ttf1,j,k,l,i) ; cAA (z3, za) ;
		 if (0)
		   {
		     za = t4E (ttf2,j,k,l,i) ; cAA (z3, za) ;
		     za = t4E (ttf3,i,j,k,l) ; cAA (z3, za) ;
		     za = t4E (ttf4,i,j,k,l) ; cAA (z3, za) ;
		     za = t4E (ttf5,i,j,k,l) ; cAA (z3, za) ;
		     za = t4E (ttf6,i,j,k,l) ; cAA (z3, za) ;
		   }
		 z = t4E(ZHHKK, i,j,k,l) ;
		 za = t4E(ZHHKKA, i,j,k,l) ; cAA (z, za) ;
		 z->x *= -1.0 ; z->y *= -1.0 ;
		 za = t4E(ZHHKKB, i,j,k,l) ; cAA (z, za) ;
	       }
       
       if (1)
	 {
	   printf ("Contre terme a 4 scalaires ZHHKKA\n") ;
	   if (1) t4Show (ZHHKKA) ;
	   printf ("Contre terme a 4 scalaires ZHHKKB\n") ;
	   if (1) t4Show (ZHHKKB) ;
	   printf ("Contre terme a 4 scalaires ttf1\n") ;
	   if (1) t4Show (ttf1) ;
	   printf ("Contre terme a 4 scalaires ZHHKKf\n") ;
	   if (1) t4Show (ZHHKKf) ;
	   printf ("Contre terme a 4 scalaires ZH2K2\n") ;
	   if (1) t4Show (ZHHKK) ;
	 }
   }

   /* masse du top */
   if (0)
     {
       double mw, mz, mt, mt2, s2 ;
       printf("Masse du top\n") ;
       printf ("MW=80.403+-0.029; MZ=91.187+-0.002, 1-s2(theta)=m2w/m2Z->4 cases\n") ;
       mw = 80.403+0.029 ; mz=91.187+0.002;
       s2 = 1 - (mw/mz)*(mw/mz) ; printf ("\n%g", s2) ;
       mt = 4.0 * sqrt((1.0-(s2 - .25)/(s2 - .5))/3) * mw ;printf ("\t%g", mt) ;
       mt2 = 4.0 * sqrt((4*s2)/3) * mw ;printf ("\tmt2=%g", mt2) ;
       mw = 80.403-0.029 ; mz=91.187+0.002;
       s2 = 1 - (mw/mz)*(mw/mz) ; printf ("\n%g", s2) ;
       mt = 4.0 * sqrt((1.0-(s2 - .25)/(s2 - .5))/3) * mw ;printf ("\t%g", mt) ;
       mt2 = 4.0 * sqrt((4*s2)/3) * mw ;printf ("\tmt2=%g", mt2) ;
       mw = 80.403+0.029 ; mz=91.187-0.002;
       s2 = 1 - (mw/mz)*(mw/mz) ; printf ("\n%g", s2) ;
       mt = 4.0 * sqrt((1.0-(s2 - .25)/(s2 - .5))/3) * mw ;printf ("\t%g", mt) ;
       mt2 = 4.0 * sqrt((4*s2)/3) * mw ;printf ("\tmt2=%g", mt2) ;
       mw = 80.403-0.029 ; mz=91.187-0.002;
       s2 = 1 - (mw/mz)*(mw/mz) ; printf ("\n%g", s2) ;
       mt = 4.0 * sqrt((1.0-(s2 - .25)/(s2 - .5))/3) * mw ; printf ("\t%g", mt) ;
       mt2 = 4.0 * sqrt((4*s2)/3) * mw ;printf ("\tmt2=%g", mt2) ;
     }
   printf("\n\n") ;
   ac_free (h) ;
  return 0 ;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


