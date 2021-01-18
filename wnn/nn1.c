#include <ac.h>
#include <matrix.h>

int main (int argc, char *argv[])
{
  AC_HANDLE h = ac_new_handle () ;

  MX a = mxCreate (h, "a", MX_INT, 4, 2, 0) ;
  MX at = mxCreate (h, "at", MX_INT, 2, 4, 0) ;
  MX b = mxCreate (h, "b", MX_INT, 4, 2, 0) ;
  MX c = mxCreate (h, "c", MX_INT, 4, 2, 0) ;
  int ii[] = {1,2,3,4,5,6,7,8} ;
  int jj[] = {1,1,1,1,2,2,2,2} ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  mxSet (b, ii) ;
  mxSet (c, jj) ;

  mxAdd (a, b, c, 0) ;

  mxShow (b) ;
  mxShow (c) ;
  mxShow (a) ;
  mxTranspose (at, a, 0,1, 0) ;
  mxShow (at) ;
  mxShow (mxTranspose (0,b,0,1,0)) ;
  return 0 ;
}
