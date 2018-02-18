/*  Last edited: Jan 28 16:28 1994 (srk) */
/*  $Id: next.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */
#ifdef NEXT

int _Xdebug = 0  ;
_XQEvent * _qfree = NULL ;
Display *_XHeadOfDisplayList = 0 ;

 int (*_XErrorFunction)() ;

 int (*_XExtensionErrorFunction)() ;
 int (*_XIOErrorFunction)() ;

#endif
