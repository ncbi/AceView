/*  Last edited: Mar 24 12:46 1993 (mieg) */

/* $Id: graphdos_.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */
/* device specific definitions - must include key definitions */

#ifdef __MSDOS__        /* TurboC specific - uses Borland graphics */
#define __COLORS      /* prevents TC definitions of BLACK, WHITE, RED etc */
#include <graphics.h>
  struct DevStruct
    { int         left,top,right,bottom ;
      unsigned int       size ;
      int         isRetained ;
      void        *map ;
    } ;
#endif


#define DEV_DEFINED
