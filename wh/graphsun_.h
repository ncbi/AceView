/* graphsun_.h
   version of graph.h for inclusion within graphics package in
     files that use SUNVIEW calls.
   contains device dependent private information and then 
     includes graph_.h, which declares device independent private 
     information and includes graph.h (public information).
*/

/* $Id: graphsun_.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

#include <suntool/sunview.h>
#include <suntool/canvas.h>
#include <sunwindow/pixwin.h>
#include <sunwindow/notify.h>
#include <pixrect/pr_line.h>
  
struct DevStruct
  { Frame     frame ;
    Canvas    canvas ;
    Pixwin    *pixwin ;	
    Frame     message ;
  } ;

#define DEV_DEFINED

#include "graph_.h"
 
/****** end of file ******/
