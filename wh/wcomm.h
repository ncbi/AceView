/*  File: wcomm.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:  definitions for worm communications
 *              NOTE: "WComm_v1" flags the old way of doing things (objects are ints).
 *              just for temporary use; the default (WComm_v1 not defined) is the
 *              new way: objects are strings
 *
 *              Sometime all the "NeedFunctionPrototypes" need removing,
 *              after all, if we are not doing ANSI C it's all hopeless.
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 13 14:23 1998 (edgrif)
 * * Sep 30 16:32 1998 (edgrif): Added function prototypes for init/finish
 *              wcomm stuff and added this standard format header.
 * * Aug 03 00:59 1994 (rd)
 * * Jul 15 16:32 1991 (?????): t.friedman 7/15/91 - first major overhaul
 * Created: Wed Sep 30 16:30:16 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: wcomm.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */


#include <X11/Xlib.h>


#if defined(__cplusplus) || defined(c_plusplus)
#define C_CALL	extern "C"	/* c++ linkage information */
#else
#define C_CALL		/* nothing needed for normal C */
#endif /* c++ external definitions */

typedef enum {
		ImHere,			/* maintenance function, no message */
		ImLeaving,		/* ditto */	
			    /* room for other functions */
		MaintainLink=5,		/* maintenance function, no message */
		MakeAnnotation,		/* to wcs: create an annotation */
		ShowAnnotation,		/* to wcs: display one */
		ShowGenMap,		/* to ace: display genetic map */
		ShowObject,		/* to ace: display tree */
		ShowPhysMap,		/* to ace: display physical map */
		ShowSummary,		/* to ace: not implemented yet */
		SearchString		/* to wcs: search on a string */
	} WMsgType;

typedef enum {unset=0, readonly, readwrite, writeonly} WCommType;

typedef struct {
#ifdef WComm_v1
	long	obj;	/* old internal reference */
#else
	char	*obj;	/* co-program internal reference */
#endif
	char	*note;	/* various uses */
	    /* either vbl can be void, depending on message type */
    } WHandle;

typedef struct wcb {
		Atom	WC_atom;		/* gets X atom */
		Window	WC_window; 
		Display	*WC_display;
		char	*WC_name;
		int	WC_name_len;
		WCommType  WC_comm_type;
		struct wcb  *WC_nextblock;
		void	(*WC_handler)();	/* routine to handle data */
       } WCommBlock;
typedef	WCommBlock *	WCommLink;

/* #define	GOT_WCOMM	1		temporary, while batch testing */
/*#define	WCDefaultLink  (WCommLink) 0  until we build up */
#define	WCNoLink  (WCommLink) -1  	/*  no usable external program */

/* subroutine return type definitions */

C_CALL void wcommInit (Display *d, Window w) ;
C_CALL void wcommFinish (void) ;


C_CALL void	WCleanup ();
C_CALL void    WFree (
#if NeedFunctionPrototypes
    WHandle     *h_buffer
#endif
);

C_CALL Bool WPropMessage (
#if NeedFunctionPrototypes
    XEvent *xevt
#endif
);

#ifdef WCOMM_INTERN
C_CALL Bool    WGetMsg (
#if NeedFunctionPrototypes
    XPropertyEvent *myPevent,   /* redefine myevent structure */
    WCommLink   *who_from,
    WMsgType    *msg_type,
    int         *num_handles,
    WHandle     **h_buffer
#endif
);
#endif	/*WCOMM_INTERN*/

C_CALL Bool    WInitComm (
#if NeedFunctionPrototypes
    char    *my_name,       /* my program name */
    Display *zdisplay, /* if non-zero, use the Disp & Window given */
    Window  zwindow,
    WCommType mycomm,	/* what comm my program will do */
    void    (* handler)()   /* default routine to handle incoming data */
#endif
);

C_CALL int WSendMsg (
#if NeedFunctionPrototypes
    WCommLink   who_to,     /* ignored for now */
    WMsgType    msg_type,
    int         num_handles,
    WHandle     *h_buff
#endif
);

C_CALL int  W_show_external_reference (
#if NeedFunctionPrototypes
    WCommLink whoto,   /* ignored for now */
#ifdef WComm_v1 
    int ext_obj,
#else 
    char *ext_obj,
#endif 
    char *name
#endif
);

C_CALL WCommLink  WOpenCommLink (
#if NeedFunctionPrototypes
    char * whoto,   /* moniker of external program */
    WCommType comm_type,
    void    (* handler)()   /* routine to handle incoming data */
#endif
);

  int	WEndComm  (WCommLink whowith) ;
  void WCleanup () ;
/*
  void wcsAccept (WMsgType, int, WHandle*, WCommLink) ;
*/

