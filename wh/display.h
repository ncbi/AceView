/*  File: display.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: header for display.c, 
 *              a module that drives acedb windowing displays
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 09:49 1998 (fw)
 * * Sep 17 16:34 1998 (edgrif): Removed displayCreate(int d) definition,
 *              should be redundant now.
 * Created: Thu Sep 17 16:34:52 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: display.h,v 1.6 2013/10/08 02:37:03 mieg Exp $ */

#ifndef _DISPLAY_H
#define _DISPLAY_H

#include "acedb.h"

#include "../whooks/disptype.h"
#include "menu.h"		/* new libfree menu package */

/************************************/
typedef struct acedbDisplayStuct  AcedbDisplay ;
struct acedbDisplayStuct 
  { 
    char title[42], help[32] ;
    float x, y, width, height ;
    int type ;  /* from GraphType */
    } ;
/************************************/

extern DisplayFunc displayFunc[] ;
extern ParseFunc parseFunc[] ;
extern DumpFunc dumpFunc[] ;

#ifndef XLIB_DISPLAY
BOOL	display (KEY key, KEY from, char *displayName) ;
#endif

Graph   displayCreate(char *displayName) ;

void	displayPreserve (void) ;

BOOL	displayBlock (BlockFunc func, char * message) ;	/* call func on next picked key */
BOOL displayMaybeBlock(BlockFunc func, char *message, 
		       float height, BOOL *enabled);
void    displayUnBlock (void) ;
void    displayMaybeUnblock(void);
void    displayRepeatBlock (void) ;
BOOL	isDisplayBlocked (void) ;

void*   keySetNewDisplay (KEYSET ks, const char *title) ; /* display keyset in new window */

void*	keySetShow (KEYSET set, void* handle) ; /* h can be zero, or previous return */
void*   keySetMessageShow (KEYSET set, void* handle, char *message) ;
char*   keySetCurrent (void* handle, KEY k) ;
void    newKeySet (char *title); /* create a new keyset list display */

void pickDefaultDisplays (void);
FREEOPT *pickGetDisplayMenu(void);
void pickGetDisplayTypes(void);
void pickSetGraphTypes (void);
char *pickMainTitle(void);
void pickSetDisplaySize (const char *cp, float x, float y, float w, float h);

/************ from action.c **********/
/***** should they be here ???????? */

void externalDisplay (KEY key) ;
void acedbMailer (KEY key, KEYSET ks, Stack sText);

#endif /* _DISPLAY_H */
