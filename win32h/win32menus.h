/* ACEMenu.h : header file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:  Contains WIN32 specific ACEDB menu data type declarations
 * HISTORY:
 * Last edited: Jun 11 00:40 1996 (rbrusk): normal looking freemenus
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
/* $Id: win32menus.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */

// ACEDB Menu Types
enum BoxMenuType { NO_MENU, NORMAL_MENU, FREE_MENU, NEWBOX_MENU  } ;

extern CString& MenuType( BoxMenuType boxMenuType ) ;

typedef CTypedPtrMap<CMapWordToPtr, WORD, void*> MENUMAP ;

typedef struct normalMenuBits 
{
	CMenu *menuObj;
	MENUOPT *opts;
	// menuMap maps WIN32 menu item id's to callbacks
	MENUMAP *menuMap ;
} NORMALMENUBITS ;

typedef struct freemenubits 
{
	CMenu *menuObj;
	FREEOPT *opts;
	FreeMenuFunction proc;
	MENUMAP *menuMap ;
} FREEMENUBITS ;

typedef struct newmenubits 
{
	CMenu *menuObj ;
	MENU m ;
	MENUMAP *menuMap ;
	// !isaMainMenu => isaSubmenu
	BOOL isaMainMenu ;
} NEWMENUBITS ;

 
 
 
