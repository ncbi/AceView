/*  File: tree.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for treedisp.c
 *              included in graphical acedb compilations
 * HISTORY:
 * Last edited: Nov 24 10:33 1998 (fw)
 * Created: Tue Nov 24 09:55:12 1998 (fw)
 *-------------------------------------------------------------------
 *
 * $Id: tree.h,v 1.3 2007/01/26 05:22:08 mieg Exp $
 */


#ifndef _TREE_H
#define _TREE_H

#include "acedb.h"

/************************************************************/

void treeUpdate (void) ;
BOOL treeDisplay (KEY key, KEY from, BOOL isOldGraph);
BOOL treeChooseTagFromModel(int *type, int *targetClass, int classe, 
			    KEY *tagp, Stack s, int continuation);

/************************************************************/

#endif /* _TREE_H */

/**************************** eof ***************************/
