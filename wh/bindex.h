/*  File: bindex.h
 *  Author: Richard Durbin (Durbin@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: bindex.h,v 1.5 2015/09/18 22:13:50 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 10 15:14 2001 (edgrif)
 * * May  9 16:27 2000 (edgrif): Added symbolic type for call to bIndexInit()
 * Created: Thu Aug  7 15:11:29 1997 (Durbin)
 *-------------------------------------------------------------------
 */
#ifndef ACEDB_BINDEX_H
#define ACEDB_BINDEX_H

#include "bs.h"					    /* for OBJ type */

/* Used to specify type of reindexing by bIndexInit():                       */
/*                                                                           */
/* BINDEX_AUTO - used by sesssion.c, indicates automatic reindexing, happens */
/*               only if ACEDB_NO_INDEX not set, indices are out of date,    */
/*               or a reinitialised db and user wants to reindex..phew.      */
/* BINDEX_FORCE - forces reindexing, only used by bindex.c currently.           */
/* BINDEX_AFTER_READMODELS - used by mainpick.c & command.c, this is user    */
/*                           requested full reindexing.                      */
/*                                                                           */
typedef enum _BindexInitType {BINDEX_AUTO, 
			      BINDEX_FORCE = 2, BINDEX_AFTER_READMODELS = 3} BindexInitType ;

/* bIndexFind() can return one of three results, it is NOT SAFE to assume    */
/* that the tag is present if non-zero is returned.                          */
/* NOTE, I explicitly define the values because almost certainly there is    */
/* code in acedb which relies on these numerical values.                     */
typedef enum _BindexFindResult {BINDEX_TAG_ABSENT = 0, BINDEX_TAG_UNCLEAR = 1,
				BINDEX_TAG_PRESENT = 2} BindexFindResult ;

BindexFindResult bIndexFind(KEY key, KEY tag) ;		    /* the one the users want */
int bIndexVersion (int parent_version) ;		    /* used in bootstrap */
void bIndexNewModels(void) ;				    /* when you get new models */
void bIndexInit(BindexInitType ask_reindex) ;		    /* Redo indexes, called from sessionInit
							       and others. */
BOOL bIndexSave(void) ;					    /* from sessionInit */
int bIndexStatus (int *nTablep, int *nKbp) ;		    /* used by status report */

void bIndexObject (KEY key, OBJ obj) ;
unsigned long int bIndexTagCount (int cl, KEY tag) ;
void bIndexShutDown (void) ;

/*
 * These functions find tags etc. making optimal use of the index,
 * only opening objects to get data if strictly necessary.
 */
BOOL bIndexTag(KEY key, KEY tag) ;			    /* Is a tag in an object ? */
BOOL bIndexGetKey(KEY key, KEY tag, KEY *key_out) ;	    /* Retrieve key following a tag if
							       present. */
BOOL bIndexGetTag2Key(KEY key, KEY tag, KEY *key_out) ;	    /* Retrieve key following tag but
							       where tag is part of a tag2 system. */
#endif /* !ACEDB_BINDEX_H */

/************************** eof ************************************/
 
 
