/*  File: bs_.h
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: fused from 3 headers used in all the same files
 		private BS stuff, plus interaction between BS and cache
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 16 11:34 2003 (edgrif)
 * Created: Tue Mar 10 03:40:07 1992 (mieg)
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_BS__H
#define ACEDB_BS__H

/* $Id: bs_.h,v 1.9 2016/01/21 22:39:14 mieg Exp $ */

#include "acedb.h"

       /* CAVEAT : If sizeof(BSdata) changes data on disk
                   can no longer be read , see b_.h */

/* model data keys are all tags, so use table bits for flags */
#define UNIQUE_BIT KEYMAKE(1,0)
#define DELETE_BIT KEYMAKE(2,0)
#define CHECK_BIT KEYMAKE(4,0)
#define COORD_BIT KEYMAKE(8,0)


/* CAVEAT : If sizeof(BSdata) changes, data on disk can no longer be read , see b_.h */
typedef KEY TEXTKEY ;

/*
* Note: This union is copied directly into the disk file, so you can't change
* it without breaking existing databases.
*/
typedef union  { int i;
                 unsigned int u;
#ifndef ACEDB4
                 long L; /* in ACEDB4, we suppress the longs to gain space
			  *  but can t do that before else we cannot read the data
                          *  anymore
			  */
                 unsigned long U;
#endif
                 float f;
		 mytime_t time ;
                 TEXTKEY text; /* to dereference a text entry */
		 KEY key ;    /* used in models */
               } BSdata;

typedef struct bshow *BS;
typedef struct btext *BT;

/*
* bshow is a single tag or data item that is stored in an object.  They
* are arranged in a tree.
*
*/
        struct bshow {
		      BS  up, 
				/* XQX
				* up can be one of two things:
				* - points at the node before us at this level
				* - points at the node "above" us at the next
				*   higher level
				* Apparently, the code does not distinguish
				* between these two uses.  That is, as you
				* follow "up", you first pass over some sibling
				* nodes, then into a parent node, then into
				* siblibings of the parent, etc.
				*/
			  down,
				/* XQX
				* down is another tag or data item at the
				* same (!) level of the object tree
				*/
			  right;
				/* XQX
				* right is a more deeply nested tag or data
				* item in the object tree
				*/
                      BT bt ;
                      BSdata n ;
				/* XQX
				* If key is one of the magic tags that
				* indicates a data type (_Int, etc), then n
				* contains the actual data.  Otherwise, n
				* is unused.
				*
				* n is copied directly to the disk - if you
				* change it's size, the database will be
				* incompatible.
				*/
                      KEY key ;
				/* XQX
				* key is the numeric representation of this
				* tag.  It can also be one of the special tags
				* (_Int, _Float) that indicates that n contains
				* a data element.
				*/
                      int size ;  /* used in BSstore */
		      KEY timeStamp ;
 		      int tagCount ;
				/* XQX
				* locked implies that this object cannot be 
				* modified (new feature under development)
				*/
                     } ;

        struct btext { char *cp ; /* holds a copy of text, do free */
		       BS   bsm ; /* only used during treedisp update */
                     } ;

/* flag fields for overloading of bs->size as bs->flag */

#define MODEL_FLAG	0x0001 /* part of model (painted cyan) */
#define UNIQUE_FLAG	0x0002 /* unique part of original tree */
#define ADD_DATA_FLAG	0x0004 /* new node - text entry */
#define ADD_COMMENT_FLAG 0x0008 /* new node - text entry */
#define ADD_KEY_FLAG	0x0010 /* new node = text entry and displayBlock */
#define EDIT_NODE_FLAG  0x0020 /* old node = text entry */
#define SUBTYPE_FLAG    0x0040 /* adds subtype to displayed model */
#define ON_FLAG		0x0080 /* if on dont display branch */
#define ATTACH_FLAG	0x0100 /* this tag is attached from another object */
#define NEW_FLAG	0x0200 /* new this session */
                               /* all waiting operations */
#define WAITING_FLAGS   (ADD_DATA_FLAG | ADD_COMMENT_FLAG | \
                         ADD_KEY_FLAG | EDIT_NODE_FLAG)

/********* cache stuff ********/

#ifndef DEF_CACHE_HANDLE
typedef BS CACHE_HANDLE ;
#define DEF_CACHE_HANDLE
#endif

#ifndef DEF_CACHE
#define DEF_CACHE
  typedef void* CACHE ;
#endif

#ifndef DEFINE_OBJ
#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#endif

struct sobj { KEY key ;
	      int magic ;
              CACHE x ;
	      KEY localModel ;
	      BS
		root,
		localRoot,
		curr, 	/* XQX a cursor pointing into the object data */
		modCurr;/* XQX a cursor pointing into the model */
	      Stack xref ;
              unsigned int flag ;
	    } ;

#define OBJFLAG_TOUCHED		0x0001

/*******************************/

#define bsText(bs) ((bs) && (bs)->bt && (bs)->bt->cp ? (bs)->bt->cp : 0)

char* getBStext (BS bs, BOOL isModel, BOOL dontExpand);

BS BSalloc (void) ;
void BSfree (BS p) ;

BT BTalloc (void) ;
void BTfree (BT p) ;

#include "bstree.h"

#define bsIsTag(bs)	((bs) && !class((bs)->key) && (bs)->key > _LastN)

#if !defined(NEW_MODELS)
#define bsIsType(bs)	((bs) && class((bs)->key) && \
			 (bs)->key == KEYMAKE(class((bs)->key), 1))
#else
#define bsIsType(bs)    ((bs) && class((bs)->key) == _VModel)
#endif

#define bsIsComment(bs)	((bs) && (class((bs)->key) == _VComment || \
			 class((bs)->key) == _VUserSession ))


void bsFuseModel (OBJ obj, KEY tag1, int justKnownTags) ;/* fuses in model for tree display updater */
void bsMakePaths (int table) ;
BOOL bsModelMatch (BS bs, BS *bsmp) ;
BS bsModelRight (BS bsm) ;
BS bsModelRoot (OBJ obj) ;
BS bsTransferBranch (OBJ newObj, OBJ obj, BOOL whole) ;

void cacheMark (CACHE cache) ; /* Call after updating to ensure it will go to disk */
BOOL isCacheModified (CACHE cache) ;
BOOL isCacheLocked (CACHE cache) ;

#endif
/******************************************************************/
