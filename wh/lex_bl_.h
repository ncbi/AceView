/*  File: lex_bl_.h
 *  Author: Jean Thierry-Mieg
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for common structures 
 *              shared between lexer and block subs
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 23 18:03 1998 (fw)
 * Created: Mon Nov 23 18:01:16 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: lex_bl_.h,v 1.7 2015/08/17 16:12:51 mieg Exp $ */

 
/**********************************************************************/

#ifndef _LEX_BL__H
#define _LEX_BL__H
 
void lexSetDisk(KEY key, DISK dk) ; 

#include "disk__.h"   /* defines BP and DISK */

typedef struct alloc *ALLOCP;    /* control structure of the cache */
struct alloc        { BP p;
                      int ispinned;
                      int ismodified;
		      KEY oneofthekeys ;
                      ALLOCP up,next,right ;
                    }
                      ; 
 
typedef unsigned int LEXOFFST;    /*offset in the vocabulary */
 
typedef struct lexalpha 
        {KEY key;       /*a provision if we want to shift*/
                        /*some of the vocabulary to disk */
        } LEXA, *LEXAP;
 
typedef union { DISK dk ; int lx ; KEY alias ;} DISK_LEX2 ;

/*
* There is a LEXI1 struct for every object in the database.  Each class
* has an array of LEXI1 structs (see Lexi1 in "w5/lexsubs4.c").  That
* array is indexed by the key of the object with the class number masked
* off.  The result is the LEXI1 for a particular object.
*
*/
typedef struct lexique1   
        {
	 KEY  nameoffset;  /*name of the object*/
		/*
		* class type A,B
		* nameoffset is the offset into the name string.  See
		* Voc in "w5/lexsubs4.c".  It is really an integer, not
		* a KEY.
		*
		* class D
		* nameoffset is the parent key
		*/
	 DISK_LEX2 dlx ;        /* its disk address, also the alias key and the lex2 index */
	 unsigned char lock;   /* flow control */
	 unsigned char isMask ;   /* class hierarchy */
	 unsigned char is2 ;   /* the DIK_LEX union is a lexi2 */	 
        } LEXI1, *LEXP1 ;
 
#define lexiFormat "kkccc"

typedef struct lexique2    
       { DISK_LEX2 dlx ;        /* its disk address, also the alias key and the lex2 index */
	 ALLOCP addr;      /* addr->p is its memory address*/
	 void * cache;       /* address in the secondary cache */
         KEY key ;      /* the owning key */
       } LEXI2, *LEXP;   /* LEXP as used in everything but in lexsubs.c */

typedef struct {
  KEY created;
  KEY updated;
} TIMESTAMP;

#define timeStampFormat "kk"

#define SIZELEX sizeof(struct lexique)
 
 
DISK lexDisk(KEY key) ;
void lexSetDisk(KEY key, DISK dk);
LEXP KEY2LEX(KEY kk); /* returns q if(iskey(kk)), 0 otherwise */

#endif /* _LEX_BL__H */
 


